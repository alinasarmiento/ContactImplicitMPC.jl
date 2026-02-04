# PREAMBLE

# PKG_SETUP

# ## Setup
 
using ContactImplicitMPC
using RoboDojo
import ContactImplicitMPC: simulate!
using LinearAlgebra
using Infiltrator
using DelimitedFiles
using MeshCat
using YAML
using Sockets

function get_yaml()
    ctrl_params = YAML.load_file(joinpath(@__DIR__,"block_costs.yaml"))
    return ctrl_params
end

# ## Simulation
s = get_simulation("block_1D", "flat_2D_lc", "flat", model_variable_name="block_system_1D");
model = s.model
env = s.env

# ## Reference Trajectory
ctrl_params = get_yaml()
h = ctrl_params["ctrl_dt"]
H = 100 #100
ref_traj = contact_trajectory(model, env, H, h)
ref_traj.h

ee_init = deepcopy(ctrl_params["ee_init"])
ee_des = deepcopy(ctrl_params["ee_des"])
block_init = deepcopy(ctrl_params["block_init"])
block_des = deepcopy(ctrl_params["block_des"])

qref = [ee_init[1];       # ee [x]
        block_init[1]; block_init[2]; 0.0;] # block [x,z,th]

uref = [0];
f_Nee = -uref[1]
f_Ng = (model.m_block*9.81)
normal_ref = [f_Nee, f_Ng/2, f_Ng/2];
fric_ref = [0, 0, -(f_Ng/2)*model.μ_world, 0, -(f_Ng/2)*model.μ_world,0]

# ur = zeros(model.nu) #ones(model.nu)
# γr = zeros(model.nc)
# br = zeros(model.nc * friction_dim(env))
ur = ones(model.nu).*uref*h
γr = ones(model.nc).*normal_ref*h
br = ones(model.nc*friction_dim(env)).*fric_ref*h
ψr = zeros(model.nc)
ηr = zeros(model.nc * friction_dim(env))
wr = zeros(model.nw)

# ## Set Reference
function create_block_push_ref(dist, vel, q_init, dt, T)
    """
    dist: distance between EE and block when in contact
    vel: velocity for EE/block
    T: number of steps
    """
    q_traj = []
    u_traj = []
    push!(q_traj, q_init)
    push!(u_traj, [0])
    q_t = q_init
    for i=1:T-1
        bc = (q_t[2] - q_t[1] <= dist)
        push!(q_traj, q_t + [dt*vel; bc*dt*vel; 0; 0])
        mb = model.m_block*bc
        ux = [(block_xvel*(mb+model.m)/(i*dt)) + (model.μ_world*mb*9.81)] * dt
        push!(u_traj, ux)
        q_t = q_traj[end]
    end
    return q_traj, u_traj
end

block_xvel = ctrl_params["block_vel_des"];
q_t = [ee_init[1]; block_init[1]; block_des[2]; 0.0;]
q_goal = [ee_init[1]; block_des[1]; block_des[2]; 0.0;]
qref_traj, uref_traj= create_block_push_ref(((model.xlen_block/2)+model.r), block_xvel, q_t, h, H+2)
    
for t = 1:H
    ref_traj.z[t] = pack_z(model, env, qref, γr, br, ψr, ηr)
    ref_traj.θ[t] = pack_θ(model, qref, qref, ur, wr, model.μ_world, ref_traj.h)
    ref_traj.q[t] = q_goal #qref_traj[t] #qref + [t*h*block_xvel; t*h*block_xvel; 0; 0]
    ref_traj.u[t] = [0] #uref_traj[t]  #ur
    ref_traj.γ[t] = [0, γr[2], γr[3]] #[-uref_traj[t][1], γr[2], γr[3]] #γr
    ref_traj.b[t] = br
end
ref_traj.q[H+1] = qref_traj[H+1]
ref_traj.q[H+2] = qref_traj[H+2]
update_friction_coefficient!(ref_traj, model, env)

# ## Initial conditions

q1 = [ee_init[1];       # ee [x,z]
        block_init[1]; block_init[2]; 0.0;] # block [x,z,th]
v1 = [0.0;
      0.0; 0.0; 0.0;]

# ## Simulator
sim = simulator(s, H, h=0.001)

# ## Simulate -- simulates entire trajectory?
status = simulate!(sim, q1, v1, verbose=true)

##########################
#Visualizer
# vis = ContactImplicitMPC.Visualizer()
# ContactImplicitMPC.render(vis)

# ## Visualize
# vis_traj = contact_trajectory(s.model, s.env, H, h)
# anim = visualize_robot!(vis, model, sim.traj, sample = 1, h=h)
# # u_mat = mapreduce(permutedims, vcat, sim.traj.u)
# v_mat = mapreduce(permutedims, vcat, sim.traj.v)
# q_mat = mapreduce(permutedims, vcat, sim.traj.q)
# force_mat = mapreduce(permutedims, vcat, sim.traj.γ)
# tanf_mat = mapreduce(permutedims, vcat, sim.traj.b)

# # @infiltrate
##########################


function obj_policy(model, env, ctrl_params, N_sample)
    H_mpc = ctrl_params["H_mpc"]
    h_sim = ctrl_params["ctrl_dt"] / N_sample
    H_sim = 1000
    κ_mpc = ctrl_params["kappa"] #0.7e-3
    ## Cost
    q_scale = deepcopy(ctrl_params["q_scale"])
    q_vec = q_scale .* ctrl_params["q_vec"]
    v_scale = deepcopy(ctrl_params["v_scale"])
    v_vec = v_scale .* ctrl_params["v_vec"]
    u_scale = deepcopy(ctrl_params["u_scale"])
    u_vec = u_scale .* ctrl_params["u_vec"]
    qlim_scale = deepcopy(ctrl_params["qlim_scale"])
    qlim_vec = qlim_scale .* ctrl_params["qlim_vec"]
    ulim_scale = deepcopy(ctrl_params["ulim_scale"])
    ulim_vec = ulim_scale .* ctrl_params["ulim_vec"]
    print("creating objective\n")
    q_obj_vec = [(t/H_mpc)^2 for t=1:H_mpc-0]
    q_obj_vec[H_mpc] = 1
    println(q_obj_vec)
    q_obj = [Diagonal(q_vec) .* q_obj_vec[t] for t=1:H_mpc-0]
    println(q_obj)
    obj = TrackingVelocityObjective(model, env, H_mpc,
                                    q = q_obj, #[Diagonal(q_vec) .* (t/H_mpc) for t = 1:H_mpc-0],
                                    v = [Diagonal(v_vec) .* (t/H_mpc) for t = 1:H_mpc-0],
                                    # q = [Diagonal(q_vec) for t = 1:H_mpc-0],
                                    # v = [Diagonal(v_vec) for t = 1:H_mpc-0],
                                    u = [Diagonal(u_vec) for t = 1:H_mpc-0],
                                    γ = [Diagonal(1.0e-100 * ones(model.nc)) for t = 1:H_mpc-0],
                                    b = [Diagonal(1.0e-100 * ones(model.nc * friction_dim(env))) for t = 1:H_mpc],
                                    qlim = [Diagonal(qlim_vec) for t = 1:H_mpc-0],
                                    ulim = [Diagonal(ulim_vec) for t = 1:H_mpc-0] );
    # ## Policy
    print("policy\n")
    p = ci_mpc_policy(ref_traj, s, obj,
        H_mpc = H_mpc,
        N_sample = N_sample,
        κ_mpc = κ_mpc,
        # ip_opts = InteriorPointOptions(
        #                       undercut = 1.0,
        #                       κ_tol = κ_mpc,
        #                       r_tol = 1.0e-8,
        #                       diff_sol = true,
        #                       solver = :empty_solver,
        #                       max_time = 1e5,),
        n_opts = NewtonOptions(
                    r_tol = 3e-4,
                    max_iter = 10,
                    max_time = ref_traj.h/2, # HARD REAL TIME
                    ),
                      mpc_opts = CIMPCOptions());
    return obj, p
end


# ## MPC setup 
N_sample = 5
obj, p = obj_policy(model, env, ctrl_params, N_sample)
κ_sim = ctrl_params["kappa_sim"] #0.7e-3
h_sim = ctrl_params["ctrl_dt"] / N_sample
H_sim = 1000
# ## Initial Conditions
q1_sim = q1
v1_sim = v1

print("simulation\n")

# ## Simulator
sim_ip_opts = InteriorPointOptions(
    undercut = 1.0,
    κ_tol = κ_sim,
    r_tol = 1.0e-8,
    diff_sol = true,
    max_time = 1e5,)

sim = simulator(s, H_sim, h=h_sim, policy=p)
# sim = simulator(s, H_sim, h=h_sim, policy=p, solver_opts=sim_ip_opts) #, dist=d)

# ## Simulate
status = simulate!(sim, q1_sim, v1_sim, verbose=true)

# ## Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.render(vis)

# ## Visualize
vis_traj = contact_trajectory(s.model, s.env, H_sim, h_sim)
anim = visualize_robot!(vis, model, sim.traj, sample = 1, h=h_sim)
# @infiltrate

using StaticArrays
function find_JTf(mod, envi, q, f, b, idx)
    cf = []
    for i in 1:3
        cf = push!(cf, b[idx, 2*i-1]+b[idx, 2*i])
        cf = push!(cf, f[idx,i])
    end
    ff = SVector{6}(cf)
    return transpose(ContactImplicitMPC.J_func(mod, envi, q[idx-1,1:end])) * ff
end

# ## Timing result

# Julia is [JIT-ed](https://en.wikipedia.org/wiki/Just-in-time_compilation) so re-run the MPC setup through Simulate for correct timing results.
process!(sim.stats, N_sample) # Time budget
H_sim * h_sim / sum(sim.stats.policy_time) # Speed ratio

u_mat = mapreduce(permutedims, vcat, sim.traj.u)
v_mat = mapreduce(permutedims, vcat, sim.traj.v)
q_mat = mapreduce(permutedims, vcat, sim.traj.q)
force_mat = mapreduce(permutedims, vcat, sim.traj.γ)
tanf_mat = mapreduce(permutedims, vcat, sim.traj.b)

# writedlm("usol_fcimpc.csv",u_mat/h_sim,',')
