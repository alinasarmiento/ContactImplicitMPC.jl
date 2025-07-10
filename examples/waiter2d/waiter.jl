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
using Sockets

# ## Simulation
s = get_simulation("waiter_2D", "flat_2D_lc", "flat");
model = s.model
env = s.env

# ## Reference Trajectory
h = 0.01
H = 100
ref_traj = contact_trajectory(model, env, H, h)
ref_traj.h

# qref = [0.5; 0.485;
#         0.5; 0.5; 0.0;]
qref = [0.5; 0.42;
        0.65; 0.4831; 0.0;]

ur = zeros(model.nu)
γr = zeros(model.nc)
br = zeros(model.nc * friction_dim(env))
ψr = zeros(model.nc)
ηr = zeros(model.nc * friction_dim(env))
wr = zeros(model.nw)

# ## Set Reference
for t = 1:H
    ref_traj.z[t] = pack_z(model, env, qref, γr, br, ψr, ηr)
    ref_traj.θ[t] = pack_θ(model, qref, qref, ur, wr, model.μ_world, ref_traj.h)
end

# ## Initial conditions
# q0 = ContactImplicitMPC.SVector{2}([0.0 * π, 0.0])
# for instantiation BEFORE controller created

q1 = [0.5; 0.42;
      0.65; 0.4831; 0.0;] #0.483
# q1 = [0.5; 0.6;
#         0.5; 0.616; 0.0;]
v1 = [0.0; 2.0;
      0.0; 0.0; 0.0;]

# ## Simulator
sim = simulator(s, H, h=h)

# ## Simulate -- simulates entire trajectory?
status = simulate!(sim, q1, v1, verbose=true)

##########################
# Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.render(vis)

# ## Visualize
vis_traj = contact_trajectory(s.model, s.env, H, h)
anim = visualize_robot!(vis, model, sim.traj, sample = 1, h=h)
@infiltrate
##########################

# ## MPC setup 
N_sample = 2
H_mpc = 40
h_sim = h / N_sample
H_sim = 200
κ_mpc = 2.0e-4

## Cost
q_scale = 1e-1
q_vec = q_scale .* [1., 1., 1., 1., 1.] # x_ee, z_ee, x_tray, z_tray, θ_tray

v_scale = 1e-1
v_vec = v_scale .* [1., 1., 1., 1., 1.] # x_ee, z_ee, x_tray, z_tray, θ_tray

u_scale = 10
u_vec = [1, .01]

print("creating objective\n")
obj = TrackingVelocityObjective(model, env, H_mpc,
	q = [Diagonal(q_vec .* (t/H_mpc)^2) for t = 1:H_mpc-0],
	v = [Diagonal(v_vec ./ (h^2.0)) for t = 1:H_mpc-0],
	u = [Diagonal(u_vec) for t = 1:H_mpc-0],
	γ = [Diagonal(1.0e-100 * ones(model.nc)) for t = 1:H_mpc-0],
	b = [Diagonal(1.0e-100 * ones(model.nc * friction_dim(env))) for t = 1:H_mpc]);

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


# ## Initial Conditions
q1_sim = q1
v1_sim = v1

print("simulation\n")
# ## Simulator
# sim_ip_opts = InteriorPointOptions(
#     undercut = 1.0,
#     κ_tol = κ_mpc,
#     r_tol = 1.0e-8,
#     diff_sol = true,
#     max_time = 1e5,)

sim = simulator(s, H_sim, h=h_sim, policy=p) #, solver_opts=sim_ip_opts) #, dist=d)

# ## Simulate
status = simulate!(sim, q1_sim, v1_sim, verbose=true)

# ## Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.render(vis)

# ## Visualize
vis_traj = contact_trajectory(s.model, s.env, H_sim, h_sim)
anim = visualize_robot!(vis, model, sim.traj, sample = 1, h=h_sim)
@infiltrate

# ## Timing result

# Julia is [JIT-ed](https://en.wikipedia.org/wiki/Just-in-time_compilation) so re-run the MPC setup through Simulate for correct timing results.
process!(sim.stats, N_sample) # Time budget
H_sim * h_sim / sum(sim.stats.policy_time) # Speed ratio

u_mat = mapreduce(permutedims, vcat, sim.traj.u)
v_mat = mapreduce(permutedims, vcat, sim.traj.v)
q_mat = mapreduce(permutedims, vcat, sim.traj.q)
force_mat = mapreduce(permutedims, vcat, sim.traj.γ)
tanf_mat = mapreduce(permutedims, vcat, sim.traj.b)

