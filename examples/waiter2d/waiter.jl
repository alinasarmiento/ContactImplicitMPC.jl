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
h = 0.005
H = 100
ref_traj = contact_trajectory(model, env, H, h)
ref_traj.h

qref = [0.5; 0.484;
        0.5; 0.5; 0.0;]
ur = [0.0; 0.37*9.81] #zeros(model.nu)
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
q1 = [0.5; 0.43;
      0.65; 0.485; 0.0;]
# q1 = [0.5; 0.484;
#       0.5; 0.5; 0.0;]
v1 = [0.0; 0.0;
      0.0; 0.0; 0.0;]

# ## Simulator
sim = simulator(s, H, h=h)

# ## Simulate -- simulates entire trajectory?
status = simulate!(sim, q1, v1)
# ## Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.render(vis)

# ## Visualize
vis_traj = contact_trajectory(s.model, s.env, H, h)
anim = visualize_robot!(vis, model, sim.traj, sample = 1, h=h)
@infiltrate

# ## MPC setup 
N_sample = 2
H_mpc = 40
h_sim = h / N_sample
H_sim = 1000
κ_mpc = 1.0e-4

## Cost
q_scale = 20.0
q_vec = q_scale .* [1., 1., 1., 1., 1.] # x_ee, z_ee, x_tray, z_tray, θ_tray

v_scale = 10.0
v_vec = v_scale .* [1., 1., 1., 1., 0.] # x_ee, z_ee, x_tray, z_tray, θ_tray

u_scale = 1.0
u_vec = [1., 1.]

print("creating objective\n")
obj = TrackingVelocityObjective(model, env, H_mpc,
	q = [Diagonal(q_vec .* ones(5) .* (t/H_mpc)^2) for t = 1:H_mpc-0],
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
    n_opts = NewtonOptions(
		r_tol = 3e-3,
		max_iter = 10,
		max_time = ref_traj.h/2, # HARD REAL TIME
		),
                  mpc_opts = CIMPCOptions());

# ## Disturbances
# idx_d1 = 20
# idx_d2 = idx_d1 + 200
# idx_d3 = idx_d2 + 80
# idx_d4 = idx_d3 + 200
# idx_d5 = idx_d4 + 30
# idx = [idx_d1, idx_d2, idx_d3, idx_d4, idx_d5]
# impulses = [[0.0; 0.0], [0.0; 0.0], [0.0; 0.0], [0.0; 0.0], [0.0; 0.0]]
# d = impulse_disturbances(impulses, idx);

# ## Initial Conditions
q1_sim = q1
v1_sim = v1

print("simulation\n")
# ## Simulator
sim = simulator(s, H_sim, h=h_sim, policy=p) #, dist=d)

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
@infiltrate
