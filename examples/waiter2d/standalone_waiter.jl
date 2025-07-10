# PREAMBLE

# PKG_SETUP

# ## Setup
 
using ContactImplicitMPC
using RoboDojo
import ContactImplicitMPC: simulate!, newton_solve!
using LinearAlgebra
using LCMCore
import ContactImplicitMPC: callback_sim, lcmt_robot_output, lcmt_robot_input, debug_callback, callback_sim_py
using PyCall
using Infiltrator

# ## Simulation
s = get_simulation("waiter_2D", "flat_2D_lc", "flat");
model = s.model
env = s.env
# @infiltrate

# ## Reference Trajectory
h = 0.01
H = 100
ref_traj = contact_trajectory(model, env, H, h)
ref_traj.h
qref = [0.5; 0.42;
        0.65; 0.4831; 0.0;]
ur = zeros(model.nu)
γr = zeros(model.nc)
br = zeros(model.nc * friction_dim(env))
ψr = zeros(model.nc)
ηr = zeros(model.nc * friction_dim(env))
wr = zeros(model.nw)

## Set Reference
for t = 1:H
	ref_traj.z[t] = pack_z(model, env, qref, γr, br, ψr, ηr)
	ref_traj.θ[t] = pack_θ(model, qref, qref, ur, wr, model.μ_world, ref_traj.h)
end

## Initial conditions
q1 = [0.5; 0.42;
        0.65; 0.4831; 0.0;]
v1 = [0.0; 0.0;
      0.0; 0.0; 0.0;]

## Simulator
sim = simulator(s, H, h=h)

## Simulate -- set initial conditions
status = simulate!(sim, q1, v1)

# ## MPC setup 
N_sample = 2
H_mpc = 40
h_sim = h / N_sample
H_sim = 200
κ_mpc = 1.0e-5

## Cost
q_scale = 1.0
q_vec = q_scale .* [1., 10., 1., 1., 1.] # x_ee, z_ee, x_tray, z_tray, θ_tray

v_scale = 1.0
v_vec = v_scale .* [1., 1., 1., 1., 1.] # x_ee, z_ee, x_tray, z_tray, θ_tray

u_scale = 1
u_vec = [1e-8, 1e-8]

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
    ip_opts = InteriorPointOptions(
                          undercut = 1.0,
                          κ_tol = κ_mpc,
                          r_tol = 1.0e-8,
                          diff_sol = true,
                          solver = :empty_solver,
                          max_time = 1e5,),
    n_opts = NewtonOptions(
                r_tol = 3e-4,
                max_iter = 10,
                max_time = ref_traj.h/2, # HARD REAL TIME
                ),
                  mpc_opts = CIMPCOptions());


# ## Initial Conditions
q1_sim = q1
v1_sim = v1

## Simulator
sim_ip_opts = InteriorPointOptions(
                          undercut = 1.0,
                          κ_tol = κ_mpc,
                          r_tol = 1.0e-8,
                          diff_sol = true,
                          # solver = :empty_solver,
                          max_time = 1e5,)
sim = simulator(s, H_sim, h=h_sim, policy=p, solver_opts=sim_ip_opts) #, dist=d)

# set up for warm start
newton_solve!(sim.policy.newton, sim.policy.s, sim.policy.q0, q1_sim,
              sim.policy.im_traj, sim.policy.traj, warm_start = false)

x_lcm_channel = "ALL_STATE_SIMULATION"
u_lcm_channel = "ROBOT_INPUT"

## LCM + Drake loop in Python
sys = pyimport("sys")
pushfirst!(sys."path","")
lcm = pyimport("lcm")
lcm_py_callback = pyimport("lcmtypes.lcm_py_callback")
lc = lcm.LCM()
subscription = lc.subscribe(x_lcm_channel, lcm_py_callback.py_handler(lc,
                                                                      sim,
                                                                      model,
                                                                      env,
                                                                      u_lcm_channel,
                                                                      q0_sim=q1,
                                                                      hp=0.01))
print("\n LCM ready.")

while true
    lc.handle()
end

## Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.render(vis)

## Visualize
vis_traj = contact_trajectory(s.model, s.env, H_sim, h_sim)
anim = visualize_robot!(vis, model, vis_traj, sample = 1)
pθ_right = generate_pusher_traj(d, vis_traj, side=:right)
pθ_left  = generate_pusher_traj(d, vis_traj, side=:left)
visualize_disturbance!(vis, model, pθ_right, anim=anim, sample=1, offset=0.05, name=:PusherRight);
visualize_disturbance!(vis, model, pθ_left,  anim=anim, sample=1, offset=0.05, name=:PusherLeft);

## Timing result
# Julia is [JIT-ed](https://en.wikipedia.org/wiki/Just-in-time_compilation) so re-run the MPC setup through Simulate for correct timing results.
process!(sim.stats, N_sample) # Time budget
H_sim * h_sim / sum(sim.stats.policy_time) # Speed ratio
