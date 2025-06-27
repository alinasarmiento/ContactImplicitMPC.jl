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
s = get_simulation("pushbot", "flat_2D_lc", "flat");
model = s.model
env = s.env
# @infiltrate

# ## Reference Trajectory
h = 0.05 #0.04
H = 1000
ref_traj = contact_trajectory(model, env, H, h)
ref_traj.h
qref = [0.0; 0.0]
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
q1 = [0.0,0.0] 
v1 = [0.0; 0.0]

## Simulator
sim = simulator(s, H, h=h)

## Simulate -- set initial conditions
status = simulate!(sim, q1, v1)

## MPC setup 
N_sample = 5
H_mpc = 40
h_sim = 0.01 #h / N_sample
H_sim = 10000
κ_mpc = 1.0e-4

# Slow Recovery
obj = TrackingVelocityObjective(model, env, H_mpc,
	q = [Diagonal([20*(t/H_mpc)^2; 4.0*(t/H_mpc)^4]) for t = 1:H_mpc-0],
	v = [Diagonal([1; 0.01] ./ (h^2.0)) for t = 1:H_mpc-0],
	u = [Diagonal([100; 5.0]) for t = 1:H_mpc-0],
	γ = [Diagonal(1.0e-100 * ones(model.nc)) for t = 1:H_mpc-0],
	b = [Diagonal(1.0e-100 * ones(model.nc * friction_dim(env))) for t = 1:H_mpc]);

## mass == 0.001
# obj = TrackingVelocityObjective(model, env, H_mpc,
#     q = [Diagonal([40*(t/H_mpc)^2; 5*(t/H_mpc)])^2 for t = 1:H_mpc-0],
#     v = [Diagonal([1.0; 0.5] ./ (h^2.0)) for t = 1:H_mpc-0],
#     u = [Diagonal([10; .001]) for t = 1:H_mpc-0],
#     γ = [Diagonal(1.0e-100 * ones(model.nc)) for t = 1:H_mpc-0],
#     b = [Diagonal(1.0e-100 * ones(model.nc * friction_dim(env))) for t = 1:H_mpc]);

## mass == 1.0
# obj = TrackingVelocityObjective(model, env, H_mpc,
#     q = [Diagonal([4000*(t/H_mpc)^2; 5*(t/H_mpc)])^2 for t = 1:H_mpc-0],
#     v = [Diagonal([1.0; 0.5] ./ (h^2.0)) for t = 1:H_mpc-0],
#     u = [Diagonal([1; .0001]) for t = 1:H_mpc-0],
#     γ = [Diagonal(1.0e-10 * ones(model.nc)) for t = 1:H_mpc-0],
#     b = [Diagonal(1.0e-100 * ones(model.nc * friction_dim(env))) for t = 1:H_mpc]);


## Policy
p = ci_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
		r_tol = 3e-4,
		max_iter = 10,
		max_time = ref_traj.h/2, # HARD REAL TIME
		),
                  mpc_opts = CIMPCOptions());

## Initial Conditions
q1_sim = [0.1, 0.0]
v1_sim = [0.0; 0.0]

## Simulator
sim = simulator(s, H_sim, h=h_sim, policy=p) #, dist=d)

# set up for warm start
newton_solve!(sim.policy.newton, sim.policy.s, sim.policy.q0, q1_sim,
              sim.policy.im_traj, sim.policy.traj, warm_start = false)

x_lcm_channel = "PUSHBOT_STATE_SIMULATION"
u_lcm_channel = "PUSHBOT_INPUT"

## LCM + Drake loop in Python
sys = pyimport("sys")
pushfirst!(sys."path","")
lcm = pyimport("lcm")
lcm_py_callback = pyimport("lcmtypes.lcm_py_callback")
lc = lcm.LCM()
subscription = lc.subscribe(x_lcm_channel, lcm_py_callback.py_handler(lc, sim, u_lcm_channel, q1_sim))
print("\n LCM ready.")

while true
    lc.handle()
end
