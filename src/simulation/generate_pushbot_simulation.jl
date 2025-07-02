################################################################################
# PushBot
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/pushbot")
dir_sim   = joinpath(module_dir(), "src/simulation/pushbot")
dir_sim_env   = joinpath(dir_sim, "flat")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
model = deepcopy(pushbot)
env = deepcopy(flat_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "flat/residual.jld2")
path_jac = joinpath(dir_sim, "flat/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env, jacobians=:approx)
save_expressions(expr_res, path_res, overwrite=true)
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac, jacobians=:approx)

################################################################################
# Waiter2D
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/waiter_2D")
dir_sim   = joinpath(module_dir(), "src/simulation/waiter_2D")
model = deepcopy(waiter_2D)
env = deepcopy(flat_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "flat/residual.jld2")
path_jac = joinpath(dir_sim, "flat/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)
