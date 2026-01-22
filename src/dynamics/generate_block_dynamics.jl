################################################################################
# Block
################################################################################
dir = joinpath(@__DIR__, "block")
model = deepcopy(block_system)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")

expr_base = generate_base_expressions(model, M_analytical=true, C_analytical=false)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)


################################################################################
# Block 1D
################################################################################
dir = joinpath(@__DIR__, "block_1D")
model = deepcopy(block_system_1D)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")

expr_base = generate_base_expressions(model, M_analytical=true, C_analytical=false)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)
