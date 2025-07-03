abstract type Objective end

mutable struct TrackingObjective{Q,U,C,B} <: Objective
    q::Vector{Q}
    u::Vector{U}
    γ::Vector{C}
    b::Vector{B}
end

function TrackingObjective(model, env, H::Int;
    q = [Diagonal(zeros(SizedVector{model.nq})) for t = 1:H],
    u = [Diagonal(zeros(SizedVector{model.nu})) for t = 1:H],
    γ = [Diagonal(zeros(SizedVector{model.nc})) for t = 1:H],
    b = [Diagonal(zeros(SizedVector{model.nc * friction_dim(env)})) for t = 1:H])
    return TrackingObjective(q, u, γ, b)
end

mutable struct TrackingVelocityObjective{Q,V,U,C,B} <: Objective
    q::Vector{Q}
    v::Vector{V}
    u::Vector{U}
    γ::Vector{C}
    b::Vector{B}
end

function TrackingVelocityObjective(model, env, H::Int;
    q = [Diagonal(zeros(SizedVector{model.nq})) for t = 1:H],
    v = [Diagonal(zeros(SizedVector{model.nq})) for t = 1:H],
    u = [Diagonal(zeros(SizedVector{model.nu})) for t = 1:H],
    γ = [Diagonal(zeros(SizedVector{model.nc})) for t = 1:H],
    b = [Diagonal(zeros(SizedVector{model.nc * friction_dim(env)})) for t = 1:H])
    return TrackingVelocityObjective(q, v, u, γ, b)
end

function eval_obj(model, env, im_traj::ImplicitTrajectory, H::Int, obj::Objective, t::Int)
    z = im_traj.lin[t].z
    q2, γ1, b1, ψ1, s1, η1, s2 = unpack_z(model, env, z)
    print(z)
    # qobj = transpose(q2) * obj.q * q2
end

    
