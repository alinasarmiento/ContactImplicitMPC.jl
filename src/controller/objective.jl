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
    qlim::Vector{QLim}
    ulim::Vector{ULim}
end

function TrackingVelocityObjective(model, env, H::Int;
    q = [Diagonal(zeros(SizedVector{model.nq})) for t = 1:H],
    v = [Diagonal(zeros(SizedVector{model.nq})) for t = 1:H],
    u = [Diagonal(zeros(SizedVector{model.nu})) for t = 1:H],
    γ = [Diagonal(zeros(SizedVector{model.nc})) for t = 1:H],
    b = [Diagonal(zeros(SizedVector{model.nc * friction_dim(env)})) for t = 1:H],
    qlim = [Diagonal(zeros(SizedVector{model.nq})) for t = 1:H],
    ulim = [Diagonal(zeros(SizedVector{model.nu})) for t = 1:H])
    return TrackingVelocityObjective(q, v, u, γ, b, qlim, ulim)
end

function eval_obj(model, env, z, q1, p, obj::Objective)
    q2, γ1, b1, ψ1, s1, η1, s2 = unpack_z(model, env, z)
    # v = (q2-q1)/0.01
    
    # qobj = transpose(q2) * obj.q * q2
    # vobj = transpose(v) * obj.v * v
    # uobj = transpose(u) * obj.u * u
    # gam_obj = transpose(γ1) * obj.γ * γ1
    # bobj = transpose(b1) * obj.b * b1
    return γ1
end

    
