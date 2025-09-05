abstract type Objective end

mutable struct TrackingObjective{Q,U} <: Objective
    q::Vector{Q}
    u::Vector{U}
end

function TrackingObjective(model, env, H::Int;
    q = [Diagonal(zeros(SizedVector{model.nq})) for t = 1:H],
    u = [Diagonal(zeros(SizedVector{model.nu})) for t = 1:H])
    return TrackingObjective(q, u)
end

function gradient!(res, obj::TrackingObjective, core, traj, ref_traj)
    for t = 1:traj.H
        # Cost function
        delta!(core.Δq[t], traj.q[t+2], ref_traj.q[t+2])
        delta!(core.Δu[t], traj.u[t], ref_traj.u[t])

        res.q2[t] .+= obj.q[t] * core.Δq[t]
        res.u1[t] .+= obj.u[t] * core.Δu[t]
    end
end

function hessian!(hess, obj::TrackingObjective)
    for t = 1:length(obj.u)
        # Cost function
        hess.obj_q2[t] .+= obj.q[t]
        hess.obj_u1[t] .+= obj.u[t]
    end
end

mutable struct TrackingVelocityObjective{Q,V,U} <: Objective
    q::Vector{Q}
    v::Vector{V}
    u::Vector{U}
end

function TrackingVelocityObjective(model, env, H::Int;
    q = [Diagonal(zeros(SizedVector{model.nq})) for t = 1:H],
    v = [Diagonal(zeros(SizedVector{model.nq})) for t = 1:H],
    u = [Diagonal(zeros(SizedVector{model.nu})) for t = 1:H])
    return TrackingVelocityObjective(q, v, u)
end

function gradient!(res, obj::TrackingVelocityObjective, core, traj, ref_traj)
    for t = 1:traj.H
        # Cost function
        delta!(core.Δq[t], traj.q[t+2], ref_traj.q[t+2])
        delta!(core.Δu[t], traj.u[t], ref_traj.u[t])

        res.q2[t] .+= obj.q[t] * core.Δq[t]
        res.u1[t] .+= obj.u[t] * core.Δu[t]

        # limits (new) : qlim / ulim in format [[min], [max]] Vector{Vector{T}}
        nq = size(traj.qlim[1])
        nu = size(traj.ulim[1])
        q_vio_min = max.(zeros(nq), traj.qlim[1] - traj.q[t+2])
        q_vio_max = max.(zeros(nq), traj.q[t+2] - traj.qlim[2])
        u_vio_min = max.(zeros(nu), traj.ulim[1] - traj.u[t])
        u_vio_max = max.(zeros(nu), traj.u[t] - traj.ulim[2])

        res.q2[t] .+= obj.qlim[t] * (q_vio_min + q_vio_max)
        res.u1[t] .+= obj.ulim[t] * (u_vio_min + u_vio_max)

        # velocity
        res.q2[t] .+= obj.v[t] * (traj.q[t+2] - traj.q[t+1])
        t == 1 && continue
        res.q2[t-1] .-= obj.v[t] * (traj.q[t+2] - traj.q[t+1])
    end
end

function hessian!(hess, obj::TrackingVelocityObjective)
    for t = 1:length(obj.u)
        # Cost function
        hess.obj_q2[t] .+= obj.q[t]
        hess.obj_u1[t] .+= obj.u[t]

        # velocity
        hess.obj_q2[t] .+= obj.v[t]
        t == 1 && continue
        hess.obj_q2[t-1] .+= obj.v[t]
        hess.obj_q1q2[t-1] .-= obj.v[t]
        hess.obj_q2q1[t-1] .-= obj.v[t]
    end
end
