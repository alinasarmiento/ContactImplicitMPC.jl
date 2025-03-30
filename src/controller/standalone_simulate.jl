using LCMCore, StaticArrays
using Infiltrator

# 1. receive LCM state
# 2. Set Julia state from LCM message
# 3. run controller and calculate u
# 4. extract control action and transmit LCM message

mutable struct state_vector_t <: LCMType
    timestamp::Float64
    position::SVector{2, Float64}
    velocity::SVector{2, Float64}
end

mutable struct u_vector_t <: LCMType
    timestamp::Float64
    input::SVector{2, Float64}
end

@lcmtypesetup(state_vector_t)
@lcmtypesetup(u_vector_t)

function callback_sim(lcm, sim, u_lcm_channel)
    return function(channel::String, msg)
        @show channel
        @show msg
        print(msg)
        
        p = sim.policy
        traj = sim.traj
        q1 = msg.position

        newton_solve!(p.newton, p.s, p.q0, q1,
                       p.im_traj, p.traj, warm_start=true)
        update!(p.im_traj, p.traj, p.s, p.altitude, p.κ[1], p.traj.H)

        rot_n_stride!(p.traj, p.traj_cache, p.stride)
        p.q0 .= q1

        # scale control
        if p.newton_mode == :direct
            p.u .= p.newton.traj.u[1] 
            p.u ./= p.N_sample
        elseif p.newton_mode == :structure
            p.u .= p.newton.u[1] 
            p.u ./= p.N_sample
        else
            println("newton mode specified not available")
        end

        # lcm broadcast p.u
        u_lcm = u_vector_t(msg.timestamp, p.u)
        publish(lcm, u_lcm_channel, u_lcm)
        
    end
end
