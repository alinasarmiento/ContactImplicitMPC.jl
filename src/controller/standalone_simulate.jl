using LCMCore, StaticArrays
import LCMCore: encode, decode
# using Infiltrator

# 1. receive LCM state
# 2. Set Julia state from LCM message
# 3. run controller and calculate u
# 4. extract control action and transmit LCM message

mutable struct lcmt_robot_output <: LCMType
    utime::Int64
    num_positions::Int32
    num_velocities::Int32
    num_efforts::Int32

    position_names::Vector{String}
    position::Vector{Float64}
    
    velocity_names::Vector{String}
    velocity::Vector{Float64}
    
    effort_names::Vector{String}
    effort::Vector{Float64}

    imu_accel::SVector{3, Float64}
end

mutable struct lcmt_robot_input <: LCMType
    utime::Int64
    num_efforts::Int32

    effort_names::Vector{String}
    efforts::Vector{Float64}
end

@lcmtypesetup(lcmt_robot_output,
              position => (num_positions,),
              position_names => (num_positions,),
              velocity => (num_velocities,),
              velocity_names => (num_velocities,),
              effort => (num_efforts,),
              effort_names => (num_efforts,)
              )
@lcmtypesetup(lcmt_robot_input,
              effort_names => (num_efforts,),
              efforts => (num_efforts,)
              )

function debug_callback(channel::String, msg)
    print(msg)
    # @infiltrate
end

    
function callback_sim(lcm, sim, u_lcm_channel)
    return function(channel::String, msg)
        # print("\n decode \n")
        print("\n received")
        #msg = decode(msg, lcmt_robot_output)
        msg = lcmt_robot_output(0, 2, 2, 2, ["base_joint", "push_joint"], [0.7731647863882035, 0.20706555274702418], ["base_jointdot", "push_jointdot"], [4.551114023405554e-10, 6.184129494732598e-10], ["base_motor", "push_motor"], [0.0, 0.0], [0.0, 0.0, 0.0])
        print("\n x: ")
        print(msg)

        p = sim.policy
        # print(p)
        traj = sim.traj
        # q1 = msg.position
        q1 = [0.0, 0.0]
        print("\n sim policy properties\n")
        print(propertynames(p))
        print("\n q0\n")
        print(p.q0)

        # newton_solve!(p.newton, p.s, p.q0, q1,
        #                     p.im_traj, p.traj, warm_start=true)
        # update!(p.im_traj, p.traj, p.s, p.altitude, p.κ[1], p.traj.H)
        # rot_n_stride!(p.traj, p.traj_cache, p.stride)
        # p.q0 .= q1

        # # scale control
        # if p.newton_mode == :direct
        #     p.u .= p.newton.traj.u[1] 
        #     p.u ./= p.N_sample
        # elseif p.newton_mode == :structure
        #     p.u .= p.newton.u[1] 
        #     p.u ./= p.N_sample
        # else
        #     println("newton mode specified not available")
        # end

        # lcm broadcast p.u
        # u = lcmt_robot_input(msg.utime, msg.num_efforts, msg.effort_names, p.u)
        # print(u)
        # u_lcm = encode(u)
        # publish(lcm, u_lcm_channel, u_lcm)        
    end
end
