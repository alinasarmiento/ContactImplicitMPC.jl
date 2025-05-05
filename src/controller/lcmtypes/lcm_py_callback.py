from julia.api import Julia
jl = Julia(compiled_modules=False)
from julia import Main
cimpc_path = "/home/grey/research/ContactImplicitMPC.jl/"
Main.eval('using Pkg; Pkg.add(path="%s");' % cimpc_path)
Main.eval('using ContactImplicitMPC')
Main.eval('using StaticArrays')
# Main.include(cimpc_path + "src/controller/standalone_simulate.jl")
import julia.ContactImplicitMPC as cimpc

import sys
import os
sys.path.append(os.environ['LCMT_PATH'])
import lcm
from dairlib import lcmt_robot_input, lcmt_robot_output
# from IPython import embed; embed()

def py_callback_sim(sim, u_lcm_channel):
    def handler(channel, msg):
        print("\n received")
        msg = lcmt_robot_output.decode(msg)
        print("\n x:", msg)
        p = sim.policy
        traj = sim.traj
        q1 = msg.position

        cimpc.newton_solve_b(p.newton, p.s, p.q0, q1,
                            p.im_traj, p.traj, warm_start=true)
        cimpc.update_b(p.im_traj, p.traj, p.s, p.altitude, p.κ[1], p.traj.H)
        cimpc.rot_n_stride_b(p.traj, p.traj_cache, p.stride)
        p.q0 = q1

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

    return handler
