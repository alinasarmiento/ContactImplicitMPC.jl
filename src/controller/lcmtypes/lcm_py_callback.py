from julia.api import Julia
jl = Julia(compiled_modules=False)
from julia import Main
from julia import convert as jlconvert
Main.eval('using ContactImplicitMPC')
Main.eval('using StaticArrays')
import julia.ContactImplicitMPC as cimpc

import sys
import os
sys.path.append(os.environ['LCMT_PATH'])
import lcm
from dairlib import lcmt_robot_input, lcmt_robot_output

def py_handler(lc, sim, u_lcm_channel):
    def handler(channel, msg):
        # print("\n received")
        msg = lcmt_robot_output.decode(msg)
        p = sim.policy
        traj = sim.traj
        q1 = msg.position
        print("\n pos:",q1)
        q1 = jlconvert(Main.Vector, list(q1))

        cimpc.newton_solve_b(p.newton, p.s, p.q0, q1,
                            p.im_traj, p.traj, warm_start=True)
        cimpc.update_b(p.im_traj, p.traj, p.s, p.altitude, p.κ[0], p.traj.H)
        cimpc.rot_n_stride_b(p.traj, p.traj_cache, p.stride)
        cimpc.update_q0_u_b(p, q1)

        # lcm broadcast p.u
        u = lcmt_robot_input()
        u.utime = msg.utime
        u.num_efforts = msg.num_efforts
        u.effort_names = msg.effort_names
        u.efforts = p.u
        print("\n u:",u.efforts)
        lc.publish(u_lcm_channel, u.encode())

    return handler
