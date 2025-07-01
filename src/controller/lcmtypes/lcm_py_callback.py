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
import numpy as np
from dairlib import lcmt_robot_input, lcmt_robot_output

def py_handler(lc, sim, u_lcm_channel, q0_sim=[0,0]):
    def handler(channel, msg):
        # print("\n received")
        msg = lcmt_robot_output.decode(msg)
        p = sim.policy
        traj = sim.traj
        q1 = msg.position
        # print("\n pos:",q1)
        q1 = jlconvert(Main.Vector, list(q1))

        # check if utime is next h
        t_now = msg.utime/1e6
        if t_now%0.05 <= 0.0101:
            # if t_now == 0:
            #     cimpc.set_initial_q0_b(p,q0_sim)
            cimpc.newton_solve_b(p.newton, p.s, p.q0, q1,
                                p.im_traj, p.traj, warm_start=True)
            cimpc.update_b(p.im_traj, p.traj, p.s, p.altitude, p.κ[0], p.traj.H)
            cimpc.rot_n_stride_b(p.traj, p.traj_cache, p.stride)
            cimpc.update_q0_u_b(p, q1)

        # apply B to u
        # B = sim.model.base.B([0,0]) # q doesnt matter
        # B = jlconvert(Main.Matrix, B)
        # B = np.array(B)
        # u = np.array(p.u)
        # u = B @ u

        # lcm broadcast p.u
        u_lcm = lcmt_robot_input()
        u_lcm.utime = msg.utime
        u_lcm.num_efforts = msg.num_efforts
        u_lcm.effort_names = msg.effort_names
        u_lcm.efforts = p.u / 0.01
        print("u:",u_lcm.efforts,"h:", sim.h, "t:",t_now)
        lc.publish(u_lcm_channel, u_lcm.encode())

    return handler
