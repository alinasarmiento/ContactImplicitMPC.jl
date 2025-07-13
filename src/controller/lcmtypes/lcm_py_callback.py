from julia.api import Julia
jl = Julia(compiled_modules=False)
from julia import Main
from julia import convert as jlconvert
Main.eval('using ContactImplicitMPC')
Main.eval('using RoboDojo')
Main.eval('using StaticArrays')
import julia.ContactImplicitMPC as cimpc
import julia.RoboDojo as rodo

import sys
import os
sys.path.append(os.environ['LCMT_PATH'])
import lcm
import numpy as np
from dairlib import lcmt_robot_input, lcmt_robot_output, lcmt_object_state
from IPython import embed
import math
 
def euler_from_quaternion(w, x, y, z):
        """
        Convert a quaternion into euler angles (roll, pitch, yaw)
        roll is rotation around x in radians (counterclockwise)
        pitch is rotation around y in radians (counterclockwise)
        yaw is rotation around z in radians (counterclockwise)
        """
        t0 = +2.0 * (w * x + y * z)
        t1 = +1.0 - 2.0 * (x * x + y * y)
        roll_x = math.atan2(t0, t1)
     
        t2 = +2.0 * (w * y - z * x)
        t2 = +1.0 if t2 > +1.0 else t2
        t2 = -1.0 if t2 < -1.0 else t2
        pitch_y = math.asin(t2)
     
        t3 = +2.0 * (w * z + x * y)
        t4 = +1.0 - 2.0 * (y * y + z * z)
        yaw_z = math.atan2(t3, t4)
     
        return roll_x, pitch_y, yaw_z # in radians

def py_handler_pushbot(lc, sim, model, env, u_lcm_channel, q0_sim=[0,0], hp=0.05):
    def handler(channel, msg):
        msg = lcmt_robot_output.decode(msg)
        p = sim.policy
        q1 = msg.position
        q1 = jlconvert(Main.Vector, list(q1))

        # check if utime is next h
        t_now = msg.utime/1e6
        t_ctrl = int(t_now / hp)
        if t_now%hp <= 0.0101:
            # if t_now == 0:
            #     cimpc.set_initial_q0_b(p,q0_sim)
            cimpc.newton_solve_b(p.newton, p.s, p.q0, q1,
                                 p.im_traj, p.traj, warm_start=t_ctrl>0)
            cimpc.update_b(p.im_traj, p.traj, p.s, p.altitude, p.κ[0], p.traj.H)
            cimpc.rot_n_stride_b(p.traj, p.traj_cache, p.stride)
            cimpc.update_q0_u_b(p, q1)
        
        # lcm broadcast p.u
        u_lcm = lcmt_robot_input()
        u_lcm.utime = msg.utime
        ## pushbot
        u_lcm.num_efforts = msg.num_efforts
        u_lcm.effort_names = msg.effort_names

        # embed()
        u_lcm.efforts = p.u / 0.01
        if np.abs(u_lcm.efforts[0]) > 0.3:
            u_lcm.efforts[0] = 0.3*np.sign(u_lcm.efforts[0])
            print('#################### caught. p.u:', p.u)
        print("u:",u_lcm.efforts,"h:", sim.h, "t:", t_now%hp)
    
        lc.publish(u_lcm_channel, u_lcm.encode())

    return handler

def py_handler_waiter(lc, sim, model, env, u_lcm_channel, q0_sim=[0,0], hp=0.005):
    def handler(channel, msg):
        msg = lcmt_object_state.decode(msg)
        p = sim.policy
        q_ee = list(msg.position[0:2])
        t_xz = [msg.position[-3], msg.position[-1]]
        t_euler = euler_from_quaternion(*msg.position[2:6])
        t_th = [t_euler[2]]

        q1 = q_ee + t_xz + t_th
        q1 = jlconvert(Main.Vector, list(q1))

        # check if utime is next h
        t_now = msg.utime/1e6
        t_ctrl = int(t_now / hp)
        print(t_now)
        if t_now%hp <= hp*1.01:
                
            cimpc.newton_solve_b(p.newton, p.s, p.q0, q1,
                                 p.im_traj, p.traj, warm_start=t_ctrl>0)
            cimpc.update_b(p.im_traj, p.traj, p.s, p.altitude, p.κ[0], p.traj.H)
            cimpc.rot_n_stride_b(p.traj, p.traj_cache, p.stride)
            cimpc.update_q0_u_b(p, q1)
        
        # lcm broadcast p.u
        u_lcm = lcmt_robot_input()
        u_lcm.utime = msg.utime
        u_lcm.num_efforts = 2
        u_lcm.effort_names = ["x_motor", "z_motor"]

        u_lcm.efforts = p.u / (hp/p.N_sample)
        if any(u_lcm.efforts == np.nan):
            u_lcm.efforts = [0.0,0.0]
        print("u:",u_lcm.efforts,"h:", sim.h, "t:", t_now%hp)
    
        lc.publish(u_lcm_channel, u_lcm.encode())

    return handler
