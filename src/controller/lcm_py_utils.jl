function update_q0_u!(p, q1)
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
end

function set_initial_q0!(p, q0)
    p.q0 .= q0
    set_trajectory!(p.traj, p.ref_traj)
    set_implicit_trajectory!(p.im_traj, p.im_traj_cache)
end
