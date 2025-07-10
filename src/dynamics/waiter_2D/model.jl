using YAML
using Rotations

"""
    - 2D plate subject to contact forces

    - configuration: q = (x, z) ∈ R²
    - impact force (magnitude): γ ∈ R₊
    - friction force: β ∈ R²₊
        - friction coefficient: μ ∈ R₊

    Discrete Mechanics and Variational Integrators
        pg. 363
"""
mutable struct Waiter2D{T} <: Model{T}
    nq::Int ## EE (x,z), tray (x,z)
    nu::Int
    nw::Int
    nc::Int
    m::T # mass
    g::T # gravity
    m_tray::T # tray mass
    μ_world::T # friction coefficient
    μ_joint::T
    r::T # radius of EE
    d::T # depth of EE
    r_tray::T # radius of tray
    d_tray::T # depth of EE
    
    supp_1::SVector # pose of left (back) support point
    supp_2::SVector # pose of right (front) support point

    base::BaseMethods
    dyn::DynamicsMethods

    joint_friction::SVector
    u_min::SVector
    u_max::SVector
    u_vio_weight::T
end

function lagrangian(model::Waiter2D, q, q̇)
    L = 0.0

    # EE
    L += 0.5 * model.m * transpose(q̇[1:2]) * q̇[1:2]
    L -= model.m * model.g * q[2]

    # tray
    L += 0.5 * model.m_tray * transpose(q̇[3:4]) * q̇[3:4]
    L -= model.m_tray * model.g * q[4]

    return L
end

function kinematics(model::Waiter2D, q; mode=:contacts)
    # supposed to return pose of each contact point
    # (why also defined in visuals.jl??)
    if mode == :contacts
        ee1 = SVector{2}([q[1]-model.r, q[2]+(model.d/2)])
        ee2 = SVector{2}([q[1]+model.r, q[2]+(model.d/2)])
        
        ee3 = SVector{2}([q[1]-model.r, q[2]-(model.d/2)])
        ee4 = SVector{2}([q[1]+model.r, q[2]-(model.d/2)])          
        return SVector{12}([ee1; ee2; model.supp_1; model.supp_2;
                           ee3; ee4;])
    elseif mode == :ee
        return q[1:2]
    elseif mode == :tray
        return q[3:5]
    else
        @error "incorrect mode"
        return
    end
end

# mass matrix
function M_func(model::Waiter2D, q)
    m = model.m
    mt = model.m_tray
    h_t = model.d_tray
    w_t = 2*model.r_tray
    I_tray = (1/12)*mt*(h_t^2 + w_t^2)

    Diagonal(@SVector [m, m, mt, mt, I_tray])
end

# gravity
function C_func(model::Waiter2D, q, q̇)
    m = model.m
    mt = model.m_tray
    g = model.g

    @SVector [0.0, m*g, 0.0, mt*g, 0.0]
end

function dist_tray(model::Waiter2D, p, pt)
    # p: [x, z], pt: [xtray, ztray, θtray]

    diff = p-pt[1:2]
    # beta = atan(diff[1]/diff[2]) + pt[3] # angle between vector and tray-vertical
    R = [cos(pt[3]) -sin(pt[3]); sin(pt[3]) cos(pt[3])];
    xdiff,zdiff = R*(diff)
    xdiff = abs(xdiff)
    zdiff = abs(zdiff)
        
    # zdiff = abs(norm(diff)*cos(beta))
    zdist = zdiff-(model.d_tray/2)

    # xdiff = abs(norm(diff)*sin(beta))
    # xdist = xdiff-model.r_tray

    # dist_neg = min(0, max(xdist, zdist))
    
    # zdist = max(0, zdist)
    # xdist = max(0, xdist)
    
    # return norm([xdist, zdist]) + dist_neg
    return zdist
end

# signed distance function
function ϕ_func(model::Waiter2D, env::Environment, q)
    # ee_back-tray, ee_front-tray, tray-supp_back, tray-supp_front
    cp = kinematics(model, q, mode=:contacts)
    ee1 = cp[1:2]
    ee2 = cp[3:4]
    tray = q[3:5] 

    ee1_dist = dist_tray(model, ee1, tray)
    ee2_dist = dist_tray(model, ee2, tray)
    supp1_dist = dist_tray(model, model.supp_1, tray)
    supp2_dist = dist_tray(model, model.supp_2, tray)
    
    ee_ground = cp[10]
    
    return SVector{6}([ee1_dist; ee2_dist; supp1_dist; supp2_dist; ee_ground; ee_ground])
end

# control Jacobian
function B_func(model::Waiter2D, q)
    B = zeros(5,2)
    B[1,1] = 1
    B[2,2] = 1
    B = SMatrix{5,2}(B)
    return B
end

# disturbance Jacobian
function A_func(model::Waiter2D, q)
    A = zeros(5,2)
    A[1,1] = 1
    A[2,2] = 1
    A = SMatrix{5,2}(A)
    return A
end

function _jacobian(model::Waiter2D, q; mode=:ee_t)
    # J'λ = \tau
    # ee_t := EE-tray contact
    # t_supp := tray-support contact

    x_ee, z_ee, x_t, z_t, th_t = q
    x_ee1 = deepcopy(x_ee) - model.r
    x_ee2 = deepcopy(x_ee) + model.r
    # if mode == :ee_t
    #     j = SMatrix{4,5}([-1.0 0.0 cos(th_t) -sin(th_t) -(x_ee1-x_t)*sin(th_t);
    #                       0.0 -1.0 sin(th_t) cos(th_t) (x_ee1 - x_t)*cos(th_t);
    #                       -1.0 0.0 cos(th_t) -sin(th_t) -(x_ee2-x_t)*sin(th_t);
    #                       0.0 -1.0 sin(th_t) cos(th_t) (x_ee2 - x_t)*cos(th_t)])
    #     return j
        
    # elseif mode == :t_supp
    #     j = SMatrix{4,5}([0.0 0.0 cos(th_t) -sin(th_t) -(x_t-model.supp_1[1])*sin(th_t);
    #                       0.0 0.0 sin(th_t) cos(th_t) (x_t -model.supp_1[1])*cos(th_t);
    #                       0.0 0.0 cos(th_t) -sin(th_t) -(x_t-model.supp_2[1])*sin(th_t);
    #                       0.0 0.0 sin(th_t) cos(th_t) (x_t -model.supp_2[1])*cos(th_t)])
    #     return j
    if mode == :ee_t
        j = SMatrix{4,5}([1.0 0.0 1.0 0.0 -(x_t-x_ee1)*tan(th_t);
                          0.0 1.0 0.0 1.0 (x_t-x_ee1);
                          1.0 0.0 1.0 0.0 -(x_t-x_ee2)*tan(th_t);
                          0.0 1.0 0.0 1.0 (x_t-x_ee2)])
        return j
        
    elseif mode == :t_supp
        j = SMatrix{4,5}([0.0 0.0 1.0 0.0 -(x_t-model.supp_1[1])*tan(th_t);
                          0.0 0.0 0.0 1.0 (x_t-model.supp_1[1]);
                          0.0 0.0 1.0 0.0 -(x_t-model.supp_2[1])*tan(th_t);
                          0.0 0.0 0.0 1.0 (x_t-model.supp_2[1])])
        return j
    elseif mode == :ground
        j = SMatrix{4,5}([1.0 0.0 0.0 0.0 0.0;
                          0.0 1.0 0.0 0.0 0.0;
                          1.0 0.0 0.0 0.0 0.0;
                          0.0 1.0 0.0 0.0 0.0])
        return j

    end
end

# contact Jacobian
function J_func(model::Waiter2D, env::Environment, q)
    return SMatrix{12, 5}([_jacobian(model, q, mode=:ee_t);
                          _jacobian(model, q, mode=:t_supp);
                          _jacobian(model, q, mode=:ground);])
end

# translates the two variables normal force (γ) and tangential forces (b) into a single vector for jacobian
function contact_forces(model::Waiter2D, env::Environment{<:World, LinearizedCone}, γ1, b1, q2, k)
    # γ1: force vector (size num contacts)
    # b1: idk but size 2*num contacts
    # k: also size 2*num contacts
    
    m = friction_mapping(env) # what is this
    SVector{12}([transpose(rotation(env, k[1:1])) * [m * b1[1:2]; γ1[1]];
                 transpose(rotation(env, k[3:3])) * [m * b1[3:4]; γ1[2]];
                 transpose(rotation(env, k[5:5])) * [m * b1[5:6]; γ1[3]];
                 transpose(rotation(env, k[7:7])) * [m * b1[7:8]; γ1[4]];
                 transpose(rotation(env, k[9:9])) * [m * b1[9:10]; γ1[5]];
                 transpose(rotation(env, k[11:11])) * [m * b1[11:12]; γ1[6]];])
end

function velocity_stack(model::Waiter2D, env::Environment{<:World, LinearizedCone}, q1, q2, k, h)
    v = J_func(model, env, q2) * (q2 - q1) / h[1]
    v1_surf = rotation(env, k[1:1]) * v[1:2]
    v2_surf = rotation(env, k[3:3]) * v[3:4]
    v3_surf = rotation(env, k[5:5]) * v[5:6]
    v4_surf = rotation(env, k[7:7]) * v[7:8]
    v5_surf = rotation(env, k[9:9]) * v[9:10]
    v6_surf = rotation(env, k[11:11]) * v[11:12]
    
    SVector{12}([transpose(friction_mapping(env)) * v1_surf[1];
                transpose(friction_mapping(env)) * v2_surf[1];
                transpose(friction_mapping(env)) * v3_surf[1];
                transpose(friction_mapping(env)) * v4_surf[1];
                transpose(friction_mapping(env)) * v5_surf[1];
                transpose(friction_mapping(env)) * v6_surf[1];])
end


# Working Parameters
params = YAML.load_file(joinpath(@__DIR__,"params.yaml"))

supp1 = deepcopy(params["supp_pos"]) # back point
supp1[1] -= params["supp_xdim"]/2
supp1[2] += params["supp_zdim"]/2

supp2 = deepcopy(params["supp_pos"]) # front point
supp2[1] -= params["supp_xdim"]/2 - 0.1
supp2[2] += params["supp_zdim"]/2

# nq, nu, nw, nc, m, g, mt, mu_world, mu_joint, r, d, r_tray, d_tray, supp1, supp2

waiter_2D = Waiter2D(5, 2, 2, 6,
                     params["m_ee"], params["gravity"], params["m_tray"],
                     params["mu_world"], params["mu_joint"],
                     params["r_ee"], params["d_ee"], params["r_tray"], params["d_tray"],
                     SVector{2}(supp1), SVector{2}(supp2),
	             BaseMethods(), DynamicsMethods(),
	             SVector{5}(zeros(5)),
                     SVector{2}([-20, -20]),
                     SVector{2}([20, 20]),
                     1.0)

function friction_coefficients(model::Waiter2D) 
	return [model.μ_world]
end

function initialize_z!(z, model::Waiter2D, idx::RoboDojo.IndicesZ, q)
    z .= 1.0
    z[idx.q] .= q
end
