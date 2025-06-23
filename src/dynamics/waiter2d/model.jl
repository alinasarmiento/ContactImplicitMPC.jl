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

function kinematics(::Waiter2D, q; mode=:ee)
    if mode == :ee
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

    Diagonal(@SVector [m, m, mt, mt, mt])
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
    beta = pi - atan(diff[2]/diff[1]) - ((pi/2)-pt[3]) # angle between vector and tray-vertical

    zdiff = norm(diff)*cos(beta)
    zdist = max(0, zdiff-model.d_tray)

    xdiff = norm(diff)*sin(beta)
    xdist = max(0, xdiff-model.r_tray)

    return norm([xdist, zdist])
end

# signed distance function
function ϕ_func(model::Waiter2D, env::Environment, q)
    # ee_back-tray, ee_front-tray, tray-supp_back, tray-supp_front
    ee1 = SVector{2}([q[1]-model.r, q[2]])
    ee2 = SVector{2}([q[1]+model.r, q[2]])
    tray = q[3:5]

    ee1_dist = dist_tray(model, ee1, tray)
    ee2_dist = dist_tray(model, ee2, tray)
    supp1_dist = dist_tray(model, model.supp_1, tray)
    supp2_dist = dist_tray(model, model.supp_2, tray)
    
    return SVector{4}([ee1_dist; ee2_dist; supp1_dist; supp2_dist])
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
	SMatrix{2, 2}([1.0 0.0;
                       0.0 1.0])
end

function _jacobian(model::Waiter2D, q; mode=:ee_t)
    # J'λ = \tau
    # ee_t := EE-tray contact
    # t_supp := tray-support contact

    x_ee, z_ee, x_t, z_t, th_t = q
    if mode == :ee_t
        j = SMatrix{2,5}([0.0 1.0 sin(th_t) cos(th_t) (x_ee - x_t)*cos(th_t);
                          1.0 0.0 cos(th_t) -sin(th_t) -(x_ee-x_t)*sin(th_t)])
        return j
        
    elseif mode == :t_supp
        j = SMatrix{2,5}([0.0 0.0 sin(th_t) cos(th_t) (x_ee - x_t)*cos(th_t);
                          0.0 0.0 cos(th_t) -sin(th_t) -(x_ee-x_t)*sin(th_t)])
        return j
    end
end

# contact Jacobian
function J_func(model::Waiter2D, env::Environment, q)
	SMatrix{8, 2}([_jacobian(model, :ee_t);
                       _jacobian(model, :ee_t);
                       _jacobian(model, :t_supp);
                       _jacobian(model, :t_supp);])
end

# idk what this does
function contact_forces(model::Waiter2D, env::Environment{<:World, LinearizedCone}, γ1, b1, q2, k)
	m = friction_mapping(env)
	SVector{2}(transpose(rotation(env, k)) * [m * b1; γ1])
end

function velocity_stack(model::Waiter2D, env::Environment{<:World, LinearizedCone}, q1, q2, k, h)
	v = J_func(model, env, q2) * (q2 - q1) / h[1]
	v1_surf = rotation(env, k) * v

	SVector{2}(transpose(friction_mapping(env)) * v1_surf[1])
end


# Working Parameters
gravity = 9.81
μ_world = 0.4
μ_joint = 0.0

m_ee = 0.37
m_tray = 1
r_ee = 0.0725
d_ee = 0.01
r_tray = 0.2286
d_tray = 0.022

supp1 = SVector{2,1}([0.6,0.447]) # back point
supp2 = SVector{2,1}([0.7,0.447]) # front point

# nq, nu, nw, nc, m, g, mt, mu_world, mu_joint, r, d, r_tray, d_tray, supp1, supp2

waiter_2D = Waiter2D(5, 2, 2, 4,
                     m_ee, gravity, m_tray,
                     μ_world, μ_joint,
                     r_ee, d_ee, r_tray, d_tray,
                     supp1, supp2,
	             BaseMethods(), DynamicsMethods(),
	             SVector{2}(zeros(2)))

function friction_coefficients(model::Waiter2D) 
	return [model.μ_world]
end

function initialize_z!(z, model::Waiter2D, idx::RoboDojo.IndicesZ, q)
    z .= 1.0
    z[idx.q] .= q
end
