# using Pkg; Pkg.activate(joinpath(@__DIR__, "../../../"))
using YAML
using Rotations
# using StaticArrays
# using RoboDojo

"""
    - 2D plate subject to contact forces

    - configuration: q = (x, z) ∈ R²
    - impact force (magnitude): γ ∈ R₊
    - friction force: β ∈ R²₊
        - friction coefficient: μ ∈ R₊

    Discrete Mechanics and Variational Integrators
        pg. 363
"""
mutable struct Block{T} <: Model{T}
    nq::Int ## EE (x,z), tray (x,z)
    nu::Int
    nw::Int
    nc::Int
    m::T # mass
    g::T # gravity
    m_block::T # block mass
    μ_world::T # friction coefficient
    μ_block::T
    r::T # radius of EE
    xlen_block::T
    ylen_block::T
    zlen_block::T
    
    base::BaseMethods
    dyn::DynamicsMethods

    joint_friction::SVector
    u_min::SVector
    u_max::SVector
    q_min::SVector
    q_max::SVector
    
end

function lagrangian(model::Block, q, q̇)
    L = 0.0

    # EE
    L += 0.5 * model.m * transpose(q̇[1:2]) * q̇[1:2]
    L -= model.m * model.g * q[2]

    # block
    L += 0.5 * model.m_block * transpose(q̇[3:4]) * q̇[3:4]
    # I_b = (1/12) * model.m_block * model.zlen_block * (model.xlen_block^3)
    I_b = (1/12) * model.m_block * (model.zlen_block^2) * (model.xlen_block^2)
    L += 0.5*I_b*q̇[5]^2
    L -= model.m_block * model.g * q[4]

    return L
end

function kinematics(model::Block, q; mode=:contacts)
    # supposed to return pose of each contact point
    # (why also defined in visuals.jl??)
    if mode == :contacts
        ee = SVector{2}([q[1], q[2]])
        block1 = SVector{2}([q[3]+(model.xlen_block/2), q[4]-(model.zlen_block/2)])
        block2 = SVector{2}([q[3]-(model.xlen_block/2), q[4]-(model.zlen_block/2)])
        return SVector{6}([ee; block1; block2;])
    elseif mode == :ee
        return q[1:2]
    elseif mode == :block
        return q[3:5]
    else
        @error "incorrect mode"
        return
    end
end

# mass matrix
function M_func(model::Block, q)
    m = model.m
    mb = model.m_block
    h_b = model.zlen_block
    w_b = model.xlen_block
    I_b = (1/12)*mt*(h_b^2 + w_b^2)

    Diagonal(@SVector [m, m, mb, mb, I_b])
end

# gravity
function C_func(model::Block, q, q̇)
    m = model.m
    mb = model.m_block
    g = model.g

    @SVector [0.0, m*g, 0.0, mb*g, 0.0]
end

function dist_block(model::Block, p, pt)
    # p: [x, z], pt: [xblock, zblock, θblock]

    differ = p-pt[1:2]
    # beta = atan(diff[1]/diff[2]) + pt[3] # angle between vector and tray-vertical
    R = [cos(pt[3]) -sin(pt[3]); sin(pt[3]) cos(pt[3])];
    xdiff,zdiff = R*(differ)

    ## just halfplane
    zdist = -zdiff - model.zlen_block/2
    
    ## full
    # xdiff = abs(xdiff)
    # zdiff = abs(zdiff)        
    # zdist = zdiff-(model.d_tray/2)
    # xdist = xdiff-model.r_tray

    # dist_neg = min(0, max(xdist, zdist))
    
    # zdist = max(0, zdist)
    # xdist = max(0, xdist)
    
    # return norm([xdist, zdist]) + dist_neg
   return zdist
end

# signed distance function
function ϕ_func(model::Block, env::Environment, q)
    # ee-block, block_front-ground, block_back-ground
    cp = kinematics(model, q, mode=:contacts)
    ee = cp[1:2]
    block1 = cp[3:4]
    block2 = cp[5:6]
    block_q = q[3:5]    

    ee_block_dist = dist_block(model, ee, block_q)
    block1_dist = block_q[2]
    block2_dist = block_q[2]
    
    return SVector{3}([ee_block_dist; block1_dist; block2_dist])
end

# control Jacobian
function B_func(model::Block, q)
    B = zeros(2,5)
    B[1,1] = 1
    B[2,2] = 1
    B = SMatrix{2,5}(B)
    return B
end

# disturbance Jacobian
function A_func(model::Block, q)
    A = zeros(5,2)
    A[1,1] = 1
    A[2,2] = 1
    A = SMatrix{5,2}(A)
    return A
end

function _jacobian(model::Block, q; mode=:ee_b)
    # J'λ = \tau
    # ee_b := EE-block contact
    # b_g := block-ground contact
    # ee_g := EE-ground contact

    x_ee, z_ee, x_b, z_b, th_b = q

    #contacts: ee-b, b-g-front, b-g-back
    #x, z for each contact
    if mode == :ee_b
        j = SMatrix{2,5}([-1.0 0.0 1.0 0.0 -(x_b-x_ee)*tan(th_b);
                          0.0 -1.0 0.0 1.0 (x_b-x_ee)])
        return j
        
    elseif mode == :b_g
        j = SMatrix{4,5}([0.0 0.0 1.0 0.0 -(x_b-(model.xlen_block/2))*tan(th_b);
                          0.0 0.0 0.0 1.0 (x_b-(model.xlen_block/2));
                          0.0 0.0 1.0 0.0 -(x_b+(model.xlen_block/2))*tan(th_b);
                          0.0 0.0 0.0 1.0 (x_b-(model.xlen_block/2))])
        return j
    # elseif mode == :ee_g
    #     j = SMatrix{4,5}([1.0 0.0 0.0 0.0 0.0;
    #                       0.0 1.0 0.0 0.0 0.0;
    #                       1.0 0.0 0.0 0.0 0.0;
    #                       0.0 1.0 0.0 0.0 0.0])
    #     return j

    end
end

# contact Jacobian
function J_func(model::Block, env::Environment, q)
    return SMatrix{6, 5}([_jacobian(model, q, mode=:ee_b);
                          _jacobian(model, q, mode=:b_g);])
                          # _jacobian(model, q, mode=:ground);])
end

# translates the two variables normal force (γ) and tangential forces (b) into a single vector for jacobian
function contact_forces(model::Block, env::Environment{<:World, LinearizedCone}, γ1, b1, q2, k)
    # γ1: force vector (size num contacts)
    # b1: idk but size 2*num contacts
    # k: also size 2*num contacts
    
    m = friction_mapping(env) # what is this
    SVector{6}([transpose(rotation(env, k[1:1])) * [m * b1[1:2]; γ1[1]];
                 transpose(rotation(env, k[3:3])) * [m * b1[3:4]; γ1[2]];
                 transpose(rotation(env, k[5:5])) * [m * b1[5:6]; γ1[3]];])
                 # transpose(rotation(env, k[7:7])) * [m * b1[7:8]; γ1[4]];
                 # transpose(rotation(env, k[9:9])) * [m * b1[9:10]; γ1[5]];
                 # transpose(rotation(env, k[11:11])) * [m * b1[11:12]; γ1[6]];])
end

function velocity_stack(model::Block, env::Environment{<:World, LinearizedCone}, q1, q2, k, h)
    v = J_func(model, env, q2) * (q2 - q1) / h[1]
    v1_surf = rotation(env, k[1:1]) * v[1:2]
    v2_surf = rotation(env, k[3:3]) * v[3:4]
    v3_surf = rotation(env, k[5:5]) * v[5:6]
    # v4_surf = rotation(env, k[7:7]) * v[7:8]
    # v5_surf = rotation(env, k[9:9]) * v[9:10]
    # v6_surf = rotation(env, k[11:11]) * v[11:12]
    
    SVector{6}([transpose(friction_mapping(env)) * v1_surf[1];
                transpose(friction_mapping(env)) * v2_surf[1];
                transpose(friction_mapping(env)) * v3_surf[1];])
                # transpose(friction_mapping(env)) * v4_surf[1];
                # transpose(friction_mapping(env)) * v5_surf[1];
                # transpose(friction_mapping(env)) * v6_surf[1];])
end


# Working Parameters
params = YAML.load_file(joinpath(@__DIR__,"params.yaml"))

# nq, nu, nw, nc, m, g, m_block, μ_world, μ_block, r_ee, xlen_block, ylen_block, zlen_block
                         
block_system = Block(5, 2, 2, 3,
                     params["m_ee"], params["gravity"], params["m_block"],
                     params["mu_ground"], params["mu_block"],
                     params["r_ee"], params["xlen_block"], params["ylen_block"], params["zlen_block"],
	             BaseMethods(), DynamicsMethods(),
	             SVector{5}(zeros(5)), # joint friction
                     SVector{2}([-10, -10]), # u min
                     SVector{2}([10, 10]),   # u max
                     SVector{5}([-1,0, 0,-5,-5]), # q min (x,y,z, xtray,ytray,thtray)
                     SVector{5}([1,1.5, 5,5,5]),   # q max
)

function friction_coefficients(model::Block) 
	return [model.μ_world]
end

function initialize_z!(z, model::Block, idx::RoboDojo.IndicesZ, q)
    z .= 1.0
    z[idx.q] .= q
end
