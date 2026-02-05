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
mutable struct Block1D{T} <: Model{T}
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

function lagrangian(model::Block1D, q, q̇)
    L = 0.0

    # EE
    L += 0.5 * model.m * transpose(q̇[1]) * q̇[1]

    # block
    L += 0.5 * model.m_block * transpose(q̇[2:3]) * q̇[2:3]
    # I_b = (1/12) * model.m_block * (model.zlen_block^2) * (model.xlen_block^2)
    # L += 0.5*I_b*(q̇[4]^2)
    L -= model.m_block * model.g * q[3]

    return L
end

function kinematics(model::Block1D, q; mode=:contacts)
    # supposed to return pose of each contact point
    # (why also defined in visuals.jl??)
    if mode == :contacts
        # option 1: EE location
        # ee = SVector{2}([q[1], q[2]-model.r])

        # option 2: project EE in -Z-world to block
        # solve for projected point from EE to block surface
        # thb = q[5]
        # d = model.zlen_block/2
        # l = (1/cos(thb)) * (q[1]-q[3]+(d*sin(thb)))
        # ee_projx = q[1] #q[3] - (d*sin(thb)) + (l*cos(thb))
        # ee_projz = q[4] + (d*cos(thb)) - (l*sin(thb))
        # ee = SVector{2}([ee_projx, ee_projz])

        # option 3: closest point to block (side) plane
        d = model.zlen_block/2
        w = model.xlen_block/2
        xe, xb, zb = q
        ze = 0.05
        # thb = -thb
        # l = (-sin(thb))*(xe-xb) + (-cos(thb))*(ze-zb) - w
        # ee = SVector{2}([xe,ze]+[(l*cos(thb)),(l*sin(thb))])

        ## block contact points
        # block1 = SVector{2}([xb,zb] + d*[sin(thb),-cos(thb)] + w*[cos(thb),sin(thb)])
        # block2 = SVector{2}([xb,zb] + d*[sin(thb),-cos(thb)] - w*[cos(thb),sin(thb)])
        # block = SVector{2}([xb,zb] + d*[sin(thb),-cos(thb)])


        ## NO ROTATION:
        ee = SVector{2}([xb-w, zb])
        block = SVector{2}([xb, zb-d])
        
        return SVector{4}([ee; block;])
    elseif mode == :ee
        return q[1]
    elseif mode == :block
        return q[2:3] #4]
    else
        @error "incorrect mode"
        return
    end
end

# mass matrix
function M_func(model::Block1D, q)
    m = model.m
    mb = model.m_block
    h_b = model.zlen_block
    w_b = model.xlen_block
    I_b = (1/12)*mt*(h_b^2 + w_b^2)

    Diagonal(@SVector [m, mb, mb]) #, I_b])
end

# gravity
function C_func(model::Block1D, q, q̇)
    m = model.m
    mb = model.m_block
    g = model.g

    @SVector [0.0, 0.0, mb*g] #, 0.0]
end

function dist_block(model::Block1D, p, pt)
    # p: [x, z], pt: [xblock, zblock, θblock]

    differ = p-pt[1:2]
    # beta = atan(diff[1]/diff[2]) + pt[3] # angle between vector and tray-vertical
    # R = [cos(pt[3]) -sin(pt[3]); sin(pt[3]) cos(pt[3])];
    # xdiff,zdiff = R*(differ)
    xdiff, zdiff = differ
    ## just halfplane
    # zdist = zdiff - model.zlen_block/2 - model.r
    
    ## full
    xdiff = abs(xdiff)
    zdiff = abs(zdiff)        
    zdist = zdiff-(model.zlen_block/2)
    xdist = xdiff-(model.xlen_block/2)

    dist_neg = min(0, max(xdist, zdist))
    
    zdist = max(0, zdist)
    xdist = max(0, xdist)
    
    return norm([xdist, zdist]) + dist_neg
   # return zdist
end

# signed distance function
function ϕ_func(model::Block1D, env::Environment, q)
    # ee-block, block_front-ground, block_back-ground
    cp = kinematics(model, q, mode=:contacts)
    ee = [q[1], 0.05]
    block = cp[3:4]
    # block2 = cp[5:6]
    block_q = q[2:3] #4]    

    ee_block_dist = dist_block(model, ee, block_q) - model.r
    # block1_dist = block1[2] #- (model.zlen_block/2)
    # block2_dist = block2[2] #- (model.zlen_block/2)
    block_dist = block[2]
    
    return SVector{2}([ee_block_dist; block_dist;])
end

# control Jacobian
function B_func(model::Block1D, q)
    B = zeros(1,model.nq)
    B[1,1] = 1
    B = SMatrix{1,model.nq}(B)
    return B
end

# disturbance Jacobian
function A_func(model::Block1D, q)
    A = zeros(model.nq,1)
    A[1,1] = 1
    A = SMatrix{model.nq,1}(A)
    return A
end

function _jacobian(model::Block1D, q, dists; mode=:ee_b)
    # J'λ = \tau
    # ee_b := EE-block contact
    # b_g := block-ground contact
    # ee_g := EE-ground contact

    x_ee, x_b, z_b = q
    z_ee = 0.05
    # rot_th = [cos(th_b) sin(th_b);
    #           -sin(th_b) cos(th_b)]
    # x_ee_b, z_ee_b = rot_th * [(x_ee-x_b); (z_ee-z_b)]

    #contacts: ee-b, b-g-front, b-g-back
    # tangent, normal for each contact
    # remember Y axis is FLIPPED! so CW is positive theta
    if mode == :ee_b
        # convention: tangent faces up (+z) and normal faces back (-x)
        # j = SMatrix{2,4}([-sin(th_b) sin(th_b) cos(th_b) 0; #(model.xlen_block/2);
        #                   -cos(th_b) cos(th_b) -sin(th_b) 0;]) #-z_ee_b;])
        j = SMatrix{2,3}([0 0 0;
                          -1 1 0;])
        return j
        
    elseif mode == :b_g
        
        j = SMatrix{2,3}([0.0 1.0 0.0;
                          0.0 0.0 1.0;])
                          # 0.0 1.0 0.0 -(model.zlen_block/2)*cos(th_b);
                          # 0.0 0.0 1.0 (model.xlen_block/2)*cos(th_b)])
        return j
    end
end

# contact Jacobian
function J_func(model::Block1D, env::Environment, q)
    dists = ϕ_func(model::Block1D, env::Environment, q)
    return SMatrix{4, 3}([_jacobian(model, q, dists, mode=:ee_b);
                          _jacobian(model, q, dists, mode=:b_g);])
end

# translates the two variables normal force (γ) and tangential forces (b) into a single vector for jacobian
function contact_forces(model::Block1D, env::Environment{<:World, LinearizedCone}, γ1, b1, q2, k)
    # γ1: normal force vector (size num contacts)
    # b1: tangent force vector
    # k: contact point coordinates from kinematics()
    
    m = friction_mapping(env) # what is this
    SVector{4}([transpose(rotation(env, k[1:1])) * [m * b1[1:2]; γ1[1]];
                  transpose(rotation(env, k[3:3])) * [m * b1[3:4]; γ1[2]];])
                 # transpose(rotation(env, k[5:5])) * [m * b1[5:6]; γ1[3]];])
end

function velocity_stack(model::Block1D, env::Environment{<:World, LinearizedCone}, q1, q2, k, h)
    v = J_func(model, env, q2) * (q2 - q1) / h[1]
    v1_surf = rotation(env, k[1:1]) * v[1:2]
    v2_surf = rotation(env, k[3:3]) * v[3:4]
    # v3_surf = rotation(env, k[5:5]) * v[5:6]

    SVector{4}([transpose(friction_mapping(env)) * v1_surf[1];
                transpose(friction_mapping(env)) * v2_surf[1];])
                # transpose(friction_mapping(env)) * v3_surf[1];])
end

function load_params()
    params_path = abspath(joinpath(@__DIR__,"params.yaml"))
    params = YAML.load_file(params_path)
    return params
end

# Working Parameters
params = load_params()

# nq, nu, nw, nc, m, g, m_block, μ_world, μ_block, r_ee, xlen_block, ylen_block, zlen_block
                         
block_system_1D = Block1D(3, 1, 1, 2, #4,
                     params["m_ee"], params["gravity"], params["m_block"],
                     params["mu_ground"], params["mu_block"],
                     params["r_ee"], params["xlen_block"], params["ylen_block"], params["zlen_block"],
	             BaseMethods(), DynamicsMethods(),
	             SVector{3}(zeros(3)), # joint friction
                     SVector{1}([-10]), # u min
                     SVector{1}([10]),   # u max
                     SVector{3}([-1, -5,-5]), # q min (x,z, xblock,zblock,thblock)
                     SVector{3}([1, 5,5]),   # q max
)

function friction_coefficients(model::Block1D) 
	return [model.μ_world]
end

function initialize_z!(z, model::Block1D, idx::RoboDojo.IndicesZ, q)
    z .= 1.0
    z[idx.q] .= q
end
