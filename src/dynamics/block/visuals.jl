using Colors
using CoordinateTransformations
using FileIO
using GeometryBasics
using MeshCat, MeshIO, Meshing
using Rotations
using YAML

params = YAML.load_file(joinpath(@__DIR__,"params.yaml"))

function plot_lines!(vis::Visualizer, model::Block, q::AbstractVector;
		r=0.05, size=10, name::Symbol=:waiter_2D, col::Bool=true)
	# Point Traj
	com_point = Vector{Point{3,Float64}}()

	for qi in q
		push!(com_point, Point(qi[1], 0.0, qi[2] + r))
	end

	# Set lines
	orange_mat, blue_mat, black_mat = get_line_material(size)
	if col
		setobject!(vis[name]["lines"]["com"], MeshCat.Line(com_point, orange_mat))
	else
		setobject!(vis[name]["lines"]["com"], MeshCat.Line(com_point, black_mat))
	end
	return nothing
end

function build_robot!(vis::Visualizer, model::Block; name::Symbol=:Block, r=params["r_ee"], α=1.0)
    nc = model.nc
    r = convert(Float32, r)
    
    body_mat = MeshPhongMaterial(color = RGBA(13/255, 152/255, 186/255, α))
    contact_mat = MeshPhongMaterial(color = RGBA(1.0, 165/255, 0.0, α))
    block_mat = MeshPhongMaterial(color = RGBA(99/255, 97/255, 93/255, 1.0))
    support_mat = MeshPhongMaterial(color = RGBA(0.7, 0.7, 0.7, 1.0))

    default_background!(vis)

    setobject!(vis[name][:robot]["ee"],
               GeometryBasics.Cylinder(GeometryBasics.Point3f0(0,0,-d/2),
                                       GeometryBasics.Point3f0(0,0,d/2),
                                       convert(Float32, r)),
	       body_mat)

    setobject!(vis[name][:object]["block"],
               GeometryBasics.Rect3f(GeometryBasics.Point3f0(0,0,0),
                                       GeometryBasics.Vec3f0(params["xlen_block"], params["ylen_block"], params["zlen_block"]),
               block_mat))
    return nothing
end

function set_robot!(vis::Visualizer, model::Block, q::AbstractVector; name::Symbol=:Block, r=params["r_ee"])
    r = convert(Float32, r)
    
    settransform!(vis[name][:robot]["ee"], Translation(q[1], 0.0, q[2]))

    block_pos = Translation(q[3], 0.0, q[4])
    block_rot = LinearMap(RotY(q[5]))
    block_tf = compose(block_pos, block_rot)
    settransform!(vis[name][:object]["block"], block_tf)
    
    return nothing
end

function contact_point(model::Block, q::AbstractVector; r=params["r_ee"])
    xlen = params["xlen_block"]
    zlen = params["zlen_block"]
               
    pt_ee = [q[1], 0.0, q[2]]

    pt_block1 = [q[3]+(xlen/2), 0.0, q[4]-(zlen/2)]
    pt_block2 = [q[3]-(xlen/2), 0.0, q[4]-(zlen/2)]

    pc = [pt_ee, pt_block1, pt_block2]
    
    return pc
end
