using Colors
using CoordinateTransformations
using FileIO
using GeometryBasics
using MeshCat, MeshIO, Meshing
using Rotations

function plot_lines!(vis::Visualizer, model::Waiter2D, q::AbstractVector;
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

function build_robot!(vis::Visualizer, model::Waiter2D; name::Symbol=:Waiter2D, d=0.01, r=0.0725, α=1.0)
    nc = model.nc
    r = convert(Float32, r)
    d = convert(Float32, d)
    support_pos = [0.8; 0.0; 0.447]
    
    body_mat = MeshPhongMaterial(color = RGBA(13/255, 152/255, 186/255, α))
    contact_mat = MeshPhongMaterial(color = RGBA(1.0, 165/255, 0.0, α))
    tray_mat = MeshPhongMaterial(color = RGBA(99/255, 97/255, 93/255, 1.0))
    support_mat = MeshPhongMaterial(color = RGBA(0.7, 0.7, 0.7, 1.0))

    default_background!(vis)

    setobject!(vis[name][:robot]["plate"],
               GeometryBasics.Cylinder(GeometryBasics.Point3f(0,0,-d/2),
                                       GeometryBasics.Point3f(0,0,d/2),
                                       convert(Float32, r)),
	       body_mat)
    setobject!(vis[name][:robot]["contact_front"], Sphere(Point3f0(0.0),0.01),contact_mat)
    setobject!(vis[name][:robot]["contact_back"], Sphere(Point3f0(0.0),0.01),contact_mat)

    setobject!(vis[name][:object]["tray"],
               GeometryBasics.Cylinder(GeometryBasics.Point3f(0,0,-0.011),
                                       GeometryBasics.Point3f(0,0,0.011),
                                       convert(Float32, 0.2286)),
               tray_mat)
    
    setobject!(vis[name][:env]["support"]["init"],
               Rect(Vec(0,0,0),Vec(0.4, 0.8, 0.05)),wall_mat)
    settransform!(vis[name][:env]["support"]["init"], Translation(support_pos))
    
    return nothing
end

function set_robot!(vis::Visualizer, model::Waiter2D, q::AbstractVector; name::Symbol=:Waiter2D, d=0.01, r=0.0725)
    r = convert(Float32, r)
    
    settransform!(vis[name][:robot]["plate"], Translation(q[1], 0.0, q[2]))
    settransform!(vis[name][:robot]["contact_front"], Translation(q[1]+r, 0.0, q[2]))
    settransform!(vis[name][:robot]["contact_back"], Translation(q[1]-r, 0.0, q[2]))

    tray_pos = Translation(q[3], 0.0, q[4])
    tray_rot = LinearMap(RotY(q[5]))
    tray_tf = compose(tray_pos, tray_rot)
    settransform!(vis[name][:object]["tray"], tray_tf)
    
    return nothing
end

function contact_point(model::Waiter2D, q::AbstractVector; r=0.0725)
    p_front = [q[1]+r, 0.0, q[2]]
    p_back = [q[1]-r, 0.0, q[2]]
    pc = [p_back, p_front]
    
    return pc
end
