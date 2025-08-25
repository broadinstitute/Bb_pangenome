import bpy
from math import acos, pi
from mathutils import Vector
import json

def hex_to_rgba(hex_color, alpha=1.0):
    """Convert hex color to RGBA tuple with values from 0 to 1."""
    hex_color = hex_color.lstrip('#')
    rgb = tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))
    return rgb + (alpha,)

def create_material(name, hex_color):
    """Create a new material with the given name and hex color."""
    material = bpy.data.materials.new(name=name)
    rgba = hex_to_rgba(hex_color)
    material.use_nodes = True
    principled = material.node_tree.nodes.get("Principled BSDF")
    if principled:
        principled.inputs["Base Color"].default_value = rgba

        # Set transparency for all materials except gray and light gray
        if name not in ["gray", "light_gray"]:
            if name == "clear":
                principled.inputs["Alpha"].default_value = 0.1
            else:
                principled.inputs["Alpha"].default_value = 0.95
            material.blend_method = 'BLEND'

    return material

def draw_network(network, edge_thickness=0.25, node_size=1, directed=False):
    """Takes assembled network/molecule data and draws to blender"""

    # Add some mesh primitives
    bpy.ops.object.select_all(action='DESELECT')
    bpy.ops.mesh.primitive_uv_sphere_add()
    sphere = bpy.context.object
    bpy.ops.mesh.primitive_cylinder_add()
    cylinder = bpy.context.object
    bpy.ops.mesh.primitive_torus_add()
    torus = bpy.context.object
    bpy.ops.mesh.primitive_cube_add()
    cube = bpy.context.object

    # Keep references to all nodes and edges
    shapes = []
    # Keep separate references to shapes to be smoothed
    shapes_to_smooth = []

    ## Create materials for edges
    #edge_material = create_material("edge", "#767676")  # Light gray
    #cylinder.active_material = edge_material
    #cone.active_material = edge_material

    # Draw nodes
    for key, node in network["nodes"].items():

        # Use the color specified in the node data
        color = node.get("color", "#0083FF")  # Default to blue if color is not specified

        # Create a new material for this node
        node_material = create_material(f"node_{key}", color)
        if node['plasmid'] is None:
            pass
        # if chromosome use cube, else if cp use torus, else if lp use cylinder, else ico
        elif 'lp' in node['plasmid']:
            # Copy mesh primitive and edit to make node
            node_cube = cube.copy()
            node_cube.data = cube.data.copy()
            node_cube.location = node["location"]
            node_cube.dimensions = Vector((0.1, 0.1, 0.75))
            node_cube.active_material = node_material
            bpy.context.collection.objects.link(node_cube)
            shapes.append(node_cube)
            shapes_to_smooth.append(node_cube)
        elif 'cp' in node['plasmid']:
            # Copy mesh primitive and edit to make node
            node_torus = torus.copy()
            node_torus.data = torus.data.copy()
            node_torus.location = node["location"]
            node_torus.dimensions = Vector((0.5,0.5,0.125))
            node_torus.active_material = node_material
            bpy.context.collection.objects.link(node_torus)
            shapes.append(node_torus)
            shapes_to_smooth.append(node_torus)
        elif 'chromosome' in node['plasmid']:
            # Copy mesh primitive and edit to make node
            node_cylinder = cylinder.copy()
            node_cylinder.data = cylinder.data.copy()
            node_cylinder.location = node["location"]
            node_cylinder.dimensions = Vector((1,1,2.5))
            node_cylinder.active_material = node_material
            bpy.context.collection.objects.link(node_cylinder)
            shapes.append(node_cylinder)
            shapes_to_smooth.append(node_cylinder)

    ## Draw edges
    #for edge in network["edges"]:
    #    source_loc = network["nodes"][f'{edge["source"]}']["location"]
    #    target_loc = network["nodes"][f'{edge["source"]}']["location"]
    #
    #    diff = [c2 - c1 for c2, c1 in zip(source_loc, target_loc)]
    #    cent = [(c2 + c1) / 2 for c2, c1 in zip(source_loc, target_loc)]
    #    mag = sum([(c2 - c1) ** 2 for c1, c2 in zip(source_loc, target_loc)]) ** 0.5
    #
    #    # Euler rotation calculation
    #    v_axis = Vector(diff).normalized()
    #    v_obj = Vector((0, 0, 1))
    #    v_rot = v_obj.cross(v_axis)
    #    angle = acos(v_obj.dot(v_axis))
    #
    #    # Copy mesh primitive to create edge
    #    edge_cylinder = cylinder.copy()
    #    edge_cylinder.data = cylinder.data.copy()
    #    edge_cylinder.dimensions = [edge_thickness] * 2 + [mag - node_size]
    #    edge_cylinder.location = cent
    #    edge_cylinder.rotation_mode = "AXIS_ANGLE"
    #    edge_cylinder.rotation_axis_angle = [angle] + list(v_rot)
    #    bpy.context.collection.objects.link(edge_cylinder)
    #    shapes.append(edge_cylinder)
    #    shapes_to_smooth.append(edge_cylinder)
    #
    #    # Copy another mesh primitive to make an arrow head
    #    if directed:
    #        arrow_cone = cone.copy()
    #        arrow_cone.data = cone.data.copy()
    #        arrow_cone.dimensions = [edge_thickness * 4.0] * 3
    #        arrow_cone.location = cent
    #        arrow_cone.rotation_mode = "AXIS_ANGLE"
    #        arrow_cone.rotation_axis_angle = [angle + pi] + list(v_rot)
    #        bpy.context.collection.objects.link(arrow_cone)
    #        shapes.append(arrow_cone)

    # Remove primitive meshes
    bpy.ops.object.select_all(action='DESELECT')
    sphere.select_set(True)
    cylinder.select_set(True)
    cube.select_set(True)
    torus.select_set(True)

    # If the starting cube is there, remove it
    if "Cube" in bpy.data.objects.keys():
        bpy.data.objects.get("Cube").select_set(True)
    bpy.ops.object.delete()
    if "Sphere" in bpy.data.objects.keys():
        bpy.data.objects.get("Sphere").select_set(True)
    bpy.ops.object.delete()
    if "Cylinder" in bpy.data.objects.keys():
        bpy.data.objects.get("Cylinder").select_set(True)
    bpy.ops.object.delete()
    if "Torus" in bpy.data.objects.keys():
        bpy.data.objects.get("Torus").select_set(True)
    bpy.ops.object.delete()

    # Smooth specified shapes
    for shape in shapes_to_smooth:
        shape.select_set(True)
    bpy.context.view_layer.objects.active = shapes_to_smooth[0]
    bpy.ops.object.shade_smooth()

    # Join shapes
    for shape in shapes:
        shape.select_set(True)
    bpy.context.view_layer.objects.active = shapes[0]
    bpy.ops.object.join()

    # Center object origin to geometry
    bpy.ops.object.origin_set(type="ORIGIN_GEOMETRY", center="MEDIAN")

    # Refresh scene
    bpy.context.view_layer.update()

# If main, load json and run
if __name__ == "__main__":
    with open("network.json") as network_file:
        network = json.load(network_file)
    draw_network(network)
