import copy
import json
import helperclasses as hc
import geometry as geom
import geometry.simple_geometry as geom_sg
import geometry.mesh as geom_mesh
import geometry.hierarchy as geom_h
import scene
import glm

# Ported from C++ by Melissa Katz
# Adapted from code by Loïc Nassif and Paul Kry


def populateVec(array: list):
    if array is None:
        return None
    if len(array) == 3:
        return glm.vec3(array[0], array[1], array[2])
    elif len(array) == 4:
        return glm.vec4(array[0], array[1], array[2], array[3])


def get_or(obj, keys, default, msg=""):
    if not isinstance(keys, list) or len(keys) == 0:
        try:
            return obj[keys]
        except KeyError:
            if msg != "":
                print(msg)
            return default

    result = []
    try:
        for key in keys:
            result.append(obj[key])
        return result
    except KeyError:
        if msg != "":
            print(msg)
        return default


def load_scene(infile):
    print("Parsing file:", infile)
    f = open(infile)
    data = json.load(f)

    # Loading camera
    cam_pos = populateVec(data["camera"]["position"])
    cam_lookat = populateVec(data["camera"]["lookAt"])
    cam_up = populateVec(data["camera"]["up"])
    cam_fov = data["camera"]["fov"]

    # Loading resolution
    width, height = get_or(data, "resolution", [1080, 720], "No resolution found, defaulting to 1080x720.")

    # Loading ambient light
    ambient = populateVec(get_or(data, "ambient", [0, 0, 0], "No ambient light defined, defaulting to [0, 0, 0]"))

    # Loading Anti-Aliasing options
    try:
        jitter = data["AA"]["jitter"]
        samples = data["AA"]["samples"]
    except KeyError:
        print("No Anti-Aliasing options found, setting to default")
        jitter = False
        samples = 1

    # loading depth of field options
    try:
        focal_length = data["DOF"]["focal_length"]
        aperture = data["DOF"]["aperture"]
        dof_samples = data["DOF"]["samples"]
    except KeyError:
        print("No Depth of Field options found, setting to default")
        focal_length = 1
        aperture = 0
        dof_samples = 1

    # loading motion blur options
    try:
        motion_time = data["motion"]["time"]
        motion_samples = data["motion"]["samples"]
        motion_final = data["motion"]["final"]
    except KeyError:
        print("No motion blur options found, setting to default")
        motion_time = 0
        motion_samples = 1
        motion_final = 0

    vc = hc.ViewportCamera() \
        .set_viewport(width, height) \
        .set_camera(cam_pos, cam_lookat, cam_up, cam_fov) \
        .set_lens(focal_length, aperture, dof_samples) \
        .set_motion(motion_time, motion_samples, motion_final)

    # Loading scene lights
    lights = []
    try:
        for light in data["lights"]:
            l_type = light["type"]
            l_name = light["name"]
            l_colour = populateVec(light["colour"])

            if l_type == "point":
                l_vector = populateVec(light["position"])
                l_power = light["power"]

            elif l_type == "directional":
                l_vector = populateVec(light["direction"])
                l_power = 1.0
            else:
                print("Unkown light type", l_type, ", skipping initialization")
                continue
            lights.append(hc.Light(l_type, l_name, l_colour, l_vector, l_power))
    except KeyError as e:
        print("Error loading lights: ", e)
        lights = []

    # Loading materials
    materials = []
    for material in data["materials"]:
        mat_name = material["name"]
        mat_id = material["ID"]

        mat_type = get_or(material, "type", "diffuse")
        mat_diffuse = populateVec(get_or(material, "diffuse", [0, 0, 0]))
        mat_specular = populateVec(get_or(material, "specular", [0, 0, 0]))
        mat_hardness = get_or(material, "hardness", 32)
        mat_tint = get_or(material, "tint", 0.0)
        mat_refr_index = get_or(material, "refr_index", 1.0)

        material = hc.Material(mat_name, mat_specular, mat_diffuse, mat_hardness, mat_id, mat_type, mat_tint)
        material.refr_index = mat_refr_index
        materials.append(material)

    # Loading geometry
    objects = []

    # Extra stuff for hierarchies
    rootNames = []
    roots = []
    for geometry in data["objects"]:
        parse_geometry(geometry, objects, rootNames, roots, materials)

    print("Parsing complete")
    sc = scene.Scene(vc, jitter, samples, ambient, lights, materials, objects)

    for obj in objects:
        obj.set_scene(sc)

    return sc


def parse_geometry(geometry, objects, rootNames, roots, materials):
    # Elements common to all objects: name, type, position, material(s)
    g_name = geometry["name"]
    g_type = geometry["type"]
    g_pos = populateVec(get_or(geometry, "position", [0, 0, 0]))
    g_mats = associate_material(materials, get_or(geometry, "materials", []))
    g_speed = populateVec(get_or(geometry, "speed", None))

    if add_basic_shape(g_name, g_type, g_pos, g_speed, g_mats, geometry, objects):
        # Non-hierarchies are straightforward
        return
    elif g_type == "node":
        g_ref = get_or(geometry, "ref", "")
        g_r = populateVec(get_or(geometry, "rotation", [0, 0, 0]))
        g_s = populateVec(get_or(geometry, "scale", [1, 1, 1]))
        g_hierarchy_type = get_or(geometry, "hierarchy_type", "union")

        if g_ref == "":
            # Brand-new hierarchy
            rootNames.append(g_name)
            node = geom_h.Hierarchy(g_name, g_type, g_mats, g_hierarchy_type, g_pos, g_r, g_s, g_speed)
            traverse_children(node, geometry["children"], materials, rootNames, roots, g_speed)
            roots.append(node)
            objects.append(node)
        else:
            # Hierarchy that depends on a previously defined one
            rid = -1
            for i in range(len(rootNames)):
                # Find hierarchy that this references
                if g_ref == rootNames[i]:
                    rid = i
                    break
            if rid != -1:
                node = copy.deepcopy(roots[rid])
                node.name = g_name
                node.materials = g_mats
                node.make_matrices(g_pos, g_r, g_s)
                objects.append(node)
            else:
                print("Node reference", g_ref, "not found, skipping creation")

    else:
        print("Unkown object type", g_type, ", skipping initialization")
        return


def add_basic_shape(g_name: str, g_type: str, g_pos: glm.vec3, g_speed, g_mats: list[hc.Material], geometry, objects: list[geom.Geometry]):
    # Function for adding non-hierarchies to a list, since there's nothing extra to do with them
    # Returns True if a shape was added, False otherwise
    if g_type == "sphere":
        g_radius = geometry["radius"]
        objects.append(geom_sg.Sphere(g_name, g_type, g_mats, g_pos, g_radius, g_speed))
    elif g_type == "plane":
        g_normal = populateVec(geometry["normal"])
        objects.append(geom_sg.Plane(g_name, g_type, g_mats, g_pos, g_normal, g_speed))
    elif g_type == "box":
        try:
            g_size = populateVec(geometry["size"])
            objects.append(geom_sg.AABB(g_name, g_type, g_mats, g_pos, g_size, g_speed))
        except KeyError:
            # Boxes can also be directly declared with a min and max position
            box = geom_sg.AABB(g_name, g_type, g_mats, g_pos, glm.vec3(0, 0, 0), g_speed)
            box.minpos = populateVec(geometry["min"])
            box.maxpos = populateVec(geometry["max"])
            objects.append(box)
    elif g_type == "mesh":
        g_path = geometry["filepath"]
        g_scale = geometry["scale"]

        g_flat_shaded = get_or(geometry, "flat_shaded", False)
        objects.append(geom_mesh.Mesh(g_name, g_type, g_mats, g_pos, g_scale, g_path, g_flat_shaded, g_speed))
    else:
        return False
    return True


def traverse_children(node: geom_h.Hierarchy, children, materials: list[hc.Material], rootNames, roots, speed):
    for geometry in children:
        # Obtain info common to all shapes like in the main body of the parser
        g_name = geometry["name"]
        g_type = geometry["type"]
        g_pos = populateVec(get_or(geometry, "position", [0, 0, 0]))
        g_mats = associate_material(materials, get_or(geometry, "materials", []))
        if speed is None:
            g_speed = None
        else:
            g_speed = speed + populateVec(get_or(geometry, "speed", [0, 0, 0]))

        if add_basic_shape(g_name, g_type, g_pos, g_speed, g_mats, geometry, node.children):
            # Nothing fancy to do for non-hierarchies
            continue
        elif g_type == "node":
            # Hierarchy within a hierarchy, recurse
            g_r = populateVec(get_or(geometry, "rotation", [0, 0, 0]))
            g_s = populateVec(get_or(geometry, "scale", [1, 1, 1]))
            g_hierarchy_type = get_or(geometry, "hierarchy_type", "union")
            inner = geom_h.Hierarchy(g_name, g_type, g_mats, g_hierarchy_type, g_pos, g_r, g_s, g_speed)
            node.children.append(inner)
            traverse_children(inner, geometry["children"], materials, rootNames, roots, speed)
        else:
            parse_geometry(geometry, node.children, rootNames, roots, materials)


def associate_material(mats: list[hc.Material], ids: list[int]):
    new_list = []
    for i in ids:
        for mat in mats:
            if i == mat.ID:
                new_list.append(mat)
    return new_list
