import re
import sys
import json
import numpy
import math

# --------------------------------------------------------------#
# variables

EPSILON = 0.01

# --------------------------------------------------------------#
# functions


def load_obj(filename):
    """function created by https://github.com/xtompok"""
    header = []
    vertices = []
    objects = {}
    with open(filename) as f:
        # Parse header
        while True:
            line = f.readline()
            if line[0] == 'v':
                break
            header.append(line)

        # Parse vertices
        while True:
            if line[0] == 'o':
                break
            if len(line) < 7 or line[0] != 'v':
                line = f.readline()
                continue
            (_, x, y, z) = line.strip().split(' ')
            vertices.append([float(x), float(y), float(z)])
            line = f.readline()

        # Parse objects
        aid = None
        faces = []
        while True:
            if line == '':
                if aid:
                    objects[aid] = faces
                break
            if len(line) < 7 or line[0] not in ('o', 'f'):
                print("Skipping line: {}".format(line))
                line = f.readline()
                continue
            if line[0] == 'o':
                if aid:
                    objects[aid] = faces
                aid = line.strip().split(' ')[1][2:-2]
                faces = []
                line = f.readline()
                continue
            preface = line.strip().split(" ")[1:]
            faces.append(list(map(int, preface)))
            line = f.readline()
    return header, vertices, objects

def create_center_groups(center_matrix):
    centers ={}
    for key_cm, items_cm in center_matrix.items():
        if items_cm !=0:
            b = key_cm
            b = re.sub('center_matrix_', '', b)
            c = re.sub('\'\)$', '', b)
            d = re.sub('\(\'','',c)
            e = re.split("\',_\'", d)
            e0= str(e[0])
            e1= str(e[1])
            pe0=e0.replace("_", "-")
            pe1=e1.replace("_", "-")

            #id_change = (f[1])
            centers[pe1] = pe0
    return centers

def remove_duplicate(duplicate_list):
    """Removes duplicates values in the input list"""
    # in list of lists removing duplicate list include
    final_list = []
    for num in duplicate_list:
        if num not in final_list:
            final_list.append(num)
    return final_list


def hight_of_object(coords):
    """Finding the high of the imput object and return it, return max z coordinates of object also"""
    zcoords_duplicate = []
    zcoords = []
    for i in coords:
        zcoords_duplicate.append(i[2])
        zcoords = remove_duplicate(zcoords_duplicate)

    zcoords.sort()
    c = abs(zcoords[1] - zcoords[0])
    # if there is 3 variables in list find first one (because gabled roof have wall up to roof)[-1]
    # but we need coord [1]
    coords_max = zcoords[1]  # return max z coords
    return c, coords_max

def parameters_of_footprint(footprintcoords):
    """In the list of lists of foodprint coordinates return side a,b and its area"""
    xcoords_duplicate = []
    for i in footprintcoords:
        xcoords_duplicate.append(i[0])
    xcoords = remove_duplicate(xcoords_duplicate)

    ycoords_duplicate = []
    for i in footprintcoords:
        ycoords_duplicate.append(i[1])
    ycoords = remove_duplicate(ycoords_duplicate)

    if len(xcoords) == 2 & len(ycoords) == 2:
        a = abs(xcoords[0] - xcoords[1])
        b = abs(ycoords[0] - ycoords[1])
        area = a * b
    elif 2 < len(xcoords) <= 4 or 2 < len(ycoords) <= 4:
        a1 = abs(xcoords[0] - xcoords[1])
        a2 = abs(xcoords[1] - xcoords[2])
        b1 = abs(ycoords[0] - ycoords[1])
        b2 = abs(ycoords[1] - ycoords[2])
        a = math.sqrt((a1 * a1) + (b1 * b1))
        b = math.sqrt((a2 * a2) + (b2 * b2))
        area = a * b
    else:
        print("Polygon contain more than 4 vertices, not valid")

    zcoords_duplicate = []
    for i in footprintcoords:
        zcoords_duplicate.append(i[2])
    zcoords = remove_duplicate(zcoords_duplicate)
    z = zcoords[0]
    return a, b, z, area

def body_volume(wallcoords, footprintcoords):
    """Return a volume of body building"""
    c, coords_max = hight_of_object(wallcoords)
    a, b, z, area = parameters_of_footprint(footprintcoords)
    volume = a * b * c
    return volume

def roof_volume(footprintcoords, rooftype, roofcoords):
    """Return a volume of the roof"""
    roofvolume = 0
    list_coords = []
    roof_up = []
    roof_down = []
    list_coords_dupl = []
    c = 0
    coords_max=0
    # for roof in rooftype:
    if rooftype == "Flat":
        roofvolume = 0
        c = 0
        coords_max =0
    elif rooftype == 'Shed':
        c, coords_max = hight_of_object(roofcoords)
        a, b, z, area = parameters_of_footprint(footprintcoords)
        roofvolume = (a * b * c) / 2
    elif rooftype == 'Gabled':
        c, coords_max = hight_of_object(roofcoords)
        a, b, z, area = parameters_of_footprint(footprintcoords)
        roofvolume = (a * b * c) / 2
    elif rooftype == "Pyramidal":
        c, coords_max = hight_of_object(roofcoords)
        a, b, z, area = parameters_of_footprint(footprintcoords)
        roofvolume = (a * b * c) / 3
    elif rooftype == "Hipped":
        c, coords_max = hight_of_object(roofcoords)
        a, b, z, area = parameters_of_footprint(footprintcoords)
        
        for k in roofcoords:
            if k not in list_coords:
                list_coords.append(k)
        for u in list_coords:
            if u[2] == coords_max:
                roof_up.append(u)
            else:
                roof_down.append(u)  # nedavala jsem i

        au = abs(roof_up[0][0] - roof_up[1][0])
        bu = abs(roof_up[0][1] - roof_up[1][1])
        a_smaller = math.sqrt((au * au) + (bu * bu))

        ad = abs(roof_down[0][0] - roof_down[1][0])
        bd = abs(roof_down[0][1] - roof_down[1][1])
        a_longer = math.sqrt((ad * ad) + (bd * bd))

        if not a_longer == a:
            b = a
            a = a_longer
        num = (2 * a + a_smaller)
        roofvolume = (b * c * num) / 6

    else:
        print("unnamed roof type")
    return roofvolume, c,coords_max
# --------------------------------------------------------------#

solution_all = json.load(open(sys.argv[1]))

cm = "center_matrix_"

# create dictionary from solution of optimization
center_matrix = dict(filter(lambda item: cm in item[0], solution_all.items()))

centers = create_center_groups(center_matrix)

# --------------------------------------------------------------#
# beginning of my work

# substitutes integer roof type value for string value
def load_roof_type_names(roof_type_dict):
    for id, value in roof_type_dict.items():
        if value == 1:
            roof_type_dict[id] = "Pyramidal"
        if value == 2:
            roof_type_dict[id] = "Flat"
        if value == 3:
            roof_type_dict[id] = "Shed"
        if value == 4:
            roof_type_dict[id] = 'Gabled'
        if value == 5:
            roof_type_dict[id] = "Hipped"
    return roof_type_dict

# divides building into its parts, calculates volumes (body, roof)
def create_auxiliary_structures(roof_types, objects, vertices, volume = True, only_aggregates = False, parts_list = None):
    object_parts = {}
    volume_dict = {}
    for obj_id, faces in objects.items():
        if only_aggregates:
            if obj_id not in parts_list:
                continue
        wallcoords = []
        wall_faces = []
        roofcoords = []
        roof_faces = [] 
        for face in faces:
            is_footprint = True
            is_roof = True
            facecoords = []
            for point in face:
                z_coord = vertices[point-1][2]
                if z_coord == 0:
                    is_roof = False
                elif z_coord > 0:
                    is_footprint = False
                facecoords.append(vertices[point-1])
            if is_footprint == True:
                footprintcoords = facecoords
                footprint_face = face
            elif is_roof == True:
                for coords in facecoords:
                    if coords not in roofcoords:
                        roofcoords.append(coords)
                roof_faces.append(face)
            else:
                for coords in facecoords:
                    if coords not in wallcoords:
                        wallcoords.append(coords)
                wall_faces.append(face)
        object_parts[obj_id] = {"footprintcoords": footprintcoords, "footprint_face": footprint_face, "wallcoords": wallcoords, "wall_faces": wall_faces, "roofcoords": roofcoords, "roof_faces": roof_faces}

        if volume:
            volume_body = body_volume(wallcoords, footprintcoords)
            volume_roof, c, coords_max = roof_volume(footprintcoords, roof_types[obj_id], roofcoords)
            volume_dict[obj_id] = {"volume_body": volume_body, "volume_roof": volume_roof, "coords_max": coords_max}
    if volume:
        return object_parts, volume_dict
    else:
        return object_parts

# finds the point through which the axis of rotation will pass
def get_point_of_rotation(objects, vertices):
    rotation_point_coords = vertices[0][0]
    rotation_point = 1
    for faces in objects.values():
        for face in faces:
            for v in face:          
                if vertices[v-1][2] == 0:
                    x = vertices[v-1][0]
                    y = vertices[v-1][1]
                    if x < rotation_point_coords:
                        rotation_point_coords = x
                        rotation_point = v
    return rotation_point

# calculates the angle between the block of buildings and the x-axis
def get_object_orientation(object_parts):
    for parts in object_parts.values():
        footprintcoords = parts["footprintcoords"]
        a = footprintcoords[1]
        b = footprintcoords[2]
        x_axis_vector = [3, 0]
        if a[1] > b[1]:
            block_vector = [a[0] - b[0], a[1] - b[1]]
        else:
            block_vector = [b[0] - a[0], b[1] - a[1]]
        numerator = x_axis_vector[0]*block_vector[0] + x_axis_vector[1]*block_vector[1]
        nominator = math.sqrt(x_axis_vector[0]**2 + x_axis_vector[1]**2)*math.sqrt(block_vector[0]**2 + block_vector[1]**2)
        cos_phi = numerator/nominator
        angle = round(math.acos(cos_phi)/math.pi*180, 0)
        break

    if angle >= 90:
        angle -= 90
    
    return angle

# moves the block of buildings by certain vector and rotates; calculates roof orientation of altered roofs
def move_and_rotate_objects(vertices, rotation_point, angle, roof_orientation):
    fixing_vertex = vertices[rotation_point-1]
    move_vector = [10 - fixing_vertex[0], 10 - fixing_vertex[1]]
    counter = 0
    for vertex in vertices:
        x = vertex[0] + move_vector[0]
        y = vertex[1] + move_vector[1]
        vertices[counter][0] = x
        vertices[counter][1] = y
        counter += 1
    move_angle = 90 - angle
    counter = 0
    # rotation matrix
    move_matrix = numpy.array([[math.cos(math.radians(move_angle)), -math.sin(math.radians(move_angle))],[math.sin(math.radians(move_angle)), math.cos(math.radians(move_angle))]])
    for vertex in vertices:
        vector = numpy.array([vertex[0] - 10, vertex[1] - 10])
        new_position = numpy.dot(move_matrix, vector)
        x = new_position[0] + 10
        y = new_position[1] + 10
        vertices[counter][0] = x
        vertices[counter][1] = y
        counter += 1

    for obj_id, orientation in roof_orientation.items():
        new_orientation = orientation + move_angle
        if new_orientation > 179:
            new_orientation = 0
        roof_orientation[obj_id] = new_orientation
    return vertices, roof_orientation

# finds min and max values of vertices on two horizontal axes 
def get_global_extremes(vertices):
    x_min_global = vertices[0][0]
    x_max_global = x_min_global
    y_min_global = vertices[0][1]
    y_max_global = y_min_global
    for v in vertices:       
        x = v[0]
        y = v[1]
        if x < x_min_global:
            x_min_global = x
        if x > x_max_global:
            x_max_global = x
        if y < y_min_global:
            y_min_global = y
        if y > y_max_global:
            y_max_global = y
    return x_min_global, x_max_global, y_min_global, y_max_global

# aligns outer wall of the building
def align_object_boundary(object, vertices, change, x_min, x_max, y_min, y_max, x_min_global, x_max_global, y_min_global, y_max_global):
    for face in object:
        for vertex in face:
            v = vertices[vertex - 1]
            if change == "x_min":
                if v[0] == x_min:
                    vertices[vertex - 1][0] = x_min_global
            elif change == "x_max":
                if v[0] == x_max:
                    vertices[vertex - 1][0] = x_max_global
            elif change == "y_min":
                if v[1] == y_min:
                    vertices[vertex - 1][1] = y_min_global
            elif change == "y_max":
                if v[1] == y_max:
                    vertices[vertex - 1][1] = y_max_global
    return vertices

# aligns outer walls of block to the closest global extreme
def align_outer_boundary(x_min_global, x_max_global, y_min_global, y_max_global, object_parts, vertices, neighborhood, roof_types):
    for obj_id, neighbors in neighborhood.items():
        directions = list(neighbors.values())
        footprintcoords = object_parts[obj_id]["footprintcoords"]
        wall_faces = object_parts[obj_id]["wall_faces"]
        roof_faces = object_parts[obj_id]["roof_faces"]
        roof_type = roof_types[obj_id]
        x_min, x_max, y_min, y_max = get_global_extremes(footprintcoords)
        if "right" in directions and "left" in directions:
            if abs(y_min - y_min_global) < abs(y_max - y_max_global):
                change = ["y_min"]
            else:
                change = ["y_max"]
        elif "lower" in directions and "upper" in directions:
            if abs(x_min - x_min_global) < abs(x_max - x_max_global):
                change = ["x_min"]
            else:
                change = ["x_max"]
        elif "upper" in directions and "left" in directions:
            change = ["x_max", "y_min"]
        elif "lower" in directions and "left" in directions:
            change = ["x_max", "y_max"]
        elif "lower" in directions and "right" in directions:
            change = ["x_min", "y_max"]
        elif "upper" in directions and "right" in directions:
            change = ["x_min", "y_min"]
        else:
            keys = list(neighbors.keys())
            footprintcoords2 = object_parts[keys[0]]["footprintcoords"]
            x_min2, x_max2, y_min2, y_max2 = get_global_extremes(footprintcoords2)
            if "right" in directions or "left" in directions:
                if abs(y_min - y_min2) < abs(y_max - y_max2):
                    change = ["y_min"]
                else:
                    change = ["y_max"]
            elif "upper" in directions or "lower" in directions:
                if abs(x_min - x_min2) < abs(x_max - x_max2):
                    change = ["x_min"]
                else:
                    change = ["x_max"]



        for ch in change:
            vertices = align_object_boundary(wall_faces, vertices, ch, x_min, x_max, y_min, y_max, x_min_global, x_max_global, y_min_global, y_max_global)
            if roof_type == "Shed" or roof_type == "Gabled":
                vertices = align_object_boundary(roof_faces, vertices, ch, x_min, x_max, y_min, y_max, x_min_global, x_max_global, y_min_global, y_max_global)
    return vertices

# creates lists of aggregates
def create_aggregate_structures(centers):
    parts_list = []
    aggregation_parts = []
    for pe1, pe0 in centers.items():
        if pe1 != pe0:
            if pe1 not in parts_list:
                parts_list.append(pe1)
            if pe0 not in parts_list:
                parts_list.append(pe0)
            part = [pe1, pe0]
            counter = 0
            changed = False
            for p in aggregation_parts:
                if part[0] in p and part[1] not in p:              
                    aggregation_parts[counter].append(part[1])
                    changed = True
                elif part[1] in p and part[0] not in p:
                    aggregation_parts[counter].append(part[0])
                    changed = True
                counter += 1
            if changed == False:
                aggregation_parts.append(part)
    return parts_list, aggregation_parts

# gets coordinates of rooftop
def get_rooftop_coordinates(aggregation_parts_dict, ap, volume_dict):
    roof_up = []
    roofcoords = aggregation_parts_dict[ap]["roofcoords"]
    coords_max = volume_dict[ap]["coords_max"]
    for k in roofcoords:
        if k[2] == coords_max:
            roof_up.append(k)
    return roof_up

# detects the properties of the aggregate
def get_aggregate_parametres(parallel_aggregate, volume_dict, aggregation_parts_dict, neighborhood, roof_types):
    volume_body = 0
    volume_roof = 0
    new_obj_id = ""
    roof_up = []
    neighbor_relation = {}
    for ap in parallel_aggregate:
        roof_type = roof_types[ap]
        volume_body += volume_dict[ap]["volume_body"]
        volume_roof += volume_dict[ap]["volume_roof"]
        new_obj_id += "{} + ".format(ap)
        neighbors = neighborhood[ap]
        for key, value in neighbors.items():
            if key not in parallel_aggregate:
                neighbor_relation[key] = value

        roof_up_parts = get_rooftop_coordinates(aggregation_parts_dict, ap, volume_dict)
        for part in roof_up_parts:
            roof_up.append(part)

    new_obj_id = new_obj_id[0:len(new_obj_id)-3]
    return volume_body, volume_roof, new_obj_id, neighbor_relation, roof_type, roof_up

# divides the aggregate into parts with the same orientation 
def divide_aggregates(aggregation_parts, i, neighborhood):
    x_parallel_aggregate = []
    y_parallel_aggregate = []
    for ap in aggregation_parts[i]:
        neighbors = neighborhood[ap]
        for key, value in neighbors.items():
            if key in aggregation_parts[i]:
                if value == "right" or value == "left":
                    if key not in x_parallel_aggregate:
                        x_parallel_aggregate.append(key)
                elif value == "lower" or value == "upper":
                    if key not in y_parallel_aggregate:
                        y_parallel_aggregate.append(key)
    return x_parallel_aggregate, y_parallel_aggregate

# calculates coordinate extremes of the aggregate; detects the position of the aggregate in relation to horizontal axes and other buildings in the block
def get_aggregate_extremes_and_position(parallel_aggregate, object_parts, x_min_global, x_max_global, y_min_global, y_max_global):
    if len(parallel_aggregate) > 0:
        obj_id = parallel_aggregate[0]
        footprintcoords = object_parts[obj_id]["footprintcoords"]
        xmin, xmax, ymin, ymax = get_global_extremes(footprintcoords)
        obj_id = parallel_aggregate[1]
        footprintcoords2 = object_parts[obj_id]["footprintcoords"]
        xmin2, xmax2, ymin2, ymax2 = get_global_extremes(footprintcoords2)
        condition = abs(xmin + xmax - xmin2 - xmax2) < abs(ymin + ymax - ymin2 - ymax2)
        if condition:
            x_parallel = False
            if abs(xmin - x_min_global) < abs(xmax - x_max_global):
                side = "left"
            elif abs(xmin - x_min_global) > abs(xmax - x_max_global):
                side = "right"
        else:
            x_parallel = True
            if abs(ymin - y_min_global) < abs(ymax - y_max_global):
                side = "lower"
            elif abs(ymin - y_min_global) > abs(ymax - y_max_global):
                side = "upper"
        ap = parallel_aggregate[1:len(parallel_aggregate)]
        for obj_id in ap:
            footprintcoords = object_parts[obj_id]["footprintcoords"]
            xmin2, xmax2, ymin2, ymax2 = get_global_extremes(footprintcoords)
            if xmin2 < xmin:
                xmin = xmin2
            if xmax2 > xmax:
                xmax = xmax2
            if ymin2 < ymin:
                ymin = ymin2
            if ymax2 > ymax:
                ymax = ymax2
    return xmin, xmax, ymin, ymax, x_parallel, side

# calculates length of the aggregate
def calculate_aggregate_length(x_parallel, side, x_min, x_max, y_min, y_max):
    if x_parallel:
        if side == "lower":
            lower_vertex = [x_min, y_min, 0]
            higher_vertex = [x_max, y_min, 0]
        elif side == "upper":
            lower_vertex = [x_min, y_max, 0]
            higher_vertex = [x_max, y_max, 0]
        length = x_max - x_min
    else:
        if side == "left":
            lower_vertex = [x_min, y_min, 0]
            higher_vertex = [x_min, y_max, 0]
        elif side == "right":
            lower_vertex = [x_max, y_min, 0]
            higher_vertex = [x_max, y_max, 0]
        length = y_max - y_min

    return length, higher_vertex, lower_vertex

# calculates new width of the aggregate
def calculate_new_width(wallcoords, length, volume_body, neighbor_relation, object_parts, lower_vertex, x_parallel, side):
    collision = False
    height = 0
    height_2 = 0
    for coords in wallcoords:
        if coords[2] != 0:
            if coords[2] != height:
                height = coords[2]
            else:
                height_2 = coords[2]
            if height != 0 and height_2 != 0:
                break
    if height_2 != 0:
        if height_2 < height:
            height = height_2
    new_width = volume_body/height/length

    directions = list(neighbor_relation.keys())
    if x_parallel:
        if side == "lower":
            key = "upper"
            if key in directions:
                collision = True               
        elif side == "upper":
            key = "lower"
            if key in directions:
                collision = True          
    else:
        if side == "left":
            key = "right"
            if key in directions:
                collision = True               
        elif side == "right":
            key = "left"
            if key in directions:
                collision = True                      
    
    if collision:
        obj_id = neighbor_relation[key]
        footprintcoords = object_parts[obj_id]["footprintcoords"]
        xmin, xmax, ymin, ymax = get_global_extremes(footprintcoords)
        if x_parallel:
            neighbor_length = ymax - ymin
        else:
            neighbor_length = xmax - xmin

        if key == "upper":
            new_wall_position = lower_vertex[1] + new_width
            new_length = ymax - new_wall_position
            if new_length/neighbor_length < 0.7:
                new_width = 0.3*neighbor_length + ymin - lower_vertex[1]
                height = volume_body/length/new_width
        elif key == "lower":
            new_wall_position = lower_vertex[1] - new_width
            new_length = new_wall_position - ymin
            if new_length/neighbor_length < 0.7:
                new_width = 0.3*neighbor_length - ymax + lower_vertex[1]
                height = volume_body/length/new_width
        elif key == "right":
            new_wall_position = lower_vertex[0] + new_width
            new_length = xmax - new_wall_position
            if new_length/neighbor_length < 0.7:
                new_width = 0.3*neighbor_length + xmin - lower_vertex[0]
                height = volume_body/length/new_width
        elif key == "left":
            new_wall_position = lower_vertex[0] - new_width
            new_length = new_wall_position - xmin
            if new_length/neighbor_length < 0.7:
                new_width = 0.3*neighbor_length - xmax + lower_vertex[0]
                height = volume_body/length/new_width
    
    return new_width, height

# calculates the body coordinates of the building from the newly acquired dimensions 
def compute_new_body_coordinates(x_parallel, lower_vertex, higher_vertex, new_width, height, side):
    A = higher_vertex
    B = lower_vertex
    E = [higher_vertex[0], higher_vertex[1], height]
    F = [lower_vertex[0], lower_vertex[1], height]
    if x_parallel:
        if side == "lower":
            C = [lower_vertex[0], lower_vertex[1] + new_width, 0]
            D = [higher_vertex[0], higher_vertex[1] + new_width, 0]
            G = [lower_vertex[0], lower_vertex[1] + new_width, height]
            H = [higher_vertex[0], higher_vertex[1] + new_width, height]
        elif side == "upper":
            C = [lower_vertex[0], lower_vertex[1] - new_width, 0]
            D = [higher_vertex[0], higher_vertex[1] - new_width, 0]
            G = [lower_vertex[0], lower_vertex[1] - new_width, height]
            H = [higher_vertex[0], higher_vertex[1] - new_width, height]
    else:
        if side == "left":
            C = [lower_vertex[0] + new_width, lower_vertex[1], 0]
            D = [higher_vertex[0] + new_width, higher_vertex[1], 0]
            G = [lower_vertex[0] + new_width, lower_vertex[1], height]
            H = [higher_vertex[0] + new_width, higher_vertex[1], height]
        elif side == "right":
            C = [lower_vertex[0] - new_width, lower_vertex[1], 0]
            D = [higher_vertex[0] - new_width, higher_vertex[1], 0]
            G = [lower_vertex[0] - new_width, lower_vertex[1], height]
            H = [higher_vertex[0] - new_width, higher_vertex[1], height]
    return [A, B, C, D, E, F, G, H]

# calculates the roof height, then the roof coordinates
def compute_new_rooftop_coordinates(roof_type, roof_orientation, volume_roof, length, new_width, x_parallel, side, roof_up, points):
    if roof_orientation > 179:
        roof_orientation -= 180
    if roof_type == "Shed":
        z = (2*volume_roof)/(length*new_width) + points[4][2]
        x_min, x_max, y_min, y_max = get_global_extremes(roof_up)
        if x_parallel:
            if round(roof_orientation, 0) == 0:
                length_parallel = True
                if abs(y_min - points[1][1]) < abs(y_max - points[2][1]):
                    I = [points[0][0], points[0][1], z]
                    J = [points[1][0], points[1][1], z]
                else:
                    I = [points[3][0], points[3][1], z]
                    J = [points[2][0], points[2][1], z]
            elif round(roof_orientation, 0) == 90:
                length_parallel = False
                if abs(x_min - points[1][0]) < abs(x_max - points[0][0]):
                    I = [points[1][0], points[1][1], z]
                    J = [points[2][0], points[2][1], z]
                else:
                    I = [points[0][0], points[0][1], z]
                    J = [points[3][0], points[3][1], z]
        else:
            if round(roof_orientation, 0) == 90:
                length_parallel = True
                if abs(x_min - points[1][0]) < abs(x_max - points[2][0]):
                    I = [points[0][0], points[0][1], z]
                    J = [points[1][0], points[1][1], z]
                else:
                    I = [points[3][0], points[3][1], z]
                    J = [points[2][0], points[2][1], z]
            elif round(roof_orientation, 0) == 0:
                length_parallel = False
                if abs(y_min - points[1][1]) < abs(y_max - points[0][1]):
                    I = [points[1][0], points[1][1], z]
                    J = [points[2][0], points[2][1], z]
                else:
                    I = [points[0][0], points[0][1], z]
                    J = [points[3][0], points[3][1], z]

    elif roof_type == "Gabled":
        z = (2*volume_roof)/(length*new_width) + points[4][2]
        if x_parallel:
            if round(roof_orientation, 0) == 0:
                length_parallel = True
                xi = points[0][0]
                xj = points[1][0]
                if side == "lower":
                    yi = points[0][1] + 0.5*new_width
                    yj = yi
                elif side == "upper":
                    yi = points[0][1] - 0.5*new_width
                    yj = yi
            elif round(roof_orientation, 0) == 90:
                length_parallel = False
                yi = points[1][1]
                yj = points[2][1]
                xi = points[0][0] - 0.5*length
                xj = xi
    
        else:
            if round(roof_orientation, 0) == 90:
                length_parallel = True
                yi = points[0][1]
                yj = points[1][1]
                if side == "left":
                    xi = points[0][0] + 0.5*new_width
                    xj = xi
                elif side == "right":
                    xi = points[0][0] - 0.5*new_width
                    xj = xi
            elif round(roof_orientation, 0) == 0:
                length_parallel = False
                xi = points[1][0]
                xj = points[2][0]
                yi = points[0][1] - 0.5*length
                xj = xi
        I = [xi, yi, z]
        J = [xj, yj, z]
    
    elif roof_type == "Pyramidal":
        z = (3*volume_roof)/(length*new_width) + points[4][2]
        if x_parallel:
            x = points[0][0] - 0.5*length
            if side == "lower":
                y = points[0][1] + new_width*0.5
            elif side == "upper":
                y = points[0][1] - new_width*0.5      
        else:
            y = points[0][1] - 0.5*length
            if side == "left":
                x = points[0][0] + new_width*0.5
            elif side == "right":
                x = points[0][0] - new_width*0.5
        I = [x, y, z]

    elif roof_type == "Hipped":
        if round(roof_orientation, 0) == 90:
            if side == "left":
                x = points[0][0] + 0.5*new_width
                length_parallel = True
            elif side == "right":
                x = points[0][0] - 0.5*new_width
                length_parallel = True
            else:
                x = points[0][0] - 0.5* length
                length_parallel = False
            y_upper = roof_up[0][1]
            y_lower = y_upper
            for v in roof_up:
                if v[1] > y_upper:
                    y_upper = v[1]
                elif v[1] < y_lower:
                    y_lower = v[1]
            if side == "lower":
                I = [x, y_lower]
                J = [x, y_upper]
            else: 
                I = [x, y_upper]
                J = [x, y_lower]
            rooftop_len = y_upper - y_lower
        
        elif round(roof_orientation, 0) == 0:
            if side == "lower":
                y = points[0][1] + 0.5*new_width
                length_parallel = True
            elif side == "upper":
                y = points[0][1] - 0.5*new_width
                length_parallel = True
            else:
                y = points[0][1] - 0.5* length
                length_parallel = False
            x_right = roof_up[0][0]
            x_left = x_right
            for v in roof_up:
                if v[0] > x_right:
                    x_right = v[0]
                elif v[0] < x_left:
                    x_left = v[0]
            if side == "left":
                I = [x_left, y]
                J = [x_right, y]
            else:
                I = [x_right, y]
                J = [x_left, y]
            rooftop_len = x_right - x_left
        if length > new_width:
            num = 2*length + rooftop_len
            c = (6*volume_roof)/(new_width*num) + points[4][2]
        else:
            num = 2*new_width + rooftop_len
            c = (6*volume_roof)/(length*num) + points[4][2]
        I.append(c)
        J.append(c)
        
    if roof_type == "Pyramidal":
        return I
    else:
        return I, J, length_parallel

# modifies values of indices by the certain value        
def modify_indices(wall_faces, increase_index_by):
    index1 = 0
    for i in wall_faces:
        index2 = 0
        for j in i:
            new_value = j + increase_index_by
            wall_faces[index1][index2] = new_value
            index2 += 1
        index1 += 1
    return wall_faces

#adds newly calculated coordinates into final structures
def add_new_geometry(points, wall_faces, roof_faces, new_vertices, new_objects, new_obj_id):
    for p in points:
        new_vertices.append(p)
    
    all_faces = []
    for wf in wall_faces:
        all_faces.append(wf)
    for rf in roof_faces:
        all_faces.append(rf)

    new_objects[new_obj_id] = all_faces
    return new_vertices, new_objects 
 
# detects position of the adjacent building 
def get_neighbor_direction(x_sum, y_sum, a_x_sum, a_y_sum, extremes, a_extremes):
    if abs(a_x_sum - x_sum) > abs(a_y_sum - y_sum):
        if a_x_sum < x_sum:
            if extremes[0] != a_extremes[0] and extremes[1] != a_extremes[1]:
                direction = "left"
            else:
                if extremes[2] < a_extremes[2]:
                    direction = "upper"
                else:
                    direction = "lower"
        else:
            if extremes[0] != a_extremes[0] and extremes[1] != a_extremes[1]:
                direction = "right"
            else:
                if extremes[2] < a_extremes[2]:
                    direction = "upper"
                else:
                    direction = "lower"
    else:
        if a_y_sum < y_sum:
            if extremes[2] != a_extremes[2] and extremes[3] != a_extremes[3]:
                direction = "lower"
            else:
                if extremes[0] < a_extremes[0]:
                    direction = "right"
                else:
                    direction = "left" 
        else:
            if extremes[2] != a_extremes[2] and extremes[3] != a_extremes[3]:
                direction = "upper"
            else:
                if extremes[0] < a_extremes[0]:
                    direction = "right"
                else:
                    direction = "left"
    return direction

# asign adjacent buildings and their position to each building
def create_objects_neighborhood(object_parts, only_single_objects = False, parts_list = None):
    neighborhood = {}
    for obj_id, parts in object_parts.items():
        if only_single_objects:
            if obj_id in parts_list:
                continue
        footprintcoords = parts["footprintcoords"]
        extremes = get_global_extremes(footprintcoords)
        x_sum = 0
        y_sum = 0
        for coord in footprintcoords:
            x = coord[0]
            x_sum += x
            y = coord[1]
            y_sum += y

        a = 1000
        a_id = str()
        a_x_sum = float()
        a_y_sum = float()
        a_extremes = list()
        b = 1000
        b_id = str()
        b_x_sum = float()
        b_y_sum = float()
        b_extremes = list()
        for obj_id2, parts2 in object_parts.items():
            if obj_id != obj_id2:
                footprintcoords2 = parts2["footprintcoords"]
                x_min2, x_max2, y_min2, y_max2 = get_global_extremes(footprintcoords2)
                distance = 1000
                for coord in footprintcoords:
                    x = coord[0]
                    y = coord[1]
                    x2_sum = 0
                    y2_sum = 0
                    for coord2 in footprintcoords2:
                        x2 = coord2[0]
                        x2_sum += x2
                        y2 = coord2[1]
                        y2_sum += y2
                        x_diff = abs(x - x2)
                        y_diff = abs(y - y2)
                        new_distance = math.sqrt(x_diff**2 + y_diff**2)
                        if new_distance < distance:
                            distance = new_distance
                if a > b:
                    if distance < a:
                        a = distance
                        a_id = obj_id2
                        a_x_sum = x2_sum
                        a_y_sum = y2_sum
                        a_extremes = x_min2, x_max2, y_min2, y_max2
                else:
                    if distance < b:
                        b = distance
                        b_id = obj_id2
                        b_x_sum = x2_sum
                        b_y_sum = y2_sum
                        b_extremes = x_min2, x_max2, y_min2, y_max2
                           
        neighbors = {}
        if a != 1000:
            direction1 = get_neighbor_direction(x_sum, y_sum, a_x_sum, a_y_sum, extremes, a_extremes)
            neighbors[a_id] = direction1
        if b != 1000:
            direction2 = get_neighbor_direction(x_sum, y_sum, b_x_sum, b_y_sum, extremes, b_extremes)
            if a != 1000 and direction2 == direction1 and b < a:
                neighbors = {}
                neighbors[b_id] = direction2
            elif a != 1000 and direction2 == direction1 and b > a:
                pass
            else:
                neighbors[b_id] = direction2
        neighborhood[obj_id] = neighbors
    return neighborhood
            
# fills gap next to the single building by changing its geometry
def fill_gaps(obj_id, obj_id1, parts_list, new_objects, new_vertices, objects, vertices, object_parts, direction, EPSILON):
    object = objects[obj_id1]
    footprintcoords = object_parts[obj_id1]["footprintcoords"]
    xmin, xmax, ymin, ymax = get_global_extremes(footprintcoords)
    x1 = footprintcoords[0][0]
    y1 = footprintcoords[0][1]
    if obj_id in parts_list:
        for new_obj_id, new_faces in new_objects.items():
            if re.search(obj_id, new_obj_id) or obj_id == new_obj_id:
                new_footprintcoords = []
                for new_face in new_faces:
                    for new_vertex in new_face:
                        if new_vertices[new_vertex-1][2] == 0:
                            if new_vertices[new_vertex-1] not in new_footprintcoords:
                                new_footprintcoords.append(new_vertices[new_vertex-1])
                    if len(new_footprintcoords) == 4:
                        break
                
                new_x = new_footprintcoords[0][0]
                new_y = new_footprintcoords[0][1]
                for new_coord in new_footprintcoords:
                    if direction == "left":
                        if new_coord[0] > new_x:
                            new_x = new_coord[0]
                    elif direction == "right":
                        if new_coord[0] < new_x:
                            new_x = new_coord[0]
                    elif direction == "lower":
                        if new_coord[1] > new_y:
                            new_y = new_coord[1]
                    elif direction == "upper":
                        if new_coord[1] < new_y:
                            new_y = new_coord[1]

                for face1 in object:
                    for vertex1 in face1:
                        v = vertices[vertex1-1]
                        if direction == "left":
                            if abs(v[0] - xmin) < EPSILON:
                                vertices[vertex1-1][0] = new_x
                        elif direction == "right":
                            if abs(v[0] - xmax) < EPSILON:
                                vertices[vertex1-1][0] = new_x
                        elif direction == "lower":
                            if abs(v[1] - ymin) < EPSILON:
                                vertices[vertex1-1][1] = new_y
                        elif direction == "upper":
                            if abs(v[1] - ymax) < EPSILON:
                                vertices[vertex1-1][1] = new_y 

    else:         
        for coord1 in footprintcoords:
            if direction == "left":
                if coord1[0] < x1:
                    x1 = coord1[0]
            elif direction == "right":
                if coord1[0] > x1:
                    x1 = coord1[0]
            elif direction == "upper":
                if coord1[1] > y1:
                    y1 = coord1[1]
            elif direction == "lower":
                if coord1[1] < y1:
                    y1 = coord1[1]

        object2 = objects[obj_id]
        footprintcoords2 = object_parts[obj_id]["footprintcoords"]
        x2 = footprintcoords2[0][0]
        y2 = footprintcoords2[0][1]
        for coord2 in footprintcoords2:
            if direction == "left":
                if coord2[0] > x2:
                    x2 = coord2[0]
            elif direction == "right":
                if coord2[0] < x2:
                    x2 = coord2[0]
            elif direction == "upper":
                if coord2[1] < y2:
                    y2 = coord2[1]
            elif direction == "lower":
                if coord2[1] > y2:
                    y2 = coord2[1]
                             
        if direction == "left" or direction == "right":
            new_x = (x1 + x2)/2
            for face1 in object:
                for v1 in face1:
                    if abs(vertices[v1 - 1][0] - x1) < EPSILON:
                        vertices[v1 - 1][0] = new_x
            for face2 in object2:
                for v2 in face2:
                    if abs(vertices[v2 - 1][0] - x2) < EPSILON:
                        vertices[v2 - 1][0] = new_x
        else:
            new_y = (y1 + y2)/2
            for face1 in object:
                for v1 in face1:
                    if abs(vertices[v1 - 1][1] - y1) < EPSILON:
                        vertices[v1 - 1][1] = new_y
            for face2 in object2:
                for v2 in face2:
                    if abs(vertices[v2 - 1][1] - y2) < EPSILON:
                        vertices[v2 - 1][1] = new_y
    return vertices

# detects colision between aggregates
def detect_colision(object_parts, aggregation_parts):
    used = []
    collision = []
    identical_aggregate = []
    for obj, parts in object_parts.items():
        footprintcoords = parts["footprintcoords"]
        x_min, x_max, y_min, y_max = get_global_extremes(footprintcoords)
        for obj2, parts2 in object_parts.items():
            if obj != obj2 and [obj2, obj] not in used:
                used.append([obj, obj2])
                footprintcoords2 = parts2["footprintcoords"]
                x_min2, x_max2, y_min2, y_max2 = get_global_extremes(footprintcoords2)
                cond1 = x_min2 > x_min and x_min2 < x_max
                cond2 = x_max2 > x_min and x_max2 < x_max
                cond3 = y_min2 > y_min and y_min2 < y_max
                cond4 = y_max2 > y_min and y_max2 < y_max
                cond5 = cond1 or cond2
                cond6 = cond3 or cond4
                if cond5 and cond6:
                    collision.append([obj, obj2])
                else:
                    parts1 = obj.split(" + ")
                    parts2 = obj2.split(" + ")
                    for p in aggregation_parts:
                        if parts1[0] in p and parts2[0] in p:
                            identical_aggregate.append([obj, obj2])
    return collision, identical_aggregate

# removes overlapping of adjacent aggregates
def move_object(collision, x_parallel, side, vertices, object_parts):
    obj_id = collision[0]
    obj_id2 = collision[1]
    footprintcoords = object_parts[obj_id]["footprintcoords"]
    x_min, x_max, y_min, y_max = get_global_extremes(footprintcoords)
    footprintcoords2 = object_parts[obj_id2]["footprintcoords"]
    x_min2, x_max2, y_min2, y_max2 = get_global_extremes(footprintcoords2)
    wall_faces = object_parts[obj_id]["wall_faces"]
    roof_faces = object_parts[obj_id]["roof_faces"]
    moved_vertex = []
    if x_parallel:
        if side == "lower":
            distance = - y_max + y_min2
        elif side == "upper":
            distance = y_max2 - y_min
    else:
        if side == "left":
            distance = - x_max + x_min2
        elif side == "upper":
            distance = x_max2 - x_min
    for face in wall_faces:
        for vertex in face:
            if vertex not in moved_vertex:
                if x_parallel:
                    vertices[vertex - 1][1] += distance
                else:
                    vertices[vertex - 1][0] += distance
            moved_vertex.append(vertex)
    for face in roof_faces:
        for vertex in face:
            if vertex not in moved_vertex:
                if x_parallel:
                    vertices[vertex - 1][1] += distance
                else:
                    vertices[vertex - 1][0] += distance
            moved_vertex.append(vertex)
    return vertices

# removes overlapping within the aggregate    
def remove_overlap(identical_object, object_parts, vertices, x_parallel, side):
    obj_id = identical_object[0]
    obj_id2 = identical_object[1]
    footprintcoords = object_parts[obj_id]["footprintcoords"]
    x_min, x_max, y_min, y_max = get_global_extremes(footprintcoords)
    footprintcoords2 = object_parts[obj_id2]["footprintcoords"]
    x_min2, x_max2, y_min2, y_max2 = get_global_extremes(footprintcoords2)
    wall_faces = object_parts[obj_id]["wall_faces"]
    roof_faces = object_parts[obj_id]["roof_faces"]
    for face in wall_faces:
        for vertex in face:
            v = vertices[vertex - 1]
            x = v[0]
            y = v[1]
            if x_parallel:
                if side == "right":
                    if x == x_max:
                        vertices[vertex - 1][0] = x_min2
                elif side == "left":
                    if x == x_min:
                        vertices[vertex - 1][0] = x_max2
            else:
                if side == "upper":
                    if y == y_max:
                        vertices[vertex - 1][1] = y_min2
                elif side == "left":
                    if y == y_min:
                        vertices[vertex - 1][1] = y_max2
    for face in roof_faces:
        for vertex in face:
            v = vertices[vertex - 1]
            x = v[0]
            y = v[1]
            if x_parallel:
                if side == "right":
                    if x == x_max:
                        vertices[vertex - 1][0] = x_min2
                elif side == "left":
                    if x == x_min:
                        vertices[vertex - 1][0] = x_max2
            else:
                if side == "upper":
                    if y == y_max:
                        vertices[vertex - 1][1] = y_min2
                elif side == "left":
                    if y == y_min:
                        vertices[vertex - 1][1] = y_max2
    return vertices

# removes aggregates from the neighborhood list
def get_single_objects_neighborhood(neighborhood, parts_list):
    keys = neighborhood.keys()
    for obj_id in parts_list:
        if obj_id in keys:
            neighborhood.pop(obj_id)
    return neighborhood

# gets rid of gaps by changing length of the single buildings
def change_length_of_single_buildings(vertices, new_vertices, neighborhood, objects, new_objects, object_parts, parts_list):
    for obj_id1, dictionary in neighborhood.items():            
        for obj_id, direction in dictionary.items():
            vertices = fill_gaps(obj_id, obj_id1, parts_list, new_objects, new_vertices, objects, vertices, object_parts, direction, EPSILON)
    return vertices

# returns list of vertices, list of wall indices, list of roof indices of the new building
def create_geometry(roof_type, points, obj_id):
    if roof_type == "Pyramidal":
        I = compute_new_rooftop_coordinates(roof_type, roof_orientation[obj_id], volume_roof, length, new_width, x_parallel, side, roof_up, points)
        points.append(I)
        wall_faces = [[1, 2, 3, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 4, 8, 7], [4, 1, 5, 8]]
        roof_faces = [[5, 6, 9], [6, 7, 9], [7, 8, 9], [8, 5, 9]]
    elif roof_type == "Flat":
        wall_faces = [[1, 2, 3, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 4, 8, 7], [4, 1, 5, 8]]
        roof_faces = [[5, 6, 7, 8]]
    else:
        I, J, length_parallel = compute_new_rooftop_coordinates(roof_type, roof_orientation[obj_id], volume_roof, length, new_width, x_parallel, side, roof_up, points)
        points.append(I)
        points.append(J)
        if roof_type == "Hipped":
            wall_faces = [[1, 2, 3, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 4, 8, 7], [4, 1, 5, 8]]
            if length_parallel:
                roof_faces = [[5, 6, 10, 9], [6, 7, 10], [7, 8, 9, 10], [8, 5, 9]]
            else:
                roof_faces = [[5, 6, 9], [6, 7, 10, 9], [7, 8, 10], [8, 5, 9, 10]]
        elif roof_type == "Gabled":
            if length_parallel:
                wall_faces = [[1, 2, 3, 4], [1, 2, 6, 5], [3, 4, 8, 7]]
                roof_faces = [[2, 3, 7, 10, 6], [4, 1, 5, 9, 8], [5, 6, 10, 9], [7, 8, 9, 10]]
            else:
                wall_faces = [[1, 2, 3, 4], [2, 3, 7, 6], [1, 5, 8, 4]]
                roof_faces = [[1, 2, 6, 9, 5], [3, 4, 8, 10, 7], [6, 7, 10, 9], [8, 5, 9, 10]]
        elif roof_type == "Shed":
            wall_faces = [[1, 2, 3, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 4, 8, 7], [4, 1, 5, 8]]
            if length_parallel:
                roof_faces = [[5, 6, 10, 9], [6, 7, 10], [7, 8, 9, 10], [8, 5, 9]]
            else:
                roof_faces = [[5, 6, 9], [6, 7, 10, 9], [7, 8, 10], [8, 5, 9, 10], [5, 6, 9]]
    return points, wall_faces, roof_faces

# detects properties of the single object
def get_single_object_parametres(obj_id, roof_types, volume_dict, object_parts, x_min_global, x_max_global, y_min_global, y_max_global, neighborhood):
    corner = False
    roof_type = roof_types[obj_id]
    volume_body = volume_dict[obj_id]["volume_body"]
    volume_roof = volume_dict[obj_id]["volume_roof"]
    footprintcoords = object_parts[obj_id]["footprintcoords"]
    roofcoords = object_parts[obj_id]["roofcoords"]
    roof_faces = object_parts[obj_id]["roof_faces"]
    wall_faces = object_parts[obj_id]["wall_faces"]
    x_min, x_max, y_min, y_max = get_global_extremes(footprintcoords)
    condition1 = x_min == x_min_global or x_max == x_max_global
    condition2 = y_min == y_min_global or y_max == y_max_global
    if condition1 and condition2:
        corner = True
    body_height = roofcoords[0][2]
    roof_top = roofcoords[0][2]
    for coord in roofcoords:
        if coord[2] < body_height:
            body_height = coord[2]
        elif coord[2] > roof_top:
            roof_top = coord[2]
    neighbors = neighborhood[obj_id]
    for neig_id, direction in neighbors.items():
        if direction == "left" or direction == "right":
            length = x_max - x_min
            width = y_max - y_min
            x_parallel = True
        else:
            width = x_max - x_min
            length = y_max - y_min
            x_parallel = False
        break
    return corner, x_parallel, body_height, roof_top, width, length, x_min, x_max, y_min, y_max, roof_type, volume_body, volume_roof, roof_faces, wall_faces

# modifies width of the single object by its new length
def modify_width_of_single_object(volume_body, length, body_height, wall_faces, x_parallel, x_min, x_max, y_min, y_max, x_min_global, x_max_global, y_min_global, y_max_global, vertices, EPSILON):
    width = volume_body/(length*body_height)
    for face in wall_faces:
        for vertex in face:
            if x_parallel:
                new_x = None
                if abs(y_min - y_min_global) < EPSILON:
                    max_change = True
                    new_y = y_min + width
                    if abs(vertices[vertex-1][1] - y_max) < EPSILON:
                        vertices[vertex-1][1] = new_y
                elif abs(y_max - y_max_global) < EPSILON:
                    max_change = False
                    new_y = y_max - width
                    if abs(vertices[vertex-1][1] - y_min) < EPSILON:
                        vertices[vertex-1][1] = new_y
            else:
                new_y = None
                if abs(x_min - x_min_global) < EPSILON:
                    max_change = True
                    new_x = x_min + width
                    if abs(vertices[vertex-1][0] - x_max) < EPSILON:
                        vertices[vertex-1][0] = new_x
                elif abs(x_max - x_max_global) < EPSILON:
                    max_change = False
                    new_x = x_max - width
                    if abs(vertices[vertex-1][0] - x_min) < EPSILON:
                        vertices[vertex-1][0] = new_x
    return vertices, width, max_change, new_x, new_y

# updates values of coordinate extremes
def update_local_extremes(x_parallel, max_change, new_x, new_y, x_min, x_max, y_min, y_max):   
    if x_parallel:
        if max_change:
            y_max = new_y
        else:
            y_min = new_y
    else:
        if max_change:
            x_max = new_x
        else:
            x_min = new_x
    return x_min, x_max, y_min, y_max

# detects errors of outer alignment of the building block
def define_alignment(footprintcoords, direction, x_min, x_max, y_min, y_max, x_min_global, x_max_global, y_min_global, y_max_global, roof_type, vertices, wall_faces, roof_faces, EPSILON, changed = 0):
    xmin, xmax, ymin, ymax = get_global_extremes(footprintcoords)
    if direction == "right" or direction == "left":
        if abs(y_max - y_max_global) < abs(y_min - y_min_global):
            if abs(y_max - ymax) > EPSILON:
                change = "y_max"
                changed = 1
        else:
            if abs(y_min - ymin) > EPSILON:
                change = "y_min"
                changed = 1
    else:
        if abs(x_max - x_max_global) < abs(x_min - x_min_global):
            if abs(x_max - xmax) > EPSILON:
                change = "x_max"
                changed = 1
        else:
            if abs(x_min - xmin) > EPSILON:
                change = "x_min"
                changed = 1
    if changed == 1:
        vertices = align_object_boundary(wall_faces, vertices, change, x_min, x_max, y_min, y_max, xmin, xmax, ymin, ymax)
        if roof_type == "Gabled" or roof_type == "Shed":
            vertices = align_object_boundary(roof_faces, vertices, change, x_min, x_max, y_min, y_max, xmin, xmax, ymin, ymax)
        if direction == "left":
            direction = "right"
        elif direction == "right":
            direction = "left"
        elif direction == "lower":
            direction = "upper"
        elif direction == "upper":
            direction = "lower"

    return vertices, changed, direction

# reestablished alignment with the outer edge of the building block
def ensure_about_alignment(vertices, neighborhood, obj_id, parts_list, object_parts, new_object_parts, roof_type, x_min, x_max, y_min, y_max, x_min_global, x_max_global, y_min_global, y_max_global, wall_faces, roof_faces, EPSILON, aligned_objects, changed = 0, direction = None):
    neighbors = neighborhood[obj_id]
    values = list(neighbors.values())
    if direction in values:  
        for key, value in neighbors:
            if value == direction:
                if key in parts_list:
                    changed = 0
                else:
                    changed = 1
                    obj_id = key
    else:
        if changed != 0:
            changed = 2
        
    keys = list(neighbors.keys())
    if changed == 0:
        for key in keys:
            if key in parts_list:
                direction = neighbors[key]
                for new_obj_id, new_parts in new_object_parts.items():
                    if re.search(key, new_obj_id):
                        footprintcoords = new_parts["footprintcoords"]
                        vertices, change, direction = define_alignment(footprintcoords, direction, x_min, x_max, y_min, y_max, x_min_global, x_max_global, y_min_global, y_max_global, roof_type, vertices, wall_faces, roof_faces, EPSILON)
        if direction == None:
            changed = 2
            pass
        else:
            if change:
                aligned_objects.append(obj_id)
                if keys[0] in parts_list and keys[1] in parts_list:
                    return vertices
                else:
                    for key in keys:
                        if key not in parts_list:
                            ensure_about_alignment(vertices, neighborhood, key, parts_list, object_parts, new_object_parts, roof_type, x_min, x_max, y_min, y_max, x_min_global, x_max_global, y_min_global, y_max_global, wall_faces, roof_faces, EPSILON, aligned_objects, change, direction) 
    elif changed == 1:
        footprintcoords = object_parts[obj_id]["footprintcoords"]
        vertices, changed, direction = define_alignment(footprintcoords, direction, x_min, x_max, y_min, y_max, x_min_global, x_max_global, y_min_global, y_max_global, roof_type, vertices, wall_faces, roof_faces, EPSILON)
        if changed:
            aligned_objects.append(obj_id)
            for key in keys:
                if key != obj_id:
                    ensure_about_alignment(vertices, neighborhood, key, parts_list, object_parts, new_object_parts, roof_type, x_min, x_max, y_min, y_max, x_min_global, x_max_global, y_min_global, y_max_global, wall_faces, roof_faces, EPSILON, aligned_objects, changed, direction)
    else:
        pass
    
    return vertices, aligned_objects

# updates properties of the building
def update_object_parameters(obj_id, object_parts, corner, x_parallel):
    footprintcoords = object_parts[obj_id]["footprintcoords"]
    x_min, x_max, y_min, y_max = get_global_extremes(footprintcoords)
    if corner:
        width = x_max - x_min
        length = y_max - y_min
    else:
        if x_parallel:
            width = y_max - y_min
            length = x_max - x_min
        else:
            width = x_max - x_min
            length = y_max - y_min
    return width, length, x_min, x_max, y_min, y_max

# modifies height of the single building body
def modify_body_height_of_single_object(volume_body, length, width, roof_faces, corner, vertices, body_height, roof_top):
    new_z = volume_body/(length*width)
    roof_top_vertices = []
    for face in roof_faces:
        for vertex in face:
            if corner:
                if vertices[vertex-1][2] == body_height:
                    vertices[vertex-1][2] = new_z
            if vertices[vertex-1][2] == roof_top:
                if vertex not in roof_top_vertices:
                    roof_top_vertices.append(vertex)
    return vertices, roof_top_vertices, new_z

# modifies geometry of rooftop vertices of the single object
def modify_rooftop_of_single_object(obj_id, roof_type, x_min, x_max, y_min, y_max, length, width, body_height, volume_roof, roof_top_vertices, vertices, x_parallel, parts_list, EPSILON):
    if roof_type == "Pyramidal":
        new_x = (x_max + x_min)/2
        new_y = (y_max + y_min)/2
        new_z = volume_roof*3/(length*width) + body_height
        vertex = roof_top_vertices[0]
        vertices[vertex-1] = [new_x, new_y, new_z]
    
    else:
        a = vertices[roof_top_vertices[0]-1]
        b = vertices[roof_top_vertices[1]-1]
        if roof_type == "Gabled":
            if abs(a[0] - b[0]) < EPSILON:
                new_x = (x_max + x_min)/2
                vertices[roof_top_vertices[0]-1][0] = new_x
                vertices[roof_top_vertices[1]-1][0] = new_x
            else:
                new_y = (y_max + y_min)/2
                vertices[roof_top_vertices[0]-1][1] = new_y
                vertices[roof_top_vertices[1]-1][1] = new_y
            new_z = 2*volume_roof/(length*width) + body_height
            vertices[roof_top_vertices[0]-1][2] = new_z
            vertices[roof_top_vertices[1]-1][2] = new_z

        elif roof_type == "Shed":
            new_z = 2*volume_roof/(length*width) + body_height
            vertices[roof_top_vertices[0]-1][2] = new_z
            vertices[roof_top_vertices[1]-1][2] = new_z

        elif roof_type == "Hipped":
            if obj_id not in parts_list:
                if x_parallel:
                    if abs(a[1] - b[1]) < EPSILON:
                        rooftop_len = length/3
                        if length > width:
                            num = 2*length + rooftop_len
                            new_z = (6*volume_roof)/(new_width*num) + body_height
                        else:
                            num = 2*width + rooftop_len
                            new_z = (6*volume_roof)/(length*num) + body_height
                        new_y = (y_max + y_min)/2
                        new_x1 = x_min + rooftop_len
                        new_x2 = x_max - rooftop_len
                        if a[0] < b[0]:
                            vertices[roof_top_vertices[0]-1] = [new_x1, new_y, new_z]
                            vertices[roof_top_vertices[1]-1] = [new_x2, new_y, new_z]
                        else:
                            vertices[roof_top_vertices[0]-1] = [new_x2, new_y, new_z]
                            vertices[roof_top_vertices[1]-1] = [new_x1, new_y, new_z]
                    else:
                        rooftop_len = width/3
                        if length > width:
                            num = 2*length + rooftop_len
                            new_z = (6*volume_roof)/(new_width*num) + body_height
                        else:
                            num = 2*width + rooftop_len
                            new_z = (6*volume_roof)/(length*num) + body_height
                        new_x = (x_max + x_min)/2
                        new_y1 = y_min + rooftop_len
                        new_y2 = y_max - rooftop_len
                        if a[1] < b[1]:
                            vertices[roof_top_vertices[0]-1] = [new_x, new_y1, new_z]
                            vertices[roof_top_vertices[1]-1] = [new_x, new_y2, new_z]
                        else:
                            vertices[roof_top_vertices[0]-1] = [new_x, new_y2, new_z]
                            vertices[roof_top_vertices[1]-1] = [new_x, new_y1, new_z]
                        
                else:
                    if abs(a[0] - b[0]) < EPSILON:
                        rooftop_len = length/3
                        if length > width:
                            num = 2*length + rooftop_len
                            new_z = (6*volume_roof)/(new_width*num) + body_height
                        else:
                            num = 2*width + rooftop_len
                            new_z = (6*volume_roof)/(length*num) + body_height
                        new_x = (x_max + x_min)/2
                        new_y1 = y_min + rooftop_len
                        new_y2 = y_max - rooftop_len
                        if a[1] < b[1]:
                            vertices[roof_top_vertices[0]-1] = [new_x, new_y1, new_z]
                            vertices[roof_top_vertices[1]-1] = [new_x, new_y2, new_z]
                        else:
                            vertices[roof_top_vertices[0]-1] = [new_x, new_y2, new_z]
                            vertices[roof_top_vertices[1]-1] = [new_x, new_y1, new_z]
                    else:
                        rooftop_len = width/3
                        if length > width:
                            num = 2*length + rooftop_len
                            new_z = (6*volume_roof)/(new_width*num) + body_height
                        else:
                            num = 2*width + rooftop_len
                            new_z = (6*volume_roof)/(length*num) + body_height
                        new_y = (y_max + y_min)/2
                        new_x1 = x_min + rooftop_len
                        new_x2 = x_max - rooftop_len
                        if a[0] < b[0]:
                            vertices[roof_top_vertices[0]-1] = [new_x1, new_y, new_z]
                            vertices[roof_top_vertices[1]-1] = [new_x2, new_y, new_z]
                        else:
                            vertices[roof_top_vertices[0]-1] = [new_x1, new_y, new_z]
                            vertices[roof_top_vertices[1]-1] = [new_x2, new_y, new_z]
    return vertices

# adds single objects to aggregate structures 
def add_single_objects_to_aggregates_structures(objects, parts_list, vertices, new_vertices_length, new_objects, new_vertices):
    for obj_id, faces in objects.items():
        if obj_id not in parts_list:
            new_faces = []
            for face in faces:
                new_face = []
                for vertex in face:
                    if vertices[vertex-1] not in new_vertices:
                        new_vertices_length += 1
                        new_vertices.append(vertices[vertex-1])
                        new_face.append(new_vertices_length)
                    else:
                        counter = 1
                        for new_vertex in new_vertices:
                            if vertices[vertex-1] == new_vertex:
                                new_face.append(counter)
                            counter += 1
                new_faces.append(new_face)
            new_objects[obj_id] = new_faces
    return new_objects, new_vertices 

# performs alignment of single buildings into the block
def solve_single_buildings(objects, vertices, object_parts, neighborhood, ensure = True):
    aligned_objects_list = []
    if type(objects) == dict:
        objects = list(objects.keys())
    for obj_id in objects:
        condition1 = " + " in obj_id
        condition2 = obj_id not in parts_list
        if condition1 or condition2:
            corner, x_parallel, body_height, roof_top, width, length, x_min, x_max, y_min, y_max, roof_type, volume_body, volume_roof, roof_faces, wall_faces = get_single_object_parametres(obj_id, roof_types, volume_dict, object_parts, x_min_global, x_max_global, y_min_global, y_max_global, neighborhood)
            if not corner:
                vertices, width, max_change, new_x, new_y = modify_width_of_single_object(volume_body, length, body_height, wall_faces, x_parallel, x_min, x_max, y_min, y_max, x_min_global, x_max_global, y_min_global, y_max_global, vertices, EPSILON)
                x_min, x_max, y_min, y_max = update_local_extremes(x_parallel, max_change, new_x, new_y, x_min, x_max, y_min, y_max)
            
            if ensure:
                aligned_objects = []
                if obj_id not in aligned_objects_list:
                    vertices, aligned_objects = ensure_about_alignment(vertices, neighborhood, obj_id, parts_list, object_parts, new_object_parts, roof_type, x_min, x_max, y_min, y_max, x_min_global, x_max_global, y_min_global, y_max_global, wall_faces, roof_faces, EPSILON, aligned_objects)
                if len(aligned_objects) > 0:
                    for o in aligned_objects:
                        aligned_objects_list.append(o)

                if obj_id in aligned_objects_list:
                    width, length, x_min, x_max, y_min, y_max = update_object_parameters(obj_id, object_parts, corner, x_parallel)
                    aligned_objects_list.remove(obj_id)

            vertices, roof_top_vertices, body_height = modify_body_height_of_single_object(volume_body, length, width, roof_faces, corner, vertices, body_height, roof_top)
            if roof_type != "Flat":
                vertices = modify_rooftop_of_single_object(obj_id, roof_type, x_min, x_max, y_min, y_max, length, width, body_height, volume_roof, roof_top_vertices, vertices, x_parallel, parts_list, EPSILON)
    if ensure:
        return vertices, aligned_objects_list
    else:
        return vertices

# overwrites vertices of aggregates
def recreate_vertices(temp_new_objects, temp_new_vertices, new_vertices, aggregates_count):
    for faces in temp_new_objects.values():
        for face in faces:
            for vertex in face:
                new_vertices[vertex - 1] = temp_new_vertices[vertex - 1]
        aggregates_count -= 1
        if aggregates_count == 0:
            break
    return new_vertices

# creates OBJ files from altered geometry
def save_obj(path, header, vertices, objects):
    with open(path, "w", encoding="utf-8") as f:
        for l in header:
            f.write(l)
        f.write("\n")
        for v in vertices:
            f.write("v {}\n".format(" ".join(map(str, v))))
        f.write("\n")
        for obj_id, faces in objects.items():
            f.write("o ['{}']\n".format(obj_id))
            for face in faces:
                f.write("f {} \n".format(" ".join(map(str, face))))

# --------------------------------------------------------------#
# calling functions

# loading inputs
optimization_input = json.load(open(sys.argv[2]))
header, vertices, objects = load_obj(sys.argv[3])

roof_type_dict = optimization_input[4]
roof_orientation = optimization_input[7]
roof_types = load_roof_type_names(roof_type_dict)

# preparation of the data for aggregation
object_parts, volume_dict = create_auxiliary_structures(roof_types, objects, vertices)

angle = get_object_orientation(object_parts)

if angle > EPSILON:
    rotation_point = get_point_of_rotation(objects, vertices)
    vertices, roof_orientation = move_and_rotate_objects(vertices, rotation_point, angle, roof_orientation)
    object_parts = create_auxiliary_structures(roof_types, objects, vertices, False)

x_min_global, x_max_global, y_min_global, y_max_global = get_global_extremes(vertices)
#7. vstup
#x_min_global, x_max_global, y_min_global, y_max_global = -0.67, 42.32, 58.961, 89.23
#4. vstup
#x_min_global, x_max_global, y_min_global, y_max_global = 17.578, 73.308, 161.81, 192.98

neighborhood = create_objects_neighborhood(object_parts)

vertices = align_outer_boundary(x_min_global, x_max_global, y_min_global, y_max_global, object_parts, vertices, neighborhood, roof_types)

object_parts = create_auxiliary_structures(roof_types, objects, vertices, False)

parts_list, aggregation_parts = create_aggregate_structures(centers)

# aggregation
aggregation_parts_dict = create_auxiliary_structures(roof_types, objects, vertices, False, True, parts_list)
increase_index_by = 0
new_vertices = []
new_objects = {}
for i in range(0, len(aggregation_parts)):
    parallel_aggregates = divide_aggregates(aggregation_parts, i, neighborhood)
    for parallel_aggregate in parallel_aggregates:
        if len(parallel_aggregate) > 0:
            x_min, x_max, y_min, y_max, x_parallel, side = get_aggregate_extremes_and_position(parallel_aggregate, object_parts, x_min_global, x_max_global, y_min_global, y_max_global)
        
            volume_body, volume_roof, new_obj_id, neighbor_relation, roof_type, roof_up = get_aggregate_parametres(parallel_aggregate, volume_dict, aggregation_parts_dict, neighborhood, roof_types)
            roof_types[new_obj_id] = roof_type
            parts_list.append(new_obj_id)       
            
            length, higher_vertex, lower_vertex = calculate_aggregate_length(x_parallel, side, x_min, x_max, y_min, y_max)
            
            wallcoords = aggregation_parts_dict[parallel_aggregate[0]]["wallcoords"]
            new_width, height = calculate_new_width(wallcoords, length, volume_body, neighbor_relation, object_parts, lower_vertex, x_parallel, side)

            points = compute_new_body_coordinates(x_parallel, lower_vertex, higher_vertex, new_width, height, side)

            points, wall_faces, roof_faces = create_geometry(roof_type, points, parallel_aggregate[0])

            if increase_index_by != 0:
                wall_faces = modify_indices(wall_faces, increase_index_by)
                roof_faces = modify_indices(roof_faces, increase_index_by)
            
            if roof_type == "Pyramidal":
                increase_index_by += 9
            elif roof_type == "Flat":
                increase_index_by += 8
            else:
                increase_index_by += 10
            
            new_vertices, new_objects = add_new_geometry(points, wall_faces, roof_faces, new_vertices, new_objects, new_obj_id)

new_object_parts = create_auxiliary_structures(roof_types, new_objects, new_vertices, False)

collision, identical_aggregate = detect_colision(new_object_parts, aggregation_parts)

if len(collision) > 0:
    for c in collision:
        parallel_aggregate = c[0].split(" + ")
        x_min, x_max, y_min, y_max, x_parallel, side = get_aggregate_extremes_and_position(parallel_aggregate, object_parts, x_min_global, x_max_global, y_min_global, y_max_global)
        new_vertices = move_object(c, x_parallel, side, new_vertices, new_object_parts)
if len(identical_aggregate) > 0:
    for i in identical_aggregate:
        parallel_aggregate = i[0].split(" + ")
        x_min, x_max, y_min, y_max, x_parallel, side = get_aggregate_extremes_and_position(parallel_aggregate, object_parts, x_min_global, x_max_global, y_min_global, y_max_global)
        parallel_aggregate = i[1].split(" + ")
        x_min, x_max, y_min, y_max, x_parallel2, side2 = get_aggregate_extremes_and_position(parallel_aggregate, object_parts, x_min_global, x_max_global, y_min_global, y_max_global)
        new_vertices = remove_overlap(i, new_object_parts, new_vertices, x_parallel, side2)

# alignment of single buildings into the block
neighborhood = get_single_objects_neighborhood(neighborhood, parts_list)

vertices = change_length_of_single_buildings(vertices, new_vertices, neighborhood, objects, new_objects, object_parts, parts_list)

vertices, aligned_objects_list = solve_single_buildings(objects, vertices, object_parts, neighborhood)

if len(aligned_objects_list) > 0: 
    vertices, aligned_objects_list = solve_single_buildings(aligned_objects_list, vertices, object_parts, neighborhood)


# control between adjacent aggregates
new_vertices_length = len(new_vertices)
temp_new_objects, temp_vertices = add_single_objects_to_aggregates_structures(objects, parts_list, vertices, new_vertices_length, new_objects, new_vertices)

temp_new_object_parts, volume_dict = create_auxiliary_structures(roof_types, temp_new_objects, temp_vertices)
temp_new_neighborhood = create_objects_neighborhood(temp_new_object_parts)

temp_new_vertices = change_length_of_single_buildings(temp_vertices, temp_vertices, temp_new_neighborhood, temp_new_objects, temp_new_objects, temp_new_object_parts, parts_list)

temp_new_vertices = solve_single_buildings(temp_new_objects, temp_new_vertices, temp_new_object_parts, temp_new_neighborhood, ensure = False)
aggregates_count = len(new_object_parts)
new_vertices = recreate_vertices(temp_new_objects, temp_new_vertices, new_vertices, aggregates_count)       

new_vertices_length = len(new_vertices)

new_objects, new_vertices = add_single_objects_to_aggregates_structures(objects, parts_list, vertices, new_vertices_length, new_objects, new_vertices)   

save_obj(sys.argv[4], header, new_vertices, new_objects)