from setup_logger import logger
logger.info("Started logging")

import bpy
from bpy.types import GPencilFrame
import bmesh
from mathutils import Vector
import math
from math import sin, cos
import numpy as np
import perlin

# Base grease pencil art functions from https://towardsdatascience.com/blender-2-8-grease-pencil-scripting-and-generative-art-cbbfd3967590

def get_grease_pencil(gpencil_obj_name='GPencil') -> bpy.types.GreasePencil:
    """
    Return the grease-pencil object with the given name. Initialize one if not already present.
    :param gpencil_obj_name: name/key of the grease pencil object in the scene
    """

    # If not present already, create grease pencil object
    if gpencil_obj_name not in bpy.context.scene.objects:
        bpy.ops.object.gpencil_add(align='WORLD', location=(0, 0, 0), type='EMPTY')
        # rename grease pencil
        bpy.context.scene.objects[-1].name = gpencil_obj_name

    # Get grease pencil object
    gpencil = bpy.context.scene.objects[gpencil_obj_name]

    return gpencil


def get_grease_pencil_layer(gpencil: bpy.types.GreasePencil, gpencil_layer_name='GP_Layer',
                            clear_layer=False) -> bpy.types.GPencilLayer:
    """
    Return the grease-pencil layer with the given name. Create one if not already present.
    :param gpencil: grease-pencil object for the layer data
    :param gpencil_layer_name: name/key of the grease pencil layer
    :param clear_layer: whether to clear all previous layer data
    """

    # Get grease pencil layer or create one if none exists
    if gpencil.data.layers and gpencil_layer_name in gpencil.data.layers:
        gpencil_layer = gpencil.data.layers[gpencil_layer_name]
    else:
        gpencil_layer = gpencil.data.layers.new(gpencil_layer_name, set_active=True)

    if clear_layer:
        gpencil_layer.clear()  # clear all previous layer data

    # bpy.ops.gpencil.paintmode_toggle()  # need to trigger otherwise there is no frame

    return gpencil_layer


# Util for default behavior merging previous two methods
def init_grease_pencil(gpencil_obj_name='GPencil', gpencil_layer_name='GP_Layer',
                       clear_layer=True) -> bpy.types.GPencilLayer:
    gpencil = get_grease_pencil(gpencil_obj_name)
    gpencil_layer = get_grease_pencil_layer(gpencil, gpencil_layer_name, clear_layer=clear_layer)
    return gpencil_layer

def draw_line(gp_frame, p0: tuple, p1: tuple):
    # Init new stroke
    gp_stroke = gp_frame.strokes.new()
    gp_stroke.display_mode = '3DSPACE'  # allows for editing

    # Define stroke geometry
    gp_stroke.points.add(count=2)
    gp_stroke.points[0].co = p0
    gp_stroke.points[1].co = p1
    return gp_stroke

def draw_circle(gp_frame, center: tuple, radius: float, segments: int):
    # Init new stroke
    gp_stroke = gp_frame.strokes.new()
    gp_stroke.display_mode = '3DSPACE'  # allows for editing
    gp_stroke.draw_cyclic = True        # closes the stroke

    # Define stroke geometry
    angle = 2*math.pi/segments  # angle in radians
    gp_stroke.points.add(count=segments)
    for i in range(segments):
        x = center[0] + radius*math.cos(angle*i)
        y = center[1] + radius*math.sin(angle*i)
        z = center[2]
        gp_stroke.points[i].co = (x, y, z)

    return gp_stroke

def rotate_stroke(stroke, angle, axis='z'):
    # Define rotation matrix based on axis
    if axis.lower() == 'x':
        transform_matrix = np.array([[1, 0, 0],
                                     [0, cos(angle), -sin(angle)],
                                     [0, sin(angle), cos(angle)]])
    elif axis.lower() == 'y':
        transform_matrix = np.array([[cos(angle), 0, -sin(angle)],
                                     [0, 1, 0],
                                     [sin(angle), 0, cos(angle)]])
    # default on z
    else:
        transform_matrix = np.array([[cos(angle), -sin(angle), 0],
                                     [sin(angle), cos(angle), 0],
                                     [0, 0, 1]])

    # Apply rotation matrix to each point
    for i, p in enumerate(stroke.points):
        p.co = transform_matrix @ np.array(p.co).reshape(3, 1)

def draw_sphere(gp_frame, radius: int, circles: int):
    angle = math.pi / circles
    for i in range(circles):
        circle = draw_circle(gp_frame, (0, 0, 0), radius, 32)
        rotate_stroke(circle, angle*i, 'x')

# End base art functions from https://towardsdatascience.com/blender-2-8-grease-pencil-scripting-and-generative-art-cbbfd3967590

def get_cone_circles(center: tuple, start_radius: int, stop_radius: int, steps: int, step_size: int, scale_z: float) -> list:
    start_x = center[0]
    start_y = center[1]
    start_z = center[2]
    unscaled_z_values = list(range(start_z, start_z + steps * step_size, step_size))
    z_values = [z * scale_z for z in unscaled_z_values]
    radius_values = list(np.linspace(start_radius, stop_radius, len(z_values)))
    origins = []

    for index, value in enumerate(z_values):
        origins.append((start_x, start_y, value, radius_values[index]))

    return origins

def draw_cone(gp_frame, circles: tuple):
    origins, radius_values = circles
    for index, origin in enumerate(origins):
        draw_circle(gp_frame, origin, radius_values[index], 20)

def draw_tornado(gp_frame, circles: list, variance: float):
    baseNoise = perlin.BaseNoise()
    simplexNoise = perlin.SimplexNoise()
    perlin_noise = []
    for index, value in enumerate(circles):
        noise = simplexNoise.noise2(index, 0)
        perlin_noise.append(noise)

    print("perlin_noise:")
    print(perlin_noise)
    print("circles:")
    print(circles)

    for index, circle in enumerate(circles):
        x, y, z, radius = circle[0], circle[1], circle[2], circle[3]
        origin = (x, y, z)
        radius = radius + perlin_noise[index] * radius * variance
        draw_circle(gp_frame=gp_frame, center=origin, radius=radius, segments=20)


#def draw_curve_of_pursuit(top_left: tuple, starting_distance: int, steps: int, step_length: int):
#    top_right = tuple(np.add((starting_distance,0,0), top_left))
#    bottom_left = tuple(np.subtract(top_left, (0, starting_distance, 0)))
#    bottom_right = tuple(np.add(top_left, (starting_distance, -starting_distance, 0)))
#    print("Top left corner corner: " + top_left)
#    print("Top right corner: " + top_right)
#    print("Bottom left corner: " + bottom_left)
#    print("Bottom right corner: " + bottom_right)

#    while (steps > 0):
#        # Forecast along edge of current square
#        fc_top_left = tuple(np.subtract(top_left), (0, step_length, 0))
#        fc_top_right = tuple(np.subtract(top_right), (step_length, 0, 0))
#        fc_bottom_left = tuple(np.add(bottom_left), (step_length, 0, 0))
#        fc_bottom_right = tuple(np.add(bottom_right), (0, step_length, 0))

#        # For top right
#        # Calculate length of line between points
#        # Inflection points for curve are
#        # top left when y < bottom left
#        #

#        fc_top_left =
#        steps = steps - 1
#    # Forecast the points along the current square by step length
#    # Find angle between old theta and forecasted
#    # Draw a line cos(theta) = step_length / draw_length -> draw_length = step_length / cos(theta)
#    # Find theta as cos(theta) = (y_current - y_forecast) / (x_current - x_forecast)


def main():
    gp_layer = init_grease_pencil()
    gp_frame = gp_layer.frames.new(0)
    # draw_line(gp_frame, (0,0,0), (0,0,100))
    # draw_circle(gp_frame, (0,0,100), 10, 10)
    # draw_circle(gp_frame, (0,0,100), 20, 3)
    # draw_curve_of_pursuit(top_left=(0,0,0), starting_distance=100, steps=10, step_length=5)
    circles = get_cone_circles(center=(0,0,0), start_radius=1, stop_radius=80, steps=100, step_size=25, scale_z=0.05)
    # draw_cone(gp_frame, circles)
    draw_tornado(gp_frame=gp_frame, circles=circles, variance=0.3)

if __name__ == '__main__':
    main()
