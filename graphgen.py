#!/usr/bin/env python3

import numpy as np
import pygame
import sys
import itertools
import pickle

"""
creates the positioned list of vertices for an a x b x c hexagon graph
returns vertex_list, center_list
"""
def hex_vertices(a, b, c, shift_factor = (0,0)):

    factor = np.sqrt(3)/2
    # makes the height and width of the whole grid too big to start
    height = 2*c + a + b
    width = int(1.5*(a + b)+1)

    #creates the whole grid of hex points (square grid for now)
    vertex_list = []
    center_list = []
    for i in range(-width, width):
        for j in range(-height, height):
            # shifts every other row left by 1.5
            if j % 2 == 1:
                if i % 3 != 2:
                    vertex_list.append((i - 3/2, j * factor))
                else:
                    center_list.append((i - 3/2, j * factor))
            # non-shifted rows
            elif i % 3 != 2:
                vertex_list.append((i, j * factor))
            elif i % 3 == 2:
                center_list.append((i, j * factor))

    #trims the list of vertices based on the following six inequalities
    trimmed_list = []
    trimmed_center_list = []
    for vertex in vertex_list:
        if (vertex[1] <= factor * (1/1.5*vertex[0] + 2*c) and # this one controls the height on the left
            vertex[1] <= factor * (-1/1.5*vertex[0] + 2*c + 2*(a - 1) + 1) and
            vertex[0] >= -0.5 and 
            vertex[0] <= (a + b - 1)*1.5 and 
            vertex[1] >= factor * (-1/1.5*vertex[0]) and 
            vertex[1] >= factor * (1/1.5*vertex[0] + 1/3 - (2*b-1))
           ):
            trimmed_list.append(vertex)

    for vertex in center_list:
        if (vertex[1] <= factor * (1/1.5*vertex[0] + 2*c) and # this one controls the height on the left
            vertex[1] <= factor * (-1/1.5*vertex[0] + 2*c + 2*(a - 1) + 1) and
            vertex[0] >= -0.5 and 
            vertex[0] <= (a + b - 1)*1.5 and 
            vertex[1] >= factor * (-1/1.5*vertex[0]) and 
            vertex[1] >= factor * (1/1.5*vertex[0] + 1/3 - (2*b-1))
           ):
            trimmed_center_list.append(vertex)

    # creates a list of all included x- and y-values
    x_vals = []
    y_vals = []
    for vertex in trimmed_list:
        x_vals.append(vertex[0])
        y_vals.append(vertex[1])

    # shifts the set of vertices down and right so the whole graph is drawn on screen
    shifted_list = []
    shifted_center_list = []
    for vertex in trimmed_list:
        shifted_list.append((vertex[0] - min(x_vals) + 0.5 + shift_factor[0], vertex[1] - min(y_vals) + 0.5 + shift_factor[1]))
    for vertex in trimmed_center_list:
        shifted_center_list.append((vertex[0] - min(x_vals) + 0.5 + shift_factor[0], vertex[1] - min(y_vals) + 0.5 + shift_factor[1]))

    return shifted_list, shifted_center_list

"""
creates a square grid of vertices and centers
"""

def square_vertices(x, y, shift_factor = (0,0)):

    height = y
    width = x

    #creates the whole grid of hex points (square grid for now)
    vertex_list = []
    center_list = []

    vertex_list = list(itertools.product(range(width), range(height)))
    shifted_vertex_list = []
    for vertex in vertex_list:
        new_vertex = (vertex[0] + shift_factor[0], vertex[1] + shift_factor[1])
        shifted_vertex_list.append(new_vertex)
        if vertex[0] != width-1 and vertex[1] != height-1:
            
            center_list.append((new_vertex[0]+0.5, new_vertex[1] + 0.5))

    # creates a list of all included x- and y-values
    x_vals = []
    y_vals = []
    for vertex in shifted_vertex_list:
        x_vals.append(vertex[0])
        y_vals.append(vertex[1])


    return shifted_vertex_list, center_list

"""
uses the scaling factor to scale the drawn vertices
"""
def scale_vertices(vertex_list, scale_factor):
    return_vertices = []
    for i in vertex_list:
        return_vertices.append(tuple([(i[0])*scale_factor, (i[1])*scale_factor]))

    return return_vertices

""" 
takes in list of UNSCALED vertices and returns a list of edges in the form edge = [(x1,y1), (x2,y2)]
also takes in a scaling factor so that when the edges are created you can stretch them on the screen
"""
def make_box(edge):

    if max(edge[0][1], edge[1][1]) - min(edge[0][1], edge[1][1]) <= 1:
        y_val = 10
        y_min = min(edge[0][1], edge[1][1]) - 5
    else: 
        y_val = max(edge[0][1], edge[1][1]) - min(edge[0][1], edge[1][1])
        y_min = min(edge[0][1], edge[1][1])

    new_box = pygame.Rect(min(edge[0][0], edge[1][0]),
                            y_min,
                            max(edge[0][0], edge[1][0]) - min(edge[0][0], edge[1][0]),
                            y_val)
    return new_box

""" 
creates the edges given a list of vertices by checking if they are a distance of one apart, 
then scales based on the global scaling factor
"""

def make_edges(my_vertices, scale_factor):
    my_edges = []
    for i in range(len(my_vertices)):
        for j in range(i, len(my_vertices)):
            distance = (my_vertices[i][0] - my_vertices[j][0])**2 + (my_vertices[i][1] - my_vertices[j][1])**2
            if distance <= 1.01 and distance >= 0.99:
                my_edges.append([(my_vertices[i][0]*scale_factor, my_vertices[i][1]*scale_factor), (my_vertices[j][0]*scale_factor, my_vertices[j][1]*scale_factor)])


    # creates invisible rectangles around each edge
    my_boxes = []
    for edge in my_edges:
        my_boxes.append(make_box(edge))
        # checks to see if it's a horizontal edge, and if so it makes the click-box
        # 10 pixels taller for ease of clicking

    return my_edges, my_boxes

""" 
creates the edges given a list of vertices by checking if they are a distance of one apart, 
then scales based on the global scaling factor
"""

def make_edges(my_vertices, scale_factor):
    my_edges = []
    for i in range(len(my_vertices)):
        for j in range(i, len(my_vertices)):
            distance = (my_vertices[i][0] - my_vertices[j][0])**2 + (my_vertices[i][1] - my_vertices[j][1])**2
            if distance <= 1.01 and distance >= 0.99:
                my_edges.append([(my_vertices[i][0]*scale_factor, my_vertices[i][1]*scale_factor), (my_vertices[j][0]*scale_factor, my_vertices[j][1]*scale_factor)])


    # creates invisible rectangles around each edge
    my_boxes = []
    for edge in my_edges:
        my_boxes.append(make_box(edge))
        # checks to see if it's a horizontal edge, and if so it makes the click-box
        # 10 pixels taller for ease of clicking

    return my_edges, my_boxes


"""
calls the draw circle function in pygame to draw a dots at the location of each vertex in the requested color
"""
def render_vertex(coord_tuple, color, thickness = 4):
    point = (coord_tuple[0], coord_tuple[1])
    pygame.draw.circle(screen, color, point, thickness)

"""
calls the draw line function in pygame to draw a line from coord1 to coord2 in the requested color
"""
def render_edge(coord1, coord2, color, thickness = 4):
    pygame.draw.line(screen, color, coord1, coord2, thickness)

"""
draws a button in the desired position with given label
"""
def button(position, label, event, mouse_location):
    font = pygame.font.Font('freesansbold.ttf', 28)
    text = font.render(label, True, black, white)
    buttonbox = text.get_rect()
    buttonbox.midleft = position

    screen.blit(text, buttonbox)
    buttonbox.inflate_ip(20,20)
    pygame.draw.rect(screen, black, buttonbox, 1)
    buttonbox.inflate_ip(10,10)
    pygame.draw.rect(screen, (200,200,200), buttonbox, 1)

    if event.type == pygame.MOUSEBUTTONDOWN and buttonbox.collidepoint(mouse_location):
        pygame.draw.rect(screen, (255, 255, 255), buttonbox, 1)
        return True

"""
draws a button with + and - adjusters on the right and left and 
returns +1 or -1 if + or - is clicked
"""
def pmbutton(position, label, event, mouse_location):
    font = pygame.font.Font('freesansbold.ttf', 28)
    start_x = position[0]
    start_y = position[1]

    minus = button(position, "-", event, mouse_location)

    buttonbox = font.render("-", False, black, white).get_rect()
    new_x_position = buttonbox[2]
    button((start_x + new_x_position + 30, start_y), label, event, mouse_location)

    buttonbox = font.render(label, False, black, white).get_rect()
    new_x_position += buttonbox[2]
    plus = button((start_x + new_x_position + 60, start_y), "+", event, mouse_location)

    if minus:
        return -1

    if plus:
        return 1

"""
renders the square grid and centers
"""
def render_squares(x,y):
    unscaled_vertex_list, unscaled_center_list = square_vertices(x,y, shift_factor = (3, 4)) 
    edge_list, edge_rect_list = make_edges(unscaled_vertex_list, scale_factor)
    vertex_list = scale_vertices(unscaled_vertex_list, scale_factor)
    center_list = scale_vertices(unscaled_center_list, scale_factor)

    # goes through the edge list to draw the edges
    for j, i in enumerate(edge_list):
        tup_1 = i[0]
        tup_2 = i[1]
        render_edge(tup_1, tup_2, black, 6)
        #render_edge(tup_1, tup_2, rainbow(j, len(edge_list)))
        # for a rainbow over all the edges, replace the color with 
        "rainbow(j, len(edge_list))"

    # goes through the list of vertices to draw the dots on top of the edges
    for j,i in enumerate(vertex_list):
        render_vertex(i, black, 6)

    for i in center_list:
        render_vertex(i, teal, 8)

    return vertex_list, center_list, edge_list

"""
renders the hexagonal grid and centers
"""
def render_hexes(a,b,c):
    unscaled_vertex_list, unscaled_center_list = hex_vertices(a, b, c, shift_factor = (4,4)) # \ by / by |
    edge_list, edge_rect_list = make_edges(unscaled_vertex_list, scale_factor)
    vertex_list = scale_vertices(unscaled_vertex_list, scale_factor)
    center_list = scale_vertices(unscaled_center_list, scale_factor)
    
    # goes through the edge list to draw the edges
    for j, i in enumerate(edge_list):
        tup_1 = i[0]
        tup_2 = i[1]
        render_edge(tup_1, tup_2, black, 6)
        #render_edge(tup_1, tup_2, rainbow(j, len(edge_list)))
        # for a rainbow over all the edges, replace the color with 
        "rainbow(j, len(edge_list))"

    # goes through the list of vertices to draw the dots on top of the edges    
    for j,i in enumerate(vertex_list):
        render_vertex(i, black, 6)

    for i in center_list:
        render_vertex(i, teal, 8)

    return vertex_list, center_list, edge_list


##############################################    
######           PYGAME STUFF          #######
##############################################

fourk = 2

scale_factor = 60
a = 6
b = 8
c = 4

x_val = 17
y_val = 10

pygame.init()
size = width, height = 1020*fourk, 900*fourk

# generate some custom colors since full defaults are garish
black = 0, 0, 0 
white = 255, 255, 255
blue = 100, 100, 255
green = 0, 150, 0
red = 255, 50, 50
purple = 150, 50, 220
teal = 100, 200, 200
magenta = 255, 100, 255
yellow = 255, 255, 50


#renders the screen
screen = pygame.display.set_mode(size, pygame.RESIZABLE)

#sets the background to be white
screen.fill(white)



# needs to be a continual loop to keep drawing the screen
mouse_location = (0,0)
while True:

    for event in pygame.event.get():

        if event.type == pygame.MOUSEBUTTONDOWN:
            mouse_location = pygame.mouse.get_pos()

        if event.type == pygame.QUIT: 
            pygame.quit()
            sys.exit()

        square_grid = button((50,50), "square grid", event, mouse_location)
        hexagonal_grid = button((300, 50), "hexagonal grid", event, mouse_location)

        x_size = pmbutton((50, 125), "width", event, mouse_location)
        if x_size == 1 or x_size == -1:
            x_val += x_size
            if x_val < 1:
                x_val = 1

            screen.fill(white)
            vertex_list, center_list, edge_list = render_squares(x_val, y_val)

        y_size = pmbutton((300, 125), "height", event, mouse_location)
        if y_size == 1 or y_size == -1:
            y_val += y_size
            if y_val < 1:
                y_val = 1

            screen.fill(white)
            vertex_list, center_list, edge_list = render_squares(x_val, y_val)

        a_size = pmbutton((50, 200), "a", event, mouse_location)
        if a_size == 1 or a_size == -1:
            a += a_size
            if a < 1:
                a = 1

            screen.fill(white)
            vertex_list, center_list, edge_list = render_hexes(a,b,c)

        b_size = pmbutton((300, 200), "b", event, mouse_location)
        if b_size == 1 or b_size == -1:
            b += b_size
            if b < 1:
                b = 1
            screen.fill(white)
            vertex_list, center_list, edge_list = render_hexes(a,b,c)

        c_size = pmbutton((550, 200), "c", event, mouse_location)
        if c_size == 1 or c_size == -1:
            c += c_size
            if c < 1:
                c = 1
            screen.fill(white)
            render_hexes(a,b,c)

        if square_grid:

            screen.fill(white)
            vertex_list, center_list, edge_list =  render_squares(x_val, y_val)
            hexagonal_grid = False

        if hexagonal_grid:

            screen.fill(white)

            vertex_list, center_list, edge_list = render_hexes(a,b,c)
            square_grid = False
            
            
        save = button((size[0] - 100, size[1] - 100), "save", event, mouse_location)
        if save:
            with open("lattice.txt", "wb") as fh:
                pickle.dump({"vertices" : vertex_list,
                             "edges" : edge_list,
                             "centers" : center_list
                             }, fh)
                print("File Saved")
            

    # display.flip() will update only a portion of the
    # screen to updated, not full area
        pygame.display.flip()

