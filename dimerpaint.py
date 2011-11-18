#!/usr/bin/python -tt

import sys
import math
import pygame
import numpy
import os
import shutil
from collections import defaultdict

# Define some colors
black    = (   0,   0,   0)
white    = ( 255, 255, 255)
blue     = (  50,  50, 255)
green    = (   0, 255,   0)
dkgreen  = (   0, 100,   0)
red      = ( 255,   0,   0)
purple   = (0xBF,0x0F,0xB5)
brown    = (0x55,0x33,0x00)
grey     = (0x7f,0x7f,0x7f) 
 
# read in a list of vertices in format (row, column, x, y).  (row, column) just
# functions as a key.  Store in a dictionary keyed by the above tuple.
def read_vertices(filename):
    coords = {}
    f = open(filename, 'r')
    for line in f:
        (row,col,x,y) = line.strip().split(',')
        coords[int(row),int(col)] = (float(x),float(y))
    f.close()

    return coords

# Rescale the x and y coordinates isometrically so that the width of the 
# entire array is the specified size.
def rescale(coords, desired_width):
    xvalues = [v[0] for v in coords.values()]
    yvalues = [v[1] for v in coords.values()]
    xmin = min(xvalues)
    ymin = min(yvalues)
    xmax = max(xvalues)
    ymax = max(yvalues)
    scale_factor = desired_width / (xmax - xmin) 
    for k in coords.keys():
        v = coords[k]
        new_x = round((v[0] - xmin) * scale_factor)
        new_y = round((ymax - v[1]) * scale_factor)
        coords[k] = (new_x, new_y)


# Read in a list of edges.  Store in a dictionary.
def read_edges(filename):
    edges = {}
    f = open(filename, 'r')
    for line in f:
        (r1,c1,r2,c2) = line.strip().split(',')
        edges[((int(r1),int(c1)),(int(r2),int(c2)))] = 1
    f.close()
    return edges

def write_edges(edges, filename):
    f = open(filename, 'w')
    for e in edges:
        outputlist = [str(i) for i in (e[0][0], e[0][1], e[1][0], e[1][1])]
        output =  outputlist[0] + "," + outputlist[1] + "," 
        output += outputlist[2] + "," + outputlist[3] + "\n"  
        f.write(output)
    f.close()

# Read in a list of hexagons.  Store in a list.
def read_hexagons(filename):
    hexagons = []
    f = open(filename, 'r')
    for line in f:
        n = [int(s) for s in line.strip().split(',')]
        hexagons.append((
                (n[0],n[1]), (n[2],n[3]), (n[4],n[5]),
                (n[6],n[7]), (n[8],n[9]), (n[10],n[11])
            ))
    f.close()
    return hexagons

def render_hexagon_center(coords, h, xoffset, yoffset):
    x = (coords[h[0]][0] + coords[h[3]][0])/2 + xoffset
    y = (coords[h[0]][1] + coords[h[3]][1])/2 + yoffset

    return pygame.draw.circle(screen, green, (x,y), 5)

# Draw an edge in a given color
def render_edge(coords, e, xoffset, yoffset, colour, width=4):
    p0 = coords[e[0]] 
    p0 = (p0[0] + xoffset, p0[1] + yoffset)
    p1 = coords[e[1]] 
    p1 = (p1[0] + xoffset, p1[1] + yoffset)
    return pygame.draw.line(screen, colour, p0, p1, width)

def render_double_edge(coords, e, xoffset, yoffset, colours):
    p0 = coords[e[0]] 
    p1 = coords[e[1]] 
    normal_x = float(p1[1] - p0[1])
    normal_y = float(p0[0] - p1[0])
    normal_length = math.sqrt(normal_x*normal_x + normal_y*normal_y)
    normal_x *= 3/normal_length
    normal_y *= 3/normal_length

    q0 = (p0[0] + xoffset + int(normal_x), p0[1] + yoffset + int(normal_y))
    q1 = (p1[0] + xoffset + int(normal_x), p1[1] + yoffset + int(normal_y))
    bb1 = pygame.draw.line(screen, colours[0], q0, q1, 4)

    q0 = (p0[0] + xoffset - int(normal_x), p0[1] + yoffset - int(normal_y))
    q1 = (p1[0] + xoffset - int(normal_x), p1[1] + yoffset - int(normal_y))
    bb2 = pygame.draw.line(screen, colours[1], q0, q1, 4)

    return bb1.union(bb2)

    
def render_vertex(coords, v, xoffset, yoffset, colour):
    p = (coords[v][0] + xoffset, coords[v][1] + yoffset)
    pygame.draw.circle(screen, colour, p, 2)

# Draw the background.  Return bounding boxes, tagged by "A" or "B".
def render_background(coords, graph, xoffset, yoffset, which_side):
    boundingboxes = []
    for e in graph.keys():
        bb = render_edge(coords, e, xoffset, yoffset, grey, 1)        
        boundingboxes.append(("edge", e, which_side, bb))
    return boundingboxes    

# Draw a matching.
def render_matching(coords, matching, xoffset, yoffset, color):
    for e in matching.keys():
        bb = render_edge(coords, e, xoffset, yoffset, color)        

# Redraw doubled edges in a pair of matchings.
def render_doubled_edges(coords, m1, m2, xoffset, yoffset):
    for e in m1.keys():
        if(e in m2 or (e[1], e[0]) in m2):
            render_double_edge(coords, e, xoffset, yoffset, [black,red])


#=============================================================
# Draw everything.  Return a list of clickable boxes.
def render_everything(background, matchings, hexagons,window):
    screen.fill(white)

    boxes = render_background(coords, background, 0*window, 0*window, 0)
    boxes += render_background(coords, background, 2*window, 0*window, 1)

    render_matching(coords, matchings[0], 0*window, 0*window, black)
    render_matching(coords, matchings[0], 1*window, 0*window, black)
    render_matching(coords, matchings[1], 1*window, 0*window, red)
    render_matching(coords, matchings[1], 2*window, 0*window, red)

    render_doubled_edges(coords, matchings[0],matchings[1], 1*window, 0*window)

    boxes += render_active_hex_centers(coords, hexagons, matchings[0], 
            0*window, 0*window, 0)
    boxes += render_active_hex_centers(coords, hexagons, matchings[1], 
            2*window, 0*window, 1)
    pygame.display.flip()
    return boxes

#==============================================================
# make an adjacency map from a matching. 
def adjacency_map(M): 
    adj = {}
    for edge in M:
        (p0, p1) = edge
        adj[p0] = p1
        adj[p1] = p0
    return adj

#==============================================================
# Can we flip a matching around a hexagon?  We take the matching as an
# adjacency map.  To do this, count the vertices around the hexagon
# which are matched to another point on the hexagon, and see if it is 6. 
def is_active(hexagon, adj):
    acc = 0
    hexdict = {}
    for p in hexagon:
        hexdict[p] = 1
    for p in hexagon:
        if ((p in adj) and (adj[p] in hexdict)): acc += 1
    return(acc == 6)
     
#=============================================================
# Return a list of the hexagons adjacent to a given hexagon
# Some of these might not exist in the graph
#
#              H                               
#                                              
#        H   X---X   H                         
#           /     \                            
#          X   H   X                           
#           \     /                            
#        H   X---X   H                         
#                                              
#              H                               

def hex_neighbours(hexagon):
    (r,c) = hexagon
    return [(r-4,c), (r-2, c+6), (r+2,c+6), (r+4,c), (r+2,c-6), (r-2,c-6)]

#==============================================================
def render_active_hex_centers(coords, hexlist, matching, xoffset, yoffset, 
            which_side):
    adj = adjacency_map(matching)
    boxes = []
    for i in range(len(hexlist)):
        h = hexlist[i]
        if is_active(h, adj): 
            bb = render_hexagon_center(coords, h, xoffset, yoffset)
            boxes.append(("hexagon", i, which_side, bb))
    return boxes

       
             
#==============================================================
# Find a path in the superposition of two matchings, starting at
# a point in the first matching.  Might be a loop.
# The matchings should be given as adjacency maps (see adjacency_map)
# The path is returned as a list of vertices.  If the path is closed then
# the starting vertex is repeated.
def find_path(adj1, adj2, start):
    path = [start];
    p1 = start;
    try:
        p2 = adj2[p1]
        path.append(p2)
        p1 = adj1[p2]
        while(p1 != start):
            path.append(p1)
            p2 = adj2[p1]
            path.append(p2)
            p1 = adj1[p2]
        path.append(p1)
        return path
    except KeyError: 
        pass

# We fall through to here if path isn't closed.  Get the rest of the path.

    p1 = start
    try:  
        p2 = adj1[p1]
        path.insert(0,p2)
        p1 = adj2[p2]
        while(p1 != start):
            path.insert(0,p1)
            p2 = adj1[p1]
            path.insert(0,p2)
            p1 = adj2[p2]

    except KeyError: 
        pass

    return path

# Flip a path in two matchings, starting at an edge in the first.
def flip_path(matchings, m1, m2, edge):
    adj1 = adjacency_map(matchings[m1])
    adj2 = adjacency_map(matchings[m2])
    path = find_path(adj1, adj2, edge[0])

    # don't do anything if user clicked a doubled edge.
    if(len(path) == 3 and path[0] == path[2]): return 


    
    # make lists of edges corresponding to the path
    loop1 = []
    loop2 = []
    for i in range(len(path) - 1):
        for e in [(path[i], path[i+1]), (path[i+1], path[i])]:
            if e in matchings[m1]: loop1.append(e)
            if e in matchings[m2]: loop2.append(e)
    for e in loop1:
        del matchings[m1][e]
        matchings[m2][e] = 1
    for e in loop2:
        del matchings[m2][e]
        matchings[m1][e] = 1

def flip_hex(matching, hexagons, index):
    h = hexagons[index]
    edges = [(h[i], h[(i+1)%6]) for i in range(6)]
    for e in edges:
        backwards_e = (e[1], e[0])
        if(e in matching): 
            del matching[e]
        else:
            if(backwards_e in matching):
                del matching[backwards_e]
            else:
                matching[e] = 1

#====================================================================
def randomize(matching, hexlist):
    for trial in range(5000):
        # Choose a random index i
        adj = adjacency_map(matching)
        activelist = []
        for i in range(len(hexlist)):
            if is_active(hexlist[i],adj):
                activelist.append(i)

        i = numpy.random.randint(len(activelist))
        if is_active(hexlist[activelist[i]], adj):
            flip_hex(matching, hexlist, activelist[i]) 


#===================================================================
# Main program

if(len(sys.argv) != 3):
    exit("usage: dimerpaint (input) (output)")

basename = sys.argv[1]
basename += "/"
outputbasename = sys.argv[2]
# Setup output directory
if not os.path.isdir(outputbasename):
    os.mkdir(outputbasename)
outputbasename += "/"

coords = read_vertices(basename + "full.vertex")
background = read_edges(basename + "full.edge")
hexagons = read_hexagons(basename + "full.hexagon")

if(basename != outputbasename): 
    shutil.copyfile(basename + "full.vertex", outputbasename + "full.vertex")
    shutil.copyfile(basename + "full.edge", outputbasename + "full.edge")
    shutil.copyfile(basename + "full.hexagon", outputbasename + "full.hexagon")

region_width = 400
window = region_width + 40
rescale(coords, region_width)

matching_A = read_edges(basename + "A.edge")
matching_B = read_edges(basename + "B.edge")
matchings = [matching_A, matching_B]

#randomize(matching_A, hexagons)
#randomize(matching_B, hexagons)

pygame.init()
screen=pygame.display.set_mode([1350,600])
done=False #Loop until the user clicks the close button.
clock=pygame.time.Clock() # Used to manage how fast the screen updates

bounding_box_data = render_everything(background, matchings, hexagons,window)
bounding_boxes = [record[3] for record in bounding_box_data]




#=================================================
# Pygame event loop
while done==False:
    for event in pygame.event.get(): # User did something
        if event.type == pygame.QUIT: # If user clicked close
            done=True # Flag that we are done so we exit this loop
        if event.type == pygame.MOUSEBUTTONUP:
            radius = 2
            ul = (event.pos[0] - radius, event.pos[1] - radius)
            clickpoint = pygame.Rect(ul , (2*radius,2*radius))
            index = clickpoint.collidelist(bounding_boxes)
            if(index != -1):
                (objtype, index, side, box) = bounding_box_data[index]
                if objtype == "edge":
                    pass
                    #if index in matchings[side]:
                    #    del matchings[side][index]
                    #else:
                    #    matchings[side][index] = 1
                else:
                    flip_hex(matchings[side], hexagons, index)

                bounding_box_data = render_everything(background,matchings, 
                            hexagons, window)
                bounding_boxes = [record[3] for record in bounding_box_data]

    # Limit to 20 frames per second
    clock.tick(20)
 
pygame.quit()

# Do the output.

write_edges(matchings[0], outputbasename + "A.edge")
write_edges(matchings[1], outputbasename + "B.edge")
