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


# Read in a list of edges.  Store in a dictionary.  Each edge is a "frozen set" 
# which is an immutable, unordered tuple
def read_edges(filename):
    edges = {}
    f = open(filename, 'r')
    for line in f:
        (r1,c1,r2,c2) = line.strip().split(',')
        e = frozenset([(int(r1),int(c1)),(int(r2),int(c2))])
        edges[e] = 1
    f.close()
    return edges

def write_edges(edges, filename):
    f = open(filename, 'w')
    for e in edges:
        ends = [point for point in e]  
        outputlist = [str(i) for i in (ends[0][0], ends[0][1], ends[1][0], ends[1][1])]
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
    x = int((coords[h[0]][0] + coords[h[3]][0])/2 + xoffset)
    y = int((coords[h[0]][1] + coords[h[3]][1])/2 + yoffset)

    return pygame.draw.circle(screen, green, (x,y), 5)

# Draw an edge in a given color
def render_edge(coords, edge, xoffset, yoffset, colour, width=4):
    e = [endpt for endpt in edge]
    p0 = coords[e[0]] 
    p0 = (p0[0] + xoffset, p0[1] + yoffset)
    p1 = coords[e[1]] 
    p1 = (p1[0] + xoffset, p1[1] + yoffset)
    return pygame.draw.line(screen, colour, p0, p1, width)

def render_double_edge(coords, edge, xoffset, yoffset, colours):
    e = [endpt for endpt in edge]
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
        bb.inflate_ip(3,3) # make it a little bigger for ease of clicking
        #pygame.draw.rect(screen, green, bb, 1)  debug
        boundingboxes.append(("edge", e, which_side, bb))
    return boundingboxes    


# Draw a matching.
def render_matching(coords, matchings, which_side, xoffset, yoffset, color):
    boundingboxes = []
    for e in matchings[which_side].keys():
        bb = render_edge(coords, e, xoffset, yoffset, color)        
        bb.inflate_ip(3,3) # make it a little bigger for ease of clicking
        boundingboxes.append(("matchededge", e, which_side, bb))
    return boundingboxes

# Redraw doubled edges in a pair of matchings.
def render_doubled_edges(coords, m1, m2, xoffset, yoffset):
    for e in m1.keys():
        if(e in m2):
            render_double_edge(coords, e, xoffset, yoffset, [black,red])

#===================================================================
# Create a button, with given text and center, and a callback function
# for what to do when it's clicked.  Callbacks get a dict, "args", as
# their arguments.  I make no attempt to police what goes into args.
#
# For alignment reasons "pos" is the middle of the left side of the button.

def drawbutton(pos, label, font, callback, args):
    text = font.render(label, True, black, white)    
    bb = text.get_rect()
    bb.midleft = pos
    screen.blit(text, bb)
    bb.inflate_ip(4,4)
    pygame.draw.rect(screen, black, bb, 1)
    return ("button", args, callback, bb)

# Draw a row of buttons.  Return a list of their bounding boxes.
# Argumnet is a list of pairs: (name, callback, args)

def draw_button_row(x, y, spacing, font, buttondata):
    boundingboxes = []
    x_current = x
    for button in buttondata:
        pos = (x_current, y)
        bb = drawbutton(pos, button[0], font, button[1], button[2])
        x_current += bb[3].width + spacing
        boundingboxes.append(bb)
    return boundingboxes

def test_callback(args):
    print "Clicked test button on side" , args["side"]

def randomize_callback(args):
    print "Clicked randomize button on side", args["side"]
    randomize(args["matching"], args["hexagons"])

def quit_callback(args):
    print "Quit button clicked."
    pygame.event.post(pygame.event.Event(pygame.QUIT, {}))

def render_dimer_buttons(x,y,matchings, side, hexagons, font):
    args = {"side":side, "matching":matchings[side], "hexagons":hexagons}

    buttonrow = [
        ("Randomize", randomize_callback, args),
        ("Testing", test_callback, args)
        ]

    return draw_button_row(x, y, 5, font, buttonrow)

def render_center_buttons(x,y,matchings, hexagons, font):
    buttonrow = [
        ("Quit", quit_callback, {}),
        ("Testing", test_callback, {"side":2} )
        ]
    return draw_button_row(x, y, 5, font, buttonrow)

#=============================================================
# Draw everything.  Return a list of clickable boxes.
def render_everything(background, matchings, hexagons,window,font):
    screen.fill(white)
    y = 30
    
    boxes = render_background(coords, background, 0*window, y, 0)
    boxes += render_background(coords, background, 2*window, y, 1)

    render_matching(coords, matchings,0, 0*window, y, black)
    boxes += render_matching(coords, matchings,0, 1*window, y, black)
    boxes += render_matching(coords, matchings,1, 1*window, y, red)
    render_matching(coords, matchings,1, 2*window, y, red)


    render_doubled_edges(coords, matchings[0],matchings[1], 1*window, y)

    boxes += render_active_hex_centers(coords, hexagons, matchings[0], 
            0*window, y, 0)
    boxes += render_active_hex_centers(coords, hexagons, matchings[1], 
            2*window, y, 1)

    boxes += render_dimer_buttons(10+0*window, 10, matchings, 0, hexagons, font)
    boxes += render_dimer_buttons(10+2*window, 10, matchings, 1, hexagons, font)

    boxes += render_center_buttons(10+window, 10, matchings, hexagons, font)
    pygame.display.flip()
    return boxes

#==============================================================
# make an adjacency map from a matching. 
def adjacency_map(M): 
    adj = {}
    for edge in M:
        endpoints = [endpt for endpt in edge]
        adj[endpoints[0]] = endpoints[1]
        adj[endpoints[1]] = endpoints[0]
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
def flip_path(matchings, m1, m2, unordered_edge):
    edge = [endpoint for endpoint in unordered_edge]
    adj1 = adjacency_map(matchings[m1])
    adj2 = adjacency_map(matchings[m2])
    path = find_path(adj1, adj2, edge[0])

    # don't do anything if user clicked a doubled edge.
    if(len(path) == 3 and path[0] == path[2]): return 

    print "made it"
    print path
    
    # make lists of edges corresponding to the path
    loop1 = []
    loop2 = []
    for i in range(len(path) - 1):
        e = frozenset([path[i], path[i+1]])
        if e in matchings[m1]: loop1.append(e)
        if e in matchings[m2]: loop2.append(e)
    for e in loop1:
        del matchings[m1][e]
        matchings[m2][e] = 1
    for e in loop2:
        del matchings[m2][e]
        matchings[m1][e] = 1

#============================================================================
# Flip one hexagon.
def flip_hex(matching, hexagons, index):
    h = hexagons[index]
    edges = [frozenset([h[i], h[(i+1)%6]]) for i in range(6)]
    for e in edges:
        if(e in matching): 
            del matching[e]
        else:
            matching[e] = 1

#====================================================================
# Run the glauber dynamics to randomize one picture
def randomize(matching, hexlist):
    for trial in range(500):
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
pygame.font.init()
font = pygame.font.Font(None, 18)

screen=pygame.display.set_mode([1350,450])
done=False #Loop until the user clicks the close button.
clock=pygame.time.Clock() # Used to manage how fast the screen updates

bounding_box_data = render_everything(background, matchings, hexagons,window,font)
bounding_boxes = [record[3] for record in bounding_box_data]




#=================================================
# Pygame event loop
while done==False:
    for event in pygame.event.get(): # User did something
        if event.type == pygame.QUIT: # If user clicked close
            done=True # Flag that we are done so we exit this loop
        if event.type == pygame.MOUSEBUTTONUP:
            radius = 1
            ul = (event.pos[0] - radius, event.pos[1] - radius)
            clickpoint = pygame.Rect(ul , (2*radius,2*radius))
            index = clickpoint.collidelist(bounding_boxes)
            if(index != -1):
                objtype = bounding_box_data[index][0]
                if objtype == "edge":
                    (objtype, index, side, box) = bounding_box_data[index]
                    if index in matchings[side]:
                        print "deleting"
                        del matchings[side][index]
                    else:
                        print "adding"
                        matchings[side][index] = 1
                elif objtype == "matchededge":
                    (objtype, index, side, box) = bounding_box_data[index]
                    flip_path(matchings, side, 1-side, index)
                elif objtype == "button":
                    (objtype, args, callback, bb) = bounding_box_data[index]
                    callback(args)
                elif objtype == "hexagon":
                    (objtype, index, side, box) = bounding_box_data[index]
                    flip_hex(matchings[side], hexagons, index)
                else:
                    print "uh why am I here?"

                bounding_box_data = render_everything(background,matchings, 
                            hexagons, window, font)
                bounding_boxes = [record[3] for record in bounding_box_data]

    # Limit to 20 frames per second
    clock.tick(20)
 
pygame.font.quit()
pygame.quit()

# Do the output.

write_edges(matchings[0], outputbasename + "A.edge")
write_edges(matchings[1], outputbasename + "B.edge")
