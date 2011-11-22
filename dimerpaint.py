#!/usr/bin/python -tt

import pickle
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
def rescale(coords, dualcoords, desired_width):
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
    for k in dualcoords.keys():
        v = dualcoords[k]
        new_x = round((v[0] - xmin) * scale_factor)
        new_y = round((ymax - v[1]) * scale_factor)
        dualcoords[k] = (new_x, new_y)


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

# read in a list of rhombi.  Store in a list.
def read_rhombi(filename):
    rhombi = {}
    f = open(filename, 'r')
    for line in f:
        n = [int(s) for s in line.strip().split(',')]
        edge = frozenset([(n[0], n[1]), (n[2],n[3])])
        rhombus = [(n[4],n[5]),(n[6],n[7]),(n[8],n[9]),(n[10],n[11])]
        rhombi[edge] = rhombus
    f.close()
    return rhombi

def render_hexagon_center(coords, h, xoffset, yoffset):
    x = int((coords[h[0]][0] + coords[h[3]][0])/2 + xoffset)
    y = int((coords[h[0]][1] + coords[h[3]][1])/2 + yoffset)

    return pygame.draw.circle(screen, green, (x,y), 4)

# Draw an edge in a given color
def render_edge(coords, edge, xoffset, yoffset, colour, width=4):
    e = [endpt for endpt in edge]
    p0 = coords[e[0]] 
    p0 = (p0[0] + xoffset, p0[1] + yoffset)
    p1 = coords[e[1]] 
    p1 = (p1[0] + xoffset, p1[1] + yoffset)
    return pygame.draw.line(screen, colour, p0, p1, width)

def render_rhombus(dualcoords, rhomb, xoffset, yoffset, colour, width=2):
    coordslist = []
    for p in rhomb:
        p0 = dualcoords[p] 
        p0 = (p0[0] + xoffset, p0[1] + yoffset)
        coordslist.append(p0)
    pygame.draw.polygon(screen, colour, coordslist, width)

def render_double_edge(coords, edge, xoffset, yoffset, colours):
    e = [endpt for endpt in edge]
    p0 = coords[e[0]] 
    p1 = coords[e[1]] 
    normal_x = float(p1[1] - p0[1])
    normal_y = float(p0[0] - p1[0])
    normal_length = math.sqrt(normal_x*normal_x + normal_y*normal_y)
    normal_x *= 1.5/normal_length
    normal_y *= 1.5/normal_length

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
        bb.inflate_ip(2,2) # make it a little bigger for ease of clicking
        #pygame.draw.rect(screen, green, bb, 1)  debug
        boundingboxes.append(("edge", e, which_side, bb))
    return boundingboxes    


# Draw a matching.
def render_matching(coords, matchings, which_side, xoffset, yoffset, color):
    boundingboxes = []
    for e in matchings[which_side].keys():
        bb = render_edge(coords, e, xoffset, yoffset, color)        
        bb.inflate_ip(2,2) # make it a little bigger for ease of clicking
        boundingboxes.append(("matchededge", e, which_side, bb))
    return boundingboxes

# Draw the box pile corresponding to a matching.  This is basically like
# drawing the tiling and then shading the tiles.
def render_boxes(dualcoords, rhombi, matchings, which_side, xoffset, yoffset, color):
    for edge in matchings[which_side].keys():
        try:
            rhomb = rhombi[edge]
            ends = [p for p in edge]
            a = ends[0]
            b = ends[1]
            inverseslope = (a[0]-b[0]) / (a[1]-b[1])
            intensity = {-1: 0.65, 0:0.9, 1:0.75}[inverseslope]
            fillcolor_components = [0,0,0]
            for i in range(3):
                fillcolor_components[i] = 255 * (intensity)  + color[i] * (1-intensity)

            fillcolor = tuple(fillcolor_components)
            print fillcolor
            render_rhombus(dualcoords, rhomb, xoffset, yoffset, fillcolor, 0)        
        except KeyError:
            pass

# Draw the tiling corresponding to a matching.
def render_tiling(dualcoords, rhombi, matchings, which_side, xoffset, yoffset, color):
    for edge in matchings[which_side].keys():
        try:
            rhomb = rhombi[edge]

            render_rhombus(dualcoords, rhomb, xoffset, yoffset, color)        
        except KeyError:
            pass

# Draw the boundary of the matched region.  There's a cheap way to do this:
# the boundary is all the edges in the dual graph that are singly covered by
# a tile edge.
def render_boundary(dualcoords, rhombi, matchings, which_side, xoffset, yoffset, colour):
    rhomb_edges = defaultdict(int)
    for edge in matchings[which_side].keys():
        try:
            rhomb = rhombi[edge]
            for i in range(4):
                j = (i+1) % 4
                rhomb_edge = frozenset([rhomb[i], rhomb[j]])
                rhomb_edges[rhomb_edge] += 1 
        except KeyError:
            pass
    for edge in rhomb_edges.keys():
        if rhomb_edges[edge] == 1:
            render_edge(dualcoords, edge, xoffset, yoffset, colour, 2)

# Redraw doubled edges in a pair of matchings.
def render_doubled_edges(coords, matchings, xoffset, yoffset, colors):
    m1 = matchings[0]
    m2 = matchings[1]
    for e in m1.keys():
        if(e in m2):
            render_double_edge(coords, e, xoffset, yoffset, colors)

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

def maximize_callback(args):
    print "Clicked minimize button on side", args["side"]
    maximize(args["matching"], args["hexagons"])

def minimize_callback(args):
    print "Clicked minimize button on side", args["side"]
    minimize(args["matching"], args["hexagons"])

def randomize_callback(args):
    print "Clicked randomize button on side", args["side"]
    randomize(args["matching"], args["hexagons"])

def quit_callback(args):
    print "Quit button clicked."
    pygame.event.post(pygame.event.Event(pygame.QUIT, {}))

# Toggle visibility of a layer
def showhide_callback(args):
    layer = args["layer"]
    show = args["show"]
    print "Toggle ", layer
    show[layer] = not show[layer]
    
def render_dimer_buttons(x,y, side, renderables, font):
    buttons = []
    matchings = renderables["matchings"] 
    hexagons = renderables["hexagons"]
    show = renderables["show"]
    args = {"side":side, "matching":matchings[side], "hexagons":hexagons}
    if(side==0):
        prefix = "A_"
    else:
        prefix = "B_"

    buttonrow1 = [
        ("Graph", showhide_callback, {"layer":prefix+"background", "show":show}),
        ("Dimer", showhide_callback, {"layer":prefix+"matching", "show":show}),
        ("Tiling", showhide_callback, {"layer":prefix+"tiling", "show":show}),
        ("Boxes", showhide_callback, {"layer":prefix+"boxes", "show":show}),
        ("Border", showhide_callback, {"layer":prefix+"boundary", "show":show}),
        ("Centers", showhide_callback, {"layer":prefix+"centers", "show":show})
        ]

    buttons =  draw_button_row(x, y, 5, font, buttonrow1)
    buttonrow2 = [
        ("Randomize", randomize_callback, args),
        ("Minimize", minimize_callback, args),
        ("Maximize", maximize_callback, args),
        ]
    return buttons + draw_button_row(x, y+20, 5, font, buttonrow2)

def render_center_buttons(x,y,renderables, font):
    buttons = []
    show = renderables["show"]
    buttonrow = [
        ("Dimer A", showhide_callback, {"layer":"center_A_matching", "show":show} ),
        ("Border A", showhide_callback, {"layer":"center_A_boundary", "show":show} ),
        ("Background", showhide_callback, {"layer":"center_background", "show":show} ),
        ("Double edges", showhide_callback, {"layer":"center_doubled_edges", "show":show} )
        ]
    buttonrow2 = [
        ("Dimer B", showhide_callback, {"layer":"center_B_matching", "show":show} ),
        ("Border B", showhide_callback, {"layer":"center_B_boundary", "show":show} ),
        ("Quit", quit_callback, {}),
        ]
    buttons =  draw_button_row(x, y, 5, font, buttonrow)
    return buttons+draw_button_row(x, y+20, 5, font, buttonrow2)

#=============================================================
# Draw everything.  Return a list of clickable boxes.
def render_everything(renderables,font):
    background = renderables["background"] 
    matchings = renderables["matchings"] 
    hexagons = renderables["hexagons"]
    rhombi = renderables["rhombi"]
    coords = renderables["coords"]
    dualcoords = renderables["dualcoords"]
    show = renderables["show"]
    positions = renderables["positions"]
     
    y = positions["y"]
    xA =positions["xA"]
    xB = positions["xB"]
    xDouble = positions["xDouble"]


    screen.fill(white)
    
    boxes = []
    if show["A_background"]:
        boxes += render_background(coords, rhombi, xA, y, 0)
    if show["A_boxes"]:
        render_boxes(dualcoords, rhombi, matchings,0,xA, y, black)
    if show["A_tiling"]:
        render_tiling(dualcoords, rhombi, matchings,0, xA, y, black)
    if show["A_matching"]:
        render_matching(coords, matchings,0, xA, y, black)
    if show["A_boundary"]:
        render_boundary(dualcoords, rhombi, matchings, 0, xA, y, black)
    if show["A_centers"]:
        boxes += render_active_hex_centers(coords, hexagons, matchings[0], 
            xA, y, 0)

    if show["B_background"]:
        boxes += render_background(coords, rhombi, xB, y, 1)
    if show["B_boxes"]:
        render_boxes(dualcoords, rhombi, matchings,1,xB, y, red)
    if show["B_tiling"]:
        render_tiling(dualcoords, rhombi, matchings, 1, xB, y, red)
    if show["B_matching"]:
        render_matching(coords, matchings,1, xB, y, red)
    if show["B_boundary"]:
        render_boundary(dualcoords, rhombi, matchings, 1, xB, y, red)
    if show["B_centers"]:
        boxes += render_active_hex_centers(coords, hexagons, matchings[1], 
            xB, y, 1)

    if show["center_background"]:
        boxes += render_background(coords, rhombi, xDouble, y, 1)
    if show["center_A_boundary"]:
        render_boundary(dualcoords, rhombi, matchings, 0, xDouble, y, black)
    if show["center_B_boundary"]:
        render_boundary(dualcoords, rhombi, matchings, 1, xDouble, y, red)
    if show["center_A_matching"]:
        boxes += render_matching(coords, matchings,0, xDouble, y, black)
    if show["center_B_matching"]:
        boxes += render_matching(coords, matchings,1, xDouble, y, red)
    if show["center_A_matching"] and show["center_B_matching"]:
        if show["center_doubled_edges"]: 
            render_doubled_edges(coords, matchings, xDouble, y, [black, red])
        else:
            render_doubled_edges(coords, matchings, xDouble, y, [white, white])


    boxes += render_dimer_buttons(10+xA, 10, 0, renderables, font)
    boxes += render_dimer_buttons(10+xB, 10, 1, renderables, font)
    boxes += render_center_buttons(10+xDouble, 10, renderables, font)
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

#====================================================================
# Find the maximum matching, with respect to height function.
def maximize(matching, hexlist):
    activelist = []
    finished = False
    while not finished:
        finished = True
        for i in range(len(hexlist)):
            h = hexlist[i]
            e = [frozenset([h[i], h[(i+1)%6]]) for i in range(6)]
            if e[0] in matching and e[2] in matching and e[4] in matching:
                del matching[e[0]]
                del matching[e[2]]
                del matching[e[4]]
                matching[e[1]] = 1
                matching[e[3]] = 1
                matching[e[5]] = 1
                finished = False

#====================================================================
# Find the minimum matching, with respect to height function.
def minimize(matching, hexlist):
    activelist = []
    finished = False
    while not finished:
        finished = True
        for i in range(len(hexlist)):
            h = hexlist[i]
            e = [frozenset([h[i], h[(i+1)%6]]) for i in range(6)]
            if e[1] in matching and e[3] in matching and e[5] in matching:
                del matching[e[1]]
                del matching[e[3]]
                del matching[e[5]]
                matching[e[0]] = 1
                matching[e[2]] = 1
                matching[e[4]] = 1
                finished = False

#=================================================
# Compute coordinates of dual vertices.
def dualcoords(hexlist):
    dualvertexlist = []
    for h in hexlist:
        x = 0;
        y = 0;
        for p in h:
            x += p[0];
            y += p[1];
        center = (x/6, y/6);    
        dualvertexlist.append(center);
    return dualvertexlist;

#===================================================================
# Main program

if(len(sys.argv) != 3):
    exit("usage: dimerpaint (data directory) (input file)")

data_directory_name = sys.argv[1]
input_filename = sys.argv[2]

basename = data_directory_name + "/" + input_filename + "/"

# Setup output directory
if not os.path.isdir(basename):
    exit("Can't find "+basename)

coords = read_vertices(basename + "full.vertex")
background = read_edges(basename + "full.edge")
hexagons = read_hexagons(basename + "full.hexagon")
dualcoords = read_vertices(basename + "full.dualvertex");
rhombi = read_rhombi(basename + "full.rhombus");

region_width = 400
window = region_width + 40
rescale(coords, dualcoords, region_width)

matching_A = read_edges(basename + "A.edge")
matching_B = read_edges(basename + "B.edge")
matchings = [matching_A, matching_B]

#randomize(matching_A, hexagons)
#randomize(matching_B, hexagons)

pygame.init()
pygame.font.init()
font = pygame.font.Font(None, 18)

screen=pygame.display.set_mode([1350,600])
done=False #Loop until the user clicks the close button.
clock=pygame.time.Clock() # Used to manage how fast the screen updates

if os.path.isfile(basename + "show.pkl"):
    showfile = open(basename + "show.pkl", "rb")
    show = pickle.load(showfile)
    showfile.close()
else:
    show = {
        "A_background": True,
        "A_matching": True,
        "A_tiling": False,
        "A_boundary": False,
        "A_centers": True,
        "A_boxes": False,

        "B_background": True,
        "B_matching": True,
        "B_tiling": False,
        "B_boundary": False,
        "B_centers": True,
        "B_boxes": False,

        "center_background": False,
        "center_A_matching": True,
        "center_B_matching": True,
        "center_A_boundary": True,
        "center_B_boundary": True,
        "center_doubled_edges": False,
    }

positions = {
    "y": 45,
    "xA":  0,
    "xB": 2*window,
    "xDouble": window,
    }

renderables = {"background":background, 
               "matchings":matchings, 
               "hexagons":hexagons, 
               "rhombi": rhombi,
               "coords": coords,
               "dualcoords":dualcoords,
               "show":show,
               "positions":positions}

bounding_box_data = render_everything(renderables,font)
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

                bounding_box_data = render_everything(renderables, font)
                bounding_boxes = [record[3] for record in bounding_box_data]

    # Limit to 20 frames per second
    clock.tick(20)
 
pygame.font.quit()
pygame.quit()

# Do the output.
write_edges(matchings[0], basename + "A.edge")
write_edges(matchings[1], basename + "B.edge")
showfile = open(basename + "show.pkl", "wb")
pickle.dump(show, showfile)
showfile.close()

