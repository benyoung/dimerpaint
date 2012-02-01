#!/usr/bin/python -tt

import pickle
import copy
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
green    = (  75, 255,   0)
dkgreen  = (   0, 100,   0)
red      = ( 255,   0,   0)
purple   = (0xBF,0x0F,0xB5)
brown    = (0x55,0x33,0x00)
grey     = (0x7f,0x7f,0x7f) 

# This is a quick hack to make it possible to use poor-quality projectors.
# Edit these lines to alter the colors used.

hi_color = green
A_color = black
B_color = red
center_color = blue

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

def render_hexagon_center(coords, h, xoffset, yoffset, radius, eps):
    x = int((coords[h[0]][0] + coords[h[3]][0])/2 + xoffset)
    y = int((coords[h[0]][1] + coords[h[3]][1])/2 + yoffset)

# I'm supposed to render these in eps I guess?  Cant think why right now
    return pygame.draw.circle(screen, center_color, (x,y), radius)

# Draw an edge in a given color
def render_edge(coords, edge, xoffset, yoffset, colour, width, eps):
    e = [endpt for endpt in edge]
    n = normal(coords, edge, width/2.0)

    p0 = coords[e[0]] 
    p1 = coords[e[1]] 
    p2 = coords[e[1]] 
    p3 = coords[e[0]] 

    eps["coords"].extend([p0, p1]) 
    ps = "%f %f %f setrgbcolor " %(colour[0] / 255.0, colour[1] / 255.0, colour[2] / 255.0) 
    ps += "%d setlinewidth newpath %d %d moveto %d %d lineto stroke" % (width, p0[0], p0[1], p1[0], p1[1])
    eps["ps"].append(ps)

    p0 = (p0[0] + xoffset + n[0], p0[1] + yoffset + n[1])
    p1 = (p1[0] + xoffset + n[0], p1[1] + yoffset + n[1])
    p2 = (p2[0] + xoffset - n[0], p2[1] + yoffset - n[1])
    p3 = (p3[0] + xoffset - n[0], p3[1] + yoffset - n[1])
    
    #return pygame.draw.line(screen, colour, p0, p1, width)
    return pygame.draw.polygon(screen, colour, [p0,p1,p2,p3], 0)

def render_line(coords, edge, xoffset, yoffset, colour, width, eps):
    e = [endpt for endpt in edge]
    p0 = coords[e[0]] 
    p1 = coords[e[1]] 

    eps["coords"].extend([p0,p1])
    ps = "%f %f %f setrgbcolor " % (colour[0] / 255.0, colour[1] / 255.0, colour[2] / 255.0)
    ps += "%d setlinewidth newpath %d %d moveto %d %d lineto stroke" % (width, p0[0], p0[1], p1[0], p1[1])
    eps["ps"].append(ps)

    p0 = (p0[0] + xoffset, p0[1] + yoffset)
    p1 = (p1[0] + xoffset, p1[1] + yoffset)


    return pygame.draw.line(screen, colour, p0, p1, width)


def render_rhombus(dualcoords, rhomb, xoffset, yoffset, colour, width, eps):
    coordslist = []
    for p in rhomb:
        p0 = dualcoords[p] 
        coordslist.append(p0)


    eps["coords"].extend(coordslist) 
    flatlist = []
    for p in coordslist:
        flatlist.extend([p[0],p[1]])
    ps = "%f %f %f setrgbcolor " % (colour[0] / 255.0, colour[1] / 255.0, colour[2] / 255.0)
    ps += "newpath %d %d moveto %d %d lineto %d %d lineto %d %d lineto closepath" % tuple(flatlist)
    if width==0:
        ps += " fill"
    else:
        ps = ("%d setlinewidth " % width) + ps + " stroke"
    eps["ps"].append(ps)

    shiftedcoordslist = []
    for p0 in coordslist:
        p0 = (p0[0] + xoffset, p0[1] + yoffset)
        shiftedcoordslist.append(p0)
    pygame.draw.polygon(screen, colour, shiftedcoordslist, width)



# Compute an integer normal vector to an edge of a given magnitude.
# The edge can be implemented as a list, tuple or set of points, of length 2.
def normal(coords, edge, magnitude):
    e = [endpt for endpt in edge]
    p0 = coords[e[0]] 
    p1 = coords[e[1]] 
    normal_x = float(p1[1] - p0[1])
    normal_y = float(p0[0] - p1[0])
    normal_length = math.sqrt(normal_x*normal_x + normal_y*normal_y)
    normal_x *= magnitude/normal_length
    normal_y *= magnitude/normal_length
    return (normal_x, normal_y)

def render_double_edge(coords, edge, xoffset, yoffset, colours, width, eps):
    e = [endpt for endpt in edge]
    p0 = coords[e[0]] 
    p1 = coords[e[1]] 

    norm = normal(coords, edge, width*2/3) 

    pts = [  (p0[0]  + norm[0], p0[1]  + norm[1]),
             (p1[0]  + norm[0], p1[1]  + norm[1]),
             (p0[0]  - norm[0], p0[1] - norm[1]),
             (p1[0]  - norm[0], p1[1] - norm[1]), ]

    eps["coords"].extend(pts) 
    ps = "%f %f %f setrgbcolor "%(colours[0][0]/255.0,colours[0][1]/255.0,colours[0][2]/255.0) 
    ps += "%d setlinewidth newpath %d %d moveto %d %d lineto stroke" % (width,pts[0][0],
            pts[0][1],pts[1][0],pts[1][1])
    eps["ps"].append(ps)
    ps = "%f %f %f setrgbcolor "%(colours[1][0]/255.0,colours[1][1]/255.0,colours[1][2]/255.0) 
    ps += "%d setlinewidth newpath %d %d moveto %d %d lineto stroke" % (width,pts[2][0],
            pts[2][1],pts[3][0],pts[3][1])
    eps["ps"].append(ps)

    norm = normal(coords, edge, width/2) 
    draw_these = [[pts[0],pts[1]],[pts[2],pts[3]]]
    boxes = []
    for i in range(2):
        edge = draw_these[i]
        poly = [
            (edge[0][0] + xoffset + norm[0], edge[0][1] + yoffset + norm[1]),
            (edge[1][0] + xoffset + norm[0], edge[1][1] + yoffset + norm[1]),
            (edge[1][0] + xoffset - norm[0], edge[1][1] + yoffset - norm[1]),
            (edge[0][0] + xoffset - norm[0], edge[0][1] + yoffset - norm[1]), ]
    
        boxes += pygame.draw.polygon(screen, colours[i], poly, 0)


    return boxes

    
def render_vertex(coords, v, xoffset, yoffset, colour):
    p = (coords[v][0] + xoffset, coords[v][1] + yoffset)
    pygame.draw.circle(screen, colour, p, 2)

# Draw the background.  Return bounding boxes, tagged by "A" or "B".
def render_background(renderables, xoffset, yoffset, which_side, eps):
    coords = renderables["coords"] 
    # Turns out that it's nicer to render edges where all the rhombi are.
    #graph = renderables["background"]  
    graph = renderables["rhombi"]
    boundingboxes = []
    for e in graph.keys():
        bb = render_line(coords, e, xoffset, yoffset, grey, 1, eps)        
        bb.inflate_ip(4,4) # make it a little bigger for ease of clicking
        boundingboxes.append(("edge", e, which_side, bb))
    return boundingboxes    


# Draw a matching.
def render_matching(renderables, which_side, xoffset, yoffset, color, eps):
    coords = renderables["coords"] 
    matchings = renderables["matchings"]
    lengths = renderables["lengths"]
    boundingboxes = []
    for e in matchings[which_side].keys():
        bb = render_edge(coords, e, xoffset, yoffset, color, lengths["dimer_width"], eps)        
        bb.inflate_ip(2,2) # make it a little bigger for ease of clicking
        boundingboxes.append(("matchededge", e, which_side, bb))
    return boundingboxes

# Draw the box pile corresponding to a matching.  This is basically like
# drawing the tiling and then shading the tiles.
def render_boxes(renderables, which_side, xoffset, yoffset, color, eps):
    dualcoords = renderables["dualcoords"]
    rhombi = renderables["rhombi"]
    matchings = renderables["matchings"]
    try:
        global_intensity = renderables["lengths"]["shading_intensity"]
    except KeyError:
        global_intensity = 1

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
                fillcolor_components[i] = 255 * intensity  + color[i] * (1-intensity)
                fillcolor_components[i] = 255 * (1-global_intensity) + fillcolor_components[i] * global_intensity 
                if fillcolor_components[i] < 0:
                    fillcolor_components[i] = 0
            fillcolor = tuple(fillcolor_components)
            render_rhombus(dualcoords, rhomb, xoffset, yoffset, fillcolor, 0, eps)        
        except KeyError:
            pass

# Draw the tiling corresponding to a matching.
def render_tiling(renderables, which_side, xoffset, yoffset, color, eps):
    dualcoords = renderables["dualcoords"]
    rhombi = renderables["rhombi"]
    matchings = renderables["matchings"]
    try:
        width = renderables["lengths"]["tile_edge_width"]
    except KeyError:
        width = 2
        renderables["lengths"]["tile_edge_width"] = width
    for edge in matchings[which_side].keys():
        try:
            rhomb = rhombi[edge]
            render_rhombus(dualcoords, rhomb, xoffset, yoffset, color, width, eps)        
        except KeyError:
            pass

def render_highlight(renderables, which_side, xoffset, yoffset, color, eps):
    #for e in renderables["matchings"][which_side]:
    for e in renderables["highlight"][which_side]:
        rhomb = renderables["rhombi"][e]
        dualcoords = renderables["dualcoords"]
        render_rhombus(dualcoords, rhomb, xoffset, yoffset,color, 0, eps)


# Draw the boundary of the matched region.  There's a cheap way to do this:
# the boundary is all the edges in the dual graph that are singly covered by
# a tile edge.
def render_boundary(renderables, which_side, xoffset, yoffset, colour, eps):
    dualcoords =renderables["dualcoords"]
    rhombi =renderables["rhombi"]
    matchings =renderables["matchings"]
    rhomb_edges = defaultdict(int)
    try:
        width = renderables["lengths"]["tile_edge_width"]
    except KeyError:
        width = 2
        renderables["lengths"]["tile_edge_width"] = width

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
            render_line(dualcoords, edge, xoffset, yoffset, colour, width, eps)

# Draw doubled edges in a pair of matchings.
def render_doubled_edges(renderables, xoffset, yoffset, colors, eps):
    coords = renderables["coords"] 
    matchings = renderables["matchings"]
    lengths = renderables["lengths"]
    m1 = matchings[0]
    m2 = matchings[1]
    for e in m1.keys():
        if(e in m2):
            render_double_edge(coords, e, xoffset, yoffset, colors, lengths["dimer_width"],eps)

# Draw edges that are not doubled in a pair of matchings.
def render_xor_edges(renderables, xoffset, yoffset, colors,eps):
    coords = renderables["coords"]
    matchings = renderables["matchings"]
    lengths = renderables["lengths"]
    boxes = []
    m0 = matchings[0]
    m1 = matchings[1]
    for e in m0.keys():
        if e not in m1:
            bb = render_edge(coords, e, xoffset, yoffset, colors[0], lengths["dimer_width"],eps)     
            boxes.append(("matchededge", e, 0, bb))
    for e in m1.keys():
        if e not in m0:
            bb = render_edge(coords, e, xoffset, yoffset, colors[1], lengths["dimer_width"],eps)
            boxes.append(("matchededge", e, 1, bb))
    return boxes


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

# Draw a plus/minus button for adjusting lengths.
def draw_adjuster(pos, label, quantity, amount, lengths,font):
    buttons = [
        ("-"+str(amount), adjust_callback, {"quantity":quantity, "amount":-amount, "lengths":lengths}),
        (label, null_callback, {}),
        ("+"+str(amount), adjust_callback, {"quantity":quantity, "amount":amount,"lengths":lengths}),
        ]
    return draw_button_row(pos[0], pos[1], -1, font, buttons)

# Draw a row of plus/minus buttons for adjusting lengths.
# The last argument is a list of tuples: (label, quantity, amount)
def draw_adjuster_row(x, y, spacing, lengths, font, adjusterlist):
    boundingboxes = []
    x_current = x
    for adjuster in adjusterlist:
        pos = (x_current, y)
        bb = draw_adjuster(pos, adjuster[0], adjuster[1], adjuster[2], lengths, font) 
        x_current += bb[0][3].width + bb[1][3].width + bb[2][3].width -2 + spacing
        boundingboxes.extend(bb)
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
    randomize(args["matching"], args["hexagons"], args["steps"])

def adjust_callback(args):
    lengths = args["lengths"]
    quantity = args["quantity"]
    amount = args["amount"]
    lengths[quantity] += amount
    if(lengths[quantity] < 0):
        lengths[quantity] = 0
    print "%s set to %f" % (quantity, lengths[quantity])

def null_callback(args):
    pass

def fullscreen_callback(args):
    return;
    renderables = args["renderables"]

    if "old_screen_size" in renderables["lengths"]:
        pygame.display.set_mode(renderables["lengths"]["old_screen_size"], pygame.RESIZABLE)
        compute_picture_sizes(renderables)
        del renderables["lengths"]["old_screen_size"]
    else:
        modes = pygame.display.list_modes()
        if modes:
            renderables["lengths"]["old_screen_size"] = screen.get_size()
            pygame.display.set_mode(modes[0], pygame.FULLSCREEN)
            compute_picture_sizes(renderables)



def clear_highlight_callback(args):
    print "Clear all highlighting."
    args["renderables"]["highlight"] = [{},{}]

def quit_callback(args):
    print "Quit button clicked."
    pygame.event.post(pygame.event.Event(pygame.QUIT, {}))

def load_callback(args):
    myname = args["myname"]
    f = args["filenames"]
    print "Load button clicked:" + args["myname"]
    old_basename = f["basename"]
    save(old_basename, args["renderables"])

    new_basename = f["data_directory"] + "/" + myname + "/"
    f["input_file"] = myname
    f["basename"] = new_basename

    renderables = args["renderables"]
    new_renderables = load(filenames["basename"])
    for item in new_renderables.keys():
        renderables[item] = new_renderables[item]
    print "Looks like I loaded it successfully"

def eps_callback(args):
    print "Print to EPS button clicked"
    coords = args["eps"]["coords"]
    ps = args["eps"]["ps"]
    filename = args["eps"]["filename"]

    xvalues = [v[0] for v in coords]
    yvalues = [v[1] for v in coords]
    boundary = 10
    xmin = min(xvalues) - boundary
    ymin = min(yvalues) - boundary
    xmax = max(xvalues) + boundary
    ymax = max(yvalues) + boundary

    yc = (ymin + ymax) / 2


    headerlines = [
        "%!PS-Adobe-3.0 EPSF-3.0",
        "%%%%BoundingBox: %d %d %d %d" % (xmin, ymin, xmax, ymax),
        "%%%%HiResBoundingBox: %d %d %d %d" % (xmin, ymin, xmax, ymax),
        "%%EndComments",
        "1 setlinecap", 
        "0 %d translate" % yc,
        "1 -1 scale",
        "0 %d translate" % -yc]

    epsfile = open(filename, "w")
    for line in headerlines:
        epsfile.write(line + "\n")
    for line in ps:
        epsfile.write(line + "\n")
    epsfile.write("showpage\n")
    epsfile.close()
    #if(os.fork() == 0): # fork a child process
    #    os.execl('./showps', '')

# Toggle visibility of a layer
def showhide_callback(args):
    layer = args["layer"]
    show = args["show"]
    print "Toggle ", layer
    show[layer] = not show[layer]

# Toggle visibility of one of the 3 displays; resize pictures.
def showhide_picture_callback(args):
    renderables = args["renderables"]
    pictures = args["pictures"]
    show = renderables["show"]
    print "Showing some pictures:", pictures
    for picture in ["A", "Center", "B"]:
        show[picture] = (picture in pictures)
    compute_picture_sizes(renderables)
    
def render_dimer_buttons(x,y, side, renderables, font, eps):
    buttons = []
    matchings = renderables["matchings"] 
    hexagons = renderables["hexagons"]
    show = renderables["show"]
    args = {"side":side, "matching":matchings[side], "hexagons":hexagons}
    filename = eps["filename"]
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
    args["steps"]=renderables["lengths"]["randomize_steps"]
    buttonrow2 = [
        ("Randomize", randomize_callback, args),
        ("Minimize", minimize_callback, args),
        ("Maximize", maximize_callback, args),
        ("EPS", eps_callback, {"eps":eps}),
        ]
    return buttons + draw_button_row(x, y+20, 5, font, buttonrow2)

def render_center_buttons(x,y,renderables, font, eps):
    buttons = []
    show = renderables["show"]
    buttonrow = [
        ("Dimer A", showhide_callback, {"layer":"center_A_matching", "show":show} ),
        ("Border A", showhide_callback, {"layer":"center_A_boundary", "show":show} ),
        ("Background", showhide_callback, {"layer":"center_background", "show":show} ),
        ]
    buttonrow2 = [
        ("Dimer B", showhide_callback, {"layer":"center_B_matching", "show":show} ),
        ("Border B", showhide_callback, {"layer":"center_B_boundary", "show":show} ),
        ("Double edges", showhide_callback, {"layer":"center_doubled_edges", "show":show} ),
        ("EPS", eps_callback, {"eps":eps}),
        ]
    buttons =  draw_button_row(x, y, 5, font, buttonrow)
    return buttons+draw_button_row(x, y+20, 5, font, buttonrow2)

def render_os_buttons(x, y, filenames, renderables, font):
    buttons = []
    datadir = filenames["data_directory"]
    files = sorted(os.listdir(datadir))
    buttonrow = [
        ("Quit", quit_callback, {}),
    ]
    buttonrow.extend([( f, 
                        load_callback, 
                        {   "myname":f, 
                            "filenames":filenames, 
                            "renderables":renderables }) for f in files])
    return draw_button_row(x,y,5,font,buttonrow)

def render_showhide_buttons(x,y, renderables, font):
    buttonrow = [
        ("Just A", showhide_picture_callback, {"pictures":["A"], "renderables":renderables} ),
        ("A union B", showhide_picture_callback,{"pictures":["Center"], "renderables":renderables} ),
        ("Just B", showhide_picture_callback, {"pictures":["B"], "renderables":renderables} ),
        ("A beside B", showhide_picture_callback,{"pictures":["A","B"], "renderables":renderables} ),
        ("All three", showhide_picture_callback,{"pictures":["A","Center", "B"], 
                "renderables":renderables} ),
        ("Clear highlight", clear_highlight_callback,{"renderables":renderables} )

        #("Full screen", fullscreen_callback, {"renderables":renderables}),
        
    ]
    return draw_button_row(x,y,5,font,buttonrow)

def render_global_adjusters(x,y, renderables, font):
    adjusterlist = [
        ("Dimer width", "dimer_width", 1),
        ("Center radius", "hex_flipper_radius", 1),
        ("Overlay offset", "overlay_offset", 1),
        ("Tile edge", "tile_edge_width", 1),
        ("Shading intensity", "shading_intensity", 0.1),
        ("Randomize steps", "randomize_steps", 100)
    ]
    return draw_adjuster_row(x,y,5,renderables["lengths"],font,  adjusterlist)

#=============================================================
# Draw everything.  Return a list of clickable boxes.
def render_everything(renderables,filenames, font):
    background = renderables["background"] 
    matchings = renderables["matchings"] 
    hexagons = renderables["hexagons"]
    rhombi = renderables["rhombi"]
    coords = renderables["coords"]
    dualcoords = renderables["dualcoords"]
    show = renderables["show"]
    lengths = renderables["lengths"]

    current_file = filenames["input_file"]
     
    y = lengths["y"]
    xA =lengths["xA"]
    xB = lengths["xB"]
    xCenter = lengths["xCenter"]
    window = lengths["window"]

    screen.fill(white)
    
    boxes = []

    if show["A"]:
        epsA = {"coords":[], "ps":[], "filename":current_file + "_A.eps"}
        if show["A_background"]:
            boxes += render_background(renderables, xA, y, 0, epsA)
        if show["A_boxes"]:
            render_boxes(renderables,0,xA, y, A_color, epsA)
        if show["Highlight"]:
            render_highlight(renderables,0, xA, y, hi_color, epsA)
        if show["A_matching"]:
            render_matching(renderables,0, xA, y, A_color, epsA)
        if show["A_boundary"]:
            render_boundary(renderables,  0, xA, y, A_color, epsA)
        if show["A_tiling"]:
            render_tiling(renderables,0, xA, y, A_color, epsA)
        if show["A_centers"]:
            boxes += render_active_hex_centers(renderables, xA, y, 0,epsA)
        boxes += render_dimer_buttons(10+xA, 10, 0, renderables, font, epsA)
        renderables["epsA"] = epsA
        if eps_only: # hack to render eps programattically
            eps_callback({"eps":epsA})

    if show["B"]:
        epsB = {"coords":[], "ps":[], "filename": current_file + "_B.eps"}
        if show["B_background"]:
            boxes += render_background(renderables, xB, y, 1, epsB)
        if show["B_boxes"]:
            render_boxes(renderables,1,xB, y, B_color, epsB)
        if show["Highlight"]:
            render_highlight(renderables,1, xB, y, hi_color, epsB)
        if show["B_matching"]:
            render_matching(renderables,1, xB, y, B_color, epsB)
        if show["B_boundary"]:
            render_boundary(renderables, 1, xB, y, B_color, epsB)
        if show["B_tiling"]:
            render_tiling(renderables, 1, xB, y, B_color, epsB)
        if show["B_centers"]:
            boxes += render_active_hex_centers(renderables, xB, y, 1, epsB)
        boxes += render_dimer_buttons(10+xB, 10, 1, renderables, font, epsB)
        renderables["epsB"] = epsB
        if eps_only: # hack to render eps programattically
            eps_callback({"eps":epsB})
        
    if show["Center"]:
        epsCenter = {"coords":[], "ps":[], "filename":current_file + "_AB.eps"}
        if show["center_background"]:
            render_background(renderables, xCenter, y, 1, epsCenter)
        if show["Highlight"]:
            render_highlight(renderables,0, xCenter, y, hi_color, epsCenter)
            render_highlight(renderables,1, xCenter, y, hi_color, epsCenter)
        if show["center_A_boundary"]:
            render_boundary(renderables, 0, xCenter+lengths["overlay_offset"], y, A_color, epsCenter)
        if show["center_B_boundary"]:
            render_boundary(renderables, 1, xCenter-lengths["overlay_offset"], y, B_color, epsCenter)
        if show["center_A_matching"] and not show["center_B_matching"]:
            boxes += render_matching(renderables,0, xCenter, y, A_color, epsCenter)
        if show["center_B_matching"] and not show["center_A_matching"]:
            boxes += render_matching(renderables,1, xCenter, y, B_color, epsCenter)
        if show["center_A_matching"] and show["center_B_matching"]:
            boxes += render_xor_edges(renderables, xCenter, y, [A_color, B_color], epsCenter)
            if show["center_doubled_edges"]: 
                render_doubled_edges(renderables, xCenter, y, [A_color, B_color], epsCenter)
        boxes += render_center_buttons(10+xCenter, 10, renderables, font, epsCenter)
        renderables["epsCenter"] = epsCenter
        if eps_only: # hack to render eps programattically
            eps_callback({"eps":epsCenter})

    y = lengths["screen_height"] - lengths["button_height"]
    
    boxes += render_showhide_buttons(10,y, renderables, font)
    boxes += render_global_adjusters(400, y, renderables, font)
    y -= lengths["button_height"]
    boxes += render_os_buttons(10,y,filenames,renderables,font)

    if eps_only:
        bigfont = pygame.font.Font(None, 72)
        text = bigfont.render("Creating eps... please wait", True, black, grey)    
        bb = text.get_rect()
        bb.center = (lengths["screen_width"]/2, lengths["screen_height"]/2)
        screen.blit(text, bb)


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
def render_active_hex_centers(renderables, xoffset, yoffset, 
            which_side, eps):
    coords=renderables["coords"] 
    hexlist=renderables["hexagons"] 
    matching=renderables["matchings"][which_side]
    lengths=renderables["lengths"]

    adj = adjacency_map(matching)
    boxes = []
    for i in range(len(hexlist)):
        h = hexlist[i]
        if is_active(h, adj): 
            bb = render_hexagon_center(coords, h, xoffset, yoffset, lengths["hex_flipper_radius"], eps)
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
def randomize(matching, hexlist, steps):
    for trial in range(steps):
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
# Find the aspect ration of one of our pictures.  
def aspectratio(coords):
    xvalues = [v[0] for v in coords.values()]
    yvalues = [v[1] for v in coords.values()]
    xmin = min(xvalues) + 0.0
    ymin = min(yvalues) + 0.0
    xmax = max(xvalues) + 0.0
    ymax = max(yvalues) + 0.0
    return (xmax-xmin)/(ymax-ymin)

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


#====================================================================
# Save everything we could have changed in the current file. 
def save(basename, renderables):
    matchings = renderables["matchings"]
    show = renderables["show"]
    lengths = renderables["lengths"]
    write_edges(matchings[0], basename + "A.edge")
    write_edges(matchings[1], basename + "B.edge")
    showfile = open(basename + "show.pkl", "wb")
    pickle.dump(show, showfile)
    showfile.close()
    lengthsfile = open(basename + "lengths.pkl", "wb")
    pickle.dump(lengths, lengthsfile)
    lengthsfile.close()

#===================================================================
# Compute the proper scale and locations of pictures, based on the
# screen size and all the things we need to draw.
def compute_picture_sizes(renderables):
    lengths = renderables["lengths"]
    show = renderables["show"]
    picturecount = 0
    if show["A"]: picturecount += 1
    if show["B"]: picturecount += 1
    if show["Center"]: picturecount += 1

    screensize = screen.get_size()

    height = screensize[1] - 4*lengths["button_height"]
    width = screensize[0]

    window = width
    if picturecount > 0:
        window = int(screensize[0] / picturecount)
    aspect = aspectratio(renderables["dualcoords"])
    if int(height * aspect) < window:
        window = int(height * aspect)

    lengths["screen_width"] = screensize[0]
    lengths["screen_height"] = screensize[1]
    lengths["window"] = window

# Compute x coordinates at which to draw pictures
    lengths["xA"] = int((width - window * picturecount)/2)
    if(show["A"]):
        lengths["xCenter"] = lengths["xA"] + window
    else:
        lengths["xCenter"] = lengths["xA"]
    if(show["Center"]):
        lengths["xB"] = lengths["xCenter"] + window
    else:
        lengths["xB"] = lengths["xCenter"]

    renderables["coords"] = copy.deepcopy(renderables["unscaled_coords"])
    renderables["dualcoords"] = copy.deepcopy(renderables["unscaled_dualcoords"])
    rescale(renderables["coords"], renderables["dualcoords"], window)


#===================================================================
# Load a dimerpaint configuration.
def load(basename):
    if not os.path.isdir(basename):
        exit("Can't find "+basename)

    coords = read_vertices(basename + "full.vertex")
    unscaled_coords = copy.deepcopy(coords)
    background = read_edges(basename + "full.edge")
    hexagons = read_hexagons(basename + "full.hexagon")
    dualcoords = read_vertices(basename + "full.dualvertex")
    unscaled_dualcoords = copy.deepcopy(dualcoords)
    rhombi = read_rhombi(basename + "full.rhombus")

    matching_A = read_edges(basename + "A.edge")
    matching_B = read_edges(basename + "B.edge")
    matchings = [matching_A, matching_B]

    if os.path.isfile(basename + "show.pkl"):
        showfile = open(basename + "show.pkl", "rb")
        show = pickle.load(showfile)
        showfile.close()
    else:
        show = {
            "A": True,
            "B": True,
            "Center": True,
            "Highlight": True,

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
            "center_doubled_edges": True,
        }
    if os.path.isfile(basename + "lengths.pkl"):
        lengthsfile = open(basename + "lengths.pkl", "rb")
        lengths = pickle.load(lengthsfile)
        lengthsfile.close()
        if "old_screen_size" in lengths:
           del lengths["old_screen_size"]
    else:
        lengths = {}

    default_lengths = {
        "button_height": 20,
        "dimer_width":3,
        "hex_flipper_radius":4,
        "overlay_offset":0,
        "tile_edge_width":2,
        "shading_intensity":1,
        "randomize_steps":500,
        "y": 45,
        }

    for param in default_lengths.keys():
        if param not in lengths:
            lengths[param] = default_lengths[param]

    renderables = {"highlight":[{},{}], # highlighted edges on left and right
                   "background":background, 
                   "matchings":matchings, 
                   "hexagons":hexagons, 
                   "rhombi": rhombi,
                   "coords": coords,
                   "unscaled_coords": coords,
                   "dualcoords":dualcoords,
                   "unscaled_dualcoords":dualcoords,
                   "show":show,
                   "lengths":lengths}
    compute_picture_sizes(renderables)
    return renderables

#===================================================================
# Highlighting code
def highlight_edge(renderables, edge, side):
    highlight = renderables["highlight"]
    if edge in highlight[side]:
        print "edge highlight off: side ", side, edge
        del highlight[side][edge]
    else:
        print "edge highlight on: side ", side, edge
        highlight[side][edge] = 1

def highlight_hexagon(renderables, hexagon, side):
    highlight = renderables["highlight"]
    edges = []
    for i in range(6):
        j = (i+1) % 6
        edges.append(frozenset([hexagon[i], hexagon[j]]))
    all_highlighted = True

    for edge in edges:
        if not(edge in highlight[side]):
            all_highlighted = False

    if(all_highlighted):
        for edge in edges:
            del highlight[side][edge]
    else:
        for edge in edges:
            highlight[side][edge] = 1 

def highlight_path(renderables, m0, m1, unordered_edge):
    matchings = renderables["matchings"]
    edge = [endpoint for endpoint in unordered_edge]
    adj0 = adjacency_map(matchings[m0])
    adj1 = adjacency_map(matchings[m1])
    path = find_path(adj0, adj1, edge[0])

    #list of edges corresponding to the path
    loop = []
    all_highlighted = True

    for i in range(len(path) - 1):
        e = frozenset([path[i], path[i+1]])
        loop.append(e)
        if e not in renderables["highlight"][m0]: 
            all_highlighted = False
        if e not in renderables["highlight"][m1]:
            all_highlighted = False

    if all_highlighted:
        for e in loop:
            del renderables["highlight"][m0][e]
            del renderables["highlight"][m1][e]
    else:
        for e in loop:
            renderables["highlight"][m0][e] = 1
            renderables["highlight"][m1][e] = 1
 

#===================================================================
# Main program


pygame.init()
pygame.font.init()
font = pygame.font.Font(None, 18)
screen=pygame.display.set_mode((1364,690), pygame.RESIZABLE)

pygame.mouse.set_cursor(*pygame.cursors.broken_x)


try:
    data_directory = sys.argv[1]
except IndexError:
    data_directory = "dimerpaint"

try:
    input_file = sys.argv[2]
except IndexError:
    possible_start_files = os.listdir(data_directory)
    input_file = possible_start_files[0]

eps_only = False
try:
    switch = sys.argv[3]
    if switch == "eps":
        eps_only = True
except IndexError:
    pass

filenames = {   "data_directory":data_directory,
                "input_file":input_file,
                "basename":data_directory + "/" + input_file + "/" } 

renderables = load(filenames["basename"])

bounding_box_data = render_everything(renderables,filenames,font)
bounding_boxes = [record[3] for record in bounding_box_data]

done=False     #Loop until the user clicks the close button....
if eps_only:   # unless we're just rendering the eps files, which has already been done.
    done=True

clock=pygame.time.Clock() # Used to manage how fast the screen updates




#=================================================
# Pygame event loop
while done==False:
    for event in pygame.event.get(): # User did something
        if event.type == pygame.QUIT: # If user clicked close
            done=True # Flag that we are done so we exit this loop
        if event.type == pygame.VIDEORESIZE:
            print "user resized screen ", event.size
            screen=pygame.display.set_mode(event.size, pygame.RESIZABLE)
            compute_picture_sizes(renderables)
            bounding_box_data = render_everything(renderables, filenames,font)
            bounding_boxes = [record[3] for record in bounding_box_data]
        if event.type == pygame.MOUSEBUTTONUP:
            radius = 1
            ul = (event.pos[0] - radius, event.pos[1] - radius)
            clickpoint = pygame.Rect(ul , (2*radius,2*radius))
            box_index = clickpoint.collidelist(bounding_boxes)
            if(box_index != -1):
                objtype = bounding_box_data[box_index][0]
                matchings = renderables["matchings"]
                hexagons = renderables["hexagons"]
                if objtype == "edge":
                    (objtype, edge, side, box) = bounding_box_data[box_index]
                    if(event.button == 3):
                        highlight_edge(renderables, edge, side)
                    else:
                        if edge in matchings[side]:
                            print "deleting"
                            del matchings[side][edge]
                        else:
                            print "adding"
                            matchings[side][edge] = 1
                elif objtype == "matchededge":
                    (objtype, matchededge, side, box) = bounding_box_data[box_index]
                    if event.button == 3:
                        highlight_path(renderables, side, 1-side, matchededge)
                    else:
                        flip_path(matchings, side, 1-side, matchededge)
                elif objtype == "button":
                    (objtype, args, callback, bb) = bounding_box_data[box_index]
                    callback(args)
                elif objtype == "hexagon":
                    (objtype, hexagon, side, box) = bounding_box_data[box_index]
                    if event.button == 3:
                        highlight_hexagon(renderables, hexagons[hexagon], side)
                    else:
                        flip_hex(matchings[side], hexagons, hexagon)
                else:
                    print "uh why am I here?"

                bounding_box_data = render_everything(renderables, filenames, font)
                bounding_boxes = [record[3] for record in bounding_box_data]

    # Limit to 20 frames per second
    clock.tick(20)
 
pygame.font.quit()
pygame.quit()

if not eps_only:
    save(filenames["basename"], renderables)


