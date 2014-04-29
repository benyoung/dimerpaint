#!/usr/bin/python

import pickle
import copy
import sys
import math
import pygame
import numpy
import os
import shutil
from collections import defaultdict

if len(sys.argv<4):
    raise Exception("Usage: symcount.py order midline filename")

#===================================================================
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

#===================================================================
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

#===============================================================
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


#===================================================================
# Load a dimerpaint configuration.
# Lifted from Dimerpaint.  Needs modification before using.
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
    # this allows us to add new keys
    show_default_dict = defaultdict(bool)
    show_default_dict.update(show)

    renderables = {"highlight":[{},{}], # highlighted edges on left and right
                   "background":background, 
                   "matchings":matchings, 
                   "hexagons":hexagons, 
                   "rhombi": rhombi,
                   "coords": coords,
                   "unscaled_coords": coords,
                   "dualcoords":dualcoords,
                   "unscaled_dualcoords":dualcoords,
                   "show":show_default_dict,
                   "lengths":lengths}
#    compute_picture_sizes(renderables)
    return renderables

#==============================================================
# make an adjacency map from a matching. 
# Lifted from Dimerpaint.  Needs modification before using.
def adjacency_map(M): 
    adj = {}
    for edge in M:
        endpoints = [endpt for endpt in edge]
        adj[endpoints[0]] = endpoints[1]
        adj[endpoints[1]] = endpoints[0]
    return adj

#==============================================================
# Find a path in the superposition of two matchings, starting at
# a node and stopping at another node.  Matchings are represented as 
# adjacency maps.
#
# Return the INDEX into the nodes array.

def find_other_end(start, adj1, adj2, nodes):
    if start in adj1:
        adj = (adj1, adj2)
    else:
        adj = (adj2, adj1)

    parity = 1
    p = adj[0][start]
    while not(p in nodes):
        p = adj[parity][p]
        parity = (parity + 1)%2
    return nodes.index(p)

def find_perm(nodes, adj1, adj2):
    return tuple([find_other_end(i, adj1, adj2, nodes) for i in nodes])
    

#==============================================================
# Unlike dimerpaint, now we want to flip 
# unidirectionally - that is, so as only to increase the number of vertices.
# To do this we make use of an undocumented feature: the vertices of
# the hexagons are always listed in the same order:
#
#                            1---2                         
#                           /     \                         
#                          0       3                        
#                           \     /                         
#                            5---4                          
#                                                           
# Here, the "coordinates" of the vertices are (row, column) in a text file,
# literally.  So just like the above.  For example, one sample hexagon out of
# a typical data file is
# 
#      ((34, 18), (32, 20), (32, 24), (34, 26), (36, 24), (36, 20))
#
# so that one has its leftmost vertex in row 34, column 18
# So we really need adj(0)=1, adj(2)=3, adj(4) = 5. 

def is_active(hexagon, adj):
    if hexagon[0] in adj and hexagon[2] in adj and hexagon[4] in adj:
        if adj[hexagon[0]] == hexagon[1]: 
            if adj[hexagon[2]] == hexagon[3]:
                if adj[hexagon[4]] == hexagon[5]:
                    return True
    return False

def is_left_active(hexagon, adj):
    if hexagon[0] in adj and hexagon[2] in adj and hexagon[4] in adj:
        if adj[hexagon[0]] == hexagon[1]: 
            if adj[hexagon[2]] == hexagon[3]:
                if adj[hexagon[4]] == hexagon[5]:
                    if int(hexagon[0][1]) < midline - 6:
                        return True
    return False

def is_mid_active(hexagon, adj):
    if hexagon[0] in adj and hexagon[2] in adj and hexagon[4] in adj:
        if adj[hexagon[0]] == hexagon[1]: 
            if adj[hexagon[2]] == hexagon[3]:
                if adj[hexagon[4]] == hexagon[5]:
                    if hexagon == dual_hex(hexagon):
                        return True
    return False

def all_left_active(hexagons, adj):
    return [H for H in hexagons if is_left_active(H, adj)]

def all_mid_active(hexagons, adj):
    return [H for H in hexagons if is_mid_active(H, adj)]

def all_active(hexagons, adj):
    return [H for H in hexagons if is_active(H, adj)]

#============================================================================
# Flip one hexagon.
def flip_hex(matching, h):
    edges = [frozenset([h[i], h[(i+1)%6]]) for i in range(6)]
    for e in edges:
        if(e in matching): 
            matching.remove(e)
        else:
            matching.add(e)
    

#============================================================================
# Finds the reflection of a hexagon across col = midline
def dual_hex(h):
    k = []
    k = [(a,2*midline-b) for (a,b) in h]
    k = tuple([k[(3-i)%6] for i in range(6)])
    return k
 

#====================================================================
# Flips symmetric pairs of hexagons, reflected across col = midline
def flip_2hex(matching, h):
    edges = [frozenset([h[i], h[(i+1)%6]]) for i in range(6)]
    k = dual_hex(h)
    edges2 = [frozenset([k[i], k[(i+1)%6]]) for i in range(6)]
    for e in edges:
        if(e in matching):
            matching.remove(e)
        else:
            matching.add(e)
    for e in edges2:
        if(e in matching):
            matching.remove(e)
        else:
            matching.add(e)

#====================================================================
# Find all the nodes in the sense of Kenyon-Wilson - that is, the 
# vertices in the double dimer model (A,B) which are degree 1.
def find_nodes(A,B):
    nodes = defaultdict(int)
    for e in A:
        for v in e:
            nodes[v] += 1
    for e in B:
        for v in e:
            nodes[v] += 1
    return [v for v in nodes if nodes[v]==1]


renderables = load(sys.argv[3] + "/")
hexagons = renderables["hexagons"]
vertices = renderables["coords"]
A = renderables["matchings"][0]
B = renderables["matchings"][1]
map_A = adjacency_map(A)
map_B = adjacency_map(B)

nodes = find_nodes(A,B)
reference_permutation = find_perm(nodes, map_A, map_B)

# OK.  this is getting to be a mess.  
# Matchings have to be represented as frozensets in order not to double-count
# anything.  But when we alter the matchings, it's convenient to convert them 
# back to sets of edges.
#
# Sadly, dimerpaint represents matchings as neither of these, but rather as
# dictionaries keyed by edges.  Hence the following ludicrous statement: 


# Ok, so now, sadly, we cannot rely on "old" and "new" because we have to keep 
# two lines back, so we will use lists, so we can recall back one or two levels,
# and we'll just add a 2 new "available" hexagon defs which either add two 
# symmetrically or on the axis of symmetry (x = midline)

start_record = (frozenset(A.keys()), frozenset(B.keys()))

#start adding to A

s = []
s = s + [set([start_record])]
series_len = int(sys.argv[1])
midline = sys.argv[2]
for i in range(series_len):
    print "counting %d boxes in part A" % i
    new = []
    for (A,B) in s[i]:
        map_A = adjacency_map(A)
        map_B = adjacency_map(B)
        for H in all_mid_active(hexagons, map_A):
            new_A = set(A)
            flip_hex(new_A, H)
            new_permutation = find_perm(nodes, adjacency_map(new_A), map_B)
            if new_permutation == reference_permutation:
                new.append((frozenset(new_A),B))
    if i > 0:
        for (A,B) in s[i-1]:
            map_A = adjacency_map(A)
            map_B = adjacency_map(B)
            for H in all_left_active(hexagons, map_A):
                new_A = set(A)
                flip_2hex(new_A, H)
                new_permutation = find_perm(nodes, adjacency_map(new_A), map_B)
                if new_permutation == reference_permutation:
                    new.append((frozenset(new_A),B))
    s = s+[set(new)]
              
#take the list from part A and modify to get the full list

t = [] + [set(s[0])]
for i in range(series_len):
    print "counting %d boxes in part B" % i
    new = []
    for (A,B) in t[i]:
        map_A = adjacency_map(A)
        map_B = adjacency_map(B)
        for H in all_mid_active(hexagons, map_B):
            new_B = set(B)
            flip_hex(new_B, H)
            new_permutation = find_perm(nodes, map_A, adjacency_map(new_B))
            if new_permutation == reference_permutation:
                new.append((A,frozenset(new_B)))
    if i > 0:
        for (A,B) in t[i-1]:
            map_A = adjacency_map(A)
            map_B = adjacency_map(B)
            for H in all_left_active(hexagons, map_B):
                new_B = set(B)
                flip_2hex(new_B, H)
                new_permutation = find_perm(nodes, map_A, adjacency_map(new_B))
                if new_permutation == reference_permutation:
                    new.append((A,frozenset(new_B)))

    t = t + [set(new) | s[i+1]]   
counts = [len(t[i]) for i in range(series_len+1)]
print counts             
 

#counts = []
#old = set([start_record])   
#series_len = int(sys.argv[1])
#for i in range(series_len):
#    print "counting %d boxes" % i
#    counts.append(len(old))
#    new = [] 
#    for (A,B) in old:
#        map_A = adjacency_map(A)
#        map_B = adjacency_map(B)
#        for H in all_active(hexagons, map_A):
#            new_A = set(A)
#            flip_hex(new_A, H)            
#            new_permutation = find_perm(nodes, adjacency_map(new_A), map_B) 
#            if new_permutation == reference_permutation:
#                new.append((frozenset(new_A), B))
#        for H in all_active(hexagons, map_B):
#            new_B = set(B)
#            flip_hex(new_B, H)            
#            new_permutation = find_perm(nodes, map_A, adjacency_map(new_B)) 
#            if new_permutation == reference_permutation:
#                new.append((A, frozenset(new_B)))
#    old = set(new)
#counts.append(len(old))
#print counts
