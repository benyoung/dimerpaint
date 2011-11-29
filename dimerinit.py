#!/usr/bin/python
import sys
import argparse

parser = argparse.ArgumentParser(description="Make input files for mklattice / dimerpaint")

if len(sys.argv) != 5:
    sys.exit("usage: dimerinit.py [dryrun|filename] a b c")

[a,b,c] = (int(u) for u in sys.argv[2:5])


print "Generating %d x %d x %d" % (a,b,c)

#==========================================
# Strings for drawing grids of various sorts
hexagon_grid = [
    "  X---X  ",
    " /     \ ",
    "X       X",
    " \     / ",
]

hexagon_joiner = [
        "   ",
        "   ",
        "---",
        "   "
]

noedge_grid = [
    "  X   X  ",
    "         ",
    "X       X",
    "         ",
]

noedge_joiner = [
        "   ",
        "   ",
        "   ",
        "   "
]


# Return a list of strings which, when printed with newlines,
# makes a hexagon-like grid.  
def hexgrid(num_rows, num_full_cols, hexagon, joiner):
    output = []
    for row in range(num_rows):
        for line in range(4):
            output.append(joiner[line].join([hexagon[line]]*num_full_cols))
    output.append( joiner[0].join([hexagon[0]]*num_full_cols))
    return output


 
# Overwrite a portion of a given row of a list of strings.
def rewrite(grid,p,string):
    (row,col)=p
    colB = col + len(string)
    grid[row] = grid[row][:col] + string + grid[row][colB:]
    
# Overwrite strings on a grid in an affine lattice 
def rewrite_lattice(grid, origin, vec_x, vec_y, size, string):
    for i in range(size[0]):
        for j in range(size[1]):
            row = origin[0] + vec_x[0] * i + vec_y[0] * j
            col = origin[1] + vec_x[1] * i + vec_y[1] * j
            rewrite(grid, (row, col), string)

# Mark the centers of hexagons on a grid.
def hexagon_centers(grid, size):
    rewrite_lattice(grid, (2, 4), (4,0), (0,12), size, "H") 
    rewrite_lattice(grid, (4, 10), (4,0), (0,12), (size[0]-1, size[1]-1), "H") 

# Put an a by b by c dimer box upper-left justified in the grid
def box(grid, a,b,c):

    row = 4*c - 2*(b%2) + 2 # trial and error for this, to be honest
    col = 6*b + 4

    rewrite(grid, (row,col), "*") # debug
    rewrite_lattice(grid, (row+2,col-1), (2,6),  (2,-6),  (a,b), "---")
    rewrite_lattice(grid, (row-1,col-3), (2,-6), (-4,0),  (b,c), "/")
    rewrite_lattice(grid, (row-1,col+3), (-4,0), (2,6),   (c,a), "\\")


gridrows = int((2*a + 2*b + 4*c) / 4) + 1  + (a%2)
gridcols = int((6*a + 6*b) / 12) + 1 + (a + b)%2

grid = hexgrid(gridrows ,gridcols, noedge_grid, noedge_joiner)
#grid = hexgrid(a,b, hexagon_grid, hexagon_joiner)
#hexagon_centers(grid, (a, b))

box(grid, a,b,c)

print "\n".join(grid)




