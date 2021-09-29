# Cellular Free Space Routing Package
This python package intends to provide certain functionalitites for constrained free space routing in a simple, clean and efficient implementation. 

You can compute alternative paths up to homotopy (relaxed to a much simpler relation which under some constraints equals homotopy) in 2D indoor maps.

# Background
Constrained Free Space in computational geometry is usually organized as a graph of obstacles, but for navigation purposes beyond pure geometry, this is not a good ideal as many properties of relevance (speed, surface parameters, etc.) for the navigation space itself. Therefore, since quite some time we follow a data model in which a bitmap (occupied vs. non-occupied) is being used to represent constrained 2D space.

For higher-dimensional space, we could use voxels, but we feel that this is easily becoming inefficient, therefore, we take larger navigation mesh elements of arbitrary compact shape (maybe 3D boxes represented with two corner coordinates) and interconnect these navigation cells.

For this package, however, our cells are the pixels of a raster and the interconnection is given by an eight neighbor system in which each pixel is connected to all eight adjacent pixels for navigation purposes unless they represent non-walkable occupied space.

This package is given as a C++ extension to Python in which we have the following functionalities that one can easily extend to special use-cases or different cellular systems:

- Given a raster (that is uint8 numpy array), derive the relevant graph
- Given two points, compute the shortest path between both of them using Dijkstra, A*, or ALT (if landmarks are available) 
- Given a set of points and a threshold, compute landmarks at each of the points
- Given two trajectories, compute whether they are homotopic by estimating whether the (oriented) polygon spanned by them is empty
- Given two points, run the Penalty algorithm for a certain number of times possibly identifying more than one homotopy group and represent each with a simplified version (constrained DP with epsilon not changing the HG)

# Usage example

# A Counter-Example as a Warning

# Implementation Details
## Graph Representation
We are using the Boost Graph Library adjacency list representation for our graph using vectors for both vertex list and adjacency list which should be quite optimal given the fact that we don't change the graph itself a lot.

For each pixel of the input map, we check all eight neighbors and if they exist and are walkable, we add an edge to the graph taking a unit length (1) or $\sqrt{2}$ as the length.

## Shortest Path
For the first version, we just compute the shortest path using Dijkstra or A^* in the graph and return it as a sequence of pixel cells (without simplification)

## Douglas Peucker
As we need Douglas Peucker to remove the oscillating behavior in a post-processing step, we export it as well. DP can be run with a map in which case we simplify within the homotopy group or without a map in which unconstrained simplification is performed.

## Penalty Algorithm
A simple algorithm to find alternative paths is to iteratively compute shortest paths and penalize the shortest path (maybe extend each edge by multiplying with 1.1). This leads to "shortest" paths sweeping across the plane preferring non-penalized edges finally jumping into new homotopy groups. However, the downside is that the amount of penalty is high near the search start which we remedy with homotopy-constrained Douglas Peucker.

Other algorithms to find alternative paths do exist, but are not part of this collection for now.

The penalty algorithm modifies the weightmap and the modified weightmap is an interesting asset in itself as it will have a very high weight for edges crossing narrow passages. This can be exported as a bitmap again using the mean (or max) of all out-edges of each vertex.

## Approximate Homotopy Predicate
We compute whether two non-selfintersecting shortest paths are homotopic approximately by checking whether the scanline-filled polygon spanned by them is empty. Note that even for 3D alternative paths, a different definition would have to be used.

## Approximate Shortest Homotopic Path

We also provide a randomized algorithm for simplification within the same homotopy class essentially by pruning random triangles if they are empty. This is not very fast (exact algorithms for shortest path in homotopy exist and are not very complicated, yet expect a triangulated space representation), but sufficient to give a good guess on the distance to be traveled.
