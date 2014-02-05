/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 25.03.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

#ifndef DECOMPOSE_H
#define DECOMPOSE_H

// include common header with graph definition and global variables
#include "common.h"

// include standard library parts
#include <limits>

// include boost components
#include <boost/graph/connected_components.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/random_spanning_tree.hpp>

// does the graph decomposition, and calls the coloring of the subgraphs
void decompose_graph(Graph& graph, std::ostream* out, int num_trees, bool no_bipartite_check);

// get a vector of all vertices with their component id. finds connected components with DFS
void connected_components_to_subgraphs(Graph& g);

// finds biconnected components with DFS
void biconnected_components_to_subgraphs(Graph& g);

// starts at a degree>3 articulation point and walks along a path to connect it to one component
void merge_biconnected_paths(Graph& g, Vertex p, Vertex v, boost::property_map < Graph, boost::edge_component_t >::type& component, std::vector<Vertex>& art_points, int& nc);

// ear decomposition of blocks
void ramachandran_ear_decomposition(Graph& g);

// identify articulation points and add them as graph property Ak
void color_articulation_points (Graph& g);

// in an ear graph just get the degree of edges participating in this ear
int degree_in_ear (Vertex& v, Graph& g, int k);

// create subgraphs for paths between articulation points (for easier coloring afterwards)
void parts_between_articulation_points_to_subgraphs (Graph& g, int k);
// recursion for parts function
void parts_recursion (Graph& g, Vertex& start, int& k, Graph* subgptr);

#endif
