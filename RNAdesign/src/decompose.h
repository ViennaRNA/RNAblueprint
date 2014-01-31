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

// typedefs for ramachandran ear decomposition
typedef std::pair<Vertex, Vertex> edge_t;
typedef std::map<edge_t, edge_t> ear_t;
// struct to remember coloring, time, parents of a vertex
struct property {
	int color;
	int preorder;
	Vertex parent;
	Vertex low;
	edge_t ear;
};
typedef std::map<Vertex, property> ear_propertymap_t;


// lexmin implementation for pairs of any type
template <typename T>
std::pair<T, T>& lexmin(std::pair<T, T>& a, std::pair<T, T>& b) {
	if (a.first == b.first) {
		return !(b.second<a.second)?a:b;
	} else {
		return !(b.first<a.first)?a:b;
	}
}

// does the graph decomposition, and calls the coloring of the subgraphs
void decompose_graph(Graph& graph, std::ostream* out, int num_trees, bool ramachandran, bool no_bipartite_check);

// get a vector of all vertices with their component id. finds connected components with DFS
void connected_components_to_subgraphs(Graph& g);

// finds biconnected components with DFS
void biconnected_components_to_subgraphs(Graph& g);

// starts at a degree>3 articulation point and walks along a path to connect it to one component
void merge_biconnected_paths(Graph& g, Vertex p, Vertex v, boost::property_map < Graph, boost::edge_component_t >::type& component, std::vector<Vertex>& art_points, int& nc);

// do a schieber ear decomposition without any statistics
void schieber_ear_decomposition (Graph& g);

// implementation of the schieber ear decomposition
void ear_decomposition (Graph& g, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex start);

// get boost random spanning tree for schieber ear decomposition
void get_random_spanning_tree (Graph& g, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex start);

// given all the parents of a spanning_tree, find the lca of a crossedge (helper for schieber ear decomposition)
Vertex get_lca(Graph& g, std::map<Vertex, Vertex>& parents, Edge& e, Vertex& r);

// recursion function for the get_lca_distance function
Vertex lca_recursion(Graph& g, std::map<Vertex, Vertex>& parents, Edge& e, Vertex u, Vertex v, Vertex& r);

// Get the distance from one vertex to another (mostly root) in the spanning tree
int get_distance(Graph& g, std::map<Vertex, Vertex>& parents, Vertex v, Vertex r);

// return the distance of a edge to the root (smallest distance of vertices)
int get_distance(Graph& g, std::map<Vertex, Vertex>& parents, Edge e, Vertex r);

// return a iterator poiting to the parent of a certain vertex in the spannin gtree
Vertex get_parent(std::map<Vertex, Vertex>& parents, Vertex v);

// identify articulation points and add them as graph property Ak
void color_articulation_points (Graph& g);

// in an ear graph just get the degree of edges participating in this ear
int degree_in_ear (Vertex& v, Graph& g, int k);

// create subgraphs for paths between articulation points (for easier coloring afterwards)
void parts_between_articulation_points_to_subgraphs (Graph& g, int k);

// ear decomposition of blocks
void ramachandran_ear_decomposition(Graph& g);

// actual dfs for ramachandran ear decomposition (dfs algroithm)
void ear_dfs(Graph& g, Vertex v, ear_propertymap_t& p, ear_t& ear, unsigned int& counter);


#endif
