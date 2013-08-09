/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 25.03.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
* Compile with: g++ -std=c++11 -g -lboost_program_options -o struct2graph struct2graph.cc
*/

#ifndef DECOMPOSE_H
#define DECOMPOSE_H

// include common header with graph definition and global variables
#include "common.h"

// include standard library parts

// include boost components

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


// does the whole graph decomposition, uses methods below
void decompose_graph(Graph& graph, std::ostream* out);

// get a vector of all vertices with their component id. finds connected components with DFS
void connected_components_to_subgraphs(Graph& g);

// finds connected components with DFS
void biconnected_components_to_subgraphs(Graph& g);

// get max degree of a graph
int get_max_degree(Graph& g);

// do a Breadth First Search to test for bipartite property
bool is_bipartite_graph(Graph& g, Vertex startVertex, Edge& ed);

// do a schieber ear decomposition without any statistics
void schieber_ear_decomposition (Graph& g);

// implementation of the schieber ear decomposition
void ear_decomposition (Graph& g, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex start);

// get boost random spanning tree for schieber ear decomposition
void get_random_spanning_tree (Graph& g, std::mt19937& r, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex start);

// given all the parents of a spanning_tree, find the lca a crossedge (helper for schieber ear decomposition)
std::pair<Vertex, int> get_lca_distance(Graph& g, std::map<Vertex, Vertex>& parents, Edge e, Vertex r);

// do a walk in the spanning tree starting at v and ending at the root r -> return a vector with the walk  (helper for schieber ear decomposition)
std::vector<Vertex> make_tree_walk(std::map<Vertex, Vertex>& parents, Vertex v, Vertex r);

// ear decomposition of blocks
void ramachandran_ear_decomposition(Graph& g);

// actual dfs for ramachandran ear decomposition (dfs algroithm)
void ear_dfs(Graph& g, Vertex v, ear_propertymap_t& p, ear_t& ear, unsigned int& counter);


#endif
