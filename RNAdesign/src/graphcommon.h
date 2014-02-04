/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 26.08.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

#ifndef GRAPHCOMMON_H
#define GRAPHCOMMON_H

// include common header with graph definition and global variables
#include "common.h"

// include standard library parts
#include <limits>

// include boost components
#include <boost/graph/breadth_first_search.hpp>

// get max degree of a graph
std::pair <int, int> get_min_max_degree(Graph& g);

// do a Breadth First Search to test for bipartite property
bool is_bipartite_graph(Graph& g, Vertex startVertex, Edge& ed);

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

void open_ear_decomposition (Graph& g, Vertex startVertex, ear_t& ear);

void ear_dfs(Graph& g, Vertex v, ear_propertymap_t& p, ear_t& ear, unsigned int& counter);

#endif
