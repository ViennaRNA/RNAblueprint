/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 25.03.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
* Compile with: g++ -std=c++11 -g -lboost_program_options -o struct2graph struct2graph.cc
*/

#ifndef COMMON_H
#define COMMON_H

// include standard library parts
#include <string>
#include <vector>
#include <array>

// include boost components
#include <boost/program_options.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/adjacency_list.hpp>

//Global variables
extern bool verbose;
extern unsigned long seed;

// get property with Graph[Vertex].bipartite_color = int;
struct vertex_property {
	int bipartite_color;
	int color;
};

struct edge_property {
	int ear;
	int color;
};

namespace boost {
	struct edge_component_t {
		enum
		{ num = 555 };
		typedef edge_property_tag kind;
	};
}

//, graph_properties 
// boost graph template
typedef boost::adjacency_list_traits< boost::vecS, boost::vecS, boost::undirectedS > Traits;
typedef boost::subgraph< boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, 
	boost::property< boost::vertex_color_t, int, vertex_property >,
	boost::property< boost::edge_index_t, int, boost::property < boost::edge_component_t, std::size_t, edge_property> > > > Graph;
typedef Graph::edge_descriptor Edge;
typedef Graph::vertex_descriptor Vertex;

#endif
