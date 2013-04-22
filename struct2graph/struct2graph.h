/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 25.03.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
* Compile with: g++ -std=c++11 -g -lboost_program_options -o struct2graph struct2graph.cc
*/

#ifndef STRUCT2GRAPH_H
#define STRUCT2GRAPH_H



// include standard library parts
#include <string>
#include <cstring>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <array>
#include <limits>

// include boost components
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>

// get property with Graph[Graph::vertex_descriptor].bipartite_color = int;
struct vertex_property {
	int bipartite_color;
	int color;
	int search_color;
};

struct edge_property {
	// put something here
};

namespace boost {
	struct edge_component_t {
		enum
		{ num = 555 };
		typedef edge_property_tag kind;
	}
	edge_component;
}

//, graph_properties 
// boost graph template
typedef boost::adjacency_list_traits< boost::vecS, boost::vecS, boost::undirectedS > Traits;
typedef boost::subgraph< boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, 
	boost::property< boost::vertex_color_t, int, vertex_property >,
	boost::property< boost::edge_index_t, int, boost::property < boost::edge_component_t, std::size_t, edge_property> > > > Graph;

// read the input file into a string
std::vector<std::string> read_input(std::istream* in);

// parse the input string into a graph
Graph parse_graph(std::vector<std::string> structures);

// get max degree of a graph
int get_max_degree(Graph& g);

// print the graph as a gml file to the output
void print_graph(Graph& g, std::ostream* out, std::string nametag);

// print all the subgraphs as GML (iterator over subgraphs)
void print_subgraphs(Graph& g, std::ostream* out, std::string nametag);

// do a Breadth First Search to test for bipartite property
bool is_bipartite_graph(Graph& g, Graph::vertex_descriptor startVertex, Graph::edge_descriptor& ed);

// get a vector of all vertices with their component id. finds connected components with DFS
void connected_components_to_subgraphs(Graph& g);

// finds connected components with DFS
void biconnected_components_to_subgraphs(Graph& g);

// typedefs for ear decomposition
typedef std::pair<Graph::vertex_descriptor, Graph::vertex_descriptor> edge_t;
typedef std::map<edge_t, edge_t> ear_t;
// struct to remember coloring, time, parents of a vertex
struct property {
	int color;
	int preorder;
	Graph::vertex_descriptor parent;
	Graph::vertex_descriptor low;
	edge_t ear;
};
typedef std::map<Graph::vertex_descriptor, property> ear_propertymap_t;

// ear decomposition of blocks
void ear_decomposition(Graph& g, Graph::vertex_descriptor startVertex);

// actual dfs for ear decomposition
void ear_dfs(Graph& g, Graph::vertex_descriptor v, ear_propertymap_t& p, ear_t& ear, unsigned int& counter);



#endif
