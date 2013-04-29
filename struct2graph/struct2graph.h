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
#include <chrono>
#include <random>

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
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/iteration_macros.hpp>

// get property with Graph[Vertex].bipartite_color = int;
struct vertex_property {
	int bipartite_color;
	int color;
	int search_color;
};

struct edge_property {
	int ear;
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
typedef Graph::edge_descriptor Edge;
typedef Graph::vertex_descriptor Vertex;

// initialise boost command line option parser
boost::program_options::variables_map init_options(int ac, char* av[]);

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

// does the whole graph decomposition, uses methods below
void decompose_graph(Graph& graph, std::ostream* out);

// do a Breadth First Search to test for bipartite property
bool is_bipartite_graph(Graph& g, Vertex startVertex, Edge& ed);

// get a vector of all vertices with their component id. finds connected components with DFS
void connected_components_to_subgraphs(Graph& g);

// finds connected components with DFS
void biconnected_components_to_subgraphs(Graph& g);

// typedefs for ear decomposition
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

// do a ear_decomposition1 using different spanning trees
void do_ear_decompositions (Graph& g, Vertex startVertex);

// ear decomposition of blocks
void ear_decomposition1(Graph& g, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex& start);

// get spanning tree with DFS
void get_spanning_tree(Graph& g, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex& start);

// change spanning tree with sampling a random new one
void change_spanning_tree(Graph& g, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex& root);

// given all the parents of a spanning_tree, find the lca a crossedge
std::pair<Vertex, int> get_lca_distance(Graph& g, std::map<Vertex, Vertex>& parents, Edge e, Vertex r);

// do a walk in the spanning tree starting at v and ending at the root r -> return a vector with the walk
std::vector<Vertex> make_tree_walk(std::map<Vertex, Vertex>& parents, Vertex v, Vertex r);

// ear decomposition of blocks
void ear_decomposition(Graph& g, Vertex startVertex);

// actual dfs for ear decomposition
void ear_dfs(Graph& g, Vertex v, ear_propertymap_t& p, ear_t& ear, unsigned int& counter);



#endif
