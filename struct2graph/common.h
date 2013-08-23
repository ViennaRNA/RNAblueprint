/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check.
*
* Created on: 25.03.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

#ifndef COMMON_H
#define COMMON_H

// include standard library parts
#include <string>
#include <vector>
#include <array>

// include boost components
#include <boost/lexical_cast.hpp>
#include <boost/config.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

//Global variables
extern bool verbose;

// random number generator
extern std::mt19937 rand_gen;

// Encode Bases to enums
enum bases {A, U, G, C, X};

// get property with Graph[Vertex].bipartite_color = int;
struct vertex_property {
	int bipartite_color;
	int color;
	char base;
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

// graph_properties 
// boost graph template
typedef boost::adjacency_list_traits< boost::vecS, boost::vecS, boost::undirectedS > Traits;
typedef boost::subgraph< boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, 
	boost::property< boost::vertex_color_t, int, vertex_property >,
	boost::property< boost::edge_index_t, int, boost::property < boost::edge_component_t, std::size_t, edge_property> > > > Graph;
typedef Graph::edge_descriptor Edge;
typedef Graph::vertex_descriptor Vertex;

// overload << operator to print vectors with any content
template <typename T>
std::ostream& operator<< (std::ostream& os, std::vector<T>& vec) {
	int i = 0;
	for (auto elem : vec) {
		os << "(" << i++ << ") " << elem << std::endl;
	}
	return os;
}

// overload << operator to print maps with any content
template <typename U, typename V>
std::ostream& operator<< (std::ostream& os, std::map<U, V>& m) {
	for (typename std::map<U, V>::iterator it = m.begin(); it != m.end(); it++) {
        	os << it->first << "," << it->second << std::endl;
	}
	return os;
}

// matrix template
template < class T, size_t ROWS, size_t COLS > using matrix = std::array< std::array< T, COLS >, ROWS >;

// function to make enum a char again and other way round
char enum_to_char(int intletter);
int char_to_enum(char charletter);

#endif
