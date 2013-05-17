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
#include <regex>

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

// initialise boost command line option parser
boost::program_options::variables_map init_options(int ac, char* av[]);
// read the input file into a string
void read_input(std::istream* in, std::vector<std::string>& input);
// parse the input string into ??? TODO
void parse_structures(std::vector<std::string> structures);




#endif

