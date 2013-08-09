/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 25.03.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
* Compile with: g++ -std=c++11 -g -lboost_program_options -o struct2graph struct2graph.cc
*/

#ifndef PRINTGRAPH_H
#define PRINTGRAPH_H

// include common header with graph definition and global variables
#include "common.h"

// include standard library parts
#include <iostream>

// include boost components
#include <boost/graph/graphml.hpp>

// print the graph as a gml file to the output
void print_graph(Graph& g, std::ostream* out, std::string outfile, std::string nametag);

// print all the subgraphs as GML (iterator over subgraphs)
void print_subgraphs(Graph& g, std::ostream* out, std::string nametag);

#endif
