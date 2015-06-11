/* This program reads secundary RNA structures in dot-bracket and
 * builds a graph for a latter ear-decomposition and bipartitness-check
 *
 * Created on: 25.03.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

#ifndef PRINTGRAPH_H
#define PRINTGRAPH_H

// include common header with graph definition and global variables
#include "common.h"
#include "graphcommon.h"

// include standard library parts
#include <iostream>
#include <sstream>

// include boost components
#include <boost/graph/graphml.hpp>

namespace design {
  namespace detail {
    // global variables
    extern std::string outfile;

    // print the graph as a gml file to the output
    void print_graph (Graph& g, std::ostream* out, std::string nametag);

    // print all the subgraphs as GML (iterator over subgraphs)
    void print_subgraphs (Graph& g, std::ostream* out, std::string nametag);

    // print all the vertex indices of a graph/subgraph
    void print_all_vertex_names (Graph& g, std::string nametag);
  }
}
#endif
