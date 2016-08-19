/*!\file printgraph.h 
 * \brief This file holds the functions to print out a graph in xml representation.
 *
 * @date 25.03.2013
 * @author Stefan Hammer <s.hammer@univie.ac.at>
 * @copyright GPLv3
 *
 * \cond INTERNAL
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

    // print the graph as a gml file to the output
    void print_graph (Graph& g, std::ostream* out);

    // print all the subgraphs as GML (iterator over subgraphs)
    void print_subgraphs (Graph& g, std::ostream* out);
  }
}
#endif

/* 
 * \endcond INTERNAL
 */