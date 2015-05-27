/* This program reads secundary RNA structures in dot-bracket and
 * builds a graph for a latter ear-decomposition and bipartitness-check
 *
 * Created on: 13.08.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

#ifndef PATHCOLORING_H
#define PATHCOLORING_H

    // include common header with graph definition and global variables
    #include "common.h"
    #include "pairing_matrix.h"

    // include standard library parts
    #include <sstream>
    #include <unordered_map>

    // include boost components
    //#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/undirected_dfs.hpp>

namespace design {
    namespace detail {
        // typedefs
        typedef std::unordered_map< Vertex, std::unordered_map<int, unsigned long long> > nosMap;

        // class definitions

        template <typename RG>
        unsigned long long color_path_cycle_graph(Graph& g, RG* rand_ptr);
    }
}
#endif
