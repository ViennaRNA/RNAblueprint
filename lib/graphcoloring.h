/* This program reads secundary RNA structures in dot-bracket and
 * builds a graph for a latter ear-decomposition and bipartitness-check
 *
 * Created on: 13.08.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

#ifndef GRAPHCOLORING_H
#define GRAPHCOLORING_H

// include common header with graph definition and global variables
#include "common.h"
#include "pathcoloring.h"
#include "probability_matrix.h"
#include "graphcommon.h"
#include "printgraph.h"

// include standard library parts
#include <unordered_map>
#include <functional>
#include <iomanip>

// include boost components


namespace design {
  namespace detail {
    // function to do the coloring of the ear decomposition
    template <typename RG>
    void color_blocks (Graph& g, ProbabilityMatrix& pm, RG* rand_ptr);

    // color Articulation Points of this ear
    template <typename RG>
    MyKey color_articulation_points (int k, ProbabilityMatrix& pm, MyKey& colorkey, MyKey& lastkey, RG* rand_ptr);
  }
}

#endif
