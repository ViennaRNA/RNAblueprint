/*!\file pathcoloring.h 
 * \brief This file holds the functions to sample a sequence for a certain path.
 *
 * Created on: 13.08.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 * \cond INTERNAL
 */

#ifndef PATHCOLORING_H
#define PATHCOLORING_H

// include common header with graph definition and global variables
#include "common.h"
#include "graphcommon.h"
#include "printgraph.h"
#include "pairing_matrix.h"
#include "probability_matrix.h"

// include standard library parts
#include <sstream>
#include <unordered_map>

// include boost components
//#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/undirected_dfs.hpp>

namespace design {
    namespace detail {
        // typedefs
        typedef std::unordered_map< Vertex, std::unordered_map<int, SolutionSizeType> > nosMap;

        // class definitions
        
        // calculate a ProbabilityMatrix for a path
        ProbabilityMatrix get_path_pm(Graph& g);
        
        // color a path
        template <typename RG>
        SolutionSizeType color_path_graph(Graph& g, RG& rand);
        
    }
}
#endif


/* 
 * \endcond INTERNAL
 */