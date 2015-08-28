/*!\file parsestruct.h 
 * \brief This file holds the functions to parse a dot-bracket representation to a boost graph and set the sequence constraints.
 *
 * Created on: 25.03.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 * \cond INTERNAL
 */

#ifndef PARSESTRUCT_H
#define PARSESTRUCT_H

// include common header with graph definition and global variables
#include "common.h"
#include "graphcommon.h"

// include standard library parts
#include <sstream>
#include <ctype.h>

// include boost components

namespace design {
    namespace detail {
        // parse the input string into a graph
        Graph parse_structures(std::vector<std::string> structures);

        // set the sequence constraints in the graph object
        void set_constraints(Graph& g, std::string constraints);
    }
}
#endif


/* 
 * \endcond INTERNAL
 */