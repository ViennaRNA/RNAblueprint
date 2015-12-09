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
#include <cctype>
#include <iterator>

// include boost components
#include <boost/regex.hpp>

namespace design {
    namespace detail {
        
        typedef std::vector< std::pair<char, char> > BracketList;
        
        // parse the input string into a graph
        Graph parse_structures(std::vector<std::string> structures);
        // recursion for different brackets
        void parse_bracket(Graph& g, std::string& structure, BracketList::iterator bracket);
        // set the sequence constraints in the graph object
        void set_constraints(Graph& g, std::string constraints);
    }
}
#endif


/* 
 * \endcond INTERNAL
 */
