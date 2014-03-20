/* This program reads secundary RNA structures in dot-bracket and
 * builds a graph for a latter ear-decomposition and bipartitness-check
 *
 * Created on: 25.03.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

#ifndef PARSESTRUCT_H
#define PARSESTRUCT_H

// include common header with graph definition and global variables
#include "common.h"

// include standard library parts

// include boost components

// parse the input string into a graph
Graph parse_structures (std::vector<std::string> structures);

#endif
