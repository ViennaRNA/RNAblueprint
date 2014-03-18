/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 26.08.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

#ifndef GRAPHCOMMON_H
#define GRAPHCOMMON_H

// include common header with graph definition and global variables
#include "common.h"

// include standard library parts
#include <limits>

// get max degree of a graph
std::pair <int, int> get_min_max_degree(Graph& g);

#endif
