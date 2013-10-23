/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 23.10.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

#ifndef BOOSTTEST_H
#define BOOSTTEST_H

// include common header with graph definition and global variables
#include "../common.h"

// include standard library parts
#include <sstream>

// include boost components
#include <boost/test/unit_test.hpp>

// typedefs

//declare global variables
bool debug = false;
std::string outfile;

#endif
