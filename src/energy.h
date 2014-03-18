/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 23.01.2014
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

#ifndef ENERGY_H
#define ENERGY_H

// include common header with graph definition and global variables
#include "common.h"

// include standard library parts
#include <sstream>

// include boost components

// typedefs

// functions
float energy_of_structure (std::string& sequence, std::string& structure);

float fold (std::string& sequence, std::string& structure);



#endif
