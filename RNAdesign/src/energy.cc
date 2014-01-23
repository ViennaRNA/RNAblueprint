/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 23.01.2014
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

// include header
#include "energy.h"

extern "C" {
	#include "ViennaRNA/fold.h"
}

// typedefs


// global variables


// main program starts here
float energy_of_structure (std::string sequence, std::string structure) {
	float energy = energy_of_structure(sequence.c_str(), structure.c_str(), 0);
	return energy; 
}
