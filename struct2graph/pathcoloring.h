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

// include standard library parts
#include <sstream>

// include boost components

// class definitions
// Class to get Fibonacci numbers
class Fibonacci {
	public:
		Fibonacci(unsigned int length);
		unsigned int get(unsigned int n) { return numbers[n-1]; };
	private:
		std::vector< unsigned int > numbers;
};

// fills the string sequence with random bases, given the first and the last base and the intended length; returns the number of possible solutions
unsigned int generate_path_seq (std::string& sequence, char first, char last, unsigned int length);

// fills the string sequence with random bases, given the start/end-base of the cycle and the intended length; returns the number of possible solutions
unsigned int generate_cycle_seq (std::string& sequence, char first, unsigned int length);


#endif
