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

// typedefs
typedef matrix< unsigned int, 4, 4 > rnaMatrix;

// class definitions
// Class to get Fibonacci numbers
class Fibonacci {
	public:
		Fibonacci(unsigned int l);
		unsigned int get(unsigned int n) { return numbers[n-1]; };
	private:
		std::vector< unsigned int > numbers;
};

// Class to get Pairing numbers
class Pairing {
	public:
		Pairing (unsigned int length);
		unsigned int get(unsigned int l, unsigned int b1, unsigned int b2);
		unsigned int get(unsigned int l, unsigned int b1);
		unsigned int get(unsigned int l);
		//matrix< unsigned int, 4, 4 > getMatrix(unsigned int l) { return p[l]; };
	private:
		std::vector< rnaMatrix > p;
		rnaMatrix multiply(rnaMatrix a, rnaMatrix b);
		unsigned int length;
};

// fills the string sequence with random bases, given the first and the last base and the intended length; returns the number of possible solutions
unsigned int generate_path_seq (std::string& sequence, int first, int last, int length);

// same for cycles (just a wrapper)
unsigned int generate_cycle_seq (std::string& sequence, int first, int length);

#endif
