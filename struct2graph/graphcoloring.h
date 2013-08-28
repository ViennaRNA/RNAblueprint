/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 13.08.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

#ifndef GRAPHCOLORING_H
#define GRAPHCOLORING_H

// include common header with graph definition and global variables
#include "common.h"

// include standard library parts

// include boost components

// typedefs

// class definitions
// Class to get Pairing numbers
class ProbabilityMatrix {
	public:
		ProbabilityMatrix (Graph& g);
		unsigned long long get(unsigned int e, unsigned int a, unsigned int b);		// e... ear (k), a Ak# (Vertex), b... Base
		unsigned long long get(unsigned int e, unsigned int a);				// sum of all Bases of this Ak
	private:
		// vector of k (ears), map of artikulation points (Vertex number) and array of bases.
		std::vector< std::map < unsigned int, std::array<unsigned long long, A_Size > > > p;
		// biggest ear number
		unsigned int my = 0;
};

// color the root graph!
void color_graph (Graph& graph);

// function to do the coloring of the ear decomposition
void color_blocks (Graph& g);

// reset the bases of the graph to X
void reset_colors(Graph& g);

#endif
