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
#include <unordered_map>
#include <functional>

// include boost components
#include <boost/functional/hash.hpp>

// typedefs
typedef std::unordered_map < int, int > MyKey;

// Class with cusom hash function for the ProbabilityMatrix
class MyKeyHash {
	public:
		std::size_t operator() (const MyKey& k) const;
};

// class definitions
// Class to get Pairing numbers
class ProbabilityMatrix {
	public:
		ProbabilityMatrix (Graph& g);
		unsigned long long get(unsigned int k, unsigned int a, unsigned int b);		// k... ear (k), a Ak# (Vertex), b... Base
		unsigned long long get(unsigned int k, unsigned int a);				// sum of all Bases of this Ak
	private:
		// vector of k (ears), map of possibilities saved by key
		std::vector< std::unordered_map < MyKey , unsigned long long , MyKeyHash> > n;
		// structure to remember Ak (attachment vertices)
		std::map<int, std::set<Vertex> > Ak;
		// biggest ear number
		unsigned int my = 0;
		// My custom hash key used for n
		friend class MyKeyHash;
		// calculate all the key probabilities and start more from there
		void calculate_probabilities(std::set<Vertex>& ap, MyKey& k);
};

// color the root graph!
void color_graph (Graph& graph);

// function to do the coloring of the ear decomposition
void color_blocks (Graph& g);

// reset the bases of the graph to X
void reset_colors(Graph& g);

#endif
