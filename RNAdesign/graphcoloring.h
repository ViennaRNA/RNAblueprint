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
#include "pathcoloring.h"

// include standard library parts
#include <unordered_map>
#include <functional>
#include <iomanip>

// include boost components
#include <boost/functional/hash.hpp>
#include <boost/iterator/reverse_iterator.hpp>

// typedefs
typedef std::map < int, int > MyKey;

// overload << operator to print maps with any content
std::ostream& operator<< (std::ostream& os, MyKey& m);

// class definitions
// Class with cusom hash function for the ProbabilityMatrix
class MyKeyHash {
	public:
		std::size_t operator() (const MyKey& k) const;
};

// struct definitions
struct SubProbability {
		int start;
		int end;
		int length;	
	};

// Class to get Pairing numbers
class ProbabilityMatrix {
	public:
		ProbabilityMatrix (Graph& g);
		// get probability for mykey... key of Aks (12/A) (4/C) ()...
		unsigned long long get(MyKey mykey);
		// get sum of probabilities
		unsigned long long get_sum (MyKey mykey, std::vector<MyKey>& key_combinations);
		// get all the articulation points for a specific ear
		std::set<Vertex> get_Ak (unsigned int k);
		// return my of this graph
		unsigned int get_my() { return my; }
		// My custom hash key used for n
		friend class MyKeyHash;
	private:
		// map of possibilities saved by key
		std::unordered_map < MyKey , unsigned long long , MyKeyHash> n;
		// remember my as number of ears
		unsigned int my = 0;
		// remember all the Articulation points for every ear
		std::vector<std::set<Vertex> > Aks;
		// to keep track of current Ak Ai we need to update them before glueing an ear
		void updateCurrentAkAi (Graph& g, int k, std::set<Vertex>& currentAk, std::set<Vertex>& currentAi);
		// calculate all the key probabilities and start more from there
		void calculate_probabilities (std::set<Vertex>& ap, MyKey& k);
		// actually calculate the probability for the given key
		unsigned long long get_probability ( MyKey mykey, Graph& g, std::set<Vertex>& Ak, std::set<Vertex>& Ai, Pairing& p, unsigned int k);
		// recursion to get base combinations done in (sum over AUGC in 6) of (sum over AUGC in 10) of ...
		void make_sum_of_sum( Graph& g, std::set<Vertex>& Ai, MyKey& mykey, MyKey& lastkey, std::vector<SubProbability>& sub_probabilities, Pairing& p, unsigned int k, unsigned long long& max_number_of_sequences);
		// do the actual calculation of the multiplied probabilities in the sum_of_sum
		unsigned long long calculate_probability (MyKey& mykey, MyKey& lastkey, std::vector<SubProbability>& sub_probabilities, Pairing& p);
		// get the color of either mykey or lastkey.
		int get_color_from_key (MyKey& mykey, MyKey& lastkey, int vertex);
		// calculates all base combinations for current articulation points (recursion)
		void calculate_combinations (std::set<int>& Ak, MyKey& mykey, std::vector<MyKey>& key_combinations);
};

// color the root graph!
void color_graph (Graph& graph);

// function to do the coloring of the ear decomposition
void color_blocks (Graph& g);

// color Articulation Points of this ear
MyKey color_articulation_points(ProbabilityMatrix& pm, MyKey& colorkey);

// reset the bases of the graph to X
void reset_colors(Graph& g);


#endif
