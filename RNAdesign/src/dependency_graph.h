/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 18.03.2014
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

#ifndef DEPENDENCY_GRAPH_H
#define	DEPENDENCY_GRAPH_H

#include "common.h"

class DependencyGraph {
	public:
		DependencyGraph (std::vector<std::string> structures);
		unsigned long long number_of_sequences() { return nos; }
		bool is_bipartite() { return bipartite; }
		Sequence get_sequence();
		Sequence mutate (int position);
		Sequence mutate ();
		~DependencyGraph();
	private:
		Graph graph;
		bool bipartite;			// if dependency graph is bipartite and a therefore a solution exists
		unsigned long long nos = 0;	// number of sequences/solutions
		//std::list<&Graph> graph_components;
		
		
};

#endif	/* DEPENDENCY_GRAPH_H */

