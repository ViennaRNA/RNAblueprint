/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 26.08.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

// include header
#include "graphcommon.h"

std::pair <int, int> get_min_max_degree(Graph& g) {
	int max_degree = 0;
	int min_degree = std::numeric_limits<int>::max();
	
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		int current_degree = boost::out_degree(v, g);
		if (current_degree > max_degree) {
			max_degree = current_degree;
		} else if (current_degree < min_degree) {
			min_degree = current_degree;
		}
	}
	return std::make_pair(min_degree, max_degree);
}

