/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 13.08.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

// include header
#include "graphcoloring.h"
#include "graphcommon.h"
#include "pathcoloring.h"

ProbabilityMatrix::ProbabilityMatrix (Graph& g) {
	// get number of ears in this ear decomposition
	BGL_FORALL_EDGES_T(e, g, Graph) {
		if (my < g[e].ear) { my = g[e].ear; }
	}
	
	// get maximal length of an ear
	unsigned int length;
	for (unsigned int k = 0; k < my+1; k++) {
		unsigned int tmp_length = 0;
		BGL_FORALL_EDGES_T(e, g, Graph) {
			if (g[e].ear == k) { tmp_length++; }
		}
		if (length < tmp_length) { length = tmp_length; }
	}
	
	// get Pairing matrix for path, TODO only initialize once for the whole program!
	Pairing p(length+1);
	
	// structure to remember Ak (attachment vertices)
	std::map<int, std::set<Vertex> > Ak;
	// iterate over all ear decomposition iterations 
	for (int k = 0; k < my+1; k++) {
		// Ak are already stored in graph as a vertex propertys
		// write vertex property into Ak[k]
		BGL_FORALL_VERTICES_T(v, g, Graph) {
			if (g[v].Ak.find(k) != g[v].Ak.end()) {
				Ak[k].insert(v);
			}
		}
	}
	
	// now start at the outermost ear
	for (unsigned int k = 0; k < my+1; k++) {
		for (auto ap : Ak[k]) {
			// find out if previous Aks are still Ak or if they are internal now
			for (unsigned int i = 0; i < A_Size; i++) {
				// Nk[A6][A10][A1] = sum(AUGC in inner Ap = 9) P[A6][x9][3 pathlength] * P[x9][A10][1] * Nk-1 [x9][A1]
				
			}
		}
	}
}

unsigned long long ProbabilityMatrix::get(unsigned int e, unsigned int a, unsigned int b) {
	if ((e > my) || (b > A_Size-1)) {
		std::cerr << "Requested a value in probability matrix which is out of range: p[" << e << "][" << a << "][" << b << "]" << std::endl;
		exit(1);
	}
	
	unsigned long long rvalue;
	// important for map: if you request with [] an entry will be created for unexisting ones.
	if (p[e].find(a) != p[e].end()) {
		rvalue = p[e][a][b];
	} else {
		rvalue = 0;
	}
	
	return rvalue;
}

unsigned long long ProbabilityMatrix::get(unsigned int e, unsigned int a) {
	// return the sum of all probabilities of articulation point
	unsigned long long sum = 0;
	for (unsigned int i = 0; i < A_Size; i++) {
		sum += get(e, a, i);
	}
	return sum;
}

void color_graph (Graph& graph) {
	// root graph
	
	// reset basese to X
	reset_colors(graph);
	
	Graph::children_iterator cc, cc_end;
	for (boost::tie(cc, cc_end) = graph.children(); cc != cc_end; ++cc) {
		// connected components
		if (get_min_max_degree(*cc).second <= 2) {
			// color connected components here
			color_path_cycle_graph (*cc);
		} else {
			Graph::children_iterator bc, bc_end;
			for (boost::tie(bc, bc_end) = (*cc).children(); bc != bc_end; ++bc) {
				// biconnected components (color blocks first!)
				if (get_min_max_degree(*bc).second > 2) {
					// blocks
					// color blocks here
					//color_blocks(*bc);
				}
			}
			
			for (boost::tie(bc, bc_end) = (*cc).children(); bc != bc_end; ++bc) {
				if (get_min_max_degree(*bc).second >= 2) {
					// biconnected component paths
					// color paths here
					color_path_cycle_graph (*bc);
				}
			}
		}
	}
}

void color_blocks (Graph& g) {
	// start with filling the matrix
	
	
	
	// backtracing
}

void reset_colors(Graph& g) {
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		g[v].base = X;
	}
}
