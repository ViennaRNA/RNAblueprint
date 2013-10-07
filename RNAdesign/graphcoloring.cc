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
		if (my < (unsigned int) g[e].ear) { my = g[e].ear; }
	}
	my++; // my is not the biggest ear index but the total number of ears!
	
	// get maximal length of an ear
	unsigned int length;
	for (unsigned int k = 0; k < my; k++) {
		unsigned int tmp_length = 0;
		BGL_FORALL_EDGES_T(e, g, Graph) {
			if ((unsigned int) g[e].ear == k) { tmp_length++; }
		}
		if (length < tmp_length) { length = tmp_length; }
	}
	
	// get Pairing matrix for path, TODO only initialize once for the whole program!
	Pairing p(length+1);
	
	// iterate over all ear decomposition iterations 
	for (unsigned int k = 0; k < my; k++) {
		// Ak are already stored in graph as a vertex propertys
		// write vertex property into Ak[k]
		BGL_FORALL_VERTICES_T(v, g, Graph) {
			if (g[v].Ak.find(k) != g[v].Ak.end()) {
				Ak[k].insert(v);
			}
		}
	}
	
	// store internal Articulation Points
	std::set< int > Ai;
	// now start at the outermost ear
	unsigned int k = 0;
	
	Graph::children_iterator ear, ear_end;
	for (boost::tie(ear, ear_end) = (g).children(); ear != ear_end; ++ear) {
		for (auto ap : Ak[k]) {
			for (unsigned int b = 0; b < A_Size; b++) {
				// Nk[A6][A10][A1] = sum(AUGC in inner Ap = 9) P[A6][x9][3 pathlength] * P[x9][A10][1] * Nk-1 [x9][A1]
				// std::vector< std::map < unsigned int, std::array<unsigned long long, A_Size > > > n;
				
				n[k][ap][b] = 0;
				
			}
		}
		
		
		// iterate over articulation points of this ear
		// find out if this Aks are still Ak of next ear or if they are internal then (-> push into Ai)
		Ai.clear();
		for (auto v : Ak[k]) {
			if (Ak[k+1].find(v) == Ak[k+1].end()) {				// vertices are no Aps in next k
				BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {			// look at all edges
					if (g[e].ear == k+1) {				// if next ear is glued here this will be internal
						Ai.insert(v);
					} else {					// else it is still an external next k
						Ak[k+1].insert(v);
					}
				}
			}
		}
		
		// now going to next ear!
		k++;
	}
}

unsigned long long ProbabilityMatrix::get(unsigned int k, unsigned int a, unsigned int b) {
	if ((k > my) || (b > A_Size-1)) {
		std::cerr << "Requested a value in probability matrix which is out of range: p[" << k << "][" << a << "][" << b << "]" << std::endl;
		exit(1);
	}
	
	unsigned long long rvalue;
	// important for map: if you request with [] an entry will be created for unexisting ones.
	if (n[k].find(a) != n[k].end()) {
		rvalue = n[k][a][b];
	} else {
		rvalue = 0;
	}
	
	return rvalue;
}

unsigned long long ProbabilityMatrix::get(unsigned int k, unsigned int a) {
	// return the sum of all probabilities of articulation point
	unsigned long long sum = 0;
	for (unsigned int i = 0; i < A_Size; i++) {
		sum += get(k, a, i);
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
					color_blocks(*bc);
				}
			}
			
			for (boost::tie(bc, bc_end) = (*cc).children(); bc != bc_end; ++bc) {
				if (get_min_max_degree(*bc).second <= 2) {
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
