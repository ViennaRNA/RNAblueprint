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