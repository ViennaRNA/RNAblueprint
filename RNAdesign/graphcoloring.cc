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

std::size_t MyKeyHash::operator() (const MyKey& k) const {
	// Start with 0 as a hash value   .
	std::size_t hash = 0;
	
	for (auto elem : k) {
		boost::hash_combine(hash, boost::hash_value(elem.first));
		boost::hash_combine(hash, boost::hash_value(elem.second));
	}
	
	return hash;
}



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
	
	// store inner Articulation Points
	std::set< int > Ai;
	// now start at the outermost ear
	unsigned int k = 0;

	
	// start at the outermost ear and process inwards
	Graph::children_iterator ear, ear_end;
	for (boost::tie(ear, ear_end) = (g).children(); ear != ear_end; ++ear) {
		// Nk[A6][A10][A1] = sum(AUGC in inner Ap = 9) P[A6][x9][3 pathlength] * P[x9][A10][1] * Nk-1 [x9][A1]
		
		std::vector<MyKey> key_combinations;		// this is what we want to fill next
		MyKey mykey;					// helper to recursively build the posibilities
		std::set<Vertex> ap = Ak[k];			// need to send a copy of the current Artikulation Points
		
		calculate_probabilities(ap, mykey, key_combinations);
		
		for (auto thiskey : key_combinations) {
			unsigned long long probability = get_probability(thiskey, *ear, Ak[k], Ai);
			if (probability != 0)
				n[k][thiskey] = probability;
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

void ProbabilityMatrix::calculate_combinations(std::set<Vertex>& ap, MyKey& mykey, std::vector<MyKey>& key_combinations) {
	
	if (ap.size() > 0) {
		std::set<Vertex>::iterator it=ap.begin();
		Vertex v = *it;
		ap.erase(it);
		
		for ( unsigned int b = 0; b < A_Size; b++ ) {
			mykey.insert(std::make_pair<int,int>(v,b));
			// recursion starts here
			calculate_probabilities(ap, mykey);
			
			if (ap.size() == 0) {
				// remember our generated key
				key_combinations.push_back(mykey);
			}
			// remove current vertex again to make space for a new base
			mykey.erase(v);
		}
		// add current vertex again
		ap.insert(v);
	}
}

unsigned long long ProbabilityMatrix::get_probability ( MyKey mykey, Graph& g, std::set<Vertex>& ap, std::set<Vertex>& ai) {
	unsigned long long max_number_of_sequences = 0;	
	
	// Nk[A6][A10][A1] = sum(AUGC in inner Ap = 9) P[A6][x9][3 pathlength] * P[x9][A10][1] * Nk-1 [x9][A1]
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		g[v].color = 0;
	}
	
	// get Aks of this 
	std::vertex<Vertex> thisaps;
	for ( v : ap ) {
		if (g.find_vertex(v).second) {
			thisaps.push_back(v);
		}
	}
		
	int length;
	Vertex v = thisaps[0];
	while (v != thisaps[1];
		v = get_length_to_next_internal_ap(g, length, v, ai);		
	}
	
	
	return max_number_of_sequences;
}

Vertex get_length_to_next_internal_ap(Graph& g, int& length, Vertex v, std::set<Vertex> ai) {
	Vertex v;
	
	BGL_FORALL_ADJ_T(vertex, adj, g, Graph) {
		if (g[adj].color == 0) {
			g[adj].color = 1;
			length++;
			
		}
	}
	
	return v;
}

unsigned long long ProbabilityMatrix::get(unsigned int k, unsigned int a, unsigned int b) {
//	if ((k > my) || (b > A_Size-1)) {
//		std::cerr << "Requested a value in probability matrix which is out of range: p[" << k << "][" << a << "][" << b << "]" << std::endl;
//		exit(1);
//	}
	
	unsigned long long rvalue;
	// important for map: if you request with [] an entry will be created for unexisting ones.
//	if (n[k].find(a) != n[k].end()) {
//		rvalue = n[k][a][b];
//	} else {
//		rvalue = 0;
//	}
	
	return rvalue;
}

unsigned long long ProbabilityMatrix::get(unsigned int k, unsigned int a) {
	// return the sum of all probabilities of articulation point
	unsigned long long sum = 0;
//	for (unsigned int i = 0; i < A_Size; i++) {
//		sum += get(k, a, i);
//	}
	return sum;
}

void color_graph (Graph& graph) {
	// root graph
	
	// reset bases to X
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
	ProbabilityMatrix pm(g);
	
	// backtracing
}

void reset_colors(Graph& g) {
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		g[v].base = X;
	}
}
