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
	unsigned int max_length;
	for (unsigned int k = 0; k < my; k++) {
		unsigned int tmp_length = 0;
		BGL_FORALL_EDGES_T(e, g, Graph) {
			if ((unsigned int) g[e].ear == k) { tmp_length++; }
		}
		if (max_length < tmp_length) { max_length = tmp_length; }
	}
	
	// get Pairing matrix for paths, TODO only initialize once for the whole program!
	Pairing p(max_length+1);
	
	// iterate over all ear decomposition iterations 
	for (unsigned int k = 0; k < my; k++) {
		// Ak are already stored in graph as a vertex properties
		// write vertex property into Ak[k]
		BGL_FORALL_VERTICES_T(v, g, Graph) {
			if (g[v].Ak.find(k) != g[v].Ak.end()) {
				Ak[k].insert(v);
			}
		}
	}
	
	// for last ear, which is a cycle, add two times the same articulaiton point (begin and end!)
	
	
	// to store inner Articulation Points
	std::set< Vertex > Ai;
	// now start at the outermost ear
	unsigned int k = 0;

	
	// start at the outermost ear and process inwards
	Graph::children_iterator ear, ear_end;
	for (boost::tie(ear, ear_end) = (g).children(); ear != ear_end; ++ear) {
		// Nk[A6][A10][A1] = sum(AUGC in inner Ap = 9) P[A6][x9][3 pathlength] * P[x9][A10][1] * Nk-1 [x9][A1]
		
		std::vector<MyKey> key_combinations;		// this is what we want to fill next
		MyKey mykey;					// helper to recursively build the posibilities
		std::set<Vertex> ap = Ak[k];			// need to send a copy of the current Artikulation Points
		// now fill the key_combinations with all kinds of bases
		calculate_combinations(g, ap, mykey, key_combinations);
		
		for (auto thiskey : key_combinations) {
			unsigned long long probability = get_probability(thiskey, *ear, ap, Ai, p, k);
			if (probability != 0)
				n[k][thiskey] = probability;
		}
		
		// iterate over articulation points of this ear
		// find out if this Aks are still Ak of next ear or if they are internal then (-> push into Ai)
		Ai.clear();
		for (auto v : Ak[k]) {
			if (Ak[k+1].find(v) == Ak[k+1].end()) {				// vertices are no Aps in next k
				BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {			// look at all edges
					if (g[e].ear == (int) k+1) {				// if next ear is glued here this will be internal
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

void ProbabilityMatrix::calculate_combinations(Graph& g, std::set<Vertex>& ap, MyKey& mykey, std::vector<MyKey>& key_combinations) {
	
	if (ap.size() > 0) {
		std::set<Vertex>::iterator it=ap.begin();
		Vertex v = *it;
		ap.erase(it);
		
		for ( unsigned int b = 0; b < A_Size; b++ ) {
			mykey.insert(std::make_pair(boost::get(boost::vertex_color_t(), g, v), b));
			// recursion starts here
			calculate_combinations(g, ap, mykey, key_combinations);
			
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

unsigned long long ProbabilityMatrix::get_probability ( MyKey mykey, Graph& g, std::set<Vertex> ap, std::set<Vertex>& ai, Pairing& p, unsigned int k) {
	
	// v--return this--v v--recursion get_sum_of_sum--v v--let's call them sub_probabilities--v   v--prob. last_ear--v
	// Nk[A6][A10][A1]  =   sum(AUGC in inner Ap = 9)   P[A6][x9][3 pathlength] * P[x9][A10][1] * Nk-1 [x9][A1]
	
	// this is what we calculate here
	unsigned long long max_number_of_sequences = 0;
	// key to get probability from last ear glueing and to store all base combinations for internal articulation points
	MyKey lastkey;
	// container to store all sub-probabilities
	std::vector<SubProbability> sub_probabilities;
	
	// reset the color of this graph to be able to run our walk_through.
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		g[v].color = 0;
	}
	// get end points of path (ear)
	std::array<Vertex, 2> start_end;
	int i = 0;
	for ( auto v : ap ) {
		if (g.find_vertex(v).second) {
			start_end[i] = v;
			i++;
		} else {
			// these articulation points are not in this particular ear, therefore we need to look their probability up
			// from last time. so lets generate a key therefore			
			lastkey.insert(std::make_pair(boost::get(boost::vertex_color_t(), g, v), mykey[boost::get(boost::vertex_color_t(), g, v)]));
		}
	}
	
	// now building this sub_probabilities (no base combinatoric yet!)
	// function will start at one end and go through to get all lengths and end-points
	Vertex nextVertex = start_end[0];
	
	while (nextVertex != start_end[1]) {
		int length = 0;
		
		// add to sub_probabilities for this start/end and length
		sub_probabilities.push_back(SubProbability());
		SubProbability &sub_probability = sub_probabilities.back();
		// start vertex first
		sub_probability.start = nextVertex;
		
		// get length and vertex descriptor to next articulation point (inner)
		std::tie(nextVertex, length) = get_length_to_next_ap(g, nextVertex, ai);
		
		// now end vertex and length
		sub_probability.end = nextVertex;
		sub_probability.length = length;
	}
	
	// calculate sum of sum for all bases colored X (= internal aps)
	// adds base combinatoric to the sub_probabilities and to the lastkey
	make_sum_of_sum(g, ai, mykey, lastkey, sub_probabilities, p, k, max_number_of_sequences);
		
	return max_number_of_sequences;
}

void ProbabilityMatrix::make_sum_of_sum(	Graph& g,
						std::set<Vertex>& ai, 
						MyKey& mykey, MyKey& lastkey, 
						std::vector<SubProbability>& sub_probabilities, 
						Pairing& p,
						unsigned int k,
						unsigned long long& max_number_of_sequences) 
{
	
	if (ai.size() > 0) {
		std::set<Vertex>::iterator it=ai.begin();
		Vertex v = *it;
		ai.erase(it);
		
		for ( unsigned int b = 0; b < A_Size; b++ ) {
			lastkey.insert(std::make_pair(boost::get(boost::vertex_color_t(), g, v), b));
			// recursion starts here
			make_sum_of_sum(g, ai, mykey, lastkey, sub_probabilities, p, k, max_number_of_sequences);
			
			if (ai.size() == 0) {
				// do the actual calculation here!
				unsigned long long multiplied_probabilities;
				
				// calculate product of sub probabilities
				for (auto sub_probability : sub_probabilities) {
					int startBase = get_color_from_key(mykey, lastkey, sub_probability.start);
					int endBase = get_color_from_key(mykey, lastkey, sub_probability.end);
					// get probability and multiply it
					multiplied_probabilities *= p.get(startBase, endBase, sub_probability.length);
				}
	
				// now add probability for last_ear
				multiplied_probabilities *= n[k][lastkey];
				// now add all the multiplied probabilities to the total to get sum over all (AUGC) in X
				max_number_of_sequences += multiplied_probabilities;
			}
			// remove current vertex again to make space for a new base
			lastkey.erase(v);
		}
		// add current vertex again
		ai.insert(v);
	}
}

int ProbabilityMatrix::get_color_from_key(MyKey& mykey, MyKey& lastkey, Vertex v) {
	if (mykey.find(v) != mykey.end()) {
		return mykey[v];
	} else if (lastkey.find(v) != lastkey.end()) {
		return lastkey[v];
	} else {
		std::cerr << "Something went wrong! I could not find the color of the vertex in any key (neither last nor current): " 
			<< v << std::endl;
		exit(1);
	}
}

std::pair<Vertex, int> ProbabilityMatrix::get_length_to_next_ap(Graph& g, Vertex start, std::set<Vertex>& ai) {
	// go throught the path (ear) and color everything 1 to avoid getting back
	// if i went forward, count length plus one, color 1 and if the reached vertex is a internal articulation point or
	// the other end, return the length and the vertex
	g[start].color = 1;
	bool end_reached = false;
	int length = 0;
	
	while(!end_reached) {
		BGL_FORALL_ADJ_T(start, adj, g, Graph) {
			if (g[adj].color == 0) {
				g[adj].color = 1;
				end_reached = false;
				length++;
				start = adj;
				break;
			} else {
				end_reached = true;
			}
		}
		// if the current vertex is a internal articulation point, return length and this vertex
		if (ai.find(start) != ai.end()) {
			return std::make_pair(boost::get(boost::vertex_color_t(), g, start), length);
		}
	}
	// return final length and vertex
	return std::make_pair(boost::get(boost::vertex_color_t(), g, start), length);
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
	//ProbabilityMatrix pm(g);
	
	// backtracing
}

void reset_colors(Graph& g) {
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		g[v].base = X;
	}
}
