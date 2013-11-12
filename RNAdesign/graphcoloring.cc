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
	
	// structure to remember current Ak (attachment vertices)
	std::set<Vertex> currentAk;
	// to store current inner Articulation Points
	std::set<Vertex> currentAi;
	// now start at the outermost ear
	unsigned int k = 0;
	
	int max_length = 0;
	Graph::children_iterator ear, ear_end;
	for (boost::tie(ear, ear_end) = (g).children(); ear != ear_end; ++ear) {
		int length = boost::num_edges(*ear);
		if (max_length < length)
			max_length = length;
	}

	// get Pairing matrix for paths, TODO only initialize once for the whole program!
	Pairing p(max_length+1);
	
	// start at the outermost ear and process inwards
	for (boost::tie(ear, ear_end) = (g).children(); ear != ear_end; ++ear) {
		
		// before doing anything, update current Ak and Ai
		updateCurrentAkAi (*ear, k, currentAk, currentAi);
		
		if (debug) { 
			std::cerr << "Current k: " << k << std::endl
				<< "currentAk:" << std::endl << currentAk << std::endl
				<< "currentAi:" << std::endl << currentAi << std::endl;
		}
		
		// Nk[A6][A10][A1] = sum(AUGC in inner Ap = 9) P[A6][x9][3 pathlength] * P[x9][A10][1] * Nk-1 [x9][A1]
		
		std::vector<MyKey> key_combinations;		// this is what we want to fill in the recursion
		MyKey mykey;					// helper to recursively build the posibilities
		std::set<Vertex> cAk = currentAk;		// create a copy for recursion
		// now fill the key_combinations with all kinds of bases
		calculate_combinations(g, cAk, mykey, key_combinations);
		
		// calculate the probabilities for every key and if not zero, add to matrix
		for (auto thiskey : key_combinations) {
			if (debug) { std::cerr << "Calculating probablity for key: " << std::endl << thiskey; }
			unsigned long long probability = get_probability(thiskey, *ear, currentAk, currentAi, p, k);
			if (debug) { std::cerr << probability << std::endl; }
			if (probability != 0)
				n[k][thiskey] = probability;
		}
		
		// now going to next ear!
		my = ++k;
	}
}

void ProbabilityMatrix::updateCurrentAkAi (Graph& g, int k, std::set<Vertex>& currentAk, std::set<Vertex>& currentAi) {
	// Ak and Ai are already stored in graph as a vertex properties
	// write vertex property into currentAk and currentAi
	currentAi.clear();
	
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		if (g[v].Ak.find(k) != g[v].Ak.end()) {
			currentAk.insert(v);
		} else if ((k > 0) && (g[v].Ai == k)) {
			currentAi.insert(v);
			// we need to keep Ak from previous glued ears, except those that became Ai this time!
			currentAk.erase(v);
		}
	}		
}

void ProbabilityMatrix::calculate_combinations (Graph& g, std::set<Vertex>& Ak, MyKey& mykey, std::vector<MyKey>& key_combinations) {
	
	if (Ak.size() > 0) {
		std::set<Vertex>::iterator it=Ak.begin();
		Vertex v = *it;
		Ak.erase(it);
		
		for ( unsigned int b = 0; b < A_Size; b++ ) {
			mykey.insert(std::make_pair(boost::get(boost::vertex_color_t(), g, v), b));
			// recursion starts here
			calculate_combinations(g, Ak, mykey, key_combinations);
			
			if (Ak.size() == 0) {
				// remember our generated key
				key_combinations.push_back(mykey);
			}
			// remove current vertex again to make space for a new base
			mykey.erase(v);
		}
		// add current vertex again
		Ak.insert(v);
	}
}

unsigned long long ProbabilityMatrix::get_probability ( MyKey mykey, Graph& g, std::set<Vertex>& Ak, std::set<Vertex>& Ai, Pairing& p, unsigned int k) {
	
	// v--return this--v v--recursion get_sum_of_sum--v v--let's call them sub_probabilities--v   v--prob. last_ear--v
	// Nk[A6][A10][A1]  =   sum(AUGC in inner Ap = 9)   P[A6][x9][3 pathlength] * P[x9][A10][1] * Nk-1 [x9][A1]
	
	// this is what we calculate here
	unsigned long long max_number_of_sequences = 0;
	// key to get probability from last ear glueing and to store all base combinations for internal articulation points
	MyKey lastkey;
	// container to store all sub-probabilities
	std::vector<SubProbability> sub_probabilities;
	
	
	// start to build last key
	for (auto v : Ak) {
		if (!g.find_vertex(v).second) {
			// these articulation points are not in this particular ear, therefore we need to look their probability up
			// from last time. so lets generate a key therefore			
			lastkey.insert(std::make_pair(boost::get(boost::vertex_color_t(), g, v), mykey[boost::get(boost::vertex_color_t(), g, v)]));
		}
	}
	
	// now building this sub_probabilities
	Graph::children_iterator part, part_end;
	for (boost::tie(part, part_end) = (g).children(); part != part_end; ++part) {
		
		// add to sub_probabilities for this start/end and length
		sub_probabilities.push_back(SubProbability());
		SubProbability &sub_probability = sub_probabilities.back();
		// begin, end and length go into sub_probability
		int i = 0;
		BGL_FORALL_VERTICES_T(v, *part, Graph) {
			if (boost::degree(v, *part) == 1) {
				if (i == 0)
					sub_probability.start = boost::get(boost::vertex_color_t(), *part, v);
				else
					sub_probability.end = boost::get(boost::vertex_color_t(), *part, v);
				i++;
			}
		}
		sub_probability.length = boost::num_edges(*part);
	}
	
	// calculate sum of sum for all bases colored X (= internal aps)
	// adds base combinatoric to the sub_probabilities and to the lastkey
	if (k > 0) {
		make_sum_of_sum(g, Ai, mykey, lastkey, sub_probabilities, p, k, max_number_of_sequences);
	} else {
		for (auto sub_probability : sub_probabilities) {
			int startBase = get_color_from_key(mykey, lastkey, sub_probability.start);
			int endBase = get_color_from_key(mykey, lastkey, sub_probability.end);
			// get probability and multiply it
			max_number_of_sequences += p.get(startBase, endBase, sub_probability.length);
		}
	}
		
	return max_number_of_sequences;
}

void ProbabilityMatrix::make_sum_of_sum(	Graph& g,
						std::set<Vertex>& Ai, 
						MyKey& mykey, MyKey& lastkey, 
						std::vector<SubProbability>& sub_probabilities, 
						Pairing& p,
						unsigned int k,
						unsigned long long& max_number_of_sequences) 
{
	
	if (Ai.size() > 0) {
		std::set<Vertex>::iterator it=Ai.begin();
		Vertex v = *it;
		Ai.erase(it);
		
		for ( unsigned int b = 0; b < A_Size; b++ ) {
			lastkey.insert(std::make_pair(boost::get(boost::vertex_color_t(), g, v), b));
			// recursion starts here
			make_sum_of_sum(g, Ai, mykey, lastkey, sub_probabilities, p, k, max_number_of_sequences);
			
			if (Ai.size() == 0) {
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
				multiplied_probabilities *= n[k-1][lastkey];
				// now add all the multiplied probabilities to the total to get sum over all (AUGC) in X
				max_number_of_sequences += multiplied_probabilities;
			}
			// remove current vertex again to make space for a new base
			lastkey.erase(v);
		}
		// add current vertex again
		Ai.insert(v);
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

unsigned long long ProbabilityMatrix::get(unsigned int k, MyKey mykey) {
	if (k > my) {
		std::cerr << "Requested a value in probability matrix where k is out of range: " << k << std::endl;
		exit(1);
	}
	
	unsigned long long rvalue;
	// important for map: if you request with [] an entry will be created for unexisting ones.
	if (n[k].find(mykey) != n[k].end()) {
		rvalue = n[k][mykey];
	} else {
		rvalue = 0;
	}
	
	return rvalue;
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

// overload << operator to print maps with any content
std::ostream& operator<< (std::ostream& os, MyKey& m) {
	for (typename MyKey::iterator it = m.begin(); it != m.end(); it++) {
        	os << "(" << it->first << "," << it->second << ")" << std::endl;
	}
	return os;
}
