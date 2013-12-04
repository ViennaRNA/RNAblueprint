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
	
	//std::cerr << "hash is: " << hash << std::endl;
	
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
		// also get my (number of ears)
		my++;
	}
	
	if (debug) { std::cerr << "My is: " << my << std::endl; }

	// get Pairing matrix for paths, TODO only initialize once for the whole program!
	Pairing p(max_length+1);
	
	// start at the outermost ear and process inwards
	for (boost::tie(ear, ear_end) = (g).children(); ear != ear_end; ++ear) {
		
		// before doing anything, update current Ak and Ai
		updateCurrentAkAi (*ear, k, currentAk, currentAi);
		
		// best is to remember the currentAis at this point
		Aks.push_back(currentAk);
		
		if (debug) { 
			auto printpairAk = std::make_pair(*ear, currentAk);
			auto printpairAi = std::make_pair(*ear, currentAi);
			std::cerr << "===========================================" << std::endl
				<< "Current k: " << k << std::endl
				<< "currentAk:" << std::endl << printpairAk << std::endl
				<< "currentAi:" << std::endl << printpairAi << std::endl;
		}
		
		// Nk[A6][A10][A1] = sum(AUGC in inner Ap = 9) P[A6][x9][3 pathlength] * P[x9][A10][1] * Nk-1 [x9][A1]
		
		std::vector<MyKey> key_combinations;		// this is what we want to fill in the recursion
		MyKey mykey;					// helper to recursively build the posibilities
		std::set<int> cAk;				// create a version with ints for recursion
		for (auto ap : currentAk) {
			cAk.insert(boost::get(boost::vertex_color_t(), g.root(), ap));
		}
		
		// now fill the key_combinations with all kinds of bases
		calculate_combinations(cAk, mykey, key_combinations);
		
		/*if (debug) {
			std::cerr << "Calculated following combinations of keys to calculate:" << std::endl;
			for (auto thiskey : key_combinations) {
				std::cerr << thiskey << std::endl;
			}
		}*/
		
		// calculate the probabilities for every key and if not zero, add to matrix
		for (auto thiskey : key_combinations) {
			if (debug) {
				std::cerr << "-------------------------------------------" << std::endl 
					<< "Calculating probablity for key: " << thiskey << std::endl;
			}
			
			unsigned long long probability = get_probability(thiskey, *ear, currentAk, currentAi, p, k);
			
			if (debug) { std::cerr << "= " << probability << std::endl; }
			
			if (probability != 0) {
				n[thiskey] = probability;
			}
		}
		
		// now going to next ear!
		k++;
		if (debug) {
			for (auto thiskey : key_combinations) {
				std::cerr << thiskey << " = " << get(thiskey) << std::endl;
			}
		}
	}
}

void ProbabilityMatrix::updateCurrentAkAi (Graph& g, int k, std::set<Vertex>& currentAk, std::set<Vertex>& currentAi) {
	// Ak and Ai are already stored in graph as a vertex properties
	// write vertex property into currentAk and currentAi
	currentAi.clear();
	
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		if (g[v].Ak.find(k) != g[v].Ak.end()) {
			currentAk.insert(g.local_to_global(v));
		} else if ((k > 0) && (g[v].Ai == k)) {
			currentAi.insert(g.local_to_global(v));
			// we need to keep Ak from previous glued ears, except those that became Ai this time!
			currentAk.erase(g.local_to_global(v));
		}
	}
	
	// for the last cycle, Ak is zero, so we need to "define" one Ai as Ak!
	if (currentAk.size() == 0) {
		auto first = currentAi.begin();
		// push the first Ai to Ak
		currentAk.insert(*first);
		// erase it in Ai
		currentAi.erase(first);
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
	
	
	// start to build lastkey
	for (auto v : Ak) {
		bool into_lastear = false;
		if (!g.find_vertex(v).second) {
			// are not in this particular ear
			into_lastear = true;
		} else {
			// check if were Aks before
			Vertex lv = g.global_to_local(v);
			for (unsigned int i = 0; i < k; i++) {
				if (g[lv].Ak.find(i) != g[lv].Ak.end()) {
					into_lastear = true;
				}
			}
		}
		
		if (into_lastear) {
			// these articulation points are not in this particular ear or were Aks before,
			// therefore we need to look their probability up
			// from last time. so lets generate a key therefore			
			lastkey.insert(std::make_pair(boost::get(boost::vertex_color_t(), g.root(), v), mykey[boost::get(boost::vertex_color_t(), g.root(), v)]));
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
	
	if (debug) {
		for (auto sub_probability : sub_probabilities) {
			std::cerr << "(" << sub_probability.start << ") --" 
				<< sub_probability.length << "-- ("
				<< sub_probability.end << ")\t*\t";
		}
		std::cerr << "N(" << lastkey << ")" << std::endl;
	}
	
	// calculate sum of sum for all bases colored X (= internal aps)
	// adds base combinatoric to the sub_probabilities and to the lastkey
	if (Ai.size() == 0) {
		max_number_of_sequences += calculate_probability (mykey, lastkey, sub_probabilities, p);
	} else {
		if (debug) { std::cerr << "Make sum of sum: " << k << std::endl; }
		make_sum_of_sum(g, Ai, mykey, lastkey, sub_probabilities, p, k, max_number_of_sequences);
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
			lastkey.insert(std::make_pair(boost::get(boost::vertex_color_t(), g.root(), v), b));
			// recursion starts here
			make_sum_of_sum(g, Ai, mykey, lastkey, sub_probabilities, p, k, max_number_of_sequences);
			
			if (Ai.size() == 0) {
				max_number_of_sequences += calculate_probability (mykey, lastkey, sub_probabilities, p);
			}
			// remove current vertex again to make space for a new base
			lastkey.erase(v);
		}
		// add current vertex again
		Ai.insert(v);
	}
}

int ProbabilityMatrix::get_color_from_key (MyKey& mykey, MyKey& lastkey, int vertex) {
	if (mykey.find(vertex) != mykey.end()) {
		return mykey[vertex];
	} else if (lastkey.find(vertex) != lastkey.end()) {
		return lastkey[vertex];
	} else {
		std::cerr << "Something went wrong! I could not find the color of the vertex in any key (neither last nor current): " 
			<< vertex << std::endl;
		exit(1);
	}
}

unsigned long long ProbabilityMatrix::calculate_probability (	MyKey& mykey, MyKey& lastkey, 
								std::vector<SubProbability>& sub_probabilities, 
								Pairing& p)
{
	// do the actual calculation here!
	unsigned long long multiplied_probabilities = 1;
	
	// calculate product of sub probabilities
	for (auto sub_probability : sub_probabilities) {
		int startBase = get_color_from_key(mykey, lastkey, sub_probability.start);
		int endBase = get_color_from_key(mykey, lastkey, sub_probability.end);
		// get probability and multiply it
		multiplied_probabilities *= p.get(sub_probability.length, startBase, endBase);
		
		if (debug) { std::cerr << "P(" << enum_to_char(startBase) << ", " 
			<< enum_to_char(endBase) << ", " << sub_probability.length << ") = " 
			<< p.get(sub_probability.length, startBase, endBase) << "\t*\t"; }
	}
	// for k = 0 we should not look up lastkey, as there is no lastkey!
	if (lastkey.size() != 0) {
		// now add probability for last_ear
		unsigned long long last_probability = get(lastkey);
		if (debug) { std::cerr << "N{" << lastkey << "} = "
				<< last_probability; }
	
		multiplied_probabilities *= last_probability;
	}
	
	if (debug) { std::cerr << "\t=\t" << multiplied_probabilities << std::endl; }
	return multiplied_probabilities;
}

void ProbabilityMatrix::calculate_combinations (std::set<int>& Ak, MyKey& mykey, std::vector<MyKey>& key_combinations) {
	
	if (Ak.size() > 0) {
		std::set<int>::iterator it = Ak.begin();
		int vertex = *it;
		Ak.erase(it);
		
		for ( unsigned int b = 0; b < A_Size; b++ ) {
			mykey.insert(std::make_pair(vertex, b));
			// recursion starts here
			calculate_combinations(Ak, mykey, key_combinations);
			
			if (Ak.size() == 0) {
				// remember our generated key
				key_combinations.push_back(mykey);
			}
			// remove current vertex again to make space for a new base
			mykey.erase(vertex);
		}
		// add current vertex again
		Ak.insert(vertex);
	}
}

unsigned long long ProbabilityMatrix::get (MyKey mykey) {

	unsigned long long returnvalue;
	
	// important for map: if you request with [] an entry will be created for unexisting ones.
	std::unordered_map < MyKey , unsigned long long , MyKeyHash>::const_iterator found = n.find(mykey);
	
	if (found != n.end()) {
		returnvalue = found->second;
	} else {
		returnvalue = 0;
	}
	
	return returnvalue;
}

unsigned long long ProbabilityMatrix::get_sum (MyKey mykey, std::vector<MyKey>& key_combinations) {
	MyKey tempkey;				// helper to build all combinations
	std::set<int> cAk;			// to get all combinations of keys we need a list of articulation points
	for (auto elem : mykey) {
		if (elem.second == X) {
			cAk.insert(elem.first);		// store to get all key combinations
			mykey.erase(elem.first);	// delete to in the end get a list of static key-elements
		}
	}
	// now lets calculate all combinations of keys for our points
	if (cAk.size() != 0) {
		calculate_combinations (cAk, tempkey, key_combinations);
	} else {
		key_combinations.push_back(MyKey());
	}
	
	// need to add the static key-elements in mykey to our keys in the key_combinations
	for (std::vector<MyKey>::iterator it = key_combinations.begin(); it != key_combinations.end(); ++it) {
		(*it).insert(mykey.begin(), mykey.end());
		//std::cerr << *it << " = " << get(*it) << std::endl;
	}
	
	// now get the sum of all probabililties
	unsigned long long sum = 0;
	for (auto thiskey : key_combinations) {
		sum += get(thiskey);
	}
	return sum;
}

std::set<Vertex> ProbabilityMatrix::get_Ak (unsigned int k) { 
	if (k < my) {
		return Aks[k];
	} else { 
		std::cerr << "k can't be bigger than my!" << std::endl;
		exit(1); 
	}
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
	// TODO Initialize just once!
	ProbabilityMatrix pm(g);
	
	// reverse iterate again over all ears to color Aks and all vertices in between
	for (unsigned int k = pm.get_my()-1; k > 0; --k) {
		if (debug) { std::cerr << "Start Backtracing at ear " << k << std::endl; }
		// get the current Articulation Points
		std::set<Vertex> Ak = pm.get_Ak(k);
		
		// translate the Ak set into set of ints (no vertex descriptors!)
		MyKey thiskey;
		for (auto ap : Ak) {
			int vertex = boost::get(boost::vertex_color_t(), g.root(), ap);
			int color = (g.root())[ap].base;
			thiskey.insert(std::make_pair(vertex, color));
		}
		// now do the random coloring of our points
		if (debug) { std::cerr << "Try to color this key: " << thiskey << std::endl; }
		MyKey colorkey = color_articulation_points(pm, thiskey);
		if (debug) { std::cerr << "Got a colored key: " << colorkey << std::endl; }
		
		// put colors onto graph
		for (auto v : Ak) {
			Vertex lv = g.global_to_local(v);
			if (g[lv].base == X) {
				g[lv].base = colorkey[boost::get(boost::vertex_color_t(), g, lv)];
			} else if (g[lv].base != colorkey[boost::get(boost::vertex_color_t(), g, lv)]) {
				std::cerr << "ERROR: Tried to change a already assigned base of an articulation point." << std::endl;
				exit(1);
			}
			if (debug) {
				std::cerr << "v " << v << ": " << enum_to_char(g[lv].base) << std::endl;
			}
		}
	}
	
	// now let's color all the vertices on the parts children
	Graph::children_iterator ear, ear_end;
	for (boost::tie(ear, ear_end) = g.children(); ear != ear_end; ++ear) {
		Graph::children_iterator part, part_end;
		for (boost::tie(part, part_end) = (*ear).children(); part != part_end; ++part) {
			
			if (debug) { std::cerr << "Coloring a part of the ear: ";
				BGL_FORALL_VERTICES_T(v, *part, Graph) {
					auto printpair = std::make_pair(*part, (*part).local_to_global(v));
					std::cerr << printpair << " ";
				}
				std::cerr << std::endl;
			}
			color_path_cycle_graph (*part);
		}
	}
}

MyKey color_articulation_points (ProbabilityMatrix& pm, MyKey& colorkey) {

	MyKey returnkey;
	
	// delare random number distribution and get a random number
	std::uniform_real_distribution<float> dist(0, 1);
	// get a random number between 0 and 1.
	float random = dist(rand_gen);
	if (debug) { std::cerr << "Got a random number: " << random << std::endl; }
	
	// now get all kombinations of keys and the sum of all possibilities
	std::vector<MyKey> key_combinations;	// this is what we want to fill now
	unsigned long long sum_of_possibilities = pm.get_sum(colorkey, key_combinations);
	if (debug) { std::cerr << "Sum of all possibilities is: " << sum_of_possibilities << std::endl; }
	
	// stochastically take one of the posibilities
	// start at the probability of first possible key and add each other key probability 
	// as long as the random number is bigger.
	unsigned long long sum = 0;
	for (auto thiskey : key_combinations) {
		sum += pm.get(thiskey);
		// if the random number is bigger than our probability, take this base as the first base!
		if (random*sum_of_possibilities < sum) {
			returnkey = thiskey;
			break;
		}
	}
	return returnkey;
}

void reset_colors(Graph& g) {
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		g[v].base = X;
	}
}

// overload << operator to print mykeys in pair with graph
std::ostream& operator<< (std::ostream& os, MyKey& m) {
	os << "[";
	for (typename MyKey::iterator it = m.begin(); it != m.end(); it++) {
        	os << "(" << std::setfill(' ') << std::setw(1) << it->first << "," << enum_to_char(it->second) << ")";
	}
	os << "]";
	return os;
}