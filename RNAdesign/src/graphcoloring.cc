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

void color_graph (Graph& graph) {
	// INPUT: root graph
	// reset bases to X
	reset_colors(graph);
	
	Graph::children_iterator cc, cc_end;
	for (boost::tie(cc, cc_end) = graph.children(); cc != cc_end; ++cc) {
		// connected components
		if (get_min_max_degree(*cc).second <= 2) {
			// color connected components here
			print_vertices(*cc, "Coloring a connected component");
			color_path_cycle_graph (*cc);
                        
		} else {
			Graph::children_iterator bc, bc_end;
			for (boost::tie(bc, bc_end) = (*cc).children(); bc != bc_end; ++bc) {
				auto min_max_degree = get_min_max_degree(*bc);
				// biconnected components (color blocks first!)
				if (min_max_degree.second > 2) {
					// blocks
					// color blocks here
					print_vertices(*bc, "Coloring a block");
					color_blocks(*bc);
				} else if (min_max_degree.first == 2 && min_max_degree.second == 2) {
					// color biconnected cycles here
					print_vertices(*bc, "Coloring a biconnected cycle");
					color_path_cycle_graph (*bc);
				}
			}
			
			for (boost::tie(bc, bc_end) = (*cc).children(); bc != bc_end; ++bc) {
				auto min_max_degree = get_min_max_degree(*bc);
				
				if (min_max_degree.first == 1 && min_max_degree.second == 2) {
					// biconnected component paths
					// color paths here
					print_vertices(*bc, "Coloring a biconnected path");
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
	// remember the current key for the next ear.
	MyKey lastkey;
	
	// reverse iterate again over all ears to color Aks and all vertices in between
	for (int k = pm.get_my(); k > 0; k--) {
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
		MyKey colorkey = color_articulation_points(k, pm, thiskey, lastkey);
		// remember colorkey for next ear iteration
		lastkey = colorkey;
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
			
			print_vertices(*part, "Coloring a part of the ear");
			color_path_cycle_graph (*part);
		}
	}
}

MyKey color_articulation_points (int k, ProbabilityMatrix& pm, MyKey& colorkey, MyKey& lastkey) {

	MyKey returnkey;
	
	// declare random number distribution and get a random number
	std::uniform_real_distribution<float> dist(0, 1);
	// get a random number between 0 and 1.
	float random = dist(rand_gen);
	if (debug) { std::cerr << "Got a random number: " << random << std::endl; }
	
	// now get all combinations of keys and the sum of all possibilities
	std::vector<MyKey> key_combinations;	// this is what we want to fill now
	std::unordered_map < MyKey , unsigned long long , MyKeyHash> probabilities;
	unsigned long long sum_of_possibilities = pm.get_sum(k, colorkey, lastkey, key_combinations, probabilities);
	if (debug) { std::cerr << "Sum of all possibilities is: " << sum_of_possibilities << std::endl; }
	
	// stochastically take one of the posibilities
	// start at the probability of first possible key and add each other key probability 
	// as long as the random number is bigger.
	unsigned long long sum = 0;
	for (auto thiskey : key_combinations) {
		sum += probabilities[thiskey];
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
