/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 25.03.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
* Compile with: g++ -std=c++11 -g -lboost_program_options -o struct2graph struct2graph.cc
*/

// include header
#include "treestatistics.h"



void do_spanning_tree_stat (Graph& g) {
	// random generator to make spanning tree sampling
	if (seed == 0) {
		unsigned long clock_seed = std::chrono::system_clock::now().time_since_epoch().count();
		seed = clock_seed;
	}
	std::mt19937 r (seed);  // mt19937 is a standard mersenne_twister_engine
	std::cerr << "Using this seed: " << seed << std::endl;
	
	// start at all vertices of the subgraph as root of the tree
	BGL_FORALL_VERTICES_T(root, g, Graph) {	
		// generate num_trees spanning trees for statistics
		for (int i = 1; i != num_trees+1; i++) {
			// remember predescessor map and all non-tree edges
			std::map<Vertex, Vertex> parents;
			std::vector<Edge> crossedges;
		
			// get a boost random spanning tree
			get_random_spanning_tree (g, r, parents, crossedges, root);
		
			// print parents, cross-edges and root vertex
			if (verbose) {
				std::cerr << "Root vertex: " << root << std::endl;
				std::cerr << "Spanning tree (vertex, parent) and cross-edges:" << std::endl;
				for (std::map<Vertex, Vertex>::iterator it=parents.begin(); it!=parents.end(); ++it) {
					std::cerr << it->first << " => " << it->second << std::endl;
				}
				for (auto elem : crossedges) {
					std::cerr << elem << std::endl;
				}
			}
			
			// do the actual ear decomposition
			ear_decomposition (g, parents, crossedges, root);
			
			// calculate the two performance critical variables alpha and beta
			// store attachment vertices per each step of the ear decomposition in Ak
			std::map<int, std::vector<Vertex> > Ak;
			unsigned int alpha;
			unsigned int beta;
			std::tie(alpha, beta) = calculate_alpha_beta(g, crossedges, Ak);
			
			// write statistics to a file
			print_ab_stat (alpha, beta, Ak, g, root, crossedges);
		}
	}
}

std::pair< unsigned int, unsigned int > calculate_alpha_beta(Graph& g, std::vector<Edge>& crossedges, std::map<int, std::vector<Vertex> >& Ak) {
	
	// structure to remember Ak (attachment vertices)
	//std::map<int, std::vector<Vertex> > Ak;
	unsigned int alpha = 0;
	unsigned int beta = 0;
	int my = crossedges.size();
	
	// reset_color
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		g[v].color = 0;
	}
	
	// iterate over all ear decomposition iterations and inside this loop over all edges of this ear
	// this edges have source and target vertices, iterate over both
	// then iterate over the adjacent out edges of these vertices and get the hightes ear number into maxear
	for (int k = 0; k != my; k++) {
		BGL_FORALL_EDGES_T(e, g, Graph) {
			if (g[e].ear == k) {
				std::vector<Vertex> adjacent_v;
				adjacent_v.push_back(boost::source(e, g));
				adjacent_v.push_back(boost::target(e, g));
				//TODO we calculate maxear twice for many vertices...
				for (auto adja : adjacent_v) {
					int maxear = 0;
					typename Graph::out_edge_iterator ei, ei_end;
					for (boost::tie(ei, ei_end) = boost::out_edges(adja, g);  ei != ei_end; ++ei) {
						if(g[*ei].ear > maxear) { maxear = g[*ei].ear; }
					}
				
					if (maxear > k) { g[adja].color = 1; }
					else { g[adja].color = 0; }
				}
			}
		}
		
		//print_graph(g, &std::cout, "test_color");
		// remember Ak for this k
		std::vector< Vertex > thisAk;
		// write colored vertices into thisAk
		BGL_FORALL_VERTICES_T(v, g, Graph) {
			if (g[v].color == 1) {
				thisAk.push_back(v);
			}
		}
		Ak[k] = thisAk;
		if (thisAk.size() > alpha) { alpha = thisAk.size(); }
		
		if (k > 0) {
			std::vector<Vertex> Akplus1_without_Ak;
			for (auto elem : thisAk) {
				std::vector<Vertex>::iterator iter = find(Ak[k-1].begin(), Ak[k-1].end(), elem);
				if (iter == Ak[k-1].end()) {
					// elem is in thisAk but not in Ak[k-1]
					Akplus1_without_Ak.push_back(elem);
				}
			}
			unsigned int thisbeta = Ak[k-1].size() + Akplus1_without_Ak.size();
			if (beta < thisbeta) { beta = thisbeta; }
		}
	}
	return std::make_pair(alpha, beta);
}	

void print_ab_stat (unsigned int alpha, unsigned int beta, std::map<int, std::vector<Vertex> > Ak, Graph& g, Vertex root, std::vector<Edge>& crossedges) {
	
	int my = crossedges.size();
	
	// write statistic output file
	std::stringstream filename;
	if (outfile != "") {
		filename << outfile << "-statistics.txt";
	} else { filename << "statistics.txt"; }
	std::ofstream statfile(filename.str(), std::ofstream::out | std::ofstream::app);
	if (statfile.is_open()) {
		// print alpha and beta
		statfile << alpha << " " << beta << " ";
		// print Ak sets
		for (int i = 0; i != my; i++) {
			for (auto elem : Ak[i]) {
				if(elem != Ak[i].back()) {
					statfile << boost::get(boost::vertex_color_t(), g, elem) << ",";
				} else {
					statfile << boost::get(boost::vertex_color_t(), g, elem);
				}
			}
			statfile << " ";
		}
		// print root of spanning tree
		statfile << root << " ";
		// also print crossedges to reproduce trees afterwards
		for (auto elem : crossedges) {
			statfile << elem;
		}
		statfile << std::endl;
		statfile.close();
		if (verbose) { std::cerr << "Statistics written to outfile!" << std::endl; }
	} else {
			std::cerr << " Unable to create statistics file!" << std::endl;
	}
}
