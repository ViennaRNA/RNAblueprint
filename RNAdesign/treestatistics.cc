/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 25.03.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

// include header
#include "treestatistics.h"



void do_spanning_tree_stat (Graph& g, int num_trees) {
		
	// start at all vertices of the subgraph as root of the tree
	BGL_FORALL_VERTICES_T(root, g, Graph) {	
		// generate num_trees spanning trees for statistics
		for (int i = 1; i != num_trees+1; i++) {
			// remember predescessor map and all non-tree edges
			std::map<Vertex, Vertex> parents;
			std::vector<Edge> crossedges;
		
			// get a boost random spanning tree
			get_random_spanning_tree (g, parents, crossedges, root);
		
			// print parents, cross-edges and root vertex
			if (debug) {
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
			
			// detect Articulation Points and push them into the graph as vertex property Ak
			color_Ak_points (g);
			
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

void color_Ak_points (Graph& g) {
	// iterate over all vertices and compare the ear numbers of its out-edges.
	// push the right ear number into the vertex property Ak vector
	std::vector< int > earnumbers;
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		// clear previous stored Articulation Points
		g[v].Ak.clear();
		// reset earnumbers
		earnumbers.clear();
		BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
			// compare this new earnumber to all the others
			for (auto numb : earnumbers) {
				if (g[e].ear < numb) {
					g[v].Ak.insert(g[e].ear);
				} else if (g[e].ear > numb) {
					g[v].Ak.insert(numb);
					earnumbers.push_back(g[e].ear);
				}
			}
			if (earnumbers.size() == 0)			// always save the first outedge-earnumber
				earnumbers.push_back(g[e].ear);		// save the current earnumber to the others
		}
		
	}
}

std::pair< unsigned int, unsigned int > calculate_alpha_beta(Graph& g, std::vector<Edge>& crossedges, std::map<int, std::vector<Vertex> >& Ak) {
	
	// structure to remember Ak (attachment vertices)
	//std::map<int, std::vector<Vertex> > Ak;
	unsigned int alpha = 0;
	unsigned int beta = 0;
	int my = crossedges.size();
	
	// iterate over all ear decomposition iterations 
	for (int k = 0; k != my; k++) {
		// Ak are already stored in graph as a vertex propertys
		// write vertex property into Ak[k]
		BGL_FORALL_VERTICES_T(v, g, Graph) {
			if (g[v].Ak.find(k) != g[v].Ak.end()) {
				Ak[k].push_back(v);
			}
		}
		// always store biggest alpha
		if (Ak[k].size() > alpha) { alpha = Ak[k].size(); }
		// calculate beta here
		if (k > 0) {
			std::vector<Vertex> Akplus1_without_Ak;
			for (auto elem : Ak[k]) {
				std::vector<Vertex>::iterator iter = find(Ak[k-1].begin(), Ak[k-1].end(), elem);
				if (iter == Ak[k-1].end()) {
					// elem is in Ak[k] but not in Ak[k-1]
					Akplus1_without_Ak.push_back(elem);
				}
			}
			unsigned int thisbeta = Ak[k-1].size() + Akplus1_without_Ak.size();
			// alywas store biggest beta
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
		if (debug) { std::cerr << "Statistics written to outfile!" << std::endl; }
	} else {
			std::cerr << " Unable to create statistics file!" << std::endl;
	}
}
