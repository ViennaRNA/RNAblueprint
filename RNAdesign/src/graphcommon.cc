/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 26.08.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

// include header
#include "graphcommon.h"

std::pair <int, int> get_min_max_degree(Graph& g) {
	int max_degree = 0;
	int min_degree = std::numeric_limits<int>::max();
	
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		int current_degree = boost::out_degree(v, g);
		if (current_degree > max_degree) {
			max_degree = current_degree;
		} else if (current_degree < min_degree) {
			min_degree = current_degree;
		}
	}
	return std::make_pair(min_degree, max_degree);
}

bool is_bipartite_graph(Graph& g, Vertex startVertex, Edge& ed) {
	// This is a Breadth First Search which checks if the graph is bipartit. 
	// If not, returns false and the fills the conflicting edge into the edge_descriptor
	
	if (debug) { 	std::cerr << "StartVertex is: " << startVertex << std::endl; 
			std::cerr << "Number of vertices: " << boost::num_vertices(g) << std::endl; }
	
	// exit value (if bipartite = true, else false)
	bool exit = true;
	// Define A BGL visitor for the BFS algorithm
	class my_bfs_visitor : public boost::default_bfs_visitor {
		public:
		my_bfs_visitor(Edge& ed, bool& exit) : m_ed(ed), m_exit(exit) {}
		Edge& m_ed;
		bool& m_exit;
		enum { WHITE, BLACK, GRAY, RED };
		void tree_edge(Edge e, Graph g) const {
			if (debug) { std::cout << "Detecting Tree edge: " << e << std::endl; }
			Vertex u = boost::source(e, g);
			Vertex v = boost::target(e, g);
			if (g[u].bipartite_color == RED) {
				g[v].bipartite_color = BLACK;
			} else {
				g[v].bipartite_color = RED;
			}
		}
		void non_tree_edge(Edge e, Graph g) const {
			if (debug) { std::cout << "Detecting Non-Tree edge: " << e << std::endl; }
			Vertex u = boost::source(e, g);
			Vertex v = boost::target(e, g);
			if (g[u].bipartite_color == g[v].bipartite_color) {
				if (debug) { std::cerr << "u and v have the same color -> not bipartite!" << std::endl; }
				m_ed = boost::edge(u,v,g).first;
				// return false if graph is not bipartite
				m_exit = false;
			}
		}
	};
	
	my_bfs_visitor vis(ed, exit);
	// Do a BGL BFS!
	// http://www.boost.org/doc/libs/1_53_0/libs/graph/doc/breadth_first_search.html
	boost::breadth_first_search(g, startVertex, boost::visitor(vis));
	return exit;
}



