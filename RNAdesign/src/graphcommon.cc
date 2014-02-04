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

void open_ear_decomposition (Graph& g, Vertex startVertex, ear_t& ear) {
	
	
	// map of ear decomposition properties for all vertices as key
	ear_propertymap_t p;
	
	// time starts at 0
	unsigned int counter = 0;

	if (debug) { std::cout << "StartVertex is: " << startVertex << std::endl; }
	// Algorithm from Ramachandran (1992) Parallel Open Ear Decomposition with Applications, page 8/9
	ear_dfs(g, startVertex, p, ear, counter);
	
	// print out all data-structures at the end
	if (debug) { 
		std::cerr << "index\tcolor\tporder\tparent\tlow\tear" << std::endl;
		for (ear_propertymap_t::iterator it = p.begin(); it != p.end(); it++) {
			std::cerr << it->first << "\t" << 
        		it->second.color << "\t" <<
        		it->second.preorder << "\t" <<
        		it->second.parent << "\t" <<
        		it->second.low << "\t" <<
        		it->second.ear.first << "," <<
        		it->second.ear.second << std::endl;
		}
		std::cerr << "index\tear" << std::endl;
		for (ear_t::iterator it = ear.begin(); it != ear.end(); it++) {
			std::cerr << it->first.first << "," <<
			it->first.second << "\t" <<
			it->second.first << "," <<
			it->second.second << std::endl;
		}
		std::cerr << "counter: " << counter << std::endl;	
	}
}

void ear_dfs(Graph& g, Vertex v, ear_propertymap_t& p, ear_t& ear, unsigned int& counter) {
	
	enum { WHITE, BLACK, GRAY };
	if (debug) { std::cout << "v is: " << v << std::endl; }
	
	// start ear decomposition
	p[v].color = GRAY;
	p[v].preorder = counter;
	counter++;
	p[v].low = boost::num_vertices(g);
	p[v].ear = std::make_pair(boost::num_vertices(g), boost::num_vertices(g));
	
	// get neighbouring vertices
	typename Graph::out_edge_iterator ei, ei_end;
	for (boost::tie(ei, ei_end) = boost::out_edges(v, g);  ei != ei_end; ++ei)
	{
		if (debug) { std::cerr << boost::target(*ei, g) <<" is neighbour through edge: " << *ei << std::endl; }
		Vertex w = boost::target(*ei, g);
		if (debug) { std::cout << "w is: " << w << std::endl; }
		
		if (p[w].color == WHITE) {
			if (debug) { std::cout << "w is white" << std::endl; }
			p[w].parent = v;
			// start new iteration here
			ear_dfs(g, w, p, ear, counter);
			//TODO: cast low vertex to integer a good idea?
			if ((int) p[w].low >= p[w].preorder) {
				ear[std::make_pair(p[w].parent, w)] = std::make_pair(std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
			} else {
				ear[std::make_pair(p[w].parent, w)] = p[w].ear;
			}
				p[v].low = std::min((int) p[v].low, (int) p[w].low);
				p[v].ear = lexmin(p[v].ear, p[w].ear);
		} else if (p[w].color == GRAY) {
			if (debug) { std::cout << "w is gray" << std::endl; }
			if (w != p[w].parent) {
				if (debug) { std::cout << "found a crossedge: " << v << w << std::endl; }
				//TODO: casting vertex in low to integer a bad idea?
				p[v].low = boost::vertex(std::min((int) p[v].low, p[w].preorder), g);
				ear[std::make_pair(w, v)] = std::make_pair(boost::vertex(p[w].preorder, g), boost::vertex(p[v].preorder, g));
				p[v].ear = lexmin(p[v].ear, ear[std::make_pair(w, v)]);
			}
		}
	}
	if (debug) { std::cout << "finishing vertex " << v << std::endl; }
}

