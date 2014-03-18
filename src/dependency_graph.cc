/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 18.03.2014
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

// include header
#include "dependency_graph.h"

#include "common.h"
#include "parsestruct.h"
#include "printgraph.h"
#include "decompose.h"
#include "graphcoloring.h"

// include boost components
#include <boost/graph/iteration_macros.hpp>

DependencyGraph::DependencyGraph (std::vector<std::string> structures) {
	// generate graph from input vector
	graph = parse_structures(structures);
	
	// decompose the graph into its connected components, biconnected
	// components and decompose blocks via ear decomposition
	bipartite = decompose_graph(graph);
	
	// now fill the vector in right color-ordering
	// TODO
	
	// TODO initialize a fibonacci pairing matrix with the length max(num_vertices(g) of (bi?)connected components)
	// TODO initialice a probability_matrix for every block in the dependency graph and store it in a unodered_map< &Graph, pm >
	
	//TODO calculate nos (number of solutions) for the whole dependency graph
	
	// do the initial coloring
	// TODO replace color_graph with a function that uses the graph_components list.
	color_graph(graph);
	
}

Sequence DependencyGraph::get_sequence() {
	Sequence sequence(boost::num_vertices(graph), X);
	
	BGL_FORALL_VERTICES_T(v, graph, Graph) {
		sequence[boost::get(boost::vertex_color_t(), graph, v)] = graph[v].base;
	}
	return sequence;
}

Sequence DependencyGraph::mutate() {
	// TODO replace with a good mutation function
	color_graph(graph);
	
	return this->get_sequence();
}

Sequence DependencyGraph::mutate(int position) {
	// TODO replace with a good mutation function
	color_graph(graph);
	
	return this->get_sequence();
}

DependencyGraph::~DependencyGraph() {
    
}