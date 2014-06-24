/* RNAdesign
 * A program for designing RNA molecules.
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

  if (debug) {
    std::cerr << "Initializing DependencyGraph..." << std::endl;
  }

  // generate graph from input vector
  graph = parse_structures(structures);

  // decompose the graph into its connected components, biconnected
  // components and decompose blocks via ear decomposition
  bipartite = decompose_graph(graph);

  // now fill the vector in right color-ordering and self initialize pm, nos, etc.
  Graph::children_iterator cc, cc_end;
  for (boost::tie(cc, cc_end) = graph.children(); cc != cc_end; ++cc) {
    // connected components
    if (get_min_max_degree(*cc).second <= 2) {
      // found connected components here (single point or path)
      graph_components.push_back(new GraphComponent(*cc, PATH));

    } else {
      Graph::children_iterator bc, bc_end;
      for (boost::tie(bc, bc_end) = (*cc).children(); bc != bc_end; ++bc) {
        auto min_max_degree = get_min_max_degree(*bc);
        // biconnected components (color blocks first!)
        if (min_max_degree.second > 2) {
          // found a block
          graph_components.push_back(new GraphComponent(*bc, BLOCK));
        } else if (min_max_degree.first == 2 && min_max_degree.second == 2) {
          // found biconnected cycles here
          graph_components.push_back(new GraphComponent(*bc, PATH));
        }
      }

      for (boost::tie(bc, bc_end) = (*cc).children(); bc != bc_end; ++bc) {
        auto min_max_degree = get_min_max_degree(*bc);

        if (min_max_degree.first == 1 && min_max_degree.second == 2) {
          // found biconnected component path
          graph_components.push_back(new GraphComponent(*bc, PATH));
        }
      }
    }
  }

  // TODO initialize a fibonacci pairing matrix with the length max(num_vertices(g) of (bi?)connected components)
  if (debug) {
    std::cerr << "Calculate number of sequences for dependency graph..." << std::endl;
  }
  //calculate nos (number of solutions) for the whole dependency graph
  for (auto component : graph_components) {
    nos += component->number_of_sequences();
  }
}

Sequence DependencyGraph::get_sequence () {
  Sequence sequence(boost::num_vertices(graph), X);

  BGL_FORALL_VERTICES_T(v, graph, Graph) {
    sequence[boost::get(boost::vertex_color_t(), graph, v)] = graph[v].base;
  }
  return sequence;
}

std::string DependencyGraph::get_sequence_string () {
  Sequence sequence = get_sequence();

  std::stringstream stream;
  stream << sequence;

  return stream.str();
}

Sequence DependencyGraph::mutate () {
  // reset all the colors to X
  reset_colors();
  // TODO replace with a good mutation function
  for (auto component : graph_components) {
    component->mutate();
  }

  return this->get_sequence();
}

Sequence DependencyGraph::mutate (int position) {
  // TODO replace with a good mutation function	
  return this->get_sequence();
}

void DependencyGraph::reset_colors () {

  BGL_FORALL_VERTICES_T(v, graph, Graph) {
    graph[v].base = X;
  }
}

DependencyGraph::~DependencyGraph () {
  // TODO clean up all the data
}

GraphComponent::GraphComponent (Graph& g, int t)
: subgraph (g), type (t) {

  if (debug) {
    std::cerr << "Initializing GraphComponent..." << std::endl;
  }
  // initialize colors here
  if (type == PATH) {
    nos = color_path_cycle_graph(subgraph);
  } else if (type == BLOCK) {
    pm = new ProbabilityMatrix(subgraph);
    nos = pm->number_of_sequences();
    // initialize coloring
    color_blocks(subgraph, *pm);
  }

  if (debug) {
    print_all_vertex_names(subgraph, "Initialized a GraphComponent");
  }

}

void GraphComponent::mutate () {
  // initialize color function here
  if (type == PATH) {
    color_path_cycle_graph(subgraph);
  } else if (type == BLOCK) {
    color_blocks(subgraph, *pm);
  }

}

void GraphComponent::mutate (int position) {
  // TODO
}

GraphComponent::~GraphComponent () {
  // clean up
  delete pm;
}