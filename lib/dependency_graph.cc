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

namespace design {
  namespace detail {
    template <typename R>
    DependencyGraph<R>::DependencyGraph (std::vector<std::string> structures, R _rand)
    : rand(_rand) {

      if (debug) {
        std::cerr << "Initializing DependencyGraph..." << std::endl;
      }
      
      // initialize random generator
      rand_ptr = &rand;

      // generate graph from input vector
      graph = parse_structures(structures);

      // decompose the graph into its connected components, biconnected
      // components and decompose blocks via ear decomposition
      bipartite = decompose_graph(graph, rand_ptr);

      // now fill the vector in right coloring-order and self initialize pm, nos, etc.
      Graph::children_iterator cc, cc_end;
      for (boost::tie(cc, cc_end) = graph.children(); cc != cc_end; ++cc) {
        // connected components
        if (get_min_max_degree(*cc).second <= 2) {
          // found connected components here (single point or path)
          graph_components.push_back(new GraphComponent<R>(*cc, PATH, rand_ptr));

        } else {
          Graph::children_iterator bc, bc_end;
          for (boost::tie(bc, bc_end) = (*cc).children(); bc != bc_end; ++bc) {
            auto min_max_degree = get_min_max_degree(*bc);
            // biconnected components (color blocks first!)
            if (min_max_degree.second > 2) {
              // found a block
              graph_components.push_back(new GraphComponent<R>(*bc, BLOCK, rand_ptr));
            } else if (min_max_degree.first == 2 && min_max_degree.second == 2) {
              // found biconnected cycles here
              graph_components.push_back(new GraphComponent<R>(*bc, PATH, rand_ptr));
            }
          }
          // iterate again to find remaining paths
          for (boost::tie(bc, bc_end) = (*cc).children(); bc != bc_end; ++bc) {
            auto min_max_degree = get_min_max_degree(*bc);

            if (min_max_degree.first == 1 && min_max_degree.second == 2) {
              // found biconnected component path
              graph_components.push_back(new GraphComponent<R>(*bc, PATH, rand_ptr));
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
    template <typename R>
    Sequence DependencyGraph<R>::get_sequence () {
      Sequence sequence(boost::num_vertices(graph), N);

      BGL_FORALL_VERTICES_T(v, graph, Graph) {
        sequence[boost::get(boost::vertex_color_t(), graph, v)] = graph[v].base;
      }
      return sequence;
    }
    template <typename R>
    std::string DependencyGraph<R>::get_sequence_string () {
      Sequence sequence = get_sequence();

      std::stringstream stream;
      stream << sequence;

      return stream.str();
    }
    template <typename R>
    void DependencyGraph<R>::mutate () {
      // reset all the colors to N
      reset_colors();
      // TODO replace with a good mutation function
      for (auto component : graph_components) {
        component->mutate();
      }
    }
    template <typename R>
    void DependencyGraph<R>::mutate (int position) {
      // TODO replace with a good mutation function	
      this->mutate();
    }
    template <typename R>
    void DependencyGraph<R>::reset_colors () {

      BGL_FORALL_VERTICES_T(v, graph, Graph) {
        graph[v].base = N;
      }
    }
    template <typename R>
    DependencyGraph<R>::~DependencyGraph () {
      // TODO clean up all the data
      for (auto component : graph_components) {
        delete component;
      }
    }
    template <typename R>
    GraphComponent<R>::GraphComponent (Graph& g, int t, R * r)
    : subgraph (g), type (t), rand_ptr(r) {

      if (debug) {
        std::cerr << "Initializing GraphComponent..." << std::endl;
      }
      // initialize colors here
      if (type == PATH) {
        nos = color_path_cycle_graph(subgraph, rand_ptr);
      } else if (type == BLOCK) {
        pm = new ProbabilityMatrix(subgraph);
        nos = pm->number_of_sequences();
        // initialize coloring
        color_blocks(subgraph, *pm, rand_ptr);
      }

      if (debug) {
        print_all_vertex_names(subgraph, "Initialized a GraphComponent");
      }

    }
    template <typename R>
    void GraphComponent<R>::mutate () {
      // initialize color function here
      if (type == PATH) {
        color_path_cycle_graph(subgraph, rand_ptr);
      } else if (type == BLOCK) {
        color_blocks(subgraph, *pm, rand_ptr);
      }

    }
    template <typename R>
    void GraphComponent<R>::mutate (int position) {
      // TODO
    }
    template <typename R>
    GraphComponent<R>::~GraphComponent () {
      // clean up
      delete pm;
    }
    
    template class DependencyGraph<std::mt19937>;
    template class GraphComponent<std::mt19937>;
  }
}