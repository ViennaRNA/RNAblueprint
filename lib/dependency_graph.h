/* This program reads secundary RNA structures in dot-bracket and
 * builds a graph for a latter ear-decomposition and bipartitness-check
 *
 * Created on: 18.03.2014
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

#ifndef DEPENDENCY_GRAPH_H
#define	DEPENDENCY_GRAPH_H

#include "common.h"
#include "graphcoloring.h"

#include <sstream>
#include <random>

namespace design {
  namespace detail {
    
    enum GraphTypes {
      PATH, BLOCK
    };
    
    template <typename R>
    class GraphComponent {
    public:
      GraphComponent (Graph& g, int t, R * r);

      unsigned long long number_of_sequences () {
        return nos;
      }
      void mutate (int position);
      void mutate ();
      ~GraphComponent ();
    private:
      Graph& subgraph;
      unsigned long long nos = 0;
      ProbabilityMatrix * pm = NULL;
      int type;
      R * rand_ptr;
    };
    
    template <typename S, typename R>
    class DependencyGraph {
    public:
      DependencyGraph (S structures, R rand);
      unsigned long long number_of_sequences () {
        return nos;
      }

      bool is_bipartite () {
        return bipartite;
      }
      Sequence get_sequence ();
      std::string get_sequence_string ();
      void mutate (int position);
      void mutate ();
      void reset_colors ();
      ~DependencyGraph ();
    private:
      Graph graph;
      bool bipartite; // if dependency graph is bipartite and a therefore a solution exists
      unsigned long long nos = 0; // number of sequences/solutions
      std::list< GraphComponent<R>* > graph_components;
      R rand;
      R * rand_ptr;
    };
    
    // implementation
    template <typename S, typename R>
    DependencyGraph<S, R>::DependencyGraph (S structures, R _rand) {

      if (debug) {
        std::cerr << "Initializing DependencyGraph..." << std::endl;
      }
      
      // initialize random generator
      rand = _rand;
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
    template <typename S, typename R>
    Sequence DependencyGraph<S, R>::get_sequence () {
      Sequence sequence(boost::num_vertices(graph), N);

      BGL_FORALL_VERTICES_T(v, graph, Graph) {
        sequence[boost::get(boost::vertex_color_t(), graph, v)] = graph[v].base;
      }
      return sequence;
    }
    template <typename S, typename R>
    std::string DependencyGraph<S, R>::get_sequence_string () {
      Sequence sequence = get_sequence();

      std::stringstream stream;
      stream << sequence;

      return stream.str();
    }
    template <typename S, typename R>
    void DependencyGraph<S, R>::mutate () {
      // reset all the colors to N
      reset_colors();
      // TODO replace with a good mutation function
      for (auto component : graph_components) {
        component->mutate();
      }
    }
    template <typename S, typename R>
    void DependencyGraph<S, R>::mutate (int position) {
      // TODO replace with a good mutation function	
      this->mutate();
    }
    template <typename S, typename R>
    void DependencyGraph<S, R>::reset_colors () {

      BGL_FORALL_VERTICES_T(v, graph, Graph) {
        graph[v].base = N;
      }
    }
    template <typename S, typename R>
    DependencyGraph<S, R>::~DependencyGraph () {
      // TODO clean up all the data
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
  }
}
#endif	/* DEPENDENCY_GRAPH_H */

