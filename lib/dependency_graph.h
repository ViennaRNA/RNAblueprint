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
#include "pairing_matrix.h"
#include "probability_matrix.h"
#include "graphcommon.h"
#include "pathcoloring.h"

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
    
    template <typename R>
    class DependencyGraph {
    public:
      DependencyGraph (std::vector<std::string> structures, R rand);
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
      R * rand_ptr;
      ~DependencyGraph ();
    private:
      Graph graph;
      bool bipartite; // if dependency graph is bipartite and a therefore a solution exists
      unsigned long long nos = 0; // number of sequences/solutions
      std::list< GraphComponent<R>* > graph_components;
      R rand;
    };
  }
}
#endif	/* DEPENDENCY_GRAPH_H */

