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

namespace design {
  namespace detail {
    
    enum GraphTypes {
      PATH, BLOCK
    };

    class GraphComponent {
    public:
      GraphComponent (Graph& g, int t);

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
    };

    class DependencyGraph {
    public:
      DependencyGraph (std::vector<std::string> structures);

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
      std::list< GraphComponent* > graph_components;


    };
  }
}
#endif	/* DEPENDENCY_GRAPH_H */

