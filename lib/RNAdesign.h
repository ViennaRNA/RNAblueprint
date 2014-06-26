/* This program reads secundary RNA structures in dot-bracket and
 * builds a graph for a latter ear-decomposition and bipartitness-check.
 *
 * Created on: 26.06.2014
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */


#ifndef RNADESIGN_H
#define RNADESIGN_H

// include header
#include "common.h"
#include "dependency_graph.h"
/*
 * RNA design library
 * @autor Stefan Hammer
 */
namespace design {
  
  /*
   * Dependency Graph which holds all structural constraints
   * This graph is used to generate valid sequences compatible to the input structures
   */
  class DependencyGraph {
  public:
    DependencyGraph(std::vector<std::string> structures);
    ~DependencyGraph();
    bool is_bipartite();
    std::string get_sequence();
    void mutate ();
    void mutate (int position);
    unsigned long long number_of_sequences ();
  private:
    detail::DependencyGraph * g;
  };
  
}

#endif