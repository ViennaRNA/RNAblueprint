/* This program reads secundary RNA structures in dot-bracket and
 * builds a graph for a latter ear-decomposition and bipartitness-check.
 *
 * Created on: 26.06.2014
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

#include "RNAdesign.h"

namespace design {
  
  DependencyGraph::DependencyGraph(std::vector<std::string> structures) {
    g = new detail::DependencyGraph(structures);
  }
  
  DependencyGraph::~DependencyGraph() {
    delete g;
  }
    
  bool DependencyGraph::is_bipartite() {
    return g->is_bipartite ();
  }
  
  std::string DependencyGraph::get_sequence() {
    return g->get_sequence_string();
  }
  
  void DependencyGraph::mutate () {
    g->mutate ();
  }
  
  void DependencyGraph::mutate (int position) {
    g->mutate (position);
  }
  
  unsigned long long DependencyGraph::number_of_sequences () {
    return g->number_of_sequences ();
  }
}