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
   * Initialize the Library
   * Set the debug boolean to true if you want to get verbose output to std:err
   * Set the random generator for stochastic processes in the library
   */
  void initialize_library(bool debug);
  
  /*
   * Dependency Graph which holds all structural constraints
   * This graph is used to generate valid sequences compatible to the input structures
   */
  template <typename S, typename R>
  class DependencyGraph {
  public:
    DependencyGraph(S structures, R rand);
    ~DependencyGraph();
    bool is_bipartite();
    std::string get_sequence();
    void mutate ();
    void mutate (int position);
    unsigned long long number_of_sequences ();
  private:
    detail::DependencyGraph<S, R> * g;
  };
  
  // implementation
  void initialize_library(bool debug) {
    *detail::debug_ptr = debug;
  }
  
  template <typename S, typename R>
  DependencyGraph<S, R>::DependencyGraph(S structures, R rand) {
    g = new detail::DependencyGraph<S, R>(structures, rand);
  }
  template <typename S, typename R>
  DependencyGraph<S, R>::~DependencyGraph() {
    delete g;
  }
  template <typename S, typename R>
  bool DependencyGraph<S, R>::is_bipartite() {
    return g->is_bipartite ();
  }
  template <typename S, typename R>
  std::string DependencyGraph<S, R>::get_sequence() {
    return g->get_sequence_string();
  }
  template <typename S, typename R>
  void DependencyGraph<S, R>::mutate () {
    g->mutate ();
  }
  template <typename S, typename R>
  void DependencyGraph<S, R>::mutate (int position) {
    g->mutate (position);
  }
  template <typename S, typename R>
  unsigned long long DependencyGraph<S, R>::number_of_sequences () {
    return g->number_of_sequences ();
  }
}

#endif