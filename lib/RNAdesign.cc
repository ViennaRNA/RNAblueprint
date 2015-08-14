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

    void initialize_library(bool debug) {
        *detail::debug_ptr = debug;
    }

    template <typename R>
    DependencyGraph<R>::DependencyGraph(std::vector<std::string> structures, std::string constraints, R rand) {
        g = new detail::DependencyGraph<R>(structures, constraints, rand);
    }

    template <typename R>
    DependencyGraph<R>::~DependencyGraph() {
        delete g;
    }

    template <typename R>
    std::string DependencyGraph<R>::get_sequence() {
        return g->get_sequence_string();
    }

    template <typename R>
    void DependencyGraph<R>::mutate() {
        g->mutate();
    }

    template <typename R>
    void DependencyGraph<R>::mutate(int position) {
        g->mutate(position);
    }

    template <typename R>
    unsigned long long DependencyGraph<R>::number_of_sequences() {
        return g->number_of_sequences();
    }

    template class DependencyGraph<std::mt19937>;

}