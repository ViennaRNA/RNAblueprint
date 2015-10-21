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
    DependencyGraph<R>::DependencyGraph(std::vector<std::string> structures, std::string constraints, unsigned long seed) {
        // initialize mersenne twister with the given seed
        std::mt19937 rand_gen;
        rand_gen.seed(seed);
        if (*detail::debug_ptr) {
            std::cerr << "Using this seed: " << seed << std::endl;
        }
        g = new detail::DependencyGraph<R>(structures, constraints, rand_gen);
    }
    
    template <typename R>
    DependencyGraph<R>::DependencyGraph(std::vector<std::string> structures, std::string constraints) {
        // initialize mersenne twister with our seed
        unsigned long seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::mt19937 rand_gen;
        rand_gen.seed(seed);
        if (*detail::debug_ptr) {
            std::cerr << "Using this seed: " << seed << std::endl;
        }
        g = new detail::DependencyGraph<R>(structures, constraints, rand_gen);
    }
    
    template <typename R>
    DependencyGraph<R>::DependencyGraph(std::vector<std::string> structures, R rand) {
        g = new detail::DependencyGraph<R>(structures, "", rand);
    }
    
    template <typename R>
    DependencyGraph<R>::DependencyGraph(std::vector<std::string> structures) {
        // initialize mersenne twister with our seed
        unsigned long seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::mt19937 rand_gen;
        rand_gen.seed(seed);
        if (*detail::debug_ptr) {
            std::cerr << "Using this seed: " << seed << std::endl;
        }
        g = new detail::DependencyGraph<R>(structures, "", rand_gen);
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
    void DependencyGraph<R>::set_sequence(std::string sequence) {
        g->set_sequence_string(sequence);
    }

    template <typename R>
    void DependencyGraph<R>::set_sequence() {
        g->set_sequence();
    }
    
    template <typename R>
    SolutionSizeType DependencyGraph<R>::mutate_local(int min_num_pos, int max_num_pos) {
        // if local_global=-1 we will mutate only actual graphs (independent of their type)
        return g->mutate_local_global(-1, min_num_pos, max_num_pos);
    }
    
    template <typename R>
    SolutionSizeType DependencyGraph<R>::mutate_local() {
        // if local_global=-1 we will mutate only actual graphs, with 0,0 we choose any
        return g->mutate_local_global(-1, 0, 0);
    }
    
    template <typename R>
    SolutionSizeType DependencyGraph<R>::mutate_global(int min_num_pos, int max_num_pos) {
        // if local_global=1 we will mutate only connected components
        return g->mutate_local_global(1, min_num_pos, max_num_pos);
    }
    
    template <typename R>
    SolutionSizeType DependencyGraph<R>::mutate_global(int connected_component_ID) {
        // mutate the connected component with the ID
        return g->mutate_global(connected_component_ID);
    }
    
    template <typename R>
    SolutionSizeType DependencyGraph<R>::mutate_global() {
        // mutate any component with any size
        return g->mutate_local_global(1, 0, 0);
    }
    
    template <typename R>
    SolutionSizeType DependencyGraph<R>::mutate(int position) {
        return g->mutate(position);
    }
    
    template <typename R>
    SolutionSizeType DependencyGraph<R>::mutate(int start, int end) {
        return g->mutate(start, end);
    }

    template <typename R>
    SolutionSizeType DependencyGraph<R>::number_of_sequences() {
        return g->number_of_sequences();
    }
    
    template <typename R>
    SolutionSizeType DependencyGraph<R>::number_of_sequences(int connected_component_ID) {
        return g->number_of_sequences(connected_component_ID);
    }
    
    template <typename R>
    int DependencyGraph<R>::number_of_connected_components() {
        return g->number_of_connected_components();
    }
    
    template <typename R>
    std::vector<int> DependencyGraph<R>::component_vertices(int connected_component_ID) {
        return g->component_vertices(connected_component_ID);
    }
    
    template <typename R>
    std::vector< int > DependencyGraph<R>::special_vertices() {
        return g->special_vertices();
    }
    
    template class DependencyGraph<std::mt19937>;

}
