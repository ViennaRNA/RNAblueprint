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
    
    std::string structures_to_graphml(std::vector<std::string> structures, std::string constraints) {
        detail::Graph graph;
        // generate graph from input vector
        try {
            graph = detail::parse_structures(structures);
        } catch (std::exception& e) {
            std::stringstream ss;
            ss << "Error while parsing the structures: " << std::endl << e.what();
            throw std::logic_error(ss.str());
        }
            
        // set sequence constraints
        detail::set_constraints(graph, constraints);
        
        std::ostringstream stream;
        detail::print_graph(graph, dynamic_cast<std::ostream*>(&stream));
        return stream.str();
    }
    
    bool graph_is_bipartite(std::vector<std::string> structures) {
        detail::Graph graph;
        // generate graph from input vector
        try {
            graph = detail::parse_structures(structures);
        } catch (std::exception& e) {
            std::stringstream ss;
            ss << "Error while parsing the structures: " << std::endl << e.what();
            throw std::logic_error(ss.str());
        }
            
        // return if graph is bipartite
        return boost::is_bipartite(graph);
    }

    template <typename R>
    DependencyGraph<R>::DependencyGraph(std::vector<std::string> structures, std::string constraints, R rand) {
        g = new detail::DependencyGraph<R>(structures, constraints, rand);
    }
    
    template <>
    DependencyGraph<std::mt19937>::DependencyGraph(std::vector<std::string> structures, std::string constraints, unsigned long seed) {
        // initialize mersenne twister with the given seed
        g = new detail::DependencyGraph<std::mt19937>(structures, constraints, std::mt19937());
        g->set_seed(seed);
    }
    
    template<>
    DependencyGraph<std::mt19937>::DependencyGraph(std::vector<std::string> structures, std::string constraints) {
        // initialize mersenne twister and set a clock seed
        g = new detail::DependencyGraph<std::mt19937>(structures, constraints, std::mt19937());
        g->set_seed();
    }
    
    template <typename R>
    DependencyGraph<R>::DependencyGraph(std::vector<std::string> structures, R rand) {
        g = new detail::DependencyGraph<R>(structures, "", rand);
    }
    
    template<>
    DependencyGraph<std::mt19937>::DependencyGraph(std::vector<std::string> structures) {
        // initialize mersenne twister with a mersenne twister engine and set a clock seed
        g = new detail::DependencyGraph<std::mt19937>(structures, "", std::mt19937());
        g->set_seed();
    }
    
    template <typename R>
    DependencyGraph<R>::DependencyGraph(const DependencyGraph& copy) {
        g = new detail::DependencyGraph<R>(*copy.g);
        // TODO set new seed from old seed to make it reproduce able!
        g->set_seed();
    }

    template <typename R>
    DependencyGraph<R>::~DependencyGraph() {
        delete g;
    }
    
    template <typename R>
    void DependencyGraph<R>::set_history_size(int size) {
        g->set_history_size(size);
    }
    
    template <typename R>
    std::string DependencyGraph<R>::get_graphml() {
        return g->get_graphml();
    }
    
    template <typename R>
    std::string DependencyGraph<R>::get_graphml(int connected_component_ID) {
        return g->get_graphml(connected_component_ID);
    }

    template <typename R>
    std::string DependencyGraph<R>::get_sequence() {
        return g->get_sequence_string();
    }
    
    template <typename R>
    SolutionSizeType DependencyGraph<R>::set_sequence(std::string sequence) {
        return g->set_sequence_string(sequence);
    }

    template <typename R>
    SolutionSizeType DependencyGraph<R>::sample() {
        return g->sample();
    }
    
    template <typename R>
    bool DependencyGraph<R>::revert_sequence() {
        return g->revert_sequence(1);
    }
    
    template <typename R>
    bool DependencyGraph<R>::revert_sequence(unsigned int jump) {
        return g->revert_sequence((int)jump);
    }
    
    template <typename R>
    SolutionSizeType DependencyGraph<R>::sample_local(int min_num_pos, int max_num_pos) {
        // if local_global=-1 we will sample only actual graphs (independent of their type)
        return g->sample_local_global(-1, min_num_pos, max_num_pos);
    }
    
    template <typename R>
    SolutionSizeType DependencyGraph<R>::sample_local() {
        // if local_global=-1 we will sample only actual graphs, with 0,0 we choose any
        return g->sample_local_global(-1, 0, 0);
    }
    
    template <typename R>
    SolutionSizeType DependencyGraph<R>::sample_global(int min_num_pos, int max_num_pos) {
        // if local_global=1 we will sample only connected components
        return g->sample_local_global(1, min_num_pos, max_num_pos);
    }
    
    template <typename R>
    SolutionSizeType DependencyGraph<R>::sample_global(int connected_component_ID) {
        // sample the connected component with the ID
        return g->sample_global(connected_component_ID);
    }
    
    template <typename R>
    SolutionSizeType DependencyGraph<R>::sample_global() {
        // sample any component with any size
        return g->sample_local_global(1, 0, 0);
    }
    
    template <typename R>
    SolutionSizeType DependencyGraph<R>::sample(int position) {
        return g->sample(position);
    }
    
    template <typename R>
    SolutionSizeType DependencyGraph<R>::sample(int start, int end) {
        return g->sample(start, end);
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
    
    template <typename R>
    std::vector< int > DependencyGraph<R>::special_vertices(int connected_component_ID) {
        return g->special_vertices(connected_component_ID);
    }
    
    template class DependencyGraph<std::mt19937>;

}
