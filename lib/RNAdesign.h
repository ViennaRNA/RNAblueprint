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
    /*! \brief Initialize the Library
     * 
     *  Set the debug boolean to true if you want to get verbose output to std:err
     */
    void initialize_library(bool debug);

    /*
     *! \brief Dependency Graph which holds all structural constraints
     * 
     *  This graph is used to generate valid sequences compatible to the input structures
     */
    template <typename R>
    class DependencyGraph {
    public:
        /*! \brief Constructor for the Dependency graph
        *
        *  A vector of strings for structures in dot-bracket notation. This will be parsed to a dependency graph
        *  A string containing the sequence constraints in IUPAC notation.
        *  A random number generator of your choice. Supported at the moment is only std::mt19937, but as it is templated, it can easily be extended to support more generators.
        */
        DependencyGraph(std::vector<std::string> structures, std::string constraints, R rand);
        ~DependencyGraph();
        /*! \brief get_sequence() returns the current RNA sequence as a std::string
        *  
        *  This sequence is only N directly after construction. You need to call mutate() first to sample a initial sequence.
        */
        std::string get_sequence();
        /*! \brief mutate() resets all bases to N in the whole dependency graph and samples a new sequence
        *  
        *  Call get_sequence() after you sampled a new sequence.
        */
        void mutate();
        /*! \brief mutate(int pos) resets only the smalles subgraph(s) containing the vertex pos to N.
        *  
        *  This way the walk through the solution space happens in much smaller hops.
        */
        void mutate(int position);
        /*! \brief number_of_sequences() returns the amount of solutions given the dependency graph and sequence constraints
        *  
        *  Number of sequences is the total amount of sequences possible for the given structural and sequence constraints.
        */
        unsigned long long number_of_sequences();
    private:
        detail::DependencyGraph<R> * g;
    };
}

#endif