/*!\file RNAdesign.h 
 * \brief This file holds the external representation of the DependencyGraph, the main construct for designing RNA molecules.
 * 
 * The dependency graph is constructed from structures in dot-bracket notation and sequence constraints following the IUPAC notation.
 * All important functions on the graph are available as member functions.
 *
 * Created on: 26.06.2014
 * @Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

/*! \mainpage RNAdesign library written in C++
 *
 * \section intro_sec Introduction
 *
 * The RNAdesign library still needs some documentation!
 *
 * \section install_sec Installation
 * 
 * 
 * Just call these commands:
 * 
 * ./autogen.sh
 * ./configure
 * make
 * make install
 * 
 * TIP: You might want call ./configure --help for all install options!
 * 
 */


#ifndef RNADESIGN_H
#define RNADESIGN_H

// include header
#include "common.h"
#include "dependency_graph.h"

/*! \brief All classes and functions for the RNA design library are under the design namespace.
 */
namespace design {
    /*! \brief Initialize the Library.
     * 
     *  Set the debug boolean to true if you want to get verbose output to std:err
     */
    void initialize_library(bool debug);

    /*! \brief Dependency Graph which holds all structural constraints.
     * 
     *  This graph is used to generate valid sequences compatible to the input structures
     */
    template <typename R>
    class DependencyGraph {
    public:
        /*! \brief constructor for the Dependency graph.
        *
        *  A vector of strings for structures in dot-bracket notation. This will be parsed to a dependency graph
        *  A string containing the sequence constraints in IUPAC notation.
        *  A random number generator of your choice. Supported at the moment is only std::mt19937, but as it is templated, it can easily be extended to support more generators.
        */
        DependencyGraph(std::vector<std::string> structures, std::string constraints, R rand);
        /*! \brief constructor for the Dependency graph using a predefined random number generator with the given seed.
         */
        DependencyGraph(std::vector<std::string> structures, std::string constraints, unsigned long seed);
        /*! \brief constructor for the Dependency graph using a predefined random number generator with a clock generated seed.
        */
        DependencyGraph(std::vector<std::string> structures, std::string constraints);
        /*! \brief constructor for the Dependency graph without any sequence constraints. Sending an empty string leads to the same result.
        */
        DependencyGraph(std::vector<std::string> structures, R rand);
        /*! \brief constructor for the Dependency graph using a predefined random number generator with a clock generated seed and an empty string as
         * sequence constraints.
        */
        DependencyGraph(std::vector<std::string> structures);
        /*! \brief Simple copy constructor which sets a new random seed in the copy sampled from the old generator.
        */
        DependencyGraph(const DependencyGraph& copy);
        /*! \brief Simple destructor for the dependency graph object.
        */
        ~DependencyGraph();
        /*! \brief Set the maximum number of previous sampled sequences to memorize (Default: 10).
        *  
        *  A history of all sampled sequences is stored in the dependency graph object. This sets the size of the stack in number of sequences.
        */ 
        void set_history_size(int size);
        /*! \brief Returns the root graph in GraphML format as a std::string
        *  
        *  Get the dependency graph in the xml based GraphML format for further analysis or visualization.
        */        
        std::string get_graphml();
        /*! \brief Returns the connected component graph with the connected_component_ID in GraphML format as a std::string
        *  
        *  Get the connected component graph in the xml based GraphML format for further analysis or visualization.
        */
        std::string get_graphml(int connected_component_ID);
        /*! \brief Returns the current RNA sequence as a std::string
        *  
        *  This sequence is only N directly after construction. You need to call mutate() first to sample a initial sequence.
        */
        std::string get_sequence();
        /*! \brief Allows you to set a initial sequence as starting point for your optimization.
        *  
        *  Only real bases are allowed and the sequence has to fulfill all structural constraints, otherwise an error is thrown.
        */
        void set_sequence(std::string sequence);
        /*! \brief Resets all bases to N in the whole dependency graph and samples a new initial sequence randomly
        *  
        *  Call get_sequence() after you sampled a new sequence.
        */
        void set_sequence();
        /*! \brief Reverts the sequence to the previous one
        *  
        *  Call revert_sequence() to set the dependency graph to the previous sampled sequence
        */
        bool revert_sequence();
        /*! \brief Reverts the sequence to a previous one being (jump) steps in the history;
        *  
        *  revert_sequence(2) would for example revert the sequence to the one before the previous one.
        */
        bool revert_sequence(unsigned int jump);
        /*! \brief Randomly chooses one path (either top-level a connected component, or within a block, etc.) with the given length and mutates all positions.
        *  
        *  Special vertices such as cut points or articulation points will stay the same. Therefore it is guaranteed that the sampling is correct,
        *  even if we only sample a small local piece of a more complex graph object.
        *  int min_num_pos, int max_num_pos set the minimal/maximal number of mutated positions, e.g., for (3,5) only paths
        *  with minimal 3 and maximal 5 non-special vertices will be chosen for mutation. 0 for the upper boundary defines infinity. The range is inclusive!
        */
        SolutionSizeType mutate_local(int min_num_pos, int max_num_pos);
        /*! \brief Randomly chooses one path (either top-level a connected component, or within a block, etc.) and mutates all positions.
        *  
        */
        SolutionSizeType mutate_local();
        /*! \brief Randomly chooses a connected component and samples a new sequence for the whole component.
        *  
        *  This is a more global way of mutating the structure as it probably exchanges a much bigger graph object.
        *  int min_num_pos, int max_num_pos set the minimal/maximal number of mutated positions, e.g., for (3,5) only connected components
        *  with minimal 3 and maximal 5 vertices will be chosen for mutation. 0 for the upper boundary defines infinity. The range is inclusive!
        * 
        *  Returns: The number of possible sequences for this mutation.
        */
        SolutionSizeType mutate_global(int min_num_pos, int max_num_pos);
        /*! \brief Takes the connected component with the ID (int connected_component_ID) and samples a new sequence for the whole component.
        *  
        *  The connected_component_IDs range from 0 to  number_of_connected_components()
        *  
        *  Returns: The number of possible sequences for this mutation.
        */
        SolutionSizeType mutate_global(int connected_component_ID);
        /*! \brief Takes a random connected component and samples a new sequence for the whole component.
        *  
        *  Returns: The number of possible sequences for this mutation.
        */
        SolutionSizeType mutate_global();
        /*! \brief Resets only the smallest subgraph(s) possible containing the vertex at the given position in the sequence.
        *  
        *  This way you can optimize by targeted mutagenesis at the given positions. All positions dependent on the chosen one
        *  will also be mutated.
        *  If your position is a special vertex, the whole connected component will be re-sampled. Else, in case of being a non-special vertex,
        *  only the smallest path containing the vertex will be mutated. Positions range from [0,N-1]
        * 
        *  Returns: The number of possible sequences for this mutation.
        */
        SolutionSizeType mutate(int position);
        /*! \brief Resets only the smallest subgraph(s) possible containing the vertices from position start to end.
        *  
        *  This way you can optimize by targeted mutagenesis at the given positions. All positions dependent on the chosen ones
        *  will also be mutated.
        *  If your positions contain special vertices, the whole connected components will be re-sampled. Else, in case of being only non-special vertex,
        *  only the smallest paths containing the vertices will be mutated. Positions range from [0,N-1]
        * 
        *  Returns: The number of possible sequences for this mutation.
        */
        SolutionSizeType mutate(int start, int end);
        /*! \brief Returns the amount of solutions given the dependency graph and sequence constraints
        *  
        *  Number of sequences is the total amount of sequences possible for the given structural and sequence constraints.
        */
        SolutionSizeType number_of_sequences();
        /*! \brief Returns the amount of solutions for the connected component with the given ID.
        *  
        *  Number of sequences is the total amount of sequences possible for the given structural and sequence constraints.
        *  The connected_component_IDs can be retrieved by calling connected_components()
        */
        SolutionSizeType number_of_sequences(int connected_component_ID);
        /*! \brief Returns the number of connected components
        *  
        */
        int number_of_connected_components();
        /*! \brief Returns a vector with all the vertices in this component of the given connected_component_ID
        *  
        */
        std::vector<int> component_vertices(int connected_component_ID);
        /*! \brief Returns a list of vertices specified as "special".
        *  
        *  These special vertices include cut points, articulation points and also sequence constraints.
        */
        std::vector< int > special_vertices();
        
        // ??? get_special_probabilities(connected components ID)
        // Returns basically the ProbabilityMatrix for all special vertices (remove make_internal() function) of the whole connected component 
        // hash of hash: vertex -> base_color -> number_of_sequences

        // ??? get_constrained_probabilities()
        // Returns a hash of hash table of vertex -> base_color -> number_of_sequences for all positions given the base_color.
        // for special vertices we can just call the function above and get the sums, 
        // for non-specials we maybe have to do a constraint recalculation for every position/base combination?
        
        
    private:
        /*! \brief Pointer to the internal dependency graph object.
         */
        detail::DependencyGraph<R> * g;
    };
}

#endif
