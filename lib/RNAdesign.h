/*!\file RNAdesign.h 
 * \brief This file holds the external representation of the DependencyGraph, the main construct for designing RNA molecules.
 * 
 * The dependency graph is constructed from structures in dot-bracket notation and sequence constraints following the IUPAC notation.
 * All important functions on the graph are available as member functions of this object.
 *
 * Created on: 26.06.2014
 * @author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

/*! \mainpage RNAdesign library written in C++
 *
 * \section intro_sec Introduction
 *
 * The RNAdesign library still needs some documentation!
 * Write here something about the theoretical background
 * Show how to cite the software
 *
 * \section dependency_sec Dependencies
 *
 * Required:
 * - GNU Automake
 * - Boost Graph Library
 * - C++ Standard Library
 * 
 * Optional:
 * - Boost Programm Options
 * - SWIG for interfaces
 * - Python for interface
 * - Perl for interface
 * - Doxygen for documentation
 * - LaTeX for PDF documentation
 * - libGMP for multiprecision integers
 * - openmp for parallel computation
 * - Boost Unit Test Framework
 *
 * \section install_sec Installation
 * 
 * Just call these commands:
 * 
 * ./autogen.sh\n
 * ./configure\n
 * make\n
 * make install\n
 * 
 * Most important configure options are:
 * - \-\-prefix           Specify an installation path prefix
 * - \-\-with-boost       Specify the installation directory of the boost library
 * - \-\-disable-program  Disable RNAdesign program compilation
 * - \-\-disable-swig     Disable all SWIG scripting interfaces
 * - \-\-enable-libGMP    Enable the calculation of big numbers with multiprecision
 * - \-\-disable-openmp   Disable the usage of parallel computation
 * 
 * TIP: You might want call ./configure --help for all install options!
 * 
 * \section cite_sec How to cite
 * 
 * This is a early release for testing purposes only and a publication for this software package is in progress.
 * We will update this information as soon as a preprint is available online!
 * 
 */


#ifndef RNADESIGN_H
#define RNADESIGN_H

// include header
#include "common.h"
#include "dependency_graph.h"
#include "parsestruct.h"
#include "decompose.h"

// boost header
#include <chrono>
#include <boost/graph/bipartite.hpp>


/*! \brief All classes and functions for the RNA design library are under the design namespace.
 */
namespace design {
    /*! \brief Initialize the Libraries global variables.
     * 
     * \param debug \b boolean whether to print debugging information to std:err (default: false)
     */
    void initialize_library(bool debug);
    /*! \brief Initialize the Libraries global variables.
     * 
     * \param debug \b boolean whether to print debugging information to std:err (default: false)
     * \param construction_timeout \b integer specifying the dependency graph construction timeout in seconds. 0 is infinite. (default: 0)
     */
    void initialize_library(bool debug, int construction_timeout);
    /*! \brief Generate a graphml representation of structural and sequence constraints
     * 
     * This function generates a graphml representation of the dependency graph given some structural
     * and sequence constraints without constructing a DependencyGraph object, with or without decomposition of the graph into subpaths.
     * It is mainly thought to be useful for developmental purposes, analysis of hard problems and vizualisation.
     * 
     * \param structures \b vector of \b string structures in dot-bracket notation.
     * \param constraints \b string containing the sequence constraints in IUPAC notation. Can also be a empty string!
     * \param decompose \b boolean Whether to decompose the dependency graph into paths and therefore draw special vertices and ears.
     * \param seed \b unsigned long Seed for the random number generator which is used for some random parts of the decomposition.
     * \exception std::exception if input is invalid or construction/decomposition fails an exception is thrown.
     * \return \b string containing the GraphML notation of the dependency graph.
     */
    std::string structures_to_graphml(std::vector<std::string> structures, std::string constraints, bool decompose, unsigned long seed);
    /*! \brief Generate a graphml representation of structural and sequence constraints
     * 
     * This function generates a graphml representation of the dependency graph given some structural
     * and sequence constraints without constructing a DependencyGraph object, with or without decomposition of the graph into subpaths.
     * It is mainly thought to be useful for developmental purposes, analysis of hard problems and vizualisation.
     * Caution: There are random parts in the decomposition algorithms. If you want to assign a seed, use structures_to_graphml(std::vector<std::string> structures, std::string constraints, bool decompose, unsigned long seed)
     * 
     * \param structures \b vector of \b string structures in dot-bracket notation.
     * \param constraints \b string containing the sequence constraints in IUPAC notation. Can also be a empty string!
     * \param decompose \b boolean Whether to decompose the dependency graph into paths and therefore draw special vertices and ears.
     * \exception std::exception if input is invalid or construction/decomposition fails an exception is thrown.
     * \return \b string containing the GraphML notation of the dependency graph.
     */
    std::string structures_to_graphml(std::vector<std::string> structures, std::string constraints, bool decompose);
    /*! \brief Generate a graphml representation of structural and sequence constraints
     * 
     * This function generates a graphml representation of the dependency graph given some structural
     * and sequence constraints without constructing a DependencyGraph object, but with decomposition of the graph into subpaths.
     * It is mainly thought to be useful for developmental purposes, analysis of hard problems and vizualisation.
     * 
     * \param structures \b vector of \b string structures in dot-bracket notation.
     * \param constraints \b string containing the sequence constraints in IUPAC notation. Can also be a empty string!
     * \exception std::exception if input is invalid or construction/decomposition fails an exception is thrown.
     * \return \b string containing the GraphML notation of the dependency graph.
     */
    std::string structures_to_graphml(std::vector<std::string> structures, std::string constraints);
    /*! \brief Returns whether the dependency graph built from the given input structures is bipartite
     * 
     * This helper function checks, if a set of input structures will generate a bipartite dependency graph.
     * If the graph is not bipartite, it cannot be used to generate a valid RNA sequence.
     * 
     * \param structures \b vector of \b string structures in dot-bracket notation.
     * \exception std::exception if input is invalid or construction fails an exception is thrown.
     * \return \b boolean specifying whether the input structures can be constructed to a bipartite dependency graph.
     */
    bool graph_is_bipartite(std::vector<std::string> structures);
    /*! \brief Returns whether the the given sequence is compatible to all the given structures
     * 
     * This function checks, if a given sequence can form all base-pairs given by a set of input structures.
     * Or the other way around, if the given structures could have produced this sequence.
     * 
     * \param sequence \b string in IUPAC notation.
     * \param structures \b vector of \b string structures in dot-bracket notation.
     * \exception std::exception if input is invalid.
     * \return \b boolean specifying whether the input sequence is compatible to all the given structures.
     */
    bool sequence_structure_compatible(std::string sequence, std::vector<std::string> structures);
    /*! \brief Returns whether the the given sequence is compatible to all the given structures
     * 
     * This function checks, if a given sequence can fold into the given structure and returns 
     * an empty vector if this is the case. Else, it returns all positions on the sequence which
     * are incompatible with the given structural constraint.
     * E.g. incompatible_sequence_positions("ANC", "(.)") would return {(0, 2)}!
     * 
     * \param sequence \b string in IUPAC notation.
     * \param structure \b string in dot-bracket notation.
     * \exception std::exception if input is invalid.
     * \return \b map of a pair of \b integers specifying the sequence positions incompatible to the structure input.
     */
    std::vector<int> incompatible_sequence_positions(std::string sequence, std::string structure);
    /*! \brief Dependency Graph which holds all structural constraints.
     * 
     * This graph is used to generate valid sequences compatible to the input structures
     */
    template <typename R>
    class DependencyGraph {
    public:
        /*! \brief Constructor for the Dependency graph.
         *
         * \param structures \b vector of \b string structures in dot-bracket notation.
         * \param constraints \b string containing the sequence constraints in IUPAC notation. Can also be a empty string!
         * \param rand random number generator of your choice. Supported at the moment is only std::mt19937, but as it is templated, it can easily be extended to support more generators.
         * \exception std::exception if input is invalid or construction fails an exception is thrown.
        */
        DependencyGraph(std::vector<std::string> structures, std::string constraints, R rand);
        /*! \brief constructor for the Dependency graph using a predefined random number generator with the given seed.
         * 
         * \param structures \b vector of \b string structures in dot-bracket notation.
         * \param constraints \b string containing the sequence constraints in IUPAC notation. Can also be a empty string!
         * \param seed \b unsigned long to be used as seed of the random number generator.
         * \exception std::exception if input is invalid or construction fails an exception is thrown.
         */
        DependencyGraph(std::vector<std::string> structures, std::string constraints, unsigned long seed);
        /*! \brief constructor for the Dependency graph using a predefined random number generator with a clock generated seed.
         * 
         * \param structures \b vector of \b string structures in dot-bracket notation.
         * \param constraints \b string containing the sequence constraints in IUPAC notation. Can also be a empty string!
         * \exception std::exception if input is invalid or construction fails an exception is thrown.
        */
        DependencyGraph(std::vector<std::string> structures, std::string constraints);
        /*! \brief constructor for the Dependency graph without any sequence constraints.
         * 
         * \param structures \b vector of \b string structures in dot-bracket notation.
         * \param rand random number generator of your choice. Supported at the moment is only std::mt19937, but as it is templated, it can easily be extended to support more generators.
         * \exception std::exception if input is invalid or construction fails an exception is thrown.
        */
        DependencyGraph(std::vector<std::string> structures, R rand);
        /*! \brief constructor for the Dependency graph using a predefined random number generator with a clock generated seed and an empty string as
         * sequence constraints.
         * 
         * \param structures \b vector of \b string structures in dot-bracket notation.
         * \exception std::exception if input is invalid or construction fails an exception is thrown.
        */
        DependencyGraph(std::vector<std::string> structures);
        /*! \brief Copy constructor which sets a new random seed in the copy sampled from the old generator.
        */
        DependencyGraph(const DependencyGraph& copy);
        /*! \brief Destructor for the dependency graph object.
        */
        ~DependencyGraph();
        /*! \brief Set the maximum number of previous sampled sequences to memorize (Default: 10).
         *  
         * A history of all sampled sequences is stored within the dependency graph object. 
         * Use this function to set the size of the storage stack.
         * 
         * \param size \b integer to set the size of the history storage container.
        */ 
        void set_history_size(int size);
        /*! \brief Returns the root graph in GraphML format as a string
         * 
         * Get the dependency graph in the XML based GraphML format for further analysis or visualization.
         * 
         * \return \b string containing the GraphML notation of the dependency graph.
        */        
        std::string get_graphml();
        /*! \brief Returns the connected component graph with the connected_component_ID in GraphML format as a string
         *  
         * Get the graph for a specific connected component in the XML based GraphML format for further analysis or visualization.
         * 
         * \param connected_component_ID \b integer specifying the connected component with its ID [ 0, number_of_connected_components() ).
         * \sa number_of_connected_components()
         * \exception std::out_of_range if connected_component_ID is invalid.
         * \return \b string containing the GraphML notation of the connected component.
        */
        std::string get_graphml(int connected_component_ID);
        /*! \brief Get the current RNA sequence as a string
         * 
         * This sequence is only N directly after construction. You need to call either set or sample a initial sequence to avoid this behavior.
         * 
         * \sa sample(), set_sequence()
         * \return \b string representing the current RNA sequence.
        */
        std::string get_sequence();
        /*! \brief Allows you to set a initial sequence as starting point for your optimization.
         * 
         * Only real bases are allowed and the sequence has to fulfill all structural constraints, otherwise an error is thrown.
         * 
         * \param sequence \b string containing the sequence constraints in IUPAC notation. Only [AUGC] are allowed as bases.
         * \exception std::exception if input string contains invalid characters or constraints cannot be fulfilled.
         * \return \b number of possible solutions for this sampling. Is 1 all the time for set_sequence()
        */
        SolutionSizeType set_sequence(std::string sequence);
        /*! \brief Reverts the sequence to the previous one.
         * 
         * This sets the dependency graph to the previous sampled sequence state.
         * 
         * \sa set_history_size(int size), revert_sequence(unsigned int jump)
         * \return \b boolean specifying if this move was possible. If no previous sequence is stored, it returns false.
        */
        bool revert_sequence();
        /*! \brief Reverts the sequence to a previous one being (jump) steps in the history;
         * 
         * \param jump \b integer specifying the length of the time-jump. 2 would for example revert the sequence to the one before the previous one.
         * \sa set_history_size(int size), revert_sequence()
         * \return \b boolean specifying if this move was possible. If no such previous sequence is stored, it returns false.
        */
        bool revert_sequence(unsigned int jump);
        /*! \brief Resets all bases in the whole dependency graph and samples a new sequence randomly.
         * 
         * Call get_sequence() after you sampled a new sequence.
         * 
         * \sa get_sequence(), sample(int position), sample(int start, int end)
         * \return \b number of possible solutions for this sampling.
        */
        SolutionSizeType sample();
        /*! \brief Resets and samples only the smallest subgraph(s) possible containing the vertex at the given position in the sequence.
         * 
         * This way you can optimize by targeted sampling at the given positions. All positions dependent on the chosen one
         * will also be sampled.
         * If your position is a special vertex, the whole connected component will be re-sampled. Else, in case of being a non-special vertex,
         * only the smallest path containing the vertex will be sampled.
         * 
         * \param position \b integer specifying the position to re-sample [ 0, N )
         * \sa get_sequence(), sample(), sample(int start, int end)
         * \return \b number of possible solutions for this sampling.
        */
        SolutionSizeType sample(int position);
        /*! \brief Resets only the smallest subgraph(s) possible containing the vertices from position start to end.
         * 
         * This way you can optimize by targeted sampling at the given positions. All positions dependent on the chosen ones
         * will also be sampled.
         * If your positions contain special vertices, the whole connected components will be re-sampled. Else, in case of being only non-special vertex,
         * only the smallest paths containing the vertices will be sampled. Positions for start and end are inclusive [ start, end ]
         *
         * \param start \b integer specifying the first position to re-sample [ 0, N )
         * \param end \b integer specifying the last position to re-sample [ 0, N )
         * \sa get_sequence(), sample(), sample(int position)
         * \return \b number of possible solutions for this sampling.
        */
        SolutionSizeType sample(int start, int end);
        /*! \brief Randomly chooses one path (either top-level a connected component, or within a block, etc.) with the given size and samples all positions.
         * 
         * Special vertices such as cut points or articulation points will stay the same. Therefore it is guaranteed that the sampling is correct,
         * even if we only sample a small local piece of a more complex graph object.
         * min_num_pos and max_num_pos set the minimal/maximal number of sampled positions, e.g., for [ 3, 5 ] only paths
         * with minimal 3 and maximal 5 non-special vertices will be chosen for sampling.
         * 
         * \param min_num_pos \b integer specifying the minimal size of the component to re-sample [ 1, N )
         * \param max_num_pos \b integer specifying the maximal size of the component to re-sample [ 1, N ). 0 defines infinity.
         * \sa sample_local(), sample_global(), sample_global(int min_num_pos, int max_num_pos), sample_global(int connected_component_ID)
         * \return \b number of possible solutions for this sampling.
        */
        SolutionSizeType sample_local(int min_num_pos, int max_num_pos);
        /*! \brief Randomly chooses one path (either top-level a connected component, or within a block, etc.) and samples all positions.
         * 
         * \sa sample_local(int min_num_pos, int max_num_pos), sample_global(), sample_global(int min_num_pos, int max_num_pos), sample_global(int connected_component_ID)
         * \return \b number of possible solutions for this sampling.
        */
        SolutionSizeType sample_local();
        /*! \brief Randomly chooses a connected component with the given size and samples a new sequence for the whole component.
         * 
         * This is a more global way of sampling a new sequence as it exchanges a much bigger graph object.
         * min_num_pos, max_num_pos set the minimal/maximal number of sampled positions, e.g., for [ 3, 5 ] only connected components
         * with minimal 3 and maximal 5 vertices will be chosen for sampling.
         * 
         * \param min_num_pos \b integer specifying the minimal size of the component to re-sample [ 1, N )
         * \param max_num_pos \b integer specifying the maximal size of the component to re-sample [ 1, N ). 0 defines infinity.
         * \sa sample_local(), sample_local(int min_num_pos, int max_num_pos), sample_global(), sample_global(int connected_component_ID)
         * \return \b number of possible solutions for this sampling.
        */
        SolutionSizeType sample_global(int min_num_pos, int max_num_pos);
        /*! \brief Takes the connected component with the specified ID and samples a new sequence for the whole component.
         * 
         * \param connected_component_ID \b integer specifying the connected component with its ID [ 0, number_of_connected_components() ).
         * \sa number_of_connected_components(), sample_local(), sample_local(int min_num_pos, int max_num_pos), sample_global(), sample_global(int min_num_pos, int max_num_pos)
         * \exception std::out_of_range if connected_component_ID is invalid.
         * \return \b number of possible sequences for this sampling.
        */
        SolutionSizeType sample_global(int connected_component_ID);
        /*! \brief Takes a random connected component and samples a new sequence for the whole component.
         * 
         * \sa sample_local(), sample_local(int min_num_pos, int max_num_pos), sample_global(int connected_component_ID), sample_global(int min_num_pos, int max_num_pos)
         * \return \b number of possible sequences for this sampling.
        */
        SolutionSizeType sample_global();
        /*! \brief Returns the amount of solutions given the dependency graph and sequence constraints
         * 
         * Number of sequences is the total amount of sequences possible for the given structural and sequence constraints.
         * This defines the size of the solution space.
         * 
         * \sa number_of_sequences(int connected_component_ID)
         * \return \b number of possible sequences for this design problem.
        */
        SolutionSizeType number_of_sequences();
        /*! \brief Returns the amount of solutions for the connected component with the given ID.
         * 
         * Number of sequences is the total amount of sequences possible for the given connected component.
         *  
         * \param connected_component_ID \b integer specifying the connected component with its ID [ 0, number_of_connected_components() ).
         * \sa number_of_sequences(), number_of_connected_components()
         * \exception std::out_of_range if connected_component_ID is invalid.
         * \return \b number of possible sequences for this connected component.
        */
        SolutionSizeType number_of_sequences(int connected_component_ID);
        /*! \brief Returns the number of connected components into which the dependency graph was decomposed.
         * 
         * \return \b integer defining the number of connected components
        */
        int number_of_connected_components();
        /*! \brief Returns a vector with all the vertices in this component of the given ID.
         * 
         * \param connected_component_ID \b integer specifying the connected component with its ID [ 0, number_of_connected_components() ).
         * \sa number_of_connected_components()
         * \exception std::out_of_range if connected_component_ID is invalid.
         * \return \b vector of \b integer values specifying the positions/vertices contained in the connected component.
        */
        std::vector<int> component_vertices(int connected_component_ID);
        /*! \brief Returns a list of vertices specified as "special".
         * 
         * These special vertices include cut points, articulation points and also a cycle opening cuts.
         * 
         * \sa special_vertices(int connected_component_ID)
         * \return \b vector of \b integer values specifying the positions/vertices marked as special.
        */
        std::vector< int > special_vertices();
        /*! \brief Returns a list of vertices in the specified connected component specified as "special".
         * 
         * These special vertices include cut points, articulation points and also a cycle opening cuts.
         * 
         * \param connected_component_ID \b integer specifying the connected component with its ID [ 0, number_of_connected_components() ).
         * \sa number_of_connected_components(), special_vertices()
         * \exception std::out_of_range if connected_component_ID is invalid.
         * \return \b vector of \b integer values specifying the positions/vertices marked as special in the connected component.
        */
        std::vector< int > special_vertices(int connected_component_ID);
        
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
