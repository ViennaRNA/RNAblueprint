/* RNAdesign
 * A program for designing RNA molecules.
 *
 * Created on: 18.03.2014
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

// include header
#include "dependency_graph.h"

#include "parsestruct.h"
#include "printgraph.h"
#include "decompose.h"

// include boost components
#include <boost/graph/iteration_macros.hpp>

namespace design
{

    namespace detail{

        template <typename R>
        DependencyGraph<R>::DependencyGraph(std::vector<std::string> structures, R _rand)
        : rand(_rand) {

            if (debug) {
                std::cerr << "Initializing DependencyGraph..." << std::endl;
            }

            // initialize random generator
            rand_ptr = &rand;

            // generate graph from input vector
            graph = parse_structures(structures);

            // decompose the graph into its connected components, biconnected
            // components and decompose blocks via ear decomposition
            bipartite = decompose_graph(graph, rand_ptr);
            
            BGL_FORALL_VERTICES_T(v, graph, Graph) {
                // reset the color tag for the calculate_probabilities function
                graph[v].color = 0;
            }

            // now calculate all the PMs
            calculate_probabilities(graph);

        }

        template <typename R>
        void DependencyGraph<R>::calculate_probabilities(Graph g) {
            
            // Remember a temporary PM which holds the current state
            ProbabilityMatrix current;
            
            Graph::children_iterator cg, cg_end;
            for (boost::tie(cg, cg_end) = g.children(); cg != cg_end; ++cg) {
                
                ProbabilityMatrix thischild;
                
                if (boost::get_property(*cg, gpt).path) {
                    // this is a path and therefore needs to be treated separately
                    // calculate PM for Path and save to subgraph
                    //TODO boost::get_property(*cg, gpt).pm = calculate_path_pm(*cg);
                } else {
                    // recursion is here
                    calculate_probabilities(*cg);
                }
                
                thischild = boost::get_property(*cg, gpt).pm;
                // Multiply current with pm of this child
                current = current * thischild;
                
                BGL_FORALL_VERTICES_T(v, *cg, Graph) {
                    if ((*cg)[v].special) {
                        // update current degrees as status of special points
                        (*cg)[v].color += boost::degree(v, *cg);
                        
                        // check if a vertex becomes internal here
                        if ((*cg)[v].color == boost::degree((*cg).local_to_global(v), (*cg).root())) {
                            
                        }
                    }
                }
            }
            
            // save final state of PM to the main graph
            boost::get_property(g, gpt).pm = current;
        }

        template <typename R>
        Sequence DependencyGraph<R>::get_sequence() {
            Sequence sequence(boost::num_vertices(graph), N);

            BGL_FORALL_VERTICES_T(v, graph, Graph) {
                sequence[boost::get(boost::vertex_color_t(), graph, v)] = graph[v].base;
            }
            return sequence;
        }

        template <typename R>
        std::string DependencyGraph<R>::get_sequence_string() {
            Sequence sequence = get_sequence();

            std::stringstream stream;
            stream << sequence;

            return stream.str();
        }

        template <typename R>
        void DependencyGraph<R>::mutate() {
            // reset all the colors to N
            reset_colors();
            // TODO replace with a good mutation function
        }

        template <typename R>
        void DependencyGraph<R>::mutate(int position) {
            // TODO replace with a good mutation function	
            this->mutate();
        }

        template <typename R>
        void DependencyGraph<R>::reset_colors() {

            BGL_FORALL_VERTICES_T(v, graph, Graph) {
                graph[v].base = N;
            }
        }

        template <typename R>
        DependencyGraph<R>::~DependencyGraph() {
            // TODO clean up all the data
        }

        template class DependencyGraph<std::mt19937>;
    }
}