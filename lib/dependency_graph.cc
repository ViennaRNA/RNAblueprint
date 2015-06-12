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

            // now calculate all the PMs
            calculate_probabilities(graph);

        }

        template <typename R>
        void DependencyGraph<R>::calculate_probabilities(Graph g) {

            Graph::children_iterator cg, cg_end;
            for (boost::tie(cg, cg_end) = g.children(); cg != cg_end; ++cg) {




            }


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