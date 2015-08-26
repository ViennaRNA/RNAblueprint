/*!\file dependency_graph.h 
 * \brief This file holds the class definitions for the internal representation of the DependencyGraph.
 *
 * Created on: 18.03.2014
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 * \cond INTERNAL
 */

#ifndef DEPENDENCY_GRAPH_H
#define	DEPENDENCY_GRAPH_H

#include "common.h"
#include "graphcommon.h"
#include "probability_matrix.h"
#include "pathcoloring.h"

#include <sstream>
#include <random>

namespace design {
    namespace detail {

        template <typename R>
        class DependencyGraph {
        public:
            DependencyGraph(std::vector<std::string> structures, std::string constraints, R rand);

            unsigned long long number_of_sequences() {
                return nos;
            }

            bool is_bipartite() {
                return bipartite;
            }
            Sequence get_sequence();
            std::string get_sequence_string();
            void set_sequence(Sequence sequence);
            void set_sequence_string(std::string seq_str);
            void set_sequence();
            unsigned long long mutate_local(int min_num_pos, int max_num_pos);
            unsigned long long mutate_global(int min_num_pos, int max_num_pos);
            unsigned long long mutate(int position);
            unsigned long long mutate(int start, int end);
            R * rand_ptr;
            ~DependencyGraph() = default;
        private:
            Graph graph;
            bool bipartite; // if dependency graph is bipartite and a therefore a solution exists
            unsigned long long nos = 0; // number of sequences/solutions
            R rand;
            void calculate_probabilities(Graph& g);
            void sample_sequence(Graph& g);
            void reset_colors(Graph g);
        };
    }
}
#endif	/* DEPENDENCY_GRAPH_H */


/* 
 * \endcond INTERNAL
 */
