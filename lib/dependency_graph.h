/* This program reads secundary RNA structures in dot-bracket and
 * builds a graph for a latter ear-decomposition and bipartitness-check
 *
 * Created on: 18.03.2014
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
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
            void mutate(int position);
            void mutate();
            void reset_colors();
            R * rand_ptr;
            ~DependencyGraph();
        private:
            Graph graph;
            bool bipartite; // if dependency graph is bipartite and a therefore a solution exists
            unsigned long long nos = 0; // number of sequences/solutions
            R rand;
            void calculate_probabilities(Graph& g);
            void sample_sequence(Graph& g);
        };
    }
}
#endif	/* DEPENDENCY_GRAPH_H */

