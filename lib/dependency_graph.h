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
        
        typedef std::unordered_map< Graph*, ProbabilityMatrix> ProbabilityMatrixStorage;

        template <typename R>
        class DependencyGraph {
        public:
            DependencyGraph(std::vector<std::string> structures, std::string constraints, R rand);

            unsigned long long number_of_sequences();
            unsigned long long number_of_sequences(int connected_component_ID);
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
            std::map< int, std::vector<int> > connected_components();
            std::vector< int > special_vertices();
            R * rand_ptr;
        private:
            Graph graph;
            ProbabilityMatrixStorage pms;
            bool bipartite; // if dependency graph is bipartite and a therefore a solution exists
            R rand;
            void calculate_probabilities(Graph& g);
            void sample_sequence(Graph& g);
            void reset_colors(Graph& g);
        };
    }
}
#endif	/* DEPENDENCY_GRAPH_H */


/* 
 * \endcond INTERNAL
 */
