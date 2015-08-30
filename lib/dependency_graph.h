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
#include <unordered_set>

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
            // call this function to mutate a random subgraph (either a path, if graph_type=-1 or a connected component, if graph_type=1)
            unsigned long long mutate_local_global(int graph_type, int min_num_pos, int max_num_pos);
            unsigned long long mutate_global(int connected_component_ID);
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
            unsigned long long sample_sequence(Graph& g);
            void reset_colors(Graph& g);
            unsigned long long mutate(Graph& g);
            // this function fills the subgraphs set with all the sugraphs of the given type (root, cc, bc, path)
            // if int type= -1, then it returns all subgraphs which are actual paths (gp.is_path == true).
            // you can specify also a minimal and maximal size of the subgraph
            void get_subgraphs(Graph& g, std::unordered_set< Graph* >& subgraphs, int type, int min_size, int max_size);
        };
    }
}
#endif	/* DEPENDENCY_GRAPH_H */


/* 
 * \endcond INTERNAL
 */
