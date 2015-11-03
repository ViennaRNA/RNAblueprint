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
#include <unordered_set>
#include <chrono>
#include <iterator>

namespace design {
    namespace detail {
        // map for storing all probability matrixes
        typedef std::unordered_map< Graph*, ProbabilityMatrix> ProbabilityMatrixStorage;
        
        // actual dependency graph
        template <typename R>
        class DependencyGraph {
        public:
            DependencyGraph(std::vector<std::string> structures, std::string constraints, R rand);
            std::string get_graphml();
            std::string get_graphml(int connected_component_ID);
            bool is_bipartite() {
                return bipartite;
            }
            SolutionSizeType number_of_sequences();
            SolutionSizeType number_of_sequences(int connected_component_ID);
            int number_of_connected_components();
            std::vector< int > component_vertices(int connected_component_ID);
            std::vector< int > special_vertices();
            unsigned long set_seed(int seed) {
                rand.seed(seed);
                return seed;
            }
            unsigned long set_seed();
            
            Sequence get_sequence();
            std::string get_sequence_string();
            void set_sequence(Sequence sequence);
            void set_sequence_string(std::string seq_str);
            void set_sequence();
            // call this function to mutate a random subgraph (either a path, if graph_type=-1 or a connected component, if graph_type=1)
            SolutionSizeType mutate_local_global(int graph_type, int min_num_pos, int max_num_pos);
            SolutionSizeType mutate_global(int connected_component_ID);
            SolutionSizeType mutate(int position);
            SolutionSizeType mutate(int start, int end);
            
            void set_history_size(int size) {
                history_size = size;
            }
            bool revert_sequence(int jump);
        private:
            Graph graph;
            ProbabilityMatrixStorage pms;
            bool bipartite; // if dependency graph is bipartite and a therefore a solution exists
            R rand;
            std::list<Sequence> history;
            unsigned int history_size;
            void remember_sequence();
            void calculate_probabilities(Graph& g);
            SolutionSizeType sample_sequence(Graph& g);
            void reset_colors(Graph& g);
            SolutionSizeType mutate(Graph& g);
            Graph* find_path_subgraph(Vertex v_global, Graph& g);
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
