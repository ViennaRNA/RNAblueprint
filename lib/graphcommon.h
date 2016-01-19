/*!\file graphcommon.h 
 * \brief This file holds all important information for the dependency graph, its definition and often used functions.
 *
 * Created on: 26.08.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 * \cond INTERNAL
 */

#ifndef GRAPHCOMMON_H
#define GRAPHCOMMON_H

// include common header with graph definition and global variables
#include "common.h"
#include "probability_matrix.h"

// include standard library parts
#include <limits>
#include <sstream>

namespace design {
    namespace detail {

        // get property with Graph[Vertex].bipartite_color = int;

        struct vertex_property {
            int color = 0;
            int base = N;
            int constraint = N;
            bool special = false;
        };

        struct edge_property {
            int ear = 0;
            int color = 0;
        };

        struct graph_property {
            // concurrent number from 0 to # of components-1
            int id;
            /*
             0 -> root
             1 -> cc
             2 -> bc
             3 -> ear
             4 -> path
             */
            int type;
            SolutionSizeType nos;
            bool is_path = false;
        };

        // graph_properties 
        // boost graph template
        typedef boost::uninduced_subgraph<
        boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS,
        boost::property< boost::vertex_color_t, int, vertex_property >,
        boost::property< boost::edge_index_t, int, edge_property >,
        boost::property< boost::graph_name_t, graph_property> > > Graph;
        typedef Graph::edge_descriptor Edge;
        typedef Graph::vertex_descriptor Vertex;
        
        // get max degree of a graph
        inline std::pair <int, int> get_min_max_degree(Graph& g) {
            int max_degree = 0;
            int min_degree = std::numeric_limits<int>::max();
            
            BGL_FORALL_VERTICES_T(v, g, Graph) {
                int current_degree = boost::out_degree(v, g);
                if (current_degree > max_degree) {
                    max_degree = current_degree;
                }
                if (current_degree < min_degree) {
                    min_degree = current_degree;
                }
            }
            
            return std::make_pair(min_degree, max_degree);
        }
        
        // get the vertex descriptor from a vertex_color_t tag
        inline Vertex int_to_vertex(unsigned int i, Graph& g) {
            Vertex v = boost::vertex(i, g.root());
            if (i >= boost::num_vertices(g) || i < 0 || v == Graph::null_vertex()) {
                std::stringstream ss;
                ss << "Error getting vertex descriptor from integer: " << i;
                throw std::out_of_range(ss.str());
            }
            return g.global_to_local(v);
        }
        
        // get the vertex_color_t tag from a vertex descriptor
        inline int vertex_to_int(Vertex v, Graph& g) {
            return boost::get(boost::vertex_color_t(), g.root(), g.local_to_global(v));
        }
        
        // get the set of vertices for a certain (sub)graph
        inline std::vector<int> getVertexList(Graph& g) {
            std::vector<int> result;

            BGL_FORALL_VERTICES(v, g, Graph) {
                result.push_back(vertex_to_int(v, g));
            }
            return result;
        }
    }
}
#endif


/* 
 * \endcond INTERNAL
 */