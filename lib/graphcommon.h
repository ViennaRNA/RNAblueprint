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

        struct edge_component_t {

            enum {
                num = 555
            };
            typedef boost::edge_property_tag kind;
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
            bool is_path = false;
            
        };

        // graph_properties 
        // boost graph template
        typedef boost::subgraph<
        boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS,
        boost::property< boost::vertex_color_t, int, vertex_property >,
        boost::property< boost::edge_index_t, int, boost::property < edge_component_t, std::size_t, edge_property> >,
        boost::property< boost::graph_name_t, graph_property> > > Graph;
        typedef Graph::edge_descriptor Edge;
        typedef Graph::vertex_descriptor Vertex;

        // get max degree of a graph
        std::pair <int, int> get_min_max_degree(Graph& g);

        // get the vertex descriptor from a vertex_color_t tag
        Vertex int_to_vertex(int i, Graph& g);
        // get the vertex_color_t tag from a vertex descriptor
        int vertex_to_int(Vertex v, Graph& g);
        // get the set of vertices for a certain (sub)graph
        std::vector<int> getVertexList(Graph& g);
    }
}
#endif


/* 
 * \endcond INTERNAL
 */