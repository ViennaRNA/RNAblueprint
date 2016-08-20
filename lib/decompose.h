/*!\file decompose.h 
 * \brief This file holds the functions to decompose the dependency graph into its subgraphs.
 *
 * @date 25.03.2013
 * @author Stefan Hammer <s.hammer@univie.ac.at>
 * @copyright GPLv3
 *
 * \cond INTERNAL
 */

#ifndef DECOMPOSE_H
#define DECOMPOSE_H


// include common header with graph definition and global variables
#include "common.h"
#include "graphcommon.h"

// include standard library parts
#include <limits>
#include <iterator>
#include <list>

// include boost components
#include <boost/graph/bipartite.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/random_spanning_tree.hpp>
#include <boost/property_map/shared_array_property_map.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/vector_property_map.hpp>
// #include <boost/graph/ear_decomposition.hpp>
#include "ear-decomposition.hpp"

namespace design {
    namespace detail {

        // does the graph decomposition, and calls the coloring of the subgraphs
        template <typename RG>
        bool decompose_graph(Graph& graph, RG& rand);
        // does the graph decomposition recursion
        template <typename RG>
        void decompose_recursion(Graph& g, RG& rand);
        
        // get a vector of all vertices with their component id. finds connected components with DFS
        void connected_components_to_subgraphs(Graph& g);

        // finds biconnected components with DFS
        void biconnected_components_to_subgraphs(Graph& g);

        // starts at a degree>3 articulation point and walks along a path to connect it to one component
        void merge_biconnected_paths(Graph& g, Vertex p, Vertex v, std::map<Edge, int>& component, std::vector<Vertex>& art_points, int& nc);

        // ear decomposition of blocks
        template <typename RG>
        void ear_decomposition_to_subgraphs(Graph& g, RG& rand, bool optimize_decomposition);
        
        // calculate alpha and beta values as measure for the complexity of this ear_decomposition
        std::pair<int, int> get_alpha_beta(Graph& g, std::vector<Vertex> att_points, int num);

        // identify attachment points and add them as graph property Ak
        void color_attachment_points(Graph& g);

        // now create subgraphs for the parts between Ak and Iks
        void parts_between_specials_to_subgraphs(Graph& g);

        // recursion for parts function
        void parts_recursion(Graph& g, Graph * subgptr, Vertex v);
    }
}
#endif


/* 
 * \endcond INTERNAL
 */
