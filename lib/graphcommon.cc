/* RNAdesign
 * A program for designing RNA molecules.
 *
 * Created on: 26.08.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

// include header
#include "graphcommon.h"

namespace design {
    namespace detail {

        std::pair <int, int> get_min_max_degree(Graph& g) {
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
        Vertex int_to_vertex(int i, Graph g) {
            
            Graph root = g.root();
            Vertex v = boost::vertex(i, root);
            return g.global_to_local(v);
            /*
            BGL_FORALL_VERTICES_T(v, g, Graph) {
                if (boost::get(boost::vertex_color_t(), g, v) == i) {
                    return v;
                }
            }
            std::cerr << "This vertex is not present in this graph!";
            exit(1);
            */
        }
        
        // get the vertex_color_t tag from a vertex descriptor
        int vertex_to_int(Vertex v, Graph g) {
            return boost::get(boost::vertex_color_t(), g, v);
        }
    }
}
