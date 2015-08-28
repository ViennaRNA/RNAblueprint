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
        Vertex int_to_vertex(int i, Graph& g) {
            Vertex v = boost::vertex(i, g.root());
            if (i >= num_vertices(g) || i < 0 || v == Graph::null_vertex()) {
                std::stringstream ss;
                ss << "Error getting vertex descriptor from integer: " << i;
                throw std::out_of_range(ss.str());
            }
            return g.global_to_local(v);
        }
        
        // get the vertex_color_t tag from a vertex descriptor
        int vertex_to_int(Vertex v, Graph& g) {
            return boost::get(boost::vertex_color_t(), g.root(), g.local_to_global(v));
        }
        
        
        std::vector<int> getVertexList(Graph& g) {
            std::vector<int> result;

            BGL_FORALL_VERTICES(v, g, Graph) {
                result.push_back(vertex_to_int(v, g));
            }
            return result;
        }
    }
}
