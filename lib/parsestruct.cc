/* RNAdesign
 * A program for designing RNA molecules.
 *
 * Created on: 25.03.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

// include header
#include "parsestruct.h"

namespace design {
    namespace detail {

        Graph parse_structures(std::vector<std::string> structures) {

            // count the number of positions
            int num_vertices = structures[0].length();
            if (debug) {
                std::cerr << "Generating Graph with " << num_vertices << " vertices." << std::endl;
            }
            Graph g(num_vertices);

            // give the vertices names
            int vertex_name = 0;

            BGL_FORALL_VERTICES_T(v, g, Graph) {
                boost::put(boost::vertex_color_t(), g, v, vertex_name++);
            }

            // iterate over structures from input
            for (auto elem : structures) {
                std::vector<int> pair_table(structures[0].length(), 0); // remember position of the open bracket
                unsigned int open = 0; // remember how many open brackets there are
                // iterate over characters from structure
                for (unsigned int pos = 0; pos < elem.length(); pos++) {
                    if (elem[pos] == '(') {
                        pair_table[open] = pos;
                        if (debug) {
                            std::cerr << elem[pos] << ", open count: " << open;
                        }
                        open++;
                    } else if (elem[pos] == ')') {
                        open--;
                        // check if edge already exists
                        bool exists_ab = boost::edge(boost::vertex(pair_table[open], g), boost::vertex(pos, g), g).second;
                        bool exists_ba = boost::edge(boost::vertex(pos, g), boost::vertex(pair_table[open], g), g).second;
                        if (!exists_ab && !exists_ba) {
                            // add edge
                            boost::add_edge(boost::vertex(pair_table[open], g), boost::vertex(pos, g), g);
                        }
                        // reset value
                        pair_table[open] = pos;
                        if (debug) {
                            std::cerr << elem[pos] << ", open count: " << open;
                        }
                    } else if (elem[pos] != '.') {
                        throw ( std::logic_error( "Unknown character in dot bracket representation" ));
                    }
                    // error handling: there can't be more closing brackets than opening ones
                    if (open < 0) {
                        throw ( std::logic_error( "Unbalanced brackets in make_pair_table" ));
                    }
                    if (debug) {
                        std::cerr << " pos count:" << pos << std::endl;
                    }
                }
                // error handling: at the end all brackets must be closed again!
                if (open != 0) {
                    throw ( std::logic_error( "Too few closed brackets in make_pair_table" ));
                }
            }
            
            // label graph as root
            boost::get_property(g, boost::graph_name).id = "root_graph";
            
            return g;
        }
        
        void set_constraints(Graph& g, std::string constraints) {
            for (int pos = 0; pos < constraints.length(); pos++) {
                g[int_to_vertex(pos, g)].constraint = char_to_enum(constraints[pos]);
                
                // set constraints other than N to special
                if (g[int_to_vertex(pos, g)].constraint != N) {
                    g[int_to_vertex(pos, g)].special = true;
                }
            }
        }
    }
}