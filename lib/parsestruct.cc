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
            
            BracketList brackets = {
                {'(',')'},
                {'{','}'},
                {'[',']'},
                {'<','>'}
            };
            
            // iterate over structures from input
            for (auto s : structures) {
                std::size_t found_illegal = s.find_first_not_of("().[]{}<>");
                if (found_illegal != std::string::npos) {
                    std::stringstream ss;
                    ss << "Unknown character [" << s[found_illegal] << "] in structure: " << s << std::endl;
                    throw ( std::logic_error( ss.str() ));
                }
                for (BracketList::iterator br = brackets.begin(); br != brackets.end(); ++br) {
                    // if this is a valid bracket, start to add these edges to the graph too
                    parse_bracket(g, s, br);
                }
            }
            
            // label graph as root
            graph_property& gprop = boost::get_property(g, boost::graph_name);
            gprop.type = 0;
            gprop.id = 0;
            
            return g;
        }
        
        void parse_bracket(Graph& g, std::string& structure, BracketList::iterator bracket) {
            std::vector<int> pair_table(structure.length(), 0); // remember position of the open bracket
            unsigned int open = 0; // remember how many open brackets there are
            // iterate over characters from structure
            for (unsigned int pos = 0; pos < structure.length(); pos++) {
                if (structure[pos] == bracket->first) {
                    pair_table[open] = pos;
                    if (debug) {
                        std::cerr << structure[pos] << ", open count: " << open;
                    }
                    open++;
                } else if (structure[pos] == bracket->second) {
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
                        std::cerr << structure[pos] << ", open count: " << open;
                    }
                }

                // error handling: there can't be more closing brackets than opening ones
                if (open < 0) {
                    throw ( std::logic_error( "Unbalanced brackets in make_pair_table" ));
                }
            }
            // error handling: at the end all brackets must be closed again!
            if (open != 0) {
                throw ( std::logic_error( "Too few closed brackets in make_pair_table" ));
            }
        }
        
        void set_constraints(Graph& g, std::string constraints) {
            for (unsigned int pos = 0; pos < constraints.length(); pos++) {
                g[int_to_vertex(pos, g)].constraint = char_to_enum(std::toupper(constraints[pos]));
                
                // set constraints other than N to special
                if (g[int_to_vertex(pos, g)].constraint != N) {
                    g[int_to_vertex(pos, g)].special = true;
                }
            }
        }
    }
}