/* RNAblueprint
 * A program for designing RNA molecules.
 *
 * @date 25.03.2013
 * @author Stefan Hammer <s.hammer@univie.ac.at>
 * @copyright GPLv3
 *
 */

// include header
#include <vector>
#include <map>

#include "parsestruct.h"

namespace design {
    namespace detail {
        
        Graph parse_structures(std::vector<std::string> structures) {
            
            // sanity check for input
            if (structures.size() == 0) {
                throw ( std::logic_error( "Cannot initialize Dependency Graph with no structures!" ));
            }
            
            // remember cutpoints if cofold is inserted
            std::map<int, char>  cutpoints;
            unsigned int checksum = 0;
            for (auto& s : structures) {
                // remove ampersands and pluses and remember them as cut points
                std::size_t found_cut;
                while (true) {
                    found_cut = s.find_last_of("&+");
                    if (found_cut == std::string::npos)
                        break;
                    cutpoints[found_cut] = s[found_cut];
                    checksum += found_cut;
                    s.erase(found_cut, 1);
                }
            }
            // lambda calculating the sum of the cut point positions
            auto sum = [] (std::map<int, char> cp) {
                unsigned int sum = 0;
                for (auto c: cp)
                    sum += c.first;
                return sum;
            };
            // now check for invalid cut points
            if (checksum / structures.size() != sum(cutpoints))
                throw (std::logic_error("Cut points are not aligned properly or additional cut points!"));
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
            for (auto& s : structures) {
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
            gprop.cutpoints = cutpoints;
            
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
            set_constraints(g, constraints, true);
        }
        
        std::vector<int> set_constraints(Graph& g, std::string constraints, bool throwerror) {
            
            // remove ampersands and pluses and remember them as cut points
            std::size_t found_cut;
            while (true) {
                found_cut = constraints.find_last_of("&+");
                if (found_cut == std::string::npos)
                    break;
                constraints.erase(found_cut, 1);
            }
            
            for (unsigned int pos = 0; pos < constraints.length(); pos++) {                
                g[int_to_vertex(pos, g)].constraint = char_to_enum(std::toupper(constraints[pos]));

                // set constraints other than N to articulation
                if (g[int_to_vertex(pos, g)].constraint != N) {
                    g[int_to_vertex(pos, g)].articulation = true;
                }
            }
            // check if constraints are compatible with structures
            PairingMatrix * p = PairingMatrix::Instance();
            std::vector<int> incompatible;
            
            BGL_FORALL_EDGES_T(e, g, Graph) {
                Vertex v = boost::source(e, g);
                Vertex u = boost::target(e, g);
                int vc = g[v].constraint;
                int uc = g[u].constraint;
                if (p->get(1, vc, uc) == 0) {
                    incompatible.push_back(vertex_to_int(v, g));
                    incompatible.push_back(vertex_to_int(u, g));
                }
            }
            
            if (throwerror) {
                if (incompatible.size() != 0) {
                    std::stringstream ss;
                    ss << "Error while setting the constraints: These constraints are not compatible to the structures!" << std::endl;
                    ss << incompatible << std::endl;
                    throw std::logic_error(ss.str());
                }
            }
            
            return incompatible;
        }
    }
}