/* RNAblueprint
 * A program for designing RNA molecules.
 *
 * Created on: 25.03.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

// include header
#include "printgraph.h"

namespace design {
    namespace detail {

        void print_graph(Graph& g, std::ostream* out) {
            // to convert enums of bases to chars
            std::map<Vertex, char> bases, constraints;
            boost::associative_property_map< std::map<Vertex, char> > base_map(bases), constraint_map(constraints);
            // to convert list of articulation points to string
            std::map< Vertex, std::string > aps;
            boost::associative_property_map< std::map<Vertex, std::string> > ap_map(aps);

            BGL_FORALL_VERTICES_T(v, g, Graph) {
                bases.insert(std::make_pair(v, enum_to_char(g[v].base)));
                constraints.insert(std::make_pair(v, enum_to_char(g[v].constraint)));
            }

            // print vertex and edge properties from my self-defined bundled properties
            boost::dynamic_properties dp;
            
            //dp.property("id", boost::get(&graph_property::id, g));
            //dp.property("path", boost::get(&graph_property::is_path, g));
            dp.property("base", base_map);
            dp.property("constraints", constraint_map);
            dp.property("name", boost::get(boost::vertex_color_t(), g));
            dp.property("ear", boost::get(&edge_property::ear, g));
            dp.property("special", boost::get(&vertex_property::special, g));
            
            boost::write_graphml(*out, g, dp, true);
            if (debug) {
                std::cerr << "created graphml!" << std::endl;
            }
        }

        void print_subgraphs(Graph& g, std::ostream* out) {
            Graph::children_iterator ci, ci_end;
            for (boost::tie(ci, ci_end) = g.children(); ci != ci_end; ++ci) {
                print_graph(*ci, out);
            }
            if (debug) {
                std::cerr << "Printed all subgraphs." << std::endl;
            }
        }
    }
}