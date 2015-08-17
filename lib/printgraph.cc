/* RNAdesign
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

        void print_graph(Graph& g, std::ostream* out, std::string nametag) {
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

            boost::ref_property_map<Graph*, std::string> graphName(boost::get_property(g, boost::graph_name).id);
            dp.property("name", graphName);
            boost::ref_property_map < Graph*, bool> graphPath(boost::get_property(g, boost::graph_name).is_path);
            dp.property("path", graphPath);

            dp.property("base", base_map);
            dp.property("constraints", constraint_map);
            dp.property("name", boost::get(boost::vertex_color_t(), g));
            dp.property("ear", boost::get(&edge_property::ear, g));
            dp.property("special", boost::get(&vertex_property::special, g));

            /*if (outfile != "") {
              std::stringstream filename;
              filename << outfile << "-" << nametag << ".graphml";
              std::ofstream graphfile(filename.str());
              if (graphfile.is_open()) {
                boost::write_graphml(graphfile, g, dp, true);
                graphfile << std::endl;
                graphfile.close();
             *out << " done!" << std::endl;
              } else {
                std::cerr << " Unable to create graphml file!" << std::endl;
              }
            } else {*/
            *out << std::endl;
            boost::write_graphml(*out, g, dp, true);
            *out << std::endl;
            std::cerr << "created graphml!" << std::endl;
            //}
        }

        void print_subgraphs(Graph& g, std::ostream* out, std::string nametag) {
            Graph::children_iterator ci, ci_end;
            int num = 1;
            for (boost::tie(ci, ci_end) = g.children(); ci != ci_end; ++ci) {
                std::stringstream name;
                name << nametag << "-" << num++;
                *out << name.str() << ":";
                print_graph(*ci, out, name.str());
            }
            if (debug) {
                std::cerr << "Printed all subgraphs." << std::endl;
            }
        }

        void print_all_vertex_names(Graph& g, std::string nametag) {
            if (debug) {
                std::cerr << nametag << ": ";

                BGL_FORALL_VERTICES_T(v, g, Graph) {
                    std::cerr << boost::get(boost::vertex_color_t(), g, v) << "(" <<
                            enum_to_char(g[v].base) << ") ";
                }
                std::cerr << std::endl;
            }
        }
    }
}