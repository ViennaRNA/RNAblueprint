/* RNAdesign
 * A program for designing RNA molecules.
 *
 * Created on: 25.03.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

// include header
#include "decompose.h"
#include "printgraph.h"

// include boost components
#include <boost/graph/iteration_macros.hpp>

namespace design {
    namespace detail {

        template <typename RG>
        bool decompose_graph(Graph& graph, RG* rand_ptr) {
            
            if (debug) {
                std::cerr << "root graph:" << std::endl;
                // print the just created subgraphs
                print_graph(graph, &std::cerr, "root");
            }
            
            connected_components_to_subgraphs(graph); // get connected components and make subgraphs

            if (debug) {
                std::cerr << "subgraphs connected components:" << std::endl;
                // print the just created subgraphs
                print_subgraphs(graph, &std::cerr, "connected-component");
            }

            // iterate over all subgraphs (connected components)
            Graph::children_iterator cc, cc_end;
            for (boost::tie(cc, cc_end) = graph.children(); cc != cc_end; ++cc) {
                
                // in case that the graph is not bipartite return false
                if (!boost::is_bipartite(*cc)) {
                    return false;
                }
                // start a recursion
                decompose_recursion(*cc, rand_ptr);
            }
            
            // return that the dependency graph is bipartite
            return true;
        }
        
        template <typename RG>
        void decompose_recursion(Graph& g, RG* rand_ptr) {
            // calculate the max degree of this graph
            int max_degree;
            int min_degree;
            std::tie(min_degree, max_degree) = get_min_max_degree(g);
            
            if (debug) {
                std::cerr << "Max degree of subgraph is: " << max_degree << std::endl;
                std::cerr << "Min degree of subgraph is: " << min_degree << std::endl;
            }
            
            // this is either a biconnected component or a block
            if (max_degree > 2) {
                std::vector<Vertex> art_points;
                boost::articulation_points(g, std::back_inserter(art_points));
                // this is a biconnected component
                if (art_points.size() != 0) {
                    biconnected_components_to_subgraphs(g);

                    if (debug) {
                        std::cerr << "subgraphs biconnected components:" << std::endl;
                        // print the just created subgraphs
                        print_subgraphs(g, &std::cerr, "biconnected-component");
                    }
                // this is a block
                } else {
                    // do the ear decomposition and create to subgraphs
                    ear_decomposition_to_subgraphs(g, rand_ptr);

                    if (debug) {
                        std::cerr << "subgraphs ear decomposition:" << std::endl;
                        // print the just created subgraphs
                        print_subgraphs(g, &std::cerr, "decomposed-ear");
                    }
                }
            // this is a circle
            } else if (min_degree == 2 && max_degree == 2) {
                // check if there are already special vertices
                int count = 0;
                BGL_FORALL_VERTICES_T(v, g, Graph) {
                    if (g[v].special) {
                        count++;
                    }
                }
                if (count < 2) {
                    // assign any two special vertices and get paths in between
                    Vertex s = boost::vertex(0, g);
                    Vertex r = boost::vertex(boost::num_vertices(g)-1, g);
                    (g)[s].special = true;
                    (g)[r].special = true;
                }
                parts_between_specials_to_subgraphs(g);
                return;
            // this is a path or a single vertex
            } else {
                parts_between_specials_to_subgraphs(g);
                return;
            }
            
            // call recursion for all children
            Graph::children_iterator gc, gc_end;
            for (boost::tie(gc, gc_end) = g.children(); gc != gc_end; ++gc) {
                decompose_recursion(*gc, rand_ptr);
            }
        }

        void connected_components_to_subgraphs(Graph& g) {

            // get list of connected components into the component vector
            // http://www.boost.org/doc/libs/1_53_0/libs/graph/doc/connected_components.html
            // http://www.boost.org/doc/libs/1_53_0/libs/graph/example/connected_components.cpp
            std::vector<int> component(boost::num_vertices(g));
            int num = boost::connected_components(g, &component[0]);

            if (debug) {
                std::cerr << "Number of connected components: " << num << std::endl;
                std::cerr << component << std::endl;
            }

            // iterate over connected component numbers
            for (int i = 0; i != num; ++i) {
                // for each component number generate a new subgraph
                Graph& subg = g.create_subgraph();
                
                int vertex = 0;
                // iterate over elements of connected_components
                for (auto elem : component) {
                    if (i == elem) {
                        // add vertex into current subgraph
                        boost::add_vertex(int_to_vertex(vertex, g), subg);
                    }
                    vertex++;
                }
                
                graph_property& gprop = boost::get_property(subg, boost::graph_name);
                gprop.type = 1;
                gprop.id = i;
            }
        }

        void biconnected_components_to_subgraphs(Graph& g) {

            // get list of biconnected components into the component property map
            // http://www.boost.org/doc/libs/1_53_0/libs/graph/doc/biconnected_components.html
            // http://www.boost.org/doc/libs/1_38_0/libs/graph/example/biconnected_components.cpp
            edge_component_t edge_component;
            boost::property_map < Graph, edge_component_t >::type component = boost::get(edge_component, g);
            unsigned int num = boost::biconnected_components(g, component);
            if (debug) {
                std::cerr << "Number of biconnected components: " << num << std::endl;
            }

            std::vector<Vertex> art_points;
            boost::articulation_points(g, std::back_inserter(art_points));
            if (debug) {
                std::cerr << "Number of articulation points: " << art_points.size() << " ( ";
                for (auto elem : art_points) {
                    std::cerr << boost::get(boost::vertex_color_t(), g, elem) << " ";
                }
                std::cerr << ")" << std::endl;
            }

            if (debug) {
                // iterate over all graph edges to print connected components table

                BGL_FORALL_EDGES_T(e, g, Graph) {
                    std::cerr << e << "\t" << "(" << boost::get(boost::vertex_color_t(), g, boost::source(e, g)) << ","
                            << boost::get(boost::vertex_color_t(), g, boost::target(e, g)) << ")"
                            << "\tcomponent: " << component[e] << std::endl;
                }
            }

            // now need to merge biconnected components that are separated by a articulation point that has a degree == 2 !
            for (auto v : art_points) {
                if (boost::degree(v, g) > 2) {
                    // mark this vertex as special point
                    g[v].special = true;
                    
                    BGL_FORALL_ADJ_T(v, adj, g, Graph) {
                        if ((boost::degree(adj, g) == 2)
                                && (std::find(art_points.begin(), art_points.end(), adj) != art_points.end())) {
                            int nc = -1;
                            merge_biconnected_paths(g, v, adj, component, art_points, nc);
                        }
                    }
                }
            }
            int j = 0;
            // write biconnected components into subgraphs:
            for (unsigned int i = 0; i != num; i++) {
                // only create a subgraph if there is really an edge associated to this component (as we merged many components before)
                bool exists = false;

                BGL_FORALL_EDGES_T(e, g, Graph) {
                    if (i == component[e])
                        exists = true;
                }
                if (!exists)
                    continue;

                // for this bicomponent number generate a new subgraph
                Graph& subg = g.create_subgraph();
                
                // iterate over edges of graph
                BGL_FORALL_EDGES_T(e, g, Graph) {
                    if (i == component[e]) {
                        // add vertex into current subgraph if not present already
                        if (!subg.find_vertex(boost::get(boost::vertex_color_t(), g, boost::target(e, g))).second) {
                            boost::add_vertex(boost::get(boost::vertex_color_t(), g, boost::target(e, g)), subg);
                        }
                        if (!subg.find_vertex(boost::get(boost::vertex_color_t(), g, boost::source(e, g))).second) {
                            boost::add_vertex(boost::get(boost::vertex_color_t(), g, boost::source(e, g)), subg);
                        }
                    }
                }
                // add properties
                graph_property& gprop = boost::get_property(subg, boost::graph_name);
                gprop.type = 2;
                gprop.id = j++;
            }
        }

        void merge_biconnected_paths(Graph& g, Vertex p, Vertex v, boost::property_map < Graph, edge_component_t >::type& component, std::vector<Vertex>& art_points, int& nc) {
            if (debug) {
                std::cerr << "Merging biconnected paths..." << std::endl;
            }

            BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
                if (nc == -1)
                    nc = component[e];
                else
                    component[e] = nc;
            }

            BGL_FORALL_ADJ_T(v, adj, g, Graph) {
                if ((adj != p) && (boost::degree(adj, g) == 2)
                        && (std::find(art_points.begin(), art_points.end(), adj) != art_points.end()))
                    merge_biconnected_paths(g, v, adj, component, art_points, nc);
            }
        }

        template <typename RG>
        void ear_decomposition_to_subgraphs(Graph& g, RG* rand_ptr) {
            // blocks need to be decomposed into path. this can be done by Ear Decomposition

            // predecessor map
            boost::vector_property_map<Vertex> pred(boost::num_vertices(g));
            // get a random spanning tree
            boost::random_spanning_tree(g, *rand_ptr, boost::predecessor_map(pred));
            // ear map (this is what we want to fill!)
            auto em = boost::get(&edge_property::ear, g);

            // run the ear decomposition
            int num = boost::ear_decomposition(g, pred, em);
            // print the EarMap
            if (debug) {

                BGL_FORALL_EDGES_T(e, g, Graph) {
                    std::cout << "(" << boost::source(e, g) << "/" << boost::target(e, g) << "): " << g[e].ear << std::endl;
                }
            }

            // create subgraphs from decomposed ears

            for (int i = 1; i != num + 1; ++i) {
                // for each ear create a new subgraph
                Graph& subg = g.create_subgraph();

                BGL_FORALL_EDGES_T(e, g, Graph) {
                    if (g[e].ear == i) {
                        if (!subg.find_vertex(g.local_to_global(boost::source(e, g))).second) {
                            boost::add_vertex(g.local_to_global(boost::source(e, g)), subg);
                        }
                        if (!subg.find_vertex(g.local_to_global(boost::target(e, g))).second) {
                            boost::add_vertex(g.local_to_global(boost::target(e, g)), subg);
                        }
                    }
                }
                
                graph_property& gprop = boost::get_property(subg, boost::graph_name);
                gprop.type = 3;
                gprop.id = i-1;
            }
            // detect attachment points and push them into the graph as vertex property Ak
            color_attachment_points(g);
        }

        void color_attachment_points(Graph& g) {
            if (debug) {
                std::cerr << "Color attachment points as special vertices..." << std::endl;
            }

            BGL_FORALL_VERTICES_T(v, g, Graph) {
                int prev_ear = -1;

                BGL_FORALL_OUTEDGES(v, e, g, Graph) {
                    if (prev_ear == -1) {
                        prev_ear = g[e].ear;
                    } else if (prev_ear != g[e].ear) {
                        // set special to true!
                        g[v].special = true;
                        if (debug) {
                            std::cout << "Vertex " << g.local_to_global(v) << " is a attachment point!" << std::endl;
                        }
                        break;
                    }
                }
            }
        }

        void parts_between_specials_to_subgraphs(Graph& g) {
            bool split = false;
            BGL_FORALL_VERTICES_T(v, g, Graph) {
                split = split || (g[v].special && (boost::degree(v, g) > 1));
            }
            if (debug && !split) {
                std::cerr << "No need to generate a subgraph as this is already a path with specials only on ends." << std::endl;
            }
            
            if (split) {
                if (debug) {
                    std::cerr << "Paths between specials to subgraphs..." << std::endl;
                }
                // reset edge colors

                BGL_FORALL_EDGES_T(e, g, Graph) {
                    g[e].color = 0;
                    g[boost::source(e, g)].color = 0;
                    g[boost::target(e, g)].color = 0;
                }

                // start recursion at all vertices of edges
                int i = 0;
                BGL_FORALL_EDGES_T(e, g, Graph) {
                    if (g[e].color == 0) {
                        g[e].color = 1;
                        Graph* subgptr = &g.create_subgraph();
                        parts_recursion(g, subgptr, boost::source(e, g));
                        parts_recursion(g, subgptr, boost::target(e, g));
                        
                        graph_property& gprop = boost::get_property(*subgptr, boost::graph_name);
                        gprop.type = 4;
                        gprop.id = i++;
                        gprop.is_path = true;
                    }
                }
                if (debug) {
                    print_subgraphs(g, &std::cerr, "parts-between-specials");
                }
            } else {
                // mark this subgraph as a path
                boost::get_property(g, boost::graph_name).is_path = true;
            }
        }

        void parts_recursion(Graph& g, Graph * subgptr, Vertex v) {
            // add vertex to subgraph
            boost::add_vertex(boost::get(boost::vertex_color_t(), g, v), *subgptr);
            g[v].color = 1;
            
            if (!g[v].special) {
                BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
                    if (g[e].color == 0) {
                        g[e].color = 1;
                        parts_recursion(g, subgptr, boost::target(e, g));
                    }
                }
            }
        }

        template bool decompose_graph<std::mt19937> (Graph&, std::mt19937*);
        template void ear_decomposition_to_subgraphs<std::mt19937> (Graph&, std::mt19937*);
    }
}
