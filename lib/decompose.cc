/* RNAblueprint
 * A program for designing RNA molecules.
 *
 * @date 25.03.2013
 * @author Stefan Hammer <s.hammer@univie.ac.at>
 * @copyright GPLv3
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
        bool decompose_graph(Graph& graph, RG& rand) {
            
            if (debug) {
                std::cerr << "root graph:" << std::endl;
                // print the just created subgraphs
                print_graph(graph, &std::cerr);
            }
            
            connected_components_to_subgraphs(graph); // get connected components and make subgraphs

            if (debug) {
                std::cerr << "subgraphs connected components:" << std::endl;
                // print the just created subgraphs
                print_subgraphs(graph, &std::cerr);
            }

            // iterate over all subgraphs (connected components)
            Graph::children_iterator cc, cc_end;
            for (boost::tie(cc, cc_end) = graph.children(); cc != cc_end; ++cc) {
                
                // in case that the graph is not bipartite return false
                if (!boost::is_bipartite(*cc)) {
                    return false;
                }
                // start a recursion
                decompose_recursion(*cc, rand);
            }
            
            // return that the dependency graph is bipartite
            return true;
        }
        
        template <typename RG>
        void decompose_recursion(Graph& g, RG& rand) {
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
                        print_subgraphs(g, &std::cerr);
                    }
                // this is a block
                } else {
                    // do the ear decomposition and create to subgraphs
                    ear_decomposition_to_subgraphs(g, rand, true);

                    if (debug) {
                        std::cerr << "subgraphs ear decomposition:" << std::endl;
                        // print the just created subgraphs
                        print_subgraphs(g, &std::cerr);
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
                decompose_recursion(*gc, rand);
            }
        }

        void connected_components_to_subgraphs(Graph& g) {

            // get list of connected components into the component vector
            // http://www.boost.org/doc/libs/1_53_0/libs/graph/doc/connected_components.html
            // http://www.boost.org/doc/libs/1_53_0/libs/graph/example/connected_components.cpp
            
            // components map
            std::map<Vertex, int> component;
            boost::associative_property_map< std::map<Vertex, int> > component_map(component);
            int num = boost::connected_components(g, component_map);

            if (debug) {
                std::cerr << "Number of connected components: " << num << std::endl;
            }
            
            // map with all subgraphs
            std::map<int, Graph*> component_graphs;
            // add all adjacent edges of the vertices to the right subgraph
            BGL_FORALL_VERTICES_T(v, g, Graph) {
                // if subgraph does not exist yet, create it
                std::map<int, Graph*>::iterator it = component_graphs.find(component[v]);
                if (it == component_graphs.end()) {
                    // create subgraph
                    component_graphs[component[v]] = &g.create_subgraph();
                    graph_property& gprop = boost::get_property(*component_graphs[component[v]], boost::graph_name);
                    gprop.type = 1;
                    gprop.id = component[v];
                }
                // add vertex
                boost::add_vertex(g.local_to_global(v), *component_graphs[component[v]]);
                // add all adjacent edges
                BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
                    boost::add_edge(g.local_to_global(e), *component_graphs[component[v]]);
                }
            }
        }

        void biconnected_components_to_subgraphs(Graph& g) {

            // get list of biconnected components into the component property map
            // http://www.boost.org/doc/libs/1_53_0/libs/graph/doc/biconnected_components.html
            // http://www.boost.org/doc/libs/1_38_0/libs/graph/example/biconnected_components.cpp
            std::map<Edge, int> component;
            boost::associative_property_map< std::map<Edge, int> > component_map(component);
            std::vector<Vertex> art_points;
            // call the algorithm
            unsigned int num;
            std::back_insert_iterator<std::vector<Vertex> > ap_it(art_points);
            boost::tie(num, ap_it) = boost::biconnected_components(g, component_map, std::back_inserter(art_points));
            
            if (debug) {
                std::cerr << "Number of biconnected components: " << num << std::endl;
                std::cerr << "Number of articulation points: " << art_points.size() << " ( ";
                for (auto elem : art_points) {
                    std::cerr << vertex_to_int(elem, g) << " ";
                }
                std::cerr << ")" << std::endl;
            }

            if (debug) {
                // iterate over all graph edges to print connected components table
                BGL_FORALL_EDGES_T(e, g, Graph) {
                    std::cerr << e << "\t" << "(" << vertex_to_int(boost::source(e, g), g) << ","
                            << vertex_to_int(boost::target(e, g), g) << ")"
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
            
            // write biconnected components into subgraphs:
            // new bcc id
            int j = 0;
            // map with subgraphs
            std::map<int, Graph*> bicomponent_graphs;
            // TODO fix bug
            // add all edges of the biconnected component to the right subgraph
            BGL_FORALL_EDGES_T(e, g, Graph) {
                std::map<int, Graph*>::iterator it = bicomponent_graphs.find(component[e]);
                if (it == bicomponent_graphs.end()) {
                    // create subgraph
                    bicomponent_graphs[component[e]] = &g.create_subgraph();
                    graph_property& gprop = boost::get_property(*bicomponent_graphs[component[e]], boost::graph_name);
                    gprop.type = 2;
                    gprop.id = j++;
                }
                // add edge to subgraph
                boost::add_edge(g.local_to_global(e), *bicomponent_graphs[component[e]]);
            }
        }

        void merge_biconnected_paths(Graph& g, Vertex p, Vertex v, std::map<Edge, int>& component, std::vector<Vertex>& art_points, int& nc) {
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
        void ear_decomposition_to_subgraphs(Graph& g, RG& rand, bool optimize_decomposition) {
            // blocks need to be decomposed into path. this can be done by Ear Decomposition
            
            // ear map (this is what we want to fill!)
            boost::property_map<Graph, int edge_property::*>::type em = boost::get(&edge_property::ear, g);
            
            //best pred map
            boost::vector_property_map<Vertex> best_pred(boost::num_vertices(g));
            // container to store attachment points
            std::vector<Vertex> att_points;
            
            int best_alpha = std::numeric_limits<int>::max();
            int best_beta = std::numeric_limits<int>::max();
            // TODO think of a better way to get the best spanning tree. 
            // for now do a stupid optimization
            int count = 0;
            while (true) {
                // predecessor map
                boost::vector_property_map<Vertex> pred(boost::num_vertices(g));
                // get a random spanning tree
                boost::random_spanning_tree(g, rand, boost::predecessor_map(pred));
                // if no optimization should be done, exit here
                if (!optimize_decomposition) {
                    if (debug)
                        std::cerr << "No spanning tree and ear decomposition optimization!" << std::endl;
                    best_pred = pred;
                    break;
                }
                // run the ear decomposition
                int num = boost::ear_decomposition(g, pred, em, std::back_inserter(att_points));
                
                // print the EarMap
                if (debug) {
                    BGL_FORALL_EDGES_T(e, g, Graph) {
                        std::cout << "(" << boost::source(e, g) << "/" << boost::target(e, g) << "): " << g[e].ear << std::endl;
                    }
                }
                // calculate alpha and beta
                int alpha, beta;
                boost::tie(alpha, beta) = get_alpha_beta(g, att_points, num);
                // clean up att_points
                att_points.clear();
                
                if ((alpha < best_alpha && beta <= best_beta) || (beta < best_beta && alpha <= best_alpha)) {
                    best_alpha = alpha;
                    best_beta = beta;
                    best_pred = pred;
                    count = 0;
                    if (debug)
                        std::cerr << "Better Solution: " << best_alpha << "/" << best_beta << std::endl;
                } else {
                    count++;
                    if (debug)
                        std::cerr << "Optimization Count: " << count << " - " << alpha << "/" << beta << std::endl;
                    // optimization exit
                    if (count > 500 || best_beta < 7 || best_alpha < 5)
                        break;
                }
            }
            // redo best ear decomposition
            // run the ear decomposition with best_pred
            if (debug) {
                std::cerr << "Best alpha/beta: " << best_alpha << "/" << best_beta << std::endl;
            }
            
            boost::ear_decomposition(g, best_pred, em, std::back_inserter(att_points));
            
            // create subgraphs from decomposed ears
            // map with subgraphs
            std::map<int, Graph*> ear_graphs;
            // add all edges of the ear to the right subgraph
            BGL_FORALL_EDGES_T(e, g, Graph) {
                std::map<int, Graph*>::iterator it = ear_graphs.find(g[e].ear);
                if (it == ear_graphs.end()) {
                    // create subgraph
                    ear_graphs[g[e].ear] = &g.create_subgraph();
                    graph_property& gprop = boost::get_property(*ear_graphs[g[e].ear], boost::graph_name);
                    gprop.type = 3;
                    gprop.id = g[e].ear - 1;
                }
                // add edge to subgraph
                boost::add_edge(g.local_to_global(e), *ear_graphs[g[e].ear]);
            }
            
            // attachment points are special vertices
            for (Vertex v : att_points) {
                g[v].special = true;
                if (debug) {
                    std::cout << "Vertex " << vertex_to_int(v, g) << " is a attachment point!" << std::endl;
                }
            }
        }
        
        std::pair<int, int> get_alpha_beta(Graph& g, std::vector<Vertex> att_points, int num) {
            // return values
            int alpha = 0;
            int beta = 0;
            // map to store all attachment points ordered by ear index
            // A[k] = { vertex, vertex, vertex }
            std::map<int, std::vector<Vertex> > Ak;
            Ak[0] = {}; // this is important, because for the first and last ear we have no attachment points
            Ak[num] = {};
            // fill the map 
            for (Vertex v : att_points) {
                int min_ear_index = std::numeric_limits<int>::max();
                int max_ear_index = 0;
                BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
                    if (min_ear_index > g[e].ear)
                        min_ear_index = g[e].ear;
                    else if (max_ear_index < g[e].ear)
                        max_ear_index = g[e].ear;
                }
                for (int k = min_ear_index; k < max_ear_index; k++)
                    Ak[k].push_back(v);
            }
            if (debug)
                std::cerr << "Ak: " << std::endl << Ak << std::endl;
            
            for (int k = 0; k < num; ++k) {
                // for making set_difference we need a sorted vector (Ak[0] is empty and therfore already sorted)
                std::sort(Ak[k+1].begin(), Ak[k+1].end());
                // get the Ak+1 but not Ak set
                std::vector<Vertex> complement(Ak[k].size() + Ak[k+1].size());
                std::vector<Vertex>::iterator c_it = std::set_difference(Ak[k+1].begin(), Ak[k+1].end(), Ak[k].begin(), Ak[k].end(), complement.begin());
                complement.resize(c_it-complement.begin());
                
                if (debug)
                    std::cerr << "complement (k=" << k << "): " << std::endl << complement << std::endl;
                // calculate alpha = max_k(|Ak|)
                int current_alpha = Ak[k].size();
                if (current_alpha > alpha)
                    alpha = current_alpha;
                // calculate beta = max_k(|Ak| + |Ak+1 \ Ak|)
                int current_beta = current_alpha + complement.size();
                if (current_beta > beta)
                    beta = current_beta;
            }
            return std::make_pair(alpha, beta);
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
                        boost::add_edge(g.local_to_global(e), *subgptr);
                        parts_recursion(g, subgptr, boost::source(e, g));
                        parts_recursion(g, subgptr, boost::target(e, g));
                        
                        graph_property& gprop = boost::get_property(*subgptr, boost::graph_name);
                        gprop.type = 4;
                        gprop.id = i++;
                        gprop.is_path = true;
                    }
                }
                if (debug) {
                    std::cerr << "subgraphs parts between specials:" << std::endl;
                    print_subgraphs(g, &std::cerr);
                }
            } else {
                // mark this subgraph as a path
                boost::get_property(g, boost::graph_name).is_path = true;
            }
        }

        void parts_recursion(Graph& g, Graph * subgptr, Vertex v) {
            // add vertex to subgraph
            //boost::add_vertex(boost::get(boost::vertex_color_t(), g, v), *subgptr);
            g[v].color = 1;
            
            if (!g[v].special) {
                BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
                    if (g[e].color == 0) {
                        g[e].color = 1;
                        boost::add_edge(g.local_to_global(e), *subgptr);
                        parts_recursion(g, subgptr, boost::target(e, g));
                    }
                }
            }
        }

        template bool decompose_graph<std::mt19937> (Graph&, std::mt19937&);
        template void ear_decomposition_to_subgraphs<std::mt19937> (Graph&, std::mt19937&, bool);
    }
}
