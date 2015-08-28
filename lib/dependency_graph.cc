/* RNAdesign
 * A program for designing RNA molecules.
 *
 * Created on: 18.03.2014
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

// include header
#include "dependency_graph.h"

#include "parsestruct.h"
#include "printgraph.h"
#include "decompose.h"

// include boost components
#include <boost/graph/iteration_macros.hpp>

namespace design
{

    namespace detail{

        template <typename R>
        DependencyGraph<R>::DependencyGraph(std::vector<std::string> structures, std::string constraints, R _rand)
        : rand(_rand) {

            if (debug) {
                std::cerr << "Initializing DependencyGraph..." << std::endl;
            }
            
            // initialize random generator
            rand_ptr = &rand;

            // generate graph from input vector
            try {
                graph = parse_structures(structures);
            } catch (std::exception& e) {
                std::stringstream ss;
                ss << "Error while parsing the structures: " << std::endl << e.what();
                std::cerr << ss.str() << std::endl;
                throw std::logic_error(ss.str());
            }
            
            // set sequence constraints
            set_constraints(graph, constraints);

            // decompose the graph into its connected components, biconnected
            // components and decompose blocks via ear decomposition
            try {
                bipartite = decompose_graph(graph, rand_ptr);
            } catch (std::exception& e) {
                std::stringstream ss;
                ss << "Error while decomposing the dependency graph: " << std::endl << e.what();
                std::cerr << ss.str() << std::endl;
                throw std::logic_error(ss.str());
            }
            // trow an exception for now if graph is not bipartite
            if (!bipartite) {
                throw std::logic_error("Graph is not bipartite! No solution exists therefore.");
            }
            
            // now calculate all the PMs
            try {
                calculate_probabilities(graph);
            } catch (std::exception& e) {
                std::stringstream ss;
                ss << "Error while calculating the probabilities: " << std::endl << e.what();
                std::cerr << ss.str() << std::endl;
                throw std::logic_error(ss.str());
            }
        }

        template <typename R>
        void DependencyGraph<R>::calculate_probabilities(Graph& g) {
            
            // Remember a temporary PM which holds the current state
            ProbabilityMatrix current;
            std::unordered_map<Vertex, int> degree_map;
            
            Graph::children_iterator cg, cg_end;
            for (boost::tie(cg, cg_end) = g.children(); cg != cg_end; ++cg) {
                
                std::cerr << "Graph (" << boost::get_property(*cg, boost::graph_name).type << "-" << boost::get_property(*cg, boost::graph_name).nummer << "):" << std::endl;
                std::cerr << "current graph storage location: " << &(*cg) << std::endl;
                std::cerr << "current graph property storage location: " << &boost::get_property(*cg, boost::graph_name) << std::endl;
                
                
                if (debug) {
                    std::cerr << "current graph path? " << boost::get_property(*cg, boost::graph_name).is_path << std::endl;
                    print_graph(*cg, &std::cerr, "current graph");
                }
                
                if (boost::get_property(*cg, boost::graph_name).is_path) {
                    // this is a path and therefore needs to be treated separately
                    // calculate PM for Path and save to subgraph
                    try {
                        pms[&*cg] = get_path_pm(*cg);
                    } catch (std::exception& e) {
                        std::stringstream ss;
                        ss << "Could not get a ProbabilityMatrix for a path: " << std::endl << e.what();
                        std::cerr << ss.str() << std::endl;
                        throw std::logic_error( ss.str() );
                    }
                    if (debug) {
                        std::cerr << "Path PM (" << boost::get_property(*cg, boost::graph_name).type << "-" << boost::get_property(*cg, boost::graph_name).nummer << "):" << std::endl
                                << pms[&*cg] << std::endl;
                    }
                } else {
                    if (debug) {
                        std::cerr << "Recursion!" << std::endl;
                    }
                    // recursion is here
                    calculate_probabilities(*cg);
                }
                
                // Multiply current with pm of this child
                current = pms[&*cg] * current;
                if (debug) {
                    std::cerr << "current PM: " << std::endl << current << std::endl;
                }
                
                // now make nodes internal and remember the current pm at such nodes (except for cc as those are independent of each other)
                if (!boost::get_property(*cg, boost::graph_name).is_cc) {
                    // only save PM before first make_internal
                    bool alreadysaved = false;
                    
                    BGL_FORALL_VERTICES_T(v, *cg, Graph) {
                        Vertex global_v = cg->local_to_global(v);
                        
                        if ((*cg)[v].special) {
                            // update current degrees as status of special points
                            if (debug) {
                                std::cerr << "updating degree: " << degree_map[global_v] << " + " << boost::degree(v, *cg) << std::endl;
                            }
                            degree_map[global_v] += boost::degree(v, *cg);

                            // check if a vertex becomes internal here
                            if (debug) {
                                std::cerr << "v" << vertex_to_int(v, *cg) << " internal? " << degree_map[global_v] << "/" << boost::degree(global_v, cg->root()) << std::endl;
                            }
                            if (degree_map[global_v] == boost::degree(global_v, cg->root())) {
                                // remember this PM here
                                // only the first time a internal node is detected, otherwise we overwrite this
                                if (!alreadysaved) {
                                    pms[&*cg] = ProbabilityMatrix(current);
                                    alreadysaved = true;
                                    if (debug) {
                                        std::cerr << "saved PM (" << boost::get_property(*cg, boost::graph_name).type << "-" << boost::get_property(*cg, boost::graph_name).nummer << "):"  << std::endl
                                                << pms[&*cg] << std::endl;
                                    }
                                }
                                // remove internal special vertex from this PM!
                                current = make_internal(current, vertex_to_int(v, *cg));
                                //TODO is there a need to remove special flag if this vertex is now internal?
                            }
                        }
                    }
                }
                
                
                std::cerr << "Graph (" << boost::get_property(*cg, boost::graph_name).type << "-" << boost::get_property(*cg, boost::graph_name).nummer << "):" << std::endl;
                std::cerr << "Graph (" << &boost::get_property(*cg, boost::graph_name).type << "-" << &boost::get_property(*cg, boost::graph_name).nummer << "):" << std::endl;
                std::cerr << "current graph storage location: " << &*cg << std::endl;
                std::cerr << "current graph property storage location: " << &boost::get_property(*cg, boost::graph_name) << std::endl;
            }
            // save final state of PM to the main graph
            
            pms[&g] = ProbabilityMatrix(current);
            if (debug) {
                std::cerr << "final PM (" << boost::get_property(g, boost::graph_name).type << "-" << boost::get_property(g, boost::graph_name).nummer << "):" << std::endl
                        << pms[&g] << std::endl;
            }
        }
        
        template <typename R>
        void DependencyGraph<R>::sample_sequence(Graph& g) {
            // reverse iterate over children
            Graph::children_iterator cg, cg_end, current;
            boost::tie(cg, cg_end) = g.children();
            for ( current = cg_end; current != cg;) {
                --current;
                
                std::cerr << "current graph storage location: " << &(*current) << std::endl;
                std::cerr << "current graph property storage location: " << &boost::get_property(*current, boost::graph_name) << std::endl;
                if (debug) {
                    std::cerr << "Sampling from: " << boost::get_property(*current, boost::graph_name).type << "-" << boost::get_property(*current, boost::graph_name).nummer << std::endl;
                    std::cerr << "Sampling from: " << &boost::get_property(*current, boost::graph_name).type << "-" << &boost::get_property(*current, boost::graph_name).nummer << std::endl;
                    std::cerr << "With PM: " << pms[&*current] << std::endl;
                }
                //print_graph(*current, &std::cerr, "iterator");
                // build a key containing the constraints of already sampled bases
                ProbabilityKey constraints;
                std::set<int> a = pms[&*current].getSpecials();
                std::cerr << "specials from this pm: " << a << std::endl;
                for (auto s : pms[&*current].getSpecials()) {
                    std::cerr << s << std::endl;
                    //std::cerr << s << "/" << int_to_vertex(s, g.root()) << "/" << enum_to_char(g.root()[int_to_vertex(s, g.root())].base) << std::endl;
                    constraints[s] = g.root()[int_to_vertex(s, g.root())].base;
                    std::cerr << "bbb"  << std::endl;
                }
                
                // randomly sample one key from the matrix
                if (debug) {
                    std::cerr << "sampling from " << boost::get_property(*current, boost::graph_name).type << "-" << boost::get_property(*current, boost::graph_name).nummer << " with key: " << std::endl
                           << constraints << std::endl << "and PM: " << std::endl
                           << pms[&*current] << std::endl;
                }
                
                ProbabilityKey colors;
                try {
                    colors = pms[&*current].sample(constraints, rand_ptr);
                } catch (std::exception& e) {
                    std::stringstream ss;
                    ss << "Error while sampling from a ProbabilityMatrix: " << std::endl << e.what();
                    std::cerr << ss.str() << std::endl;
                    throw std::logic_error(ss.str());
                }
                // write to graph
                for (auto c : colors) {
                    g.root()[int_to_vertex(c.first, g.root())].base = c.second;
                }
                // if the graph is a path, we need to color everything in between special points as well
                // else we will start the recursion
                if (boost::get_property(*current, boost::graph_name).is_path) {
                    if (debug) {
                        std::cerr << "Path Coloring!" << std::endl;
                    }
                    try {
                        color_path_graph(*current, rand_ptr);
                    } catch (std::exception& e) {
                        std::stringstream ss;
                        ss << "Error while sampling a path sequence: " << std::endl << e.what();
                        std::cerr << ss.str() << std::endl;
                        throw std::logic_error(ss.str());
                    }
                } else {
                    if (debug) {
                        std::cerr << "Recursion!" << std::endl;
                    }
                    sample_sequence(*current);
                }
            }
        }

        template <typename R>
        Sequence DependencyGraph<R>::get_sequence() {
            Sequence sequence(boost::num_vertices(graph), N);
            
            BGL_FORALL_VERTICES_T(v, graph, Graph) {
                sequence[boost::get(boost::vertex_color_t(), graph, v)] = graph[v].base;
            }
            return sequence;
        }

        template <typename R>
        std::string DependencyGraph<R>::get_sequence_string() {
            Sequence sequence = get_sequence();

            std::stringstream stream;
            stream << sequence;

            return stream.str();
        }
        
        template <typename R>
        void DependencyGraph<R>::set_sequence_string(std::string seq_str) {
            //TODO
        }
        
        template <typename R>
        void DependencyGraph<R>::set_sequence(Sequence sequence) {
            //TODO
        }

        template <typename R>
        void DependencyGraph<R>::set_sequence() {
            // reset all the colors to N
            reset_colors(graph);
            if (debug) {
                std::cerr << "Sample Sequence on Graph: Backtracing!" << std::endl;
            }
            try {
                sample_sequence(graph);
            } catch (std::exception& e) {
                std::stringstream ss;
                ss << "Error while sampling a sequence: " << std::endl << e.what();
                std::cerr << ss.str() << std::endl;
                throw std::logic_error(ss.str());
            }
        }

        template <typename R>
        unsigned long long DependencyGraph<R>::mutate_local(int min_num_pos, int max_num_pos) {
            //TODO!
            unsigned long long nos = 0;
            return nos;
        }

        template <typename R>
        unsigned long long DependencyGraph<R>::mutate_global(int min_num_pos, int max_num_pos) {
            //TODO!
            unsigned long long nos = 0;
            return nos;
        }

        template <typename R>
        unsigned long long DependencyGraph<R>::mutate(int position) {
            //TODO!
            unsigned long long nos = 0;
            return nos;
        }

        template <typename R>
        unsigned long long DependencyGraph<R>::mutate(int start, int end) {
            //TODO!
            unsigned long long nos = 0;
            return nos;
        }

        template <typename R>
        void DependencyGraph<R>::reset_colors(Graph g) {

            BGL_FORALL_VERTICES_T(v, g, Graph) {
                g[v].base = N;
            }
        }
        
        template <typename R>
        unsigned long long DependencyGraph<R>::number_of_sequences() {
            return pms[&graph].mnos();
        }
        
        template <typename R>
        unsigned long long DependencyGraph<R>::number_of_sequences(int connected_component_ID) {
            Graph::children_iterator cc, cc_end;
            // iterate over all connected component and return pm.mnos() for the one with the right ID
            for (boost::tie(cc, cc_end) = graph.children(); cc != cc_end; ++cc) {
                
                if (boost::get_property(*cc, boost::graph_name).nummer == connected_component_ID) {
                    return pms[&*cc].mnos();
                }
            }
            throw std::out_of_range("Could not find a connected component with this ID!");
        }

        template <typename R>
        std::map< int, std::vector<int> > DependencyGraph<R>::connected_components() {
            // object to return
            std::map< int, std::vector<int> > connected_components_map;
            
            Graph::children_iterator cc, cc_end;
            // iterate over all connected component and fill up the map [id, Vector of Vertices]
            for (boost::tie(cc, cc_end) = graph.children(); cc != cc_end; ++cc) {
                connected_components_map[boost::get_property(*cc, boost::graph_name).nummer] = getVertexList(*cc);
            }
            
            return connected_components_map;
        }

        template <typename R>
        std::vector< int > DependencyGraph<R>::special_vertices() {
            //TODO
        }

        template class DependencyGraph<std::mt19937>;
    }
}