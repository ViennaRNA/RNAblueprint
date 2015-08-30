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
                // get property map of graph
                graph_property& gprop = boost::get_property(*cg, boost::graph_name);
                
                if (debug) {
                    std::cerr << "current graph path? " << gprop.is_path << std::endl;
                    std::cerr << "Graph (" << gprop.type << "-" << gprop.id << "):" << std::endl;
                    print_graph(*cg, &std::cerr, "current graph");
                }
                
                if (gprop.is_path) {
                    // this is a path and therefore needs to be treated separately
                    // calculate PM for Path and save to subgraph
                    try {
                        pms[&*cg] = get_path_pm(*cg);
                    } catch (std::exception& e) {
                        std::stringstream ss;
                        ss << "Could not get a ProbabilityMatrix for a path: " << std::endl << e.what();
                        throw std::logic_error( ss.str() );
                    }
                    if (debug) {
                        std::cerr << "Path PM (" << gprop.type << "-" << gprop.id << "):" << std::endl
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
                
                // only save PM before first make_internal
                bool alreadysaved = false;
                
                // now make nodes internal and remember the current pm at such nodes 
                // (except for cc as those are independent of each other)
                if (!(gprop.type == 1)) {
                    alreadysaved = true;
                }
                
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
                                    std::cerr << "saved PM (" << gprop.type << "-" << gprop.id << "):"  << std::endl
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
            // save final state of PM to the main graph
            
            pms[&g] = current;
            if (debug) {
                graph_property& rgprop = boost::get_property(g, boost::graph_name);
                std::cerr << "final PM (" << rgprop.type << "-" << rgprop.id << "):" << std::endl
                        << pms[&g] << std::endl;
            }
        }
        
        template <typename R>
        unsigned long long DependencyGraph<R>::sample_sequence(Graph& g) {
            
            // get graph properties
            graph_property& gprop = boost::get_property(g, boost::graph_name);
            //print_graph(*current, &std::cerr, "iterator");
            // build a key containing the constraints of already sampled bases
            ProbabilityKey constraints;

            for (auto s : pms[&g].getSpecials()) {
                //std::cerr << s << "/" << int_to_vertex(s, g.root()) << "/" << enum_to_char(g.root()[int_to_vertex(s, g.root())].base) << std::endl;
                constraints[s] = g.root()[int_to_vertex(s, g.root())].base;
            }

            // randomly sample one key from the matrix
            if (debug) {
                std::cerr << "sampling from " << gprop.type << "-" << gprop.id << " with key: " << std::endl
                       << constraints << std::endl << "and PM: " << std::endl
                       << pms[&g] << std::endl;
            }

            ProbabilityKey colors;
            unsigned long long cnos = 0;
            
            try {
                std::tie(colors, cnos) = pms[&g].sample(constraints, rand_ptr);
            } catch (std::exception& e) {
                std::stringstream ss;
                ss << "Error while sampling from a ProbabilityMatrix: " << std::endl 
                        << "constrained number of sequences: " << cnos << e.what();
                throw std::logic_error(ss.str());
            }
            // write to graph
            for (auto c : colors) {
                g.root()[int_to_vertex(c.first, g.root())].base = c.second;
            }
            // if the graph is a path, we need to color everything in between special points as well
            // else we will start the recursion
            if (gprop.is_path) {
                if (debug) {
                    std::cerr << "Path Coloring!" << std::endl;
                }
                try {
                    color_path_graph(g, rand_ptr);
                } catch (std::exception& e) {
                    std::stringstream ss;
                    ss << "Error while sampling a path sequence: " << std::endl 
                        << "constrained number of sequences: " << cnos << e.what();
                    throw std::logic_error(ss.str());
                }
            } else {
                if (debug) {
                    std::cerr << "Recursion!" << std::endl;
                }
                // reverse iterate over children
                Graph::children_iterator cg, cg_end, current;
                boost::tie(cg, cg_end) = g.children();
                for ( current = cg_end; current != cg;) {
                    --current;
                    sample_sequence(*current);
                } 
            }
            return cnos;
        }

        template <typename R>
        Sequence DependencyGraph<R>::get_sequence() {
            Sequence sequence(boost::num_vertices(graph), N);
            
            BGL_FORALL_VERTICES_T(v, graph, Graph) {
                sequence[vertex_to_int(v, graph)] = graph[v].base;
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
            // get a sequence object
            Sequence sequence(seq_str.length());
            
            for (int pos = 0; pos < seq_str.length(); pos++) {
                sequence[pos] = char_to_enum(std::toupper(seq_str[pos]));
            }
            // not set this sequence to graph
            set_sequence(sequence);
        }
        
        template <typename R>
        void DependencyGraph<R>::set_sequence(Sequence sequence) {
            // remember current sequence in case of emergency
            Sequence previous = get_sequence();
            // reset all the colors to N
            reset_colors(graph);
            // write bases to graph
            for (unsigned int s = 0; s < sequence.size(); s++) {
                if (sequence[s] < A_Size) {
                    graph[int_to_vertex(s, graph)].base = sequence[s];
                } else {
                    // reset to previous sequence
                    for (unsigned int p = 0; p < previous.size(); p++) {
                        graph[int_to_vertex(p, graph)].base = previous[p];
                    }
                    std::stringstream ss;
                    ss << "Error while setting the given sequence: " << sequence << std::endl
                            << "Resetting to previous sequence: " << previous << std::endl
                            << "Only real nucleotides allowed as a fixed base" << std::endl;
                    throw std::logic_error(ss.str());
                }
            }
            try {
                sample_sequence(graph);
            } catch (std::exception& e) {
                // reset to previous sequence
                for (unsigned int p = 0; p < previous.size(); p++) {
                    graph[int_to_vertex(p, graph)].base = previous[p];
                }
                std::stringstream ss;
                ss << "Error while setting the given sequence: " << sequence << std::endl
                        << "Resetting to previous sequence: " << previous << std::endl 
                        << "Maybe constraints are not fulfilled?" << std::endl << e.what();
                throw std::logic_error(ss.str());
            }
        }

        template <typename R>
        void DependencyGraph<R>::set_sequence() {
            // remember current sequence in case of emergency
            Sequence previous = get_sequence();
            // reset all the colors to N
            reset_colors(graph);
            
            try {
                sample_sequence(graph);
            } catch (std::exception& e) {
                // reset to previous sequence
                for (unsigned int p = 0; p < previous.size(); p++) {
                    graph[int_to_vertex(p, graph)].base = previous[p];
                }
                std::stringstream ss;
                ss << "Error while setting an initial sequence!" << std::endl
                        << "Resetting to previous sequence: " << previous << e.what();
                    std::cerr << ss.str() << std::endl;
                throw std::logic_error(ss.str());
            }
        }

        template <typename R>
        unsigned long long DependencyGraph<R>::mutate_local_global(int graph_type, int min_num_pos, int max_num_pos) {
            // get all paths which fulfil the requirements of the range
            std::unordered_set< Graph* > subgraphs;
            get_subgraphs(graph, subgraphs, graph_type, min_num_pos, max_num_pos);
            if (debug) {
                for (auto s : subgraphs) {
                    std::vector<int> vertices = getVertexList(*s);
                    std::cerr << "subgraph: " << std::endl << vertices << std::endl;
                }
            }
            
            // and multiply the count it with a random number
            std::uniform_real_distribution<float> dist(0, 1);
            unsigned long long random = dist(*rand_ptr) * subgraphs.size();
            
            unsigned long long sum = 0;
            for (auto s : subgraphs) {
                sum ++;
                // if the random number is bigger than our probability, take this base as the current base!
                if (random < sum) {
                    return mutate(*s);
                    break;
                }
            }
            throw std::out_of_range("Could not find any component which fulfils your requirements!");
        }
        
        template <typename R>
        unsigned long long DependencyGraph<R>::mutate_global(int connected_component_ID) {
            Graph::children_iterator cc, cc_end;
            for (boost::tie(cc, cc_end) = graph.children(); cc != cc_end; ++cc) {
                if (boost::get_property(*cc, boost::graph_name).id == connected_component_ID) {
                    return mutate(*cc);
                }
            }
            throw std::out_of_range("Could not find a connected component with this ID!");
        }
        
        template <typename R>
        void DependencyGraph<R>::get_subgraphs(Graph& g, std::unordered_set< Graph* >& subgraphs, int type, int min_size, int max_size) {
            // if max is 0, set it to infinite, as defined in the documentation
            if (max_size == 0) {
                max_size = std::numeric_limits<int>::max();
            }
            
            if (max_size < min_size) {
                // assume that this is a mistake
                int temp = max_size;
                max_size = min_size;
                min_size = temp;
            }
            
            // if this subgraph is from the given type, insert it into the set
            if (((type == -1 && boost::get_property(g, boost::graph_name).is_path) || (boost::get_property(g, boost::graph_name).type == type))
                    && (min_size <= boost::num_vertices(g) <= max_size)) {
                subgraphs.emplace(&g);
            }
            // and now check all children, too
            Graph::children_iterator c, c_end;
            for (boost::tie(c, c_end) = g.children(); c != c_end; ++c) {
                get_subgraphs(*c, subgraphs, type, min_size, max_size);
            }
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
        unsigned long long DependencyGraph<R>::mutate(Graph& g) {
            // get graph properties
            graph_property& gprop = boost::get_property(g, boost::graph_name);
            
            // if it is a connected component, do a sample_sequence (also if it is a path too)
            if (gprop.type == 1) {
                if (debug) {
                    std::cerr << "Mutating a connected component!" << std::endl;
                }
                // reset the whole connected component and sample everything
                try {
                    reset_colors(g);
                    return sample_sequence(g);
                } catch (std::exception& e) {
                    std::stringstream ss;
                    ss << "Error while mutating a connected component: " << std::endl << e.what();
                    throw std::logic_error(ss.str());
                }
            // in case this is a path, do a path sampling!
            } else if (gprop.is_path) {
                if (debug) {
                    std::cerr << "Mutating a path!" << std::endl;
                }
                // reset all vertices, except special ones, as those are the important ends
                // which have to stay the same
                BGL_FORALL_VERTICES_T(v, g, Graph) {
                    if (!g[v].special) {
                        g[v].base = N;
                    }
                }
                try {
                    return sample_sequence(g);
                } catch (std::exception& e) {
                    std::stringstream ss;
                    ss << "Error while mutating a path sequence: " << std::endl << e.what();
                    throw std::logic_error(ss.str());
                }
            } else {
                std::stringstream ss;
                ss << "I thing it is not allowed to mutate only this subgraph: "
                        << gprop.type << "-" << gprop.id << std::endl;
                throw std::logic_error(ss.str());
            }
        }

        template <typename R>
        void DependencyGraph<R>::reset_colors(Graph& g) {

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
            // iterate over all connected component and return pm.mnos() for the one with the right ID
            Graph::children_iterator cc, cc_end;
            for (boost::tie(cc, cc_end) = graph.children(); cc != cc_end; ++cc) {
                if (boost::get_property(*cc, boost::graph_name).id == connected_component_ID) {
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
                connected_components_map[boost::get_property(*cc, boost::graph_name).id] = getVertexList(*cc);
            }
            
            return connected_components_map;
        }

        template <typename R>
        std::vector< int > DependencyGraph<R>::special_vertices() {
            // get special vertices as vector
            std::vector< int > result;
            BGL_FORALL_VERTICES_T(v, graph, Graph) {
                if (graph[v].special && (graph[v].constraint == N)) {
                    result.push_back(vertex_to_int(v, graph));
                }
            }
            return result;
        }

        template class DependencyGraph<std::mt19937>;
    }
}