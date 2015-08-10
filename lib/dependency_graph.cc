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
            graph = parse_structures(structures);
            
            // set sequence constraints
            set_constraints(graph, constraints);

            // decompose the graph into its connected components, biconnected
            // components and decompose blocks via ear decomposition
            bipartite = decompose_graph(graph, rand_ptr);
            
            // now calculate all the PMs
            calculate_probabilities(graph);
            // remember nos
            nos = boost::get_property(graph, boost::graph_name).pm->mnos();

        }

        template <typename R>
        void DependencyGraph<R>::calculate_probabilities(Graph& g) {
            
            BGL_FORALL_VERTICES_T(v, g, Graph) {
                // reset the color tag for the calculate_probabilities function
                g[v].color = 0;
            }
            
            // Remember a temporary PM which holds the current state
            ProbabilityMatrix current;
            
            Graph::children_iterator cg, cg_end;
            for (boost::tie(cg, cg_end) = g.children(); cg != cg_end; ++cg) {
                if (debug) {
                    std::cerr << "current graph path? " << boost::get_property(*cg, boost::graph_name).is_path << std::endl;
                    print_graph(*cg, &std::cerr, "current graph");
                    //std::cerr << "current PM: " << current << std::endl;
                }
                
                if (boost::get_property(*cg, boost::graph_name).is_path) {
                    // this is a path and therefore needs to be treated separately
                    // calculate PM for Path and save to subgraph                    
                    ProbabilityMatrix path_matrix = get_path_pm(*cg);
                    boost::get_property(*cg, boost::graph_name).pm = std::make_shared<ProbabilityMatrix> (path_matrix);
                    if (debug) {
                        std::cerr << "Path PM: " << std::endl << *boost::get_property(*cg, boost::graph_name).pm << std::endl;
                    }
                } else {
                    if (debug) {
                        std::cerr << "Recursion!" << std::endl;
                    }
                    // recursion is here
                    calculate_probabilities(*cg);
                }
                
                // Multiply current with pm of this child
                current = *boost::get_property(*cg, boost::graph_name).pm * current;
                if (debug) {
                    std::cerr << "current PM: " << std::endl << current << std::endl;
                }
                bool pmsaved = false;
                BGL_FORALL_VERTICES_T(v, *cg, Graph) {
                    if ((*cg)[v].special) {
                        // update current degrees as status of special points
                        (*cg)[v].color += boost::degree(v, *cg);

                        // check if a vertex becomes internal here
                        if (debug) {
                            std::cerr << "Internal? " << (*cg)[v].color << "/" << boost::degree((*cg).local_to_global(v), (*cg).root()) << std::endl;
                        }
                        if ((*cg)[v].color == boost::degree((*cg).local_to_global(v), (*cg).root())) {
                            // remember this PM here
                            // only the first time a internal node is detected, otherwise we overwrite this
                            if (!pmsaved) {
                                pmsaved = true;
                                if (!boost::get_property(*cg, boost::graph_name).is_cc) {
                                    boost::get_property(*cg, boost::graph_name).pm = std::make_shared<ProbabilityMatrix> (current);
                                    if (debug) {
                                        std::cerr << "saved PM: " << std::endl << *boost::get_property(*cg, boost::graph_name).pm << std::endl;
                                    }
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
            boost::get_property(g, boost::graph_name).pm = std::make_shared<ProbabilityMatrix> (current);
            if (debug) {
                std::cerr << "PM of " << boost::get_property(g, boost::graph_name).id << std::endl << *boost::get_property(g, boost::graph_name).pm << std::endl;
            }
        }
        
        template <typename R>
        void DependencyGraph<R>::sample_sequence(Graph& g) {
            // reverse iterate over children
            Graph::children_iterator cg, cg_end, current;
            boost::tie(cg, cg_end) = g.children();
            for ( current = cg_end; current != cg;) {
                --current;
                print_graph(*current, &std::cerr, "iterator");
                
                // recursion is here (abort is if graph is a path!)
                if (!boost::get_property(*current, boost::graph_name).is_path) {
                    if (debug) {
                        std::cerr << "Recursion!" << std::endl;
                    }
                    sample_sequence(*current);
                }
                
                // build a key containing the constraints of already sampled bases
                ProbabilityKey constraints;
                for (auto s : boost::get_property(*current, boost::graph_name).pm->getSpecials()) {
                    constraints[s] = (*current)[int_to_vertex(s, *current)].base;
                }
                
                // randomly sample one key from the matrix
                ProbabilityKey colors = boost::get_property(*current, boost::graph_name).pm->sample(constraints, rand_ptr);
                // write to graph
                for (auto c : colors) {
                    (*current)[int_to_vertex(c.first, *current)].base = c.second;
                }
                // if the graph is a path, we need to color everything in between special points as well
                if (boost::get_property(*current, boost::graph_name).is_path) {
                    if (debug) {
                        std::cerr << "Path Coloring!" << std::endl;
                    }
                    color_path_graph(*current, rand_ptr);
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
        void DependencyGraph<R>::mutate() {
            // reset all the colors to N
            reset_colors();
            // TODO replace with a good mutation function
            if (debug) {
                std::cerr << "Sample Sequence on Graph: Backtracing!" << std::endl;
            }
            sample_sequence(graph);
        }

        template <typename R>
        void DependencyGraph<R>::mutate(int position) {
            // TODO replace with a good mutation function	
            this->mutate();
        }

        template <typename R>
        void DependencyGraph<R>::reset_colors() {

            BGL_FORALL_VERTICES_T(v, graph, Graph) {
                graph[v].base = N;
            }
        }

        template <typename R>
        DependencyGraph<R>::~DependencyGraph() {
            // TODO clean up all the data
        }

        template class DependencyGraph<std::mt19937>;
    }
}