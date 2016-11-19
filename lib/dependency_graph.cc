/* RNAblueprint
 * A program for designing RNA molecules.
 *
 * @date 18.03.2014
 * @author Stefan Hammer <s.hammer@univie.ac.at>
 * @copyright GPLv3
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

    namespace detail {
        template <typename R>
        DependencyGraph<R>::DependencyGraph(std::vector<std::string> structures, std::string constraints, R _rand)
        : rand(_rand), history_size(100) {

            if (debug) {
                std::cerr << "Initializing DependencyGraph..." << std::endl;
            }
            // start to measure time
            std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
            
            // generate graph from input vector
            try {
                graph = parse_structures(structures);
            } catch (std::exception& e) {
                std::stringstream ss;
                ss << "Error while parsing the structures: " << std::endl << e.what();
                throw std::logic_error(ss.str());
            }
            
            // set sequence constraints
            try {
                set_constraints(graph, constraints);
            } catch (std::exception& e) {
                throw std::logic_error(e.what());
            }

            // decompose the graph into its connected components, biconnected
            // components and decompose blocks via ear decomposition
            try {
                bipartite = decompose_graph(graph, rand);
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
                calculate_probabilities(graph, start_time);
            } catch (std::exception& e) {
                std::stringstream ss;
                ss << "Error while calculating the probabilities: " << std::endl << e.what();
                throw std::logic_error(ss.str());
            }
            // get a first initial random sequence
            sample();
        }

        template <typename R>
        void DependencyGraph<R>::calculate_probabilities(Graph& g, std::chrono::steady_clock::time_point& start_time) {
            
            // Remember a temporary PM which holds the current state
            ProbabilityMatrix current;
            std::unordered_map<Vertex, unsigned int> degree_map;
            
            Graph::children_iterator cg, cg_end;
            for (boost::tie(cg, cg_end) = g.children(); cg != cg_end; ++cg) {
                // check if timeout is already reached
                check_timeout(start_time);
                
                // get property map of graph
                graph_property& gprop = boost::get_property(*cg, boost::graph_name);
                
                if (debug) {
                    std::cerr << "current graph path? " << gprop.is_path << std::endl;
                    std::cerr << "Graph (" << gprop.type << "-" << gprop.id << "):" << std::endl;
                    print_graph(*cg, &std::cerr);
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
                    // save number of solutions as graph property
                    gprop.nos = pms[&*cg].mnos();
                    
                    if (debug) {
                        std::cerr << "Path PM (" << gprop.type << "-" << gprop.id << ") with nos " << gprop.nos << ":" << std::endl
                                << pms[&*cg] << std::endl;
                    }
                } else {
                    if (debug) {
                        std::cerr << "Recursion!" << std::endl;
                    }
                    // recursion is here
                    calculate_probabilities(*cg, start_time);
                }
                
                // Multiply current with pm of this child
                // Multiplication is not symmetric in terms of performance, therefore
                // we want to have the matrix with more special vertices at the front.
                if (current.getSpecials().size() > pms[&*cg].getSpecials().size())
                    current = current * pms[&*cg];
                else
                    current = pms[&*cg] * current;
                
                // update max_dimensions
                if (current.getDimensions() > max_dimensions) {
                    max_dimensions = current.getDimensions();
                }
                
                if (debug) {
                    std::cerr << "current PM: " << std::endl << current << std::endl;
                }
                
                // check if timeout is already reached
                check_timeout(start_time);
                
                // only save PM before first make_internal
                bool alreadysaved = false;
                
                // now make nodes internal and remember the current pm at such nodes 
                // (except for cc as those are independent of each other)
                if (gprop.type == 1) {
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
            graph_property& rgprop = boost::get_property(g, boost::graph_name);
            // remember maximum number of solutions for this subgraph
            rgprop.nos = current.mnos();
            
            if (debug) {
                graph_property& rgprop = boost::get_property(g, boost::graph_name);
                std::cerr << "final PM (" << rgprop.type << "-" << rgprop.id << ") with nos " << rgprop.nos << ":" << std::endl
                        << pms[&g] << std::endl;
            }
        }
        
        template <typename R>
        ProbabilityFraction DependencyGraph<R>::sample_sequence(Graph& g) {
            
            // get graph properties
            graph_property& gprop = boost::get_property(g, boost::graph_name);
            // build a key containing the constraints of already sampled bases
            ProbabilityKey constraints;
            bool onlypath = true;

            for (auto s : pms[&g].getSpecials()) {
                //std::cerr << s << "/" << int_to_vertex(s, g.root()) << "/" << enum_to_char(g.root()[int_to_vertex(s, g.root())].base) << std::endl;
                constraints[s] = g.root()[int_to_vertex(s, g.root())].base;
                if (constraints[s] >= A_Size)
                    onlypath = false;
            }

            // randomly sample one key from the matrix
            if (debug) {
                std::cerr << "sampling from " << gprop.type << "-" << gprop.id << " with key: " << std::endl
                       << constraints << std::endl << "and PM: " << std::endl
                       << pms[&g] << std::endl;
            }
            
            ProbabilityFraction pf = std::make_pair(1,0);
            ProbabilityKey colors;
            
            try {
                std::tie(colors, pf) = pms[&g].sample(constraints, rand);
            } catch (std::exception& e) {
                std::stringstream ss;
                ss << "Error while sampling from a ProbabilityMatrix: " << std::endl 
                        << "constrained number of sequences: " << pf.second << e.what();
                throw std::logic_error(ss.str());
            }
            // write to graph
            for (auto c : colors) {
                g.root()[int_to_vertex(c.first, g.root())].base = c.second;
            }
            // if the graph is a path, we need to colour everything in between special points as well
            // else we will start the recursion
            if (gprop.is_path) {
                if (debug) {
                    std::cerr << "Path Coloring!" << std::endl;
                }
                try {
                    ProbabilityFraction  path_pf = color_path_graph(g, rand);
                    if (onlypath) {
                        pf = path_pf;
                    }
                } catch (std::exception& e) {
                    std::stringstream ss;
                    ss << "Error while sampling a path sequence: " << std::endl 
                        << "constrained number of sequences: " << pf.second << e.what();
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
                    ProbabilityFraction graph_pf = sample_sequence(*current);
                    pf.first *= graph_pf.first/graph_pf.second;
                } 
            }
            
            return pf;
        }
        
        template <typename R>
        std::string DependencyGraph<R>::get_graphml() {
            std::ostringstream stream;
            print_graph(graph, dynamic_cast<std::ostream*>(&stream));
            return stream.str();
        }
        
        template <typename R>
        std::string DependencyGraph<R>::get_graphml(int connected_component_ID) {
            Graph::children_iterator cc, cc_end;
            for (boost::tie(cc, cc_end) = graph.children(); cc != cc_end; ++cc) {
                if (boost::get_property(*cc, boost::graph_name).id == connected_component_ID) {
                    std::ostringstream stream;
                    print_graph(*cc, dynamic_cast<std::ostream*>(&stream));
                    return stream.str();
                }
            }
            throw std::out_of_range("Could not find a connected component with this ID!");
        }
        
        template <typename R>
        unsigned long DependencyGraph<R>::set_seed() {
            unsigned long seed = std::chrono::system_clock::now().time_since_epoch().count();
            if (debug) {
                std::cerr << "Using this seed: " << seed << std::endl;
            }
            rand.seed(seed);
            return seed;
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
            return sequence_to_string(get_sequence());
        }
        
        template <typename R>
        std::string DependencyGraph<R>::sequence_to_string(Sequence sequence) {           
            std::stringstream stream;
            stream << sequence;

            std::string result_str = stream.str();
            
            // insert cutpoints
            for (auto& c : boost::get_property(graph, boost::graph_name).cutpoints) {
                result_str.insert(result_str.begin()+c.first, c.second);
            }
            
            return result_str;
        }
        
        template <typename R>
        SolutionSizeType DependencyGraph<R>::set_sequence_string(std::string seq_str) {
            // remove cut points
            std::size_t found_cut;
            std::map<int, char> cutpoints = boost::get_property(graph, boost::graph_name).cutpoints;
            while (true) {
                found_cut = seq_str.find_last_of("&+");
                if (found_cut == std::string::npos)
                    break;
                std::map<int, char>::const_iterator found = cutpoints.find(found_cut);
                if (found == cutpoints.end())
                    throw std::logic_error("Cut points of the new sequence are not aligned properly!");
                seq_str.erase(found_cut, 1);
            }
            
            // get a sequence object
            Sequence sequence(seq_str.length());
            
            for (unsigned int pos = 0; pos < seq_str.length(); pos++) {
                sequence[pos] = char_to_enum(std::toupper(seq_str[pos]));
            }
            // not set this sequence to graph
            return set_sequence(sequence);
        }
        
        template <typename R>
        SolutionSizeType DependencyGraph<R>::set_sequence(Sequence sequence) {
            ProbabilityFraction pf;
            // reset all the colors to N
            reset_colors(graph);
            // write bases to graph
            for (unsigned int s = 0; s < sequence.size(); s++) {
                if (sequence[s] < A_Size) {
                    graph[int_to_vertex(s, graph)].base = sequence[s];
                } else {
                    revert_sequence(0);
                    std::stringstream ss;
                    ss << "Error while setting the given sequence: " << sequence << std::endl
                            << "Resetting to previous sequence: " << get_sequence_string() << std::endl
                            << "Only real nucleotides allowed as a fixed base" << std::endl;
                    throw std::logic_error(ss.str());
                }
            }
            try {
                pf = sample_sequence(graph);
            } catch (std::exception& e) {
                // reset to previous sequence
                revert_sequence(0);
                std::stringstream ss;
                ss << "Error while setting the given sequence: " << sequence << std::endl
                        << "Resetting to previous sequence: " << get_sequence_string() << std::endl 
                        << "Maybe constraints are not fulfilled?" << std::endl << e.what();
                throw std::logic_error(ss.str());
            }
            // remember this new sequence in the history
            remember_sequence();
            return pf.second;
        }

        template <typename R>
        SolutionSizeType DependencyGraph<R>::sample() {
            ProbabilityFraction pf;
            // reset all the colors to N
            reset_colors(graph);
            try {
                 pf = sample_sequence(graph);
            } catch (std::exception& e) {
                revert_sequence(0);
                std::stringstream ss;
                ss << "Error while sampling an initial sequence!" << std::endl
                    << "Resetting to previous sequence: " << get_sequence_string() << e.what();
                throw std::logic_error(ss.str());
            }
            // remember this new sequence in the history
            remember_sequence();
            return pf.second;
        }

        template <typename R>
        SolutionSizeType DependencyGraph<R>::sample_local_global(int graph_type, int min_num_pos, int max_num_pos) {
            // get all paths which fulfill the requirements of the range
            std::unordered_set< Graph* > subgraphs;
            get_subgraphs(graph, subgraphs, graph_type, min_num_pos, max_num_pos);
            // calculate maximum number of solutions for all subgraphs
            SolutionSizeType mnos = 0;
            for (auto s : subgraphs) {
                mnos += boost::get_property(*s, boost::graph_name).nos;
                if (debug) {
                    std::vector<int> vertices = getVertexList(*s);
                    std::cerr << "subgraph: " << std::endl << vertices << std::endl;
                }
            }
            
            // and multiply the count it with a random number
            RandomDistType dist(0, mnos);
            SolutionSizeType random = dist(rand);
            
            SolutionSizeType sum = 0;
            for (auto s : subgraphs) {
                sum += boost::get_property(*s, boost::graph_name).nos;
                // if the random number is bigger than our probability, take this base as the current base!
                if (random < sum) {
                    SolutionSizeType cnos = sample(*s);
                    remember_sequence();
                    return cnos;
                    break;
                }
            }
            throw std::out_of_range("Could not find any component which fulfills your requirements!");
        }
        
        template <typename R>
        SolutionSizeType DependencyGraph<R>::sample_global(int connected_component_ID) {
            Graph::children_iterator cc, cc_end;
            for (boost::tie(cc, cc_end) = graph.children(); cc != cc_end; ++cc) {
                if (boost::get_property(*cc, boost::graph_name).id == connected_component_ID) {
                    SolutionSizeType cnos = sample(*cc);
                    remember_sequence();
                    return cnos;
                }
            }
            throw std::out_of_range("Could not find a connected component with this ID!");
        }
        
        template <typename R>
        void DependencyGraph<R>::get_subgraphs(Graph& g, std::unordered_set< Graph* >& subgraphs, int type, unsigned int min_size, unsigned int max_size) {
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
            // check the graph type
            // check the graph size
            // ignore subgraphs with only one possibility
            if (((type == -1 && boost::get_property(g, boost::graph_name).is_path) || (boost::get_property(g, boost::graph_name).type == type))
                    && ((min_size <= boost::num_vertices(g)) && (boost::num_vertices(g) <= max_size))
                    && (boost::get_property(g, boost::graph_name).nos != 1)) {
                subgraphs.emplace(&g);
            }
            // and now check all children, too
            Graph::children_iterator c, c_end;
            for (boost::tie(c, c_end) = g.children(); c != c_end; ++c) {
                get_subgraphs(*c, subgraphs, type, min_size, max_size);
            }
        }
        
        template <typename R>
        Graph* DependencyGraph<R>::find_path_subgraph(Vertex v_global, Graph& g) {
            Graph::children_iterator c, c_end;
            for (boost::tie(c, c_end) = g.children(); c != c_end; ++c) {
                // search if this vertex is present
                if ((*c).find_vertex(v_global).second) {
                    Vertex v = (*c).find_vertex(v_global).first;
                    // get graph properties
                    graph_property& gprop = boost::get_property(*c, boost::graph_name);
                    // return pointer to this path or connected component in case of a special position, or go deeper
                    if (gprop.is_path || (gprop.type == 1 && g[v].special)) {
                        if (debug)
                            print_graph(*c, &std::cerr);
                        return &*c;
                    } else {
                        return find_path_subgraph(v_global, *c);
                    }
                    break;
                }
            }
            // this should never happen
            return &g;
        }

        template <typename R>
        SolutionSizeType DependencyGraph<R>::sample(int position) {
            // first get the right vertex
            Vertex v_global = int_to_vertex(position, graph);
            if (debug)
                std::cerr << "vertex is: " << v_global << std::endl;
            // search for the lowest path subgraph and sample
            Graph* g = find_path_subgraph(v_global, graph);
            SolutionSizeType cnos = sample(*g);
            remember_sequence();
            return cnos;
        }

        template <typename R>
        SolutionSizeType DependencyGraph<R>::sample(int start, int end) {
            // nos
            SolutionSizeType nos = 1;
            // set with all subgraphs to sample
            std::set<Graph*> subgraphs;
            
            for (; start <= end; start++) {
                Vertex v_global = int_to_vertex(start, graph);
                if (debug)
                    std::cerr << "vertex is: " << v_global << std::endl;
                subgraphs.insert(find_path_subgraph(v_global, graph));
            }
            // sample collected subgraphs
            for (auto sg : subgraphs) {
                nos *= sample(*sg);
            }
            remember_sequence();
            return nos;
        }
        
        template <typename R>
        SolutionSizeType DependencyGraph<R>::sample(Graph& g) {
            // get graph properties
            graph_property& gprop = boost::get_property(g, boost::graph_name);
            
            // reset the whole subgraph for CCs
            if (gprop.type == 1) {
                if (debug) {
                    std::cerr << "Sampling a connected component!" << std::endl;
                }
                // reset the whole connected component and sample everything
                reset_colors(g);
                try {
                    return sample_sequence(g).second;
                } catch (std::exception& e) {
                    std::stringstream ss;
                    ss << "Error while sampling a connected component (" << gprop.type << "-" << gprop.id << "): " << std::endl << e.what();
                    throw std::logic_error(ss.str());
                }
                // in case this is a path, only reset without specials
            } else if (gprop.is_path) {
                if (debug) {
                    std::cerr << "Sampling a path!" << std::endl;
                }
                // reset all vertices, except special ones, as those are the important ends
                // which have to stay the same
                BGL_FORALL_VERTICES_T(v, g, Graph) {
                    if (!g[v].special) {
                        g[v].base = N;
                    }
                }
                try {
                    return sample_sequence(g).second;
                } catch (std::exception& e) {
                    std::stringstream ss;
                    ss << "Error while sampling a path (" << gprop.type << "-" << gprop.id << "): " << std::endl << e.what();
                    throw std::logic_error(ss.str());
                }
            } else {
                // This should never be reached
                std::stringstream ss;
                ss << "I think it is not allowed to sample only this subgraph: "
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
        SolutionSizeType DependencyGraph<R>::number_of_sequences() {
            //return pms[&graph].mnos();
            return boost::get_property(graph, boost::graph_name).nos;
        }
        
        template <typename R>
        SolutionSizeType DependencyGraph<R>::number_of_sequences(int connected_component_ID) {
            // iterate over all connected component and return pm.mnos() or gprop.nos for the one with the right ID
            Graph::children_iterator cc, cc_end;
            for (boost::tie(cc, cc_end) = graph.children(); cc != cc_end; ++cc) {
                graph_property& gprop = boost::get_property(*cc, boost::graph_name);
                if (gprop.id == connected_component_ID) {
                    //return pms[&*cc].mnos();
                    return gprop.nos;
                }
            }
            throw std::out_of_range("Could not find a connected component with this ID!");
        }

        template <typename R>
        int DependencyGraph<R>::number_of_connected_components() {
            // object to return
            int count = 0;
            Graph::children_iterator cc, cc_end;
            // iterate over all connected component and fill up the map [id, Vector of Vertices]
            for (boost::tie(cc, cc_end) = graph.children(); cc != cc_end; ++cc) {
                count++;
            }
            return count;
        }
        
        template <typename R>
        std::vector<int> DependencyGraph<R>::component_vertices(int connected_component_ID) {
            // iterate over all connected component and fill up the vector
            Graph::children_iterator cc, cc_end;
            for (boost::tie(cc, cc_end) = graph.children(); cc != cc_end; ++cc) {
                if (boost::get_property(*cc, boost::graph_name).id == connected_component_ID) {
                    return getVertexList(*cc);
                }
            }
            throw std::out_of_range("Could not find a connected component with this ID!");
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
        
        template <typename R>
        std::vector< int > DependencyGraph<R>::special_vertices(int connected_component_ID) {
            // get special vertices as vector
            std::vector< int > result;
            // iterate over all connected component and fill up the vector
            Graph::children_iterator cc, cc_end;
            for (boost::tie(cc, cc_end) = graph.children(); cc != cc_end; ++cc) {
                if (boost::get_property(*cc, boost::graph_name).id == connected_component_ID) {
                    BGL_FORALL_VERTICES_T(v, *cc, Graph) {
                        if (graph[v].special && (graph[v].constraint == N)) {
                            result.push_back(vertex_to_int(v, graph));
                        }
                    }
                    return result;
                }
            }
            throw std::out_of_range("Could not find a connected component with this ID!");
        }
        
        template <typename R>
        bool DependencyGraph<R>::revert_sequence(unsigned int jump) {
            // check if we already reached the beginning or do a boundary jump
            if (debug) {
                std::cerr << "Going back in time some steps: " << jump << std::endl;
                std::cerr << "History size: " << history.size() << "/" << history_size << std::endl;
            }
            if (jump < history.size()) {
                if (debug) {
                    std::cerr << "Lets do the time warp again!" << std::endl;
                }
                // go to the last position
                std::list<Sequence>::iterator current = std::prev(history.end());
                // jump back in time
                try {
                    // std::advance can set the iterator forward, but also reverse if jump is negative
                    std::advance(current, ((int) jump * -1));
                } catch (std::exception& e) {
                    std::stringstream ss;
                    ss << "Error while reverting the sequence: " << std::endl << e.what();
                    throw std::logic_error(ss.str());
                }
                // set the sequence on the graph
                for (unsigned int s = 0; s < current->size(); s++) {
                    graph[int_to_vertex(s, graph)].base = (*current)[s];
                }
                // erase everything behind the new current
                history.erase(++current, history.end());
                return true;
            } else {
                if (debug) {
                    std::cerr << "We already arrived at big bang!" << std::endl;
                }
                return false;
            }
        }
        
        template <typename R>
        void DependencyGraph<R>::remember_sequence() {
            // push current sequence to the end
            history.push_back(get_sequence());
            // if our stack is bigger than the maximum, clear the oldest one
            if (history.size() > history_size)
                history.erase(history.begin());
        }
        
        template <typename R>
        void DependencyGraph<R>::set_history_size(unsigned int size) {
            history_size = size;
            if (history.size() > history_size) {
                // get an iterator pointing to the first element that should be maintained
                std::list<Sequence>::iterator end = history.begin();
                std::advance(end, history.size()-history_size);
                // erase everything old. from beginning until only history_size elements remain
                history.erase(history.begin(), end);
            }
        }
        
        template <typename R>
        std::vector< std::string > DependencyGraph<R>::get_history() {
            std::vector< std::string > result;
            // convert history stack to string
            for (auto& h : history) {
                result.push_back(sequence_to_string(h));
            }
            return result;
        }
        
        template class DependencyGraph<std::mt19937>;
    }
}