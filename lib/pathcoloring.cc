/* RNAdesign
 * A program for designing RNA molecules.
 *
 * Created on: 13.08.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

// include header
#include "pathcoloring.h"

namespace design {
    namespace detail {
        
        // calculate a ProbabilityMatrix for a path
        ProbabilityMatrix get_path_pm(Graph& g) {
            // check if given graph is indeed a path with max_degree = 2 and two ends with degree = 1;
            int max_degree;
            int min_degree;
            std::tie(min_degree, max_degree) = get_min_max_degree(g);
            // assert path
            if (max_degree > 2) {
                throw std::logic_error("This graph is no cycle or path (max degree > 2). I can't color this!");
            } else if (min_degree > 1) {
                throw std::logic_error("cannot color circles this way.");
            }
            
            ProbabilityKey key;
            std::set< int > specials;
            
            BGL_FORALL_VERTICES_T(v, g, Graph) {
                if (g[v].special) {
                    if (boost::out_degree(v, g) <= 1) {
                        // remember vertex as int with sequence constraint
                        key[vertex_to_int(v, g)] = g[v].constraint;
                        specials.emplace(vertex_to_int(v, g));
                    } else {
                        throw std::logic_error("There is a special vertex which is no path end in get_path_pm. This is not possible!");
                    }
                }
            }
            
            if (specials.size() > 2) {
                throw std::logic_error("More than two special vertices in one path. ridiculous!");
            }
            
            std::vector<ProbabilityKey> keys = permute_key(key);
            //std::cerr << "keys: " << std::endl << keys;
            int length = boost::num_vertices(g) - 1;
            //std::cerr << "length: " << length << std::endl;
            PairingMatrix * p = PairingMatrix::Instance();
            
            ProbabilityMatrix result;
            
            for (auto k : keys) {
                switch ( specials.size() )
                {
                    case 0:
                        result.put(k, p->get(length, N, N));
                        break;
                    case 1:
                        result.put(k, p->get(length, k[*(specials.begin())], N));
                        break;
                    case 2:
                        result.put(k, p->get(length, k[*(specials.begin())], k[*(++specials.begin())]));
                        break;
                    default:
                        throw std::logic_error("More than two special vertices in one path. ridiculous!");
                }
            }
            return result;
        }
        
        template <typename RG>
        unsigned long long color_path_graph(Graph& g, RG* rand_ptr) {

            unsigned long long max_number_of_sequences = 0;

            // check if given graph is indeed a path with max_degree = 2 and two ends with degree = 1;
            int max_degree;
            int min_degree;
            std::tie(min_degree, max_degree) = get_min_max_degree(g);
            
            std::cerr << "min/max: " << min_degree << "/" << max_degree << std::endl;
            // assert path
            if (max_degree > 2) {
                throw std::logic_error("This graph is no cycle or path (max degree > 2). I can't color this!");
            } else if (min_degree > 1) {
                throw std::logic_error("cannot color circles this way.");
            }
            
            /*if (debug) {
                print_graph(g, &std::cerr, "this path will be colored");
            }*/

            // visitor declaration

            class color_dfs_visitor : public boost::default_dfs_visitor {
            public:

                color_dfs_visitor(unsigned long long& max_number_of_sequences, RG * rand_ptr, PairingMatrix * pair, std::uniform_real_distribution<float>& d,
                        nosMap& n, std::unordered_map<Vertex, int>& c, int& prev)
                : mnos(max_number_of_sequences), r_ptr(rand_ptr), p(pair), dist(d), nos_map(n), colors(c), previous(prev) {
                }
                unsigned long long& mnos;
                RG * r_ptr;
                PairingMatrix * p;
                std::uniform_real_distribution<float>& dist;
                nosMap& nos_map;
                std::unordered_map<Vertex, int>& colors;
                int& previous;

                void start_vertex(Vertex s, Graph g) const {
                    
                    if (debug) {
                        std::cerr << "Start vertex: " << s << " [" << enum_to_char(g[s].base) << "]" << std::endl;
                    }

                    mnos = 0;
                    for (auto b : base_conversion[ g[s].base ]) {
                        nos_map[s][b] = p->get(0, b, b);
                        if (debug) {
                            std::cerr << "v" << s << ": " << enum_to_char(b) << ": " << nos_map[s][b] << std::endl;
                        }
                        mnos += nos_map[s][b];
                    }
                }

                void tree_edge(Edge e, Graph g) const {
                    if (debug) {
                        std::cerr << "Tree edge: " << e << std::endl;
                    }
                    Vertex v = boost::source(e, g);
                    Vertex u = boost::target(e, g);

                    // we need to calculate
                    // the number of possibilities for each base on this node.
                    // therefore we do all the combinations between the last !N node
                    // and all the combinations on this one and remember them.
                    mnos = 0;
                    for (auto u_base : base_conversion[ g[u].base ]) {
                        nos_map[u][u_base] = 0;
                        for (auto v_base : base_conversion[ g[v].base ]) {
                            nos_map[u][u_base] += nos_map[v][v_base] * p->get(1, v_base, u_base);
                        }
                        // calculate maximal number of sequences on this vertex
                        mnos += nos_map[u][u_base];
                        if (debug) {
                            std::cerr << "v" << u << ": " << enum_to_char(u_base) << ": " << nos_map[u][u_base] << std::endl;
                        }
                    }
                }

                void back_edge(Edge e, Graph g) const {
                    if (debug) {
                        throw std::logic_error("Detecting back-edge (graph is a cycle). This can't be!");
                    }
                }

                void finish_vertex(Vertex u, Graph g) const {
                    if (debug) {
                        std::cerr << "Finishing vertex: " << u << std::endl;
                    }

                    // calculate number of sequences with respect to the chosen previous base
                    unsigned long long nos = 0;
                    for (auto b : base_conversion[ g[u].base ]) {
                        if (p->get(1, b, previous) > 0) {
                            nos += nos_map[u][b];
                        }
                    }

                    if (nos == 0) {
                        std::stringstream ss;
                        ss << "The requested sequence cannot be colored! Conflict at: "
                                << enum_to_char(g[u].base) << ", " << enum_to_char(previous) << std::endl 
                                << "Length is: " << boost::num_vertices(g) << std::endl;
                        throw std::logic_error(ss.str());
                    }

                    unsigned long long random = dist(*r_ptr) * nos;

                    // stochastically take one of the possibilities
                    // start at the probability of first possible character and add each other base probability as long long as the random number is bigger.
                    unsigned long long sum = 0;
                    for (auto b : base_conversion[ g[u].base ]) {
                        if (p->get(1, b, previous) > 0) {
                            sum += nos_map[u][b];
                            // if the random number is bigger than our probability, take this base as the current base!
                            if (random < sum) {
                                colors[u] = b;
                                previous = b;
                                // don't forget to exit the loop, otherwise will always be first = C;
                                break;
                            }
                        }
                    }

                    if (debug) {
                        std::cerr << "Vertex colored: " << vertex_to_int(u, g) << "/" << enum_to_char(colors[u]) << std::endl;
                    }
                }
            };

            PairingMatrix * p = PairingMatrix::Instance();
            std::uniform_real_distribution<float> dist(0, 1);
            nosMap nos_map;
            std::unordered_map<Vertex, int> colors;
            int prev = N;
            
            color_dfs_visitor vis(max_number_of_sequences, rand_ptr, p, dist, nos_map, colors, prev);

            // start is the 0 node in case of a circle
            Vertex start = boost::vertex(0, g); //boost::graph_traits<Graph>::null_vertex();
            // or is any path-end

            BGL_FORALL_VERTICES_T(v, g, Graph) {
                if (boost::out_degree(v, g) == 1) {
                    start = v;
                    break;
                }
            }

            // Do a BGL DFS!
            // http://www.boost.org/doc/libs/1_53_0/libs/graph/doc/depth_first_search.html
            // http://www.boost.org/doc/libs/1_53_0/libs/graph/doc/undirected_dfs.html
            //boost::depth_first_search(g, boost::visitor(vis).root_vertex(start));
            boost::undirected_dfs(g, boost::root_vertex(start).visitor(vis).edge_color_map(boost::get(&edge_property::color, g)));

            BGL_FORALL_VERTICES_T(v, g, Graph) {
                g[v].base = colors[v];
            }
            return max_number_of_sequences;
        }

        template unsigned long long color_path_graph<std::mt19937> (Graph&, std::mt19937*);
    }
}