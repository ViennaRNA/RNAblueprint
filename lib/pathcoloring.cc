/* RNAblueprint
 * A program for designing RNA molecules.
 *
 * @date 13.08.2013
 * @author Stefan Hammer <s.hammer@univie.ac.at>
 * @copyright GPLv3
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
            std::set< int > articulations;
            
            BGL_FORALL_VERTICES_T(v, g, Graph) {
                if (g[v].articulation) {
                    if (boost::out_degree(v, g) <= 1) {
                        // remember vertex as int with sequence constraint
                        key[vertex_to_int(v, g)] = g[v].constraint;
                        articulations.emplace(vertex_to_int(v, g));
                    } else {
                        throw std::logic_error("There is a articulation vertex which is no path end in get_path_pm. This is not possible!");
                    }
                }
            }
            
            if (articulations.size() > 2) {
                throw std::logic_error("More than two articulation vertices in one path. ridiculous!");
            }
            
            //std::cerr << "keys: " << std::endl << keys;
            int length = boost::num_vertices(g) - 1;
            //std::cerr << "length: " << length << std::endl;
            PairingMatrix * p = PairingMatrix::Instance();
            
            ProbabilityMatrix result;
            PermuteKeyFactory pkf(key);
            
            while (true) {
                ProbabilityKey* k = pkf.key();
                switch ( articulations.size() )
                {
                    case 0:
                        result.put(*k, p->get(length, N, N));
                        break;
                    case 1:
                        result.put(*k, p->get(length, (*k)[*(articulations.begin())], N));
                        break;
                    case 2:
                        result.put(*k, p->get(length, (*k)[*(articulations.begin())], (*k)[*(++articulations.begin())]));
                        break;
                    default:
                        throw std::logic_error("More than two articulation vertices in one path. ridiculous!");
                }
                if (!pkf.next_permutation())
                    break;
            }
            return result;
        }
        
        template <typename RG>
        ProbabilityFraction color_path_graph(Graph& g, RG& rand) {

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

            // visitor declaration

            class color_dfs_visitor : public boost::default_dfs_visitor {
            public:

                color_dfs_visitor(ProbabilityFraction& probability_fraction, RG& rand, PairingMatrix * pair,
                        nosMap& n, std::unordered_map<Vertex, int>& c, int& prev)
                : pf(probability_fraction), r(rand), p(pair), nos_map(n), colors(c), previous(prev) {
                }
                ProbabilityFraction& pf;
                RG& r;
                PairingMatrix * p;
                nosMap& nos_map;
                std::unordered_map<Vertex, int>& colors;
                int& previous;

                void start_vertex(Vertex s, Graph g) const {
                    
                    if (debug) {
                        std::cerr << "Start vertex: " << s << " [" << enum_to_char(g[s].base) << "]" << std::endl;
                    }
                    pf.first = 1;
                    pf.second = 0;
                    for (auto b : base_conversion[ g[s].base ]) {
                        nos_map[s][b] = p->get(0, b, b);
                        if (debug) {
                            std::cerr << "v" << s << ": " << enum_to_char(b) << ": " << nos_map[s][b] << std::endl;
                        }
                        pf.second += nos_map[s][b];
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
                    pf.second = 0;
                    for (auto u_base : base_conversion[ g[u].base ]) {
                        nos_map[u][u_base] = 0;
                        for (auto v_base : base_conversion[ g[v].base ]) {
                            nos_map[u][u_base] += nos_map[v][v_base] * p->get(1, v_base, u_base);
                        }
                        // calculate maximal number of sequences on this vertex
                        pf.second += nos_map[u][u_base];
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
                    SolutionSizeType nos = 0;
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
                    RandomDistType dist(0, nos);
                    SolutionSizeType random = dist(r);

                    // stochastically take one of the possibilities
                    // start at the probability of first possible character and add each other base probability as long long as the random number is bigger.
                    SolutionSizeType sum = 0;
                    for (auto b : base_conversion[ g[u].base ]) {
                        if (p->get(1, b, previous) > 0) {
                            sum += nos_map[u][b];
                            // if the random number is bigger than our probability, take this base as the current base!
                            if (random < sum) {
                                colors[u] = b;
                                previous = b;
                                // in the probability fraction remember the sum of all choices made as numerator
                                pf.first *= nos_map[u][b]/nos;
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
            // remember the probability of the chosen solution as a fraction
            // numerator is the amount of solutions for the chosen value
            // denominator is the total amount of solutions for the problem
            ProbabilityFraction pf;
            PairingMatrix * p = PairingMatrix::Instance();
            nosMap nos_map;
            std::unordered_map<Vertex, int> colors;
            int prev = N;
            
            color_dfs_visitor vis(pf, rand, p, nos_map, colors, prev);
            
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
            
            // correct numerator of pf
            pf.first *= pf.second;
            return pf;
        }

        template ProbabilityFraction color_path_graph<std::mt19937> (Graph&, std::mt19937&);
    }
}
