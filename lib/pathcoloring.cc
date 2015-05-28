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
#include "graphcommon.h"

namespace design {
    namespace detail {
    
    template <typename RG>
    unsigned long long color_path_cycle_graph (Graph& g, RG* rand_ptr) {

        unsigned long long max_number_of_sequences = 0;

        // check if given graph is indeed a path with max_degree = 2 and two ends with degree = 1;
        int max_degree;
        int min_degree;
        std::tie(min_degree, max_degree) = get_min_max_degree(g);
        
        // assert path or cycle
        if (max_degree > 2) {
          std::cerr << std::endl << "This graph is no cycle or path (max degree > 2). I can't color this!" << std::endl;
          exit(1);
        }
        
        // visitor declaration
        class color_dfs_visitor : public boost::default_dfs_visitor {
        public:
            color_dfs_visitor(unsigned long long& max_number_of_sequences, RG * rand_ptr, PairingMatrix * pair, std::uniform_real_distribution<float>& d,
                                nosMap& n, std::unordered_map<Vertex, int>& c, int& prev) 
            : mnos(max_number_of_sequences), r_ptr(rand_ptr), p(pair), dist(d), nos_map(n), colors(c), previous(prev) {}
            unsigned long long& mnos;
            RG * r_ptr;
            PairingMatrix * p;
            std::uniform_real_distribution<float>& dist;
            nosMap& nos_map;
            std::unordered_map<Vertex, int>& colors;
            int& previous;

            void start_vertex(Vertex s, Graph g) const {
                if (debug) {
                    std::cerr << "Start vertex: " << s << std::endl;
                }
                
                mnos = 0;
                for (auto b : base_conversion[ g[s].base ]) {
                    nos_map[s][b] = 1;
                    std::cerr << s << ":" << enum_to_char(b) << ":" << nos_map[s][b] << std::endl;
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
                    std::cerr << u << ":" << enum_to_char(u_base) << ":" << nos_map[u][u_base] << std::endl;
                }
            }

            void back_edge(Edge e, Graph g) const {
                if (debug) {
                    std::cerr << "Detecting back-edge (graph is a cycle): " << e << std::endl;
                }
                Vertex v = boost::source(e, g);
                Vertex u = boost::target(e, g);
                //TODO really necessary?
                previous = g[u].base;
                // re-calculate mnos for circle closure
                mnos = 0;
                for (auto b : base_conversion[ g[v].base ]) {
                    if (p->get(1, b, previous) > 0) {
                        mnos += nos_map[v][b];
                    }
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
                  std::cerr << std::endl << "The requested sequence cannot be colored! Conflict at: "
                      << enum_to_char(g[u].base) << ", " << enum_to_char(previous) << std::endl;
                  exit(1);
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
                    std::cerr << "Vertex colored: " << u << "/" << enum_to_char(colors[u]) << std::endl;
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
        
        BGL_FORALL_VERTICES_T(v, g, Graph) { g[v].base = colors[v]; }
        return max_number_of_sequences;
    }

    template unsigned long long color_path_cycle_graph<std::mt19937> (Graph&, std::mt19937*);
  }
}