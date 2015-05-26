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
    
    Fibonacci::Fibonacci (unsigned int length)
    : numbers (length) {
      // Definition: F1 = 0, F2 = 1, Fn = Fn-1 + Fn-2
      numbers[0] = 0;
      numbers[1] = 1;
      for (unsigned int n = 2; n < length; n++) {
        numbers[n] = numbers[n - 1] + numbers[n - 2];
      }
    }

    Pairing::Pairing (unsigned int l)
    : p (l + 1)
    , length (l + 1) {
      // Definition:
      // length 1: p[A][U][1] = 1, p[U][A][1] = 1, p[G][C][1] = 1, p[C][G][1] = 1, p[U][G][1] = 1, p[G][U][1] = 1
      // if length greater than 2, we don't care about the first letter any more -> x
      // length n: 	p[N][A][n] = p[N][U][n-1]
      //		p[N][C][n] = p[N][G][n-1]
      //		p[N][G][n] = p[N][U][n-1] + p[N][C][n-1]
      //		p[N][U][n] = p[N][A][n-1] + p[N][G][n-1]

      // fill standard pairing matrix (pathlength = 1)
      p[1][A][U] = 1;
      p[1][U][A] = 1;
      p[1][G][C] = 1;
      p[1][C][G] = 1;
      p[1][G][U] = 1;
      p[1][U][G] = 1;

      //fill length 0 with probability 1 for same base (important for setting last base)
      p[0][A][A] = 1;
      p[0][U][U] = 1;
      p[0][G][G] = 1;
      p[0][C][C] = 1;

      // fill pathlength up to length (can be done with matrix multiplication of the pairing matrix
      for (unsigned int i = 2; i <= l; i++) {
        p[i] = multiply(p[i - 1], p[1]);
      }
/*
      if(debug) {
        std::cerr << "Pairing constructor called and filled" << std::endl;
        for (unsigned int k = 0; k <= l; k++) {
          std::cerr << std::endl << k << ":" << std::endl;
          std::cerr  << "  " << enum_to_char(0) << " " << enum_to_char(1) << " " << enum_to_char(2) << " " << enum_to_char(3) << " " << std::endl;
          for (unsigned int i = 0; i < A_Size; i++) {
            std::cerr << enum_to_char(i) << " ";
            for (unsigned int j = 0; j < A_Size; j++) {
              std::cerr << get(k, i, j) << " ";
            }
            std::cerr << std::endl;
          }
        }
      }*/
    }

    rnaMatrix Pairing::multiply (rnaMatrix A, rnaMatrix B) {
      rnaMatrix C;
      int i, j, k;
      long long sum;
      for (i = 0; i < A_Size; i++) {
        for (j = 0; j < A_Size; j++) {
          sum = 0;
          for (k = 0; k < A_Size; k++) {
            sum += A[i][k] * B[k][j];
          }
          C[i][j] = sum;
        }
      }
      return C;
    }

    unsigned long long Pairing::get (unsigned int l, unsigned int b1, unsigned int b2) {

      // if we request a probability for an unknown (N) character at one or both ends, 
      // return the sum of the probabilities for all characters at this position
      
      //std::cerr << "b1 is " << enum_to_char(b1) << b1 << ", b2 is " << enum_to_char(b2) << b2 << std::endl;
      
      if ((b1 >= A_Size) || (b2 >= A_Size)) {
        if ((b1 >= A_Size) && (b2 >= A_Size)) {
          //std::cerr << "b1>=Abet; b2>=Abet" << std::endl;
          unsigned long long sum = 0;
          for (auto i : base_conversion[b2]) {
            sum += get(l, b1, i);
          }
          return sum;
          
        } else if (b1 >= A_Size) {
          //std::cerr << "b1>=Abet" << std::endl;
          unsigned long long sum = 0;
          for (auto i : base_conversion[b1]) {
            sum += get(l, i, b2);
          }
          return sum;
          
        } else if (b2 >= A_Size) {
          //std::cerr << "b2>=Abet" << std::endl;
          return get(l, b2, b1);
        }
      } else {
        //std::cerr << "b1<Abet; b2<Abet" << std::endl;
        if ((l > length) || (b1 >= A_Size) || (b2 >= A_Size)) {
          // check if the requested length is bigger than our initialization or that a base bigger than 3 is requested
          // -> to avoid segfaults or unknown behavior!
          std::cerr << "Requested a value in pairing matrix which is out of range: p[" << l << "][" << b1 << "][" << b2 << "]" << std::endl;
          exit(1);
        }
        
        return p[l][b1][b2];
      }
    }
    
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
            color_dfs_visitor(unsigned long long& max_number_of_sequences, RG * rand_ptr, Pairing& pair, std::uniform_real_distribution<float>& d,
                                nosMap& n, std::unordered_map<Vertex, int>& c, int& prev) 
            : mnos(max_number_of_sequences), r_ptr(rand_ptr), p(pair), dist(d), nos_map(n), colors(c), previous(prev) {}
            unsigned long long& mnos;
            RG * r_ptr;
            Pairing& p;
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
                        nos_map[u][u_base] += nos_map[v][v_base] * p.get(1, v_base, u_base);
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
                
                previous = g[u].base;
                // re-calculate mnos for circle closure
                mnos = 0;
                for (auto b : base_conversion[ g[v].base ]) {
                    if (p.get(1, b, previous) > 0) {
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
                    if (p.get(1, b, previous) > 0) {
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
                    if (p.get(1, b, previous) > 0) {
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
        
        Pairing p(1);
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