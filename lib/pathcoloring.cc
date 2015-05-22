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
    unsigned long long generate_path_seq (Sequence& sequence, int first, int last, int length, RG* rand_ptr) {

      // pairing matrix for every length
      Pairing p(length + 1); //TODO initialize only once for the whole program as static content!
      // set maximum possible number of sequences for first....last
      unsigned long long max_number_of_sequences = p.get(length, first, last);
      // declare random number distribution and get a random number
      std::uniform_real_distribution<float> dist(0, 1);
      // number of possible sequences at each possible step
      unsigned long long number_of_sequences = 0;
      // container to remember possible letters after our current letter
      std::vector< int > posibilities;

      // set the first base
      if (first < A_Size) {
        sequence.push_back(first);
        length--;
      }

      /*if (debug) {
        std::cerr << "Max Number of Sequences is: " << max_number_of_sequences << std::endl;
        std::cerr << "Length is: " << length << std::endl;
        std::cerr << "Sequence is: " << sequence << std::endl;
        std::cerr << "First is: " << enum_to_char(first) << std::endl;
        std::cerr << "Last is: " << enum_to_char(last) << std::endl;
      }*/

      while (length >= 0) {
        if (first >= A_Size) {
          number_of_sequences = p.get(length, first, last);
          first = N;
        } else {
          number_of_sequences = p.get(length + 1, first, last);
        }
        // look in paring matrix for next possible character and remember them
        posibilities.clear();
        for (auto i : base_conversion[ N ]) {
          if (p.get(1, first, i) >= 1) {
            posibilities.push_back(i); 
          }
        }

        // get a random number between 0 and 1.
        float random = dist(*rand_ptr);

        // stochastically take one of the possibilities
        // start at the probability of first possible character and add each other base probability as long long as the random number is bigger.
        unsigned long long sum = 0;
        for (auto base : posibilities) {
          sum += p.get(length, base, last);
          // if the random number is bigger than our probability, take this base as the first base!
          if (random * number_of_sequences < sum) {
            sequence.push_back(base);
            // our new begin is the chosen base.
            first = base;
            length--;
            // don't forget to exit the loop, otherwise will always be first = C;
            break;
          }
        }
        if (sum == 0) {
          std::cerr << std::endl << "The requested sequence cannot be colored! Conflict at: "
              << enum_to_char(first) << ", " << enum_to_char(last) << ", " << length << std::endl;
          exit(1);
        }

        /*if (debug) {
          std::cerr << "Number of Sequences is: " << number_of_sequences << std::endl;
          std::cerr << "Random is: " << random*number_of_sequences << std::endl;
          std::cerr << "Possibilities is: " << posibilities << std::endl;
          std::cerr << "Sequence is: " << sequence << std::endl;
          std::cerr << "--->" << std::endl;
          std::cerr << "First is: " << enum_to_char(first) << std::endl;
          std::cerr << "Last is: " << enum_to_char(last) << std::endl;
          std::cerr << "Length is: " << length << std::endl;
        }*/
      }
      return max_number_of_sequences;
    }
    template <typename RG>
    unsigned long long generate_cycle_seq (Sequence& sequence, int first, int length, RG* rand_ptr) {

      // max number of sequences to return
      unsigned long long max_number_of_sequences = 0;
      // check if length is even number
      if (length % 2 != 0) {
        std::cerr << std::endl << "Length of the cycle to color is an odd number. This can't be!" << std::endl;
        exit(1);
      }

      if (first < A_Size) {
        // return a path with same begin and end, but then remove the last character again -> cycle!
        max_number_of_sequences = generate_path_seq(sequence, first, first, length, rand_ptr);
        
      } else {
        
        Pairing p(length + 1); //TODO initialize only once for the whole program as static content!
        // declare random number distribution and get a random number
        std::uniform_real_distribution<float> dist(0, 1);
        
        for (auto i : base_conversion[ first ]) {
          max_number_of_sequences += p.get(length, i, i);
        }
        
        // get a random number between 0 and 1.
        float random = dist(*rand_ptr);
        // stochastically take one of the possibilities
        // start at the probability of first possible character and add each other base probability as long long as the random number is bigger.
        unsigned long long sum = 0;
        for (auto base : base_conversion[ first ]) {
          sum += p.get(length, base, base);
          // if the random number is bigger than our probability, take this base as the first base!
          if (random * max_number_of_sequences < sum) {
            // our new begin is the chosen base.
            first = base;
            // don't forget to exit the loop, otherwise will always be first = C;
            break;
          }
        }
        generate_path_seq(sequence, first, first, length, rand_ptr);
      }
      
      sequence.pop_back();
      return max_number_of_sequences;
    }

    template <typename RG>
    unsigned long long color_path_cycle_graph (Graph& g, RG* rand_ptr) {

      unsigned long long max_number_of_sequences = 0;
      
      // check if given graph is indeed a path with max_degree = 2 and two ends with degree = 1;
      int max_degree;
      int min_degree;
      std::tie(min_degree, max_degree) = get_min_max_degree(g);
      
      if (max_degree > 2) {
        std::cerr << std::endl << "This graph is no cycle or path (max degree > 2). I can't color this!" << std::endl;
        exit(1);
      }
      
      Vertex start;
      // start is any path-end

        BGL_FORALL_VERTICES_T(v, g, Graph) {
            if (boost::out_degree(v, g) == 1) {
                start = v;
                break;
            }
        }
        // or any node in case of a circle
        if (start == boost::graph_traits<Graph>::null_vertex()) {
            start = boost::vertex(0, g)
        }

        class color_dfs_visitor : public boost::default_dfs_visitor {
        public:

            color_dfs_visitor(unsigned long long& max_number_of_sequences) : mnos(max_number_of_sequences) {
            }
            unsigned long long& mnos;
            Pairing p(length + 1);
            std::map<Vertex, std::map<int, unsigned long long>  nos_map;
            unsigned int length;
            Vertex begin;

            void start_vertex(Vertex s, Graph g) const {
                if (debug) {
                    std::cerr << "Start vertex: " << s << std::endl;
                }
                begin = s;
            }

            void discover_vertex(Vertex u, Graph g) const {
                if (debug) {
                    std::cerr << "Detecting vertex: " << u << std::endl;
                }
                
                unsigned long long nos;
                
                if (g[v].base != N) {
                    for (auto base : base_conversion[ g[v].base ]) {
                        nos_map[v][base] = p.get(length, begin, base);
                        // check nos_map if last part was with a fixed end
                        // if so add last part + numbers for this part
                        // if not, calculate the combinatorics
                        // 
                        // for (auto last : base_conversion[ g[v].base ]) {
                        // for possibilities of begin = begin_base {
                        // nos += nos_map[begin][begin_base] * p.get(length, begin_base, last);
                        //}
                        // nos_map[v][last] = nos;
                        // nos = 0;
                        // }
                    }
                    length = 0;
                } else {
                    length++;
                }
            }

            void forward_or_cross_edge(Edge e, Graph g) const {
                if (debug) {
                    std::cerr << "Detecting back-edge: " << e << std::endl;
                }
                //Vertex u = boost::source(e, g);
                //Vertex v = boost::target(e, g);
                c.push_back(e);
            }
        };

        my_dfs_visitor vis(max_number_of_sequences);

        // Do a BGL DFS!
        // http://www.boost.org/doc/libs/1_53_0/libs/graph/doc/depth_first_search.html
        // Did not work: http://www.boost.org/doc/libs/1_53_0/libs/graph/doc/undirected_dfs.html
        // boost::undirected_dfs(g, boost::visitor(vis), vcolorMap, ecolorMap, rootVertex);
        boost::depth_first_search(g, visitor(vis).root_vertex(start));
        
        return max_number_of_sequences;
    }
/*
      // check if given graph is indeed a path with max_degree = 2 and two ends with degree = 1;
      int max_degree;
      int min_degree;
      std::tie(min_degree, max_degree) = get_min_max_degree(g);

      if (max_degree > 2) {
        std::cerr << std::endl << "This graph is no cycle or path (max degree > 2). I can't color this!" << std::endl;
        exit(1);
      }

      // find out the length of the path
      unsigned int length = boost::num_edges(g);

      // find out the degree of the first and the last base and the base of all non 'N' colored
      std::vector< Vertex > ends;
      std::vector< Vertex > colored_bases;
      Sequence sequence;

      BGL_FORALL_VERTICES_T(v, g, Graph) {
        // remember ends of the path
        if (boost::out_degree(v, g) == 1) {
          ends.push_back(v);
        } else if (g[v].base != N) {
          // remember non N bases which are no ends
          colored_bases.push_back(v);
        }

        // reset vertex.color tag to 0 -> must be done for sequencestring_to_graph
        g[v].color = 0;
      }

      if (ends.size() == 2) {
        // it is a path!
        // check if any of the non-end vertices have assigend colors
        if (colored_bases.size() > 0) {
          std::cerr << std::endl << "This path already is partly colored in between. I can't color this!" << std::endl;
          exit(1);
        }
        // call generate_path_seq and color the vertices accordingly
        max_number_of_sequences = generate_path_seq(sequence, g[ends[0]].base, g[ends[1]].base, length, rand_ptr);
        if (debug) {
          std::cerr << "Sequence is: " << sequence << std::endl;
        }
        // assign this sequence of bases to the graph
        sequencestring_to_graph(g, ends[0], sequence);

      } else if (max_degree == 2 && min_degree == 2 && ends.size() == 0) {
        // it is a cycle (all vertices degree 2)
        if (colored_bases.size() == 1) {
          // start to color at exact this vertex
          max_number_of_sequences = generate_cycle_seq(sequence, g[colored_bases[0]].base, length, rand_ptr);
          if (debug) {
            std::cerr << "Sequence is: " << sequence << std::endl;
          }
          // assign this sequence of bases to the graph
          sequencestring_to_graph(g, colored_bases[0], sequence);
        } else if (colored_bases.size() == 0) {
          // start to color at any vertex with N
          max_number_of_sequences = generate_cycle_seq(sequence, g[boost::vertex(0, g)].base, length, rand_ptr);
          if (debug) {
            std::cerr << "Sequence is: " << sequence << std::endl;
          }
          // assign this sequence of bases to the graph
          sequencestring_to_graph(g, boost::vertex(0, g), sequence);
        } else {
          //TODO in this case, split the path, cycle at the constraints and color in between with fixed ends.
          std::cerr << std::endl << "This cycle already is partly colored in between. I can't color this!" << std::endl;
          exit(1);
        }
      } else if (length == 0 && boost::num_vertices(g) == 1) {
        // its a single vertex!
        max_number_of_sequences = generate_path_seq(sequence, g[boost::vertex(0, g)].base, N, length, rand_ptr);
        if (debug) {
          std::cerr << "Sequence is: " << sequence << std::endl;
        }
        sequencestring_to_graph(g, boost::vertex(0, g), sequence);
      } else {
        // this is no path  - more than two "ends"
        std::cerr << std::endl << "This graph is no cycle or path. I can't color this!" << std::endl;
        exit(1);
      }
*/

    void sequencestring_to_graph (Graph& g, Vertex vertex, Sequence& sequence) {
      // check if sequence is empty
      if (sequence.size() == 0) {
        std::cerr << std::endl << "Tried to color a path/cycle with a base sequence of zero length. This is impossible!" << std::endl;
        exit(1);
      }

      // check if we are going to overwrite an already existing assignment with a different one.
      if (g[vertex].base != N) {
        if (g[vertex].base != sequence.front()) {
          std::cerr << "Tried to color following vertex with a base, but it is already colored with another one! "
              << boost::get(boost::vertex_color_t(), g, vertex) << "/" << enum_to_char(g[vertex].base) << ", new color: " << enum_to_char(sequence.front()) << std::endl;
          exit(1);
        }
      }

      // assign the current vertex with the first character
      g[vertex].base = sequence.front();
      // delete the first character
      sequence.pop_front();
      // mark as done (color = 1)
      g[vertex].color = 1;

      // iterate over all adjacent vertices and find the uncolored one.
      // recursively call this function on this vertex again

      BGL_FORALL_ADJ_T(vertex, adj, g, Graph) {
        if (g[adj].color == 0) {
          sequencestring_to_graph(g, adj, sequence);
        }
      }
    }
    
    template unsigned long long generate_path_seq<std::mt19937> (std::deque< int >&, int, int, int, std::mt19937*);
    template unsigned long long generate_cycle_seq<std::mt19937> (std::deque< int >&, int, int, std::mt19937*);
    template unsigned long long color_path_cycle_graph<std::mt19937> (Graph&, std::mt19937*);
  }
}