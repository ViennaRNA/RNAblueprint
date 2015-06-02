/* This program reads secondary RNA structures in dot-bracket and
 * builds a graph for a latter ear-decomposition and bipartitness-check.
 *
 * Created on: 25.03.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

#ifndef COMMON_H
#define COMMON_H


// include standard library parts
#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <set>
#include <iostream>
#include <utility>

// include boost components
#include <boost/lexical_cast.hpp>
#include <boost/config.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

namespace boost {
  struct edge_component_t {

    enum {
      num = 555
    };
    typedef edge_property_tag kind;
  };
}

namespace design {
  namespace detail {
    
    //Global variables
    extern bool debug;
    extern bool * debug_ptr;

    // define size of the alphabet
    #define A_Size 4
    // Encode Bases to enums

    /*
      A = adenine
      C = cytosine
      G = guanine
      T = thymine
      R = G A (purine)
      Y = T C (pyrimidine)
      K = G T (keto)
      M = A C (amino)
      S = G C (strong bonds)
      W = A T (weak bonds)
      B = G T C (all but A)
      D = G A T (all but C)
      H = A C T (all but G)
      V = G C A (all but T)
      N = A G C T (any)
    */
    
    // it is important that 0-3 are the basic bases
    enum bases {
      A, C, G, U, R, Y, K, M, S, W, B, D, H, V, N
    };
    
    static std::unordered_map<int, std::set<int> >  base_conversion = 
    {
        { A, { A } },
        { C, { C } },
        { G, { G } },
        { U, { U } },
        { R, { G, A } },
        { Y, { U, C } },
        { K, { G, U } },
        { M, { A, C } },
        { S, { G, C } },
        { W, { A, U } },
        { B, { G, U, C } },
        { D, { G, A, U } },
        { H, { A, C, U } },
        { V, { G, C, A } },
        { N, { A, G, C, U } }
    };

    // Typedef for sequences of enums
    typedef std::deque< int > Sequence;

    // get property with Graph[Vertex].bipartite_color = int;

    struct vertex_property {
      int bipartite_color = 0;
      int color = 0;
      std::set< int > Ak;
      int Ai = 0;
      int base = N;
    };

    struct edge_property {
      int ear = 0;
      int color = 0;
    };
    
    struct graph_property {
        std::string id;
    };
    
    enum graph_property_t {
        gpt
    };

    // graph_properties 
    // boost graph template
    typedef boost::subgraph< 
    boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS,
    boost::property< boost::vertex_color_t, int, vertex_property >,
    boost::property< boost::edge_index_t, int, boost::property < boost::edge_component_t, std::size_t, edge_property> >,
    boost::property< graph_property_t, graph_property> > > Graph;
    typedef Graph::edge_descriptor Edge;
    typedef Graph::vertex_descriptor Vertex;

    // function to make enum a char again and other way round
    char enum_to_char (int intletter);
    int char_to_enum (char charletter);

    // overload << operator to print vectors with any content

    template <typename T>
    std::ostream& operator<< (std::ostream& os, std::vector<T>& vec) {
      int i = 0;
      for (auto elem : vec) {
        os << "(" << i++ << ") " << elem << std::endl;
      }
      return os;
    }

    // overload << operator to print maps with any content

    template <typename U, typename V>
    std::ostream& operator<< (std::ostream& os, std::map<U, V>& m) {
      for (typename std::map<U, V>::iterator it = m.begin(); it != m.end(); it++) {
        os << it->first << "," << it->second << std::endl;
      }
      return os;
    }

    // overload << operator to print sets with any content

    template <typename W>
    std::ostream& operator<< (std::ostream& os, std::set<W>& s) {
      for (auto elem : s) {
        os << elem << ", ";
      }
      return os;
    }

    // overload << operator to print deques with sequence information
    std::ostream& operator<< (std::ostream& os, Sequence& sequence);

    // matrix template
    template < class T, size_t ROWS, size_t COLS > using matrix = std::array< std::array< T, COLS >, ROWS >;
  }
}

#endif