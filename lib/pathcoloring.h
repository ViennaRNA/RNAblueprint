/* This program reads secundary RNA structures in dot-bracket and
 * builds a graph for a latter ear-decomposition and bipartitness-check
 *
 * Created on: 13.08.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

#ifndef PATHCOLORING_H
#define PATHCOLORING_H

    // include common header with graph definition and global variables
    #include "common.h"

    // include standard library parts
    #include <sstream>
    #include <unordered_map>

    // include boost components
    //#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/undirected_dfs.hpp>

namespace design {
    namespace detail {
        // typedefs
        typedef matrix< unsigned long long, A_Size, A_Size > rnaMatrix;
        typedef std::unordered_map< Vertex, std::unordered_map<int, unsigned long long> > nosMap;

        // class definitions
        // Class to get Fibonacci numbers

        class Fibonacci {
        public:
            Fibonacci(unsigned int l);

            unsigned int get(unsigned int n) {
                return numbers[n - 1];
            };
        private:
            std::vector< unsigned int > numbers;
        };

        // Class to get Pairing numbers implemented as a Singleton to be unique
        class Pairing {
        public:
            static Pairing* Instance(unsigned int length);
            // interface
            unsigned long long get(unsigned int l, unsigned int b1, unsigned int b2);
        protected:
            Pairing(unsigned int length);
            ~Pairing();
        private:
            static Pairing * _instance;
            std::vector< rnaMatrix > p;
            rnaMatrix multiply(rnaMatrix a, rnaMatrix b);
            unsigned int length;
        };

        template <typename RG>
        unsigned long long color_path_cycle_graph(Graph& g, RG* rand_ptr);
    }
}
#endif
