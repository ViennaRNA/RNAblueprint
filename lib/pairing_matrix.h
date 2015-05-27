/* This program reads secundary RNA structures in dot-bracket and
 * builds a graph for a latter ear-decomposition and bipartitness-check
 *
 * Created on: 27.05.2015
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

#ifndef PAIRING_MATRIX_H
#define	PAIRING_MATRIX_H

// include common header with graph definition and global variables
#include "common.h"

// include standard library parts
#include <sstream>
#include <unordered_map>

namespace design {
    namespace detail {
        //typedefs
        typedef matrix< unsigned long long, A_Size, A_Size > rnaMatrix;

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
            static Pairing* Instance();
            // interface
            unsigned long long get(unsigned int l, unsigned int b1, unsigned int b2);
        protected:
            Pairing();
            ~Pairing();
        private:
            static Pairing * _instance;
            std::vector< rnaMatrix > p;
            void extend(int newlength);
            rnaMatrix multiply(rnaMatrix a, rnaMatrix b);
            unsigned int length;
        };
    }
}


#endif	/* PAIRING_MATRIX_H */
