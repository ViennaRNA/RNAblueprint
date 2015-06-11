/* This program reads secundary RNA structures in dot-bracket and
 * builds a graph for a latter ear-decomposition and bipartitness-check
 *
 * Created on: 18.03.2014
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

#ifndef PROBABILITY_MATRIX_H
#define	PROBABILITY_MATRIX_H


// include common header with graph definition and global variables
#include "common.h"
#include "pairing_matrix.h"

// include standard library parts
#include <unordered_map>
#include <functional>
#include <iomanip>

// include boost components
#include <boost/functional/hash.hpp>


namespace design {
    namespace detail {
        // class definitions
        
        // Class with custom hash function for the ProbabilityMatrix
        typedef std::map < int, int > ProbabilityKey;
        
        class ProbabilityKeyHash {
        public:
            std::size_t operator()(const ProbabilityKey& k) const;
        };
        
        // Class to get Pairing numbers
        typedef std::unordered_map< ProbabilityKey, unsigned long long, ProbabilityKeyHash> ProbabilityMap;
        
        class ProbabilityMatrix {
        public:
            ProbabilityMatrix();
            // get probability for ProbabilityKeys... key of Aks (12/A) (4/C) ()...
            unsigned long long get(ProbabilityKey pk);
            // fill a nos for a certain key
            void put(ProbabilityKey pk, unsigned long long nos);
            // get maximal number of sequences for the whole matrix/subgraph
            unsigned long long nos();
            
            
            // My custom hash key used for n
            friend class ProbabilityKeyHash;
            ~ProbabilityMatrix();
        private:
            // map of possibilities saved by key
            ProbabilityMap pm;
            // A function to permute the keys
            std::vector<ProbabilityKey> permute_key(ProbabilityKey pk);
            void permute_impl(ProbabilityKey::iterator start, ProbabilityKey::iterator end, std::vector<ProbabilityKey>& result);
            // TODO a function to "multiply" probability matrixes
        };
        
        // overload << operator to print ProbabilityKeys with any content
        std::ostream& operator<<(std::ostream& os, ProbabilityKey& m);
    }
}

#endif	/* PROBABILITY_MATRIX_H */