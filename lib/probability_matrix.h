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
#include <algorithm>
#include <iterator>

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
            //unsigned long long get(ProbabilityKey pk);
            unsigned long long operator[] (ProbabilityKey& pk);
            // fill a nos for a certain key
            void put(ProbabilityKey& pk, unsigned long long nos);            
            // get maximal number of sequences for the whole matrix/subgraph
            unsigned long long mnos();
            // get set of special vertices
            std::set< int > getSpecials() { return specials; };
            // My custom hash key used for n
            friend class ProbabilityKeyHash;
            friend ProbabilityMatrix operator* (ProbabilityMatrix& x, ProbabilityMatrix& y);
            ~ProbabilityMatrix();
        private:
            // map of possibilities saved by key
            ProbabilityMap pm;
            // remember all special points
            std::set< int > specials;
            // TODO a function to "multiply" probability matrixes
        };
        
        
        // multiply operator overloaded which calculates new pm
        ProbabilityMatrix operator* (ProbabilityMatrix& x, ProbabilityMatrix& y);
        
        // a function which creates a new ProbabilityMatrix where the given key is internal
        ProbabilityMatrix make_internal(ProbabilityMatrix& pm, int v);
        
        // A function to permute the keys if the key contains Letters bigger than the alphabet (eg N, S, Y,...)
        std::vector<ProbabilityKey> permute_key(ProbabilityKey pk);
        // helper for permute_key
        void permute_impl(ProbabilityKey::iterator start, ProbabilityKey::iterator end, std::vector<ProbabilityKey>& result, ProbabilityKey current);
        
        // overload << operator to print ProbabilityKeys
        std::ostream& operator<<(std::ostream& os, ProbabilityKey& m);
        
        // overload << operator to print ProbabilityMatrix
        std::ostream& operator<<(std::ostream& os, ProbabilityMatrix& m);
    }
}

#endif	/* PROBABILITY_MATRIX_H */