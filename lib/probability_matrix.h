/*!\file probability_matrix.h 
 * \brief This file holds the class definitions for the Probability Matrix.
 *  
 * The ProbabilityMatrix stores the number of solutions for a particular subgraph as value
 * for a certaion ProbabilityKey.
 * The ProbabilityKey contains a list of pairs, the vertex (position) and its current base-color.
 *
 * Created on: 18.03.2014
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 * 
 * \cond INTERNAL
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
        typedef std::unordered_map< ProbabilityKey, boost::multiprecision::mpz_int, ProbabilityKeyHash> ProbabilityMap;

        class ProbabilityMatrix {
        public:
            ProbabilityMatrix();
            // default copy constructor
            //ProbabilityMatrix( const ProbabilityMatrix &pm ) = default;
            // get probability for ProbabilityKeys... key of Aks (12/A) (4/C) ()...
            //boost::multiprecision::mpz_int get(ProbabilityKey pk);
            boost::multiprecision::mpz_int operator[](ProbabilityKey& pk);
            // fill a nos for a certain key
            void put(ProbabilityKey& pk, boost::multiprecision::mpz_int nos);
            // get maximal number of sequences for the whole matrix/subgraph
            boost::multiprecision::mpz_int mnos();
            // return if this is a initialized PM or not

            bool is_initialized() {
                return initialized;
            }
            // get set of special vertices

            std::set< int > getSpecials() {
                return specials;
            };
            // sample one combination randomly given a ProbabilityKey with the constraints and a random number generator
            // return a pair with the chosen ProbabilityKey and the number_of_sequences for the given input constraints
            template <typename R>
            std::pair<ProbabilityKey, boost::multiprecision::mpz_int> sample(ProbabilityKey pk, R* rand_ptr);
            // special case without probability key (and therefore constraints)
            template <typename R>
            std::pair<ProbabilityKey, boost::multiprecision::mpz_int> sample(R* rand_ptr);
            // My custom hash key used for n
            friend class ProbabilityKeyHash;
            friend ProbabilityMatrix operator*(ProbabilityMatrix& x, ProbabilityMatrix& y);
            ~ProbabilityMatrix() = default;
        private:
            // map of possibilities saved by key
            ProbabilityMap pmap;
            // remember if initialized
            bool initialized = false;
            // remember all special points
            std::set< int > specials = {};
            // TODO a function to "multiply" probability matrixes
        };

        class PermuteKeyFactory {
        public:
            ProbabilityKey* PermuteKeyFactory(ProbabilityKey pk);
            bool next_permutation();
            bool previous_permutation();
        private:
            std::map<int, std::list<int> > storage;
            std::map<int, std::list<int>::iterator> state;
            std::map<int, std::list<int>::iterator>::iterator state_it;
            std::map<int, int&> current;
        };

        // multiply operator overloaded which calculates new pm
        ProbabilityMatrix operator*(ProbabilityMatrix& x, ProbabilityMatrix& y);

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

/* 
 * \endcond INTERNAL
 */