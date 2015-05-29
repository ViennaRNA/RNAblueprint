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
            ProbabilityMatrix(Graph& g);
            // get probability for ProbabilityKeys... key of Aks (12/A) (4/C) ()...
            unsigned long long get(ProbabilityKey pk);
            // My custom hash key used for n
            friend class ProbabilityKeyHash;
            ~ProbabilityMatrix();
        private:
            // map of possibilities saved by key
            ProbabilityMap pm;

        };
        
        // overload << operator to print ProbabilityKeys with any content
        std::ostream& operator<<(std::ostream& os, ProbabilityKey& m);


        /*
        
        struct SubProbability {
            int start;
            int end;
            int length;
        };
         * 
         inside class probabilityMAtrix
         * 
         public:
         * 
        // get sum of probabilities
            unsigned long long get_sum(int k, ProbabilityKey mykey, ProbabilityKey lastkey, std::vector<ProbabilityKey>& key_combinations, ProbabilityMap& probabilities);
            // get all the articulation points for a specific ear
            std::set<Vertex> get_Ak(unsigned int k);
            // return my of this graph

            unsigned int get_my() {
                return my;
            }
            // get number of sequences/solutions

            unsigned long long number_of_sequences() {
                return nos;
            }
         * 
         private:
         * 
        // remember my as number of ears
            unsigned int my = 0;
            // remember all the attachment points for every ear
            std::vector<std::set<Vertex> > Aks;
            // remember parts of ears (between attachment points)
            std::vector<std::vector<SubProbability> > parts;
            // max number of sequences/solutions in this pairing matrix
            unsigned long long nos = 0;
            // also remember pairing matrix
            PairingMatrix *p = nullptr;
            // to keep track of current Ak Ai we need to update them before glueing an ear
            void updateCurrentAkAi(Graph& g, int k, std::set<Vertex>& currentAk, std::set<Vertex>& currentAi);
            // calculate all the key probabilities and start more from there
            void calculate_probabilities(std::set<Vertex>& ap, ProbabilityKey& k);
            // actually calculate the probability for the given key
            unsigned long long get_probability(ProbabilityKey mykey, Graph& g, std::set<Vertex>& Ak, std::set<Vertex>& Ai, unsigned int k);
            // recursion to get base combinations done in (sum over AUGC in 6) of (sum over AUGC in 10) of ...
            void make_sum_of_sum(Graph& g, std::set<Vertex>& Ai, ProbabilityKey& mykey, ProbabilityKey& lastkey, std::vector<SubProbability>& sub_probabilities, unsigned int k, unsigned long long& max_number_of_sequences);
            // do the actual calculation of the multiplied probabilities in the sum_of_sum
            unsigned long long calculate_probability(ProbabilityKey& mykey, ProbabilityKey& lastkey, std::vector<SubProbability>& sub_probabilities);
            // get the color of either mykey or lastkey.
            int get_color_from_key(ProbabilityKey& mykey, ProbabilityKey& lastkey, int vertex);
            // calculates all base combinations for current attachment points (recursion)
            void calculate_combinations(std::set<int>& Ak, ProbabilityKey& mykey, std::vector<ProbabilityKey>& key_combinations);
         */
        
        
#endif	/* PROBABILITY_MATRIX_H */
    }
}