/* RNAdesign
 * A program for designing RNA molecules.
 *
 * Created on: 18.03.2014
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

// include header
#include "probability_matrix.h"

namespace design {
    namespace detail {

        std::size_t ProbabilityKeyHash::operator()(const ProbabilityKey& k) const {
            // Start with 0 as a hash value   .
            std::size_t hash = 0;

            for (auto elem : k) {
                boost::hash_combine(hash, boost::hash_value(elem.first));
                boost::hash_combine(hash, boost::hash_value(elem.second));
            }

            //std::cerr << "hash is: " << hash << std::endl;

            return hash;
        }

        ProbabilityMatrix::ProbabilityMatrix() {

            if (debug) {
                std::cerr << "Initializing ProbabilityMatrix..." << std::endl;
            }
            
        }

        unsigned long long ProbabilityMatrix::get(ProbabilityKey pk) {

            unsigned long long returnvalue = 0;
            
            for (auto elem : pk) {
                if (elem.second >= A_Size) {
                    
                }
            }
            

            // important for map: if you request with [] an entry will be created for unexisting ones.
            ProbabilityMap::const_iterator found = pm.find(pk);

            if (found != pm.end()) {
                returnvalue = found->second;
            } else {
                returnvalue = 0;
            }
            
            return returnvalue;
        }
        
        std::vector<ProbabilityKey> ProbabilityMatrix::permute_key(ProbabilityKey pk) {
            std::vector<ProbabilityKey> result;
            result.push_back(ProbabilityKey());
            
            permute_impl(pk.begin(), pk.end(), result);
            
            return result;
        }
    /*
        Function Hypercube(int dimensions, int current, string partialCoords)
        {
          for i=0, i<=steps, i++
          {
            if(current==dimensions)
              print partialCoords + ", " + i + ")/n";
            else if current==0
              Hypercube(dimensions, current+1, "( "+i);
            else
              Hypercube(dimensions, current+1, partialCoords+", "+i);
          }

        }
        */
        void ProbabilityMatrix::permute_impl(ProbabilityKey::iterator start, ProbabilityKey::iterator end, std::vector<ProbabilityKey>& result) {
            for (auto b : base_conversion[start->second]) {
                
                result.back()[start->first] = b;
                
                if (start == end) {
                    result.push_back(ProbabilityKey());
                } else {
                    permute_impl(start++, end, result);
                }
            }
        }
        
        void ProbabilityMatrix::put(ProbabilityKey pk, unsigned long long nos) {
            
        }
        
        unsigned long long ProbabilityMatrix::nos() {
            unsigned long long nos = 0;
            // todo: calculate nos!
            return nos;
        }
        
        
        ProbabilityMatrix::~ProbabilityMatrix() {
            // cleanup?
        }
        
        // overload << operator to print ProbabilityKeys in pair with graph
        std::ostream& operator<<(std::ostream& os, ProbabilityKey& m) {
            os << "[";
            for (typename ProbabilityKey::iterator it = m.begin(); it != m.end(); it++) {
                os << "(" << std::setfill(' ') << std::setw(1) << it->first << "," << enum_to_char(it->second) << ")";
            }
            os << "]";
            return os;
        }

    }
}