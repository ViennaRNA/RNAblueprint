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

        unsigned long long ProbabilityMatrix::operator[] (ProbabilityKey& pk) {
            // sanity check if the vertices requested are indeed stored in this pm
            for (auto pair : pk) {
                if ( specials.find(pair.first) == specials.end() ) {
                    std::cerr << "Tried to get a not allowed Vertex from ProbabilityMatrix." << std::endl;
                    exit(1);
                }
            }
            
            unsigned long long returnvalue = 0;
            std::vector<ProbabilityKey> allkeys = permute_key(pk);
            
            for (auto key : allkeys) {
                // important for map: if you request with [] an entry will be created for unexisting ones.
                ProbabilityMap::const_iterator found = pm.find(key);

                if (found != pm.end()) {
                    returnvalue += found->second;
                }
            }
            return returnvalue;
        }
        
        void ProbabilityMatrix::put (ProbabilityKey& pk, unsigned long long nos) {
            
            if (pm.size() == 0) {
                for (auto pair : pk) {
                    specials.insert(pair.first);
                }
            }
            
            for (auto pair : pk) {
                // sanity check if the bases requested are within our alphabet size
                if (pair.second >= A_Size) {
                    std::cerr << "Tried to write a base outside of the alphabet size into ProbabilityMatrix." << std::endl;
                    exit(1);
                // sanity check if the vertices requested are indeed stored in this pm
                } else if ( specials.find(pair.first) == specials.end() ) {
                    std::cerr << "Tried to write a not allowed Vertex into ProbabilityMatrix." << std::endl;
                    exit(1);
                }
            }
            // only write if nos is != 0 as we have a sparse implementation
            if (nos != 0) {
                pm[pk] = nos;
            }
        }
        
        unsigned long long ProbabilityMatrix::mnos() {
            unsigned long long mnos = 0;
            for (auto elem : pm) {
                mnos += elem.second;
            }
            return mnos;
        }
        
        
        ProbabilityMatrix::~ProbabilityMatrix() {
            // cleanup?
        }
        
        ProbabilityMatrix operator* (ProbabilityMatrix& x, ProbabilityMatrix& y) {
            ProbabilityMatrix z;
            // get all vertices that are present in both PMs
            std::set< int > intersection;
            std::set< int > xSpecial;
            std::set< int > ySpecial;
            std::insert_iterator< std::set<int> > insert_it (intersection, intersection.begin());
            std::set_intersection (xSpecial.begin(), xSpecial.end(), ySpecial.begin(), ySpecial.end(), insert_it);
            
            // remove special vertices that are present in intersection
            
            // build all keys needed for new PM
            
            
            return z;
        }
        
        std::vector<ProbabilityKey> permute_key(ProbabilityKey pk) {
            std::vector<ProbabilityKey> result;
            
            ProbabilityKey current;
            permute_impl(pk.begin(), pk.end(), result, current);
            return result;
        }
        
        void permute_impl(ProbabilityKey::iterator start, ProbabilityKey::iterator end, std::vector<ProbabilityKey>& result, ProbabilityKey current) {
            if (start == end) {
                    result.push_back(current);
            } else {
                for (auto b : base_conversion[start->second]) {
                    //std::cerr << "Base " << enum_to_char(b) << " start " << start->first << " current " << current << std::endl;
                    current[start->first] = b;
                    permute_impl(++start, end, result, current);
                    start--;
                }
            }
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