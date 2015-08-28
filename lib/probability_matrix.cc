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
                    throw( new std::logic_error( "Tried to get a not allowed Vertex from ProbabilityMatrix." ));
                }
            }
            
            unsigned long long returnvalue = 0;
            std::vector<ProbabilityKey> allkeys = permute_key(pk);
            
            for (auto key : allkeys) {
                // important for map: if you request with [] an entry will be created for unexisting ones.
                ProbabilityMap::const_iterator found = pmap.find(key);

                if (found != pmap.end()) {
                    returnvalue += found->second;
                }
            }
            return returnvalue;
        }
        
        void ProbabilityMatrix::put (ProbabilityKey& pk, unsigned long long nos) {
            if (!initialized) {
                for (auto pair : pk) {
                    specials.insert(pair.first);
                }
                initialized = true;
            }
            
            for (auto pair : pk) {
                // sanity check if the bases requested are within our alphabet size
                if (pair.second >= A_Size) {
                    throw new std::out_of_range( "Tried to write a base outside of the alphabet size into ProbabilityMatrix." );
                // sanity check if the vertices requested are indeed stored in this pm
                } else if ( specials.find(pair.first) == specials.end() ) {
                    throw new std::logic_error( "Tried to write a not allowed Vertex into ProbabilityMatrix." );
                }
            }
            // only write if nos is != 0 as we have a sparse implementation
            if (nos != 0) {
                pmap[pk] = nos;
            }
        }
        
        unsigned long long ProbabilityMatrix::mnos() {
            unsigned long long mnos = 0;
            for (auto elem : pmap) {
                mnos += elem.second;
            }
            return mnos;
        }
        
        template <typename R>
        ProbabilityKey ProbabilityMatrix::sample(R* rand_ptr) {
            ProbabilityKey pk;
            
            for (auto s : getSpecials()) {
                pk[s] = N; //TODO maybe we could check sequence constraints here?
            }
            
            return sample(pk, rand_ptr);
        }
        
        template <typename R>
        ProbabilityKey ProbabilityMatrix::sample(ProbabilityKey pk, R* rand_ptr) {
            ProbabilityKey result;
            
            // get all possible keys for the constraints set in pk
            std::vector<ProbabilityKey> possible_keys = permute_key(pk);
            unsigned long long constrained_mnos = 0;
            // get the maximal number of sequences for the input constraints set in pk
            for (auto k: possible_keys) {
                constrained_mnos += (*this)[k];
            }
            
            if (constrained_mnos == 0) {
                throw std::logic_error( "Cannot fulfill constraints while sampling a key!" );
            }
            std::uniform_real_distribution<float> dist(0, 1);
            unsigned long long random = dist(*rand_ptr) * constrained_mnos;
            // stochastically take one of the possibilities
            // start at the probability of first possible character and add each other base probability as long long as the random number is bigger.
            unsigned long long sum = 0;
            for (auto k : possible_keys) {
                sum += (*this)[k];
                // if the random number is bigger than our probability, take this base as the current base!
                if (random < sum) {
                    result = k;
                    // don't forget to exit the loop, otherwise will always be last entry
                    break;
                }
            }

            if (debug) {
                std::cerr << "Key Sampled: " << result << std::endl;
            }
            return result;
        }
        
        ProbabilityMatrix operator* (ProbabilityMatrix& x, ProbabilityMatrix& y) {
            if (!x.is_initialized()) {
                return y;
            } else if (!y.is_initialized()) {
                return x;
            } else {
                ProbabilityMatrix z;
                // get all vertices that are present in both PMs
                std::set< int > xSpecials = x.getSpecials();
                std::set< int > ySpecials = y.getSpecials();
                //std::cerr << "xSpecials: " << xSpecials << std::endl;
                //std::cerr << "ySpecials: " << ySpecials << std::endl;
                std::set< int > zSpecials;
                std::insert_iterator< std::set<int> > insert_it (zSpecials, zSpecials.begin());
                std::set_union(xSpecials.begin(), xSpecials.end(), ySpecials.begin(), ySpecials.end(), insert_it);
                //std::cerr << "zSpecials: " << zSpecials << std::endl;

                // build all keys needed for new PM
                ProbabilityKey newkey;
                for (auto s : zSpecials) {
                    newkey[s] = N; //TODO maybe we could check sequence constraints here?
                }
                std::vector<ProbabilityKey> zkeys = permute_key(newkey);

                //std::cerr << "zkeys: " << zkeys << std::endl;

                // lookup keys from previous and multiply them to insert into new
                for (auto zkey : zkeys) {
                    // generate keys for both old PMs
                    ProbabilityKey xkey;
                    for (auto s : xSpecials) {
                        xkey[s] = zkey[s];
                    }
                    ProbabilityKey ykey;
                    for (auto s : ySpecials) {
                        ykey[s] = zkey[s];
                    }
                    // read probability for this keys and multiply them
                    // insert this new probability into the new pm z
                    z.put(zkey, x[xkey]*y[ykey]);

                    //std::cerr << "x: " << xkey << ": " << x[xkey] << std::endl;
                    //std::cerr << "y: " << ykey << ": " << y[ykey] << std::endl;
                    //std::cerr << "z: " << zkey << ": " << z[zkey] << std::endl;
                }
                return z;
            }
        }
        
        ProbabilityMatrix make_internal(ProbabilityMatrix& pm, int v) {
            ProbabilityMatrix result;
            
            std::set< int > specials = pm.getSpecials();
            // find v in specials
            std::set< int >::iterator v_it = specials.find(v);
            if (v_it != specials.end()) {
                specials.erase(v_it);
                // create a key for the new pm
                ProbabilityKey newkey;
                for (auto s : specials) {
                    newkey[s] = N; //TODO maybe we could check sequence constraints here?
                }
                
                std::vector<ProbabilityKey> newkeys = permute_key(newkey);
                for (auto k : newkeys) {
                    ProbabilityKey pmkey = k;
                    pmkey[v] = N;
                    // now put the new key with the sum over the internal vertex into the new matrix
                    result.put(k, pm[pmkey]);
                }
            } else {
                // in case the vertex is not present in specials, it is already "internal"
                result = pm;
            }
            return result;
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
        
        // overload << operator to print ProbabilityKeys
        std::ostream& operator<<(std::ostream& os, ProbabilityKey& m) {
            os << "[";
            for (typename ProbabilityKey::iterator it = m.begin(); it != m.end(); it++) {
                os << "(" << std::setfill(' ') << std::setw(1) << it->first << "," << enum_to_char(it->second) << ")";
            }
            os << "]";
            return os;
        }
        
        // overload << operator to print ProbabilityMatrix
        std::ostream& operator<<(std::ostream& os, ProbabilityMatrix& m) {
            ProbabilityKey key;
            std::set<int> specials = m.getSpecials();
            for (auto s : specials) {
                key[s] = N;
            }
            std::vector<ProbabilityKey> keys = permute_key(key);
            for (auto k : keys) {
                if (m[k] != 0) {
                    os << k << ": " << m[k] << std::endl;
                }
            }
            return os;
        }
        
        template ProbabilityKey ProbabilityMatrix::sample<std::mt19937> (std::mt19937*);
        template ProbabilityKey ProbabilityMatrix::sample<std::mt19937> (ProbabilityKey, std::mt19937*);
    }
}