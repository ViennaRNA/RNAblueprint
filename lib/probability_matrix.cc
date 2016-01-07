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

        SolutionSizeType ProbabilityMatrix::operator[] (ProbabilityKey& pk) {
            // sanity check if the vertices requested are indeed stored in this pm
            for (auto pair : pk) {
                if ( specials.find(pair.first) == specials.end() ) {
                    throw( new std::logic_error( "Tried to get a not allowed Vertex from ProbabilityMatrix." ));
                }
            }
            
            SolutionSizeType returnvalue = 0;
            PermuteKeyFactory pkf(pk);
            
            while (true) {
                // important for map: if you request with [] an entry will be created for unexisting ones.
                ProbabilityMap::const_iterator found = pmap.find(*pkf.key());

                if (found != pmap.end()) {
                    returnvalue += found->second;
                }
                
                if (!pkf.next_permutation())
                    break;
            }
            return returnvalue;
        }
        
        void ProbabilityMatrix::put (ProbabilityKey& pk, SolutionSizeType nos) {
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
        
        SolutionSizeType ProbabilityMatrix::mnos() {
            SolutionSizeType mnos = 0;
            for (auto elem : pmap) {
                mnos += elem.second;
            }
            return mnos;
        }
        
        template <typename R>
        std::pair<ProbabilityKey, SolutionSizeType> ProbabilityMatrix::sample(R& rand) {
            ProbabilityKey pk;
            
            for (auto s : getSpecials()) {
                pk[s] = N; //TODO maybe we could check sequence constraints here?
            }
            
            return sample(pk, rand);
        }
        
        template <typename R>
        std::pair<ProbabilityKey, SolutionSizeType> ProbabilityMatrix::sample(ProbabilityKey pk, R& rand) {
            ProbabilityKey result;
            
            // get all possible keys for the constraints set in pk
            PermuteKeyFactory pkf(pk);
            SolutionSizeType constrained_mnos = 0;
            // get the maximal number of sequences for the input constraints set in pk
            while (true) {
                constrained_mnos += (*this)[*pkf.key()];
                if (!pkf.next_permutation())
                    break;
            }
            
            if (constrained_mnos == 0) {
                throw std::logic_error( "Cannot fulfill constraints while sampling a key!" );
            }
            RandomDistType dist(0, constrained_mnos);
            SolutionSizeType random = dist(rand);
            // stochastically take one of the possibilities
            // start at the probability of first possible character and add each other base probability as long long as the random number is bigger.
            SolutionSizeType sum = 0;
            pkf.reset();
            while (true) {
                sum += (*this)[*pkf.key()];
                // if the random number is bigger than our probability, take this base as the current base!
                if (random < sum) {
                    result = *pkf.key();
                    // don't forget to exit the loop, otherwise will always be last entry
                    break;
                }
                if (!pkf.next_permutation())
                    break;
            }

            if (debug) {
                std::cerr << "Key Sampled: " << result << " with mnos: " << constrained_mnos << std::endl;
            }
            return std::make_pair(result, constrained_mnos);
        }
        
        ProbabilityMatrix ProbabilityMatrix::operator* (ProbabilityMatrix& y) {
            if (!this->is_initialized()) {
                return y;
            } else if (!y.is_initialized()) {
                return *this;
            } else {
                ProbabilityMatrix z;
                // get all vertices that are present in both PMs
                std::set< int > xSpecials = this->getSpecials();
                std::set< int > ySpecials = y.getSpecials();
                //std::cerr << "xSpecials: " << xSpecials << std::endl;
                //std::cerr << "ySpecials: " << ySpecials << std::endl;
                std::set< int > zSpecials;
                std::insert_iterator< std::set<int> > insert_it (zSpecials, zSpecials.begin());
                std::set_union(xSpecials.begin(), xSpecials.end(), ySpecials.begin(), ySpecials.end(), insert_it);
                //std::cerr << "zSpecials: " << zSpecials << std::endl;
                
                ProbabilityMap::iterator pmap_it = this->pmap.begin();
                for (pmap_it; pmap_it != this->pmap.end(); ++pmap_it) {
                    // first get the current key and insert Ns to those vertices not present already
                    ProbabilityKey newkey = pmap_it->first;
                    for (auto s : zSpecials) {
                        if (newkey.find(s) == newkey.end())
                            newkey[s] = N;
                    }
                    
                    // now generate all combinations
                    PermuteKeyFactory pkf(newkey);
                    while (true) {
                        // now access the second value
                        ProbabilityKey ykey;
                        for (auto s : ySpecials) {
                            ykey[s] = (*pkf.key())[s];
                        }
                        // read probability for this keys and multiply them
                        // insert this new probability into the new pm z
                        z.put(*pkf.key(), pmap_it->second * y[ykey]);
                        if (!pkf.next_permutation())
                            break;
                    }
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

                PermuteKeyFactory pkf(newkey);
                while (true) {
                    ProbabilityKey pmkey(*pkf.key());
                    pmkey[v] = N;
                    // now put the new key with the sum over the internal vertex into the new matrix
                    result.put(*pkf.key(), pm[pmkey]);
                    if (!pkf.next_permutation())
                        break;
                }
            } else {
                // in case the vertex is not present in specials, it is already "internal"
                result = pm;
            }
            return result;
        }

        PermuteKeyFactory::PermuteKeyFactory(ProbabilityKey pk) {
            //std::map<int, std::list<int> > storage;
            //std::map<int, std::list<int>::iterator> state;
            
            // fill storage container with all possibilities
            for (auto k : pk) {
                for (int b : base_conversion[k.second]) {
                    storage[k.first].push_back(b);
                }
            }
            // fill state container
            reset();
            // fill current
            copy_current();
        }
        
        ProbabilityKey* PermuteKeyFactory::key() {
            return &current;
        }
        
        void PermuteKeyFactory::reset() {
            // reset to begin
            std::map<int, std::list<int> >::iterator s_it;
            for (s_it = storage.begin(); s_it != storage.end(); s_it++) {
                state[s_it->first] = s_it->second.begin();
            }
        }

        bool PermuteKeyFactory::next_permutation() {
            // move to next step
            bool result = make_next_step(state.begin());
            // update current
            copy_current();
            return result;
        }

        bool PermuteKeyFactory::make_next_step(std::map<int, std::list<int>::iterator>::iterator state_it) {
            if (state.empty())
                return false;
            else {
                state_it->second++;
                if (state_it->second == storage[state_it->first].end()) {
                    state_it->second = storage[state_it->first].begin();
                    state_it++;
                    if (state_it != state.end())
                        return make_next_step(state_it);
                    else
                        return false;
                } else {
                    return true;
                }
            }
        }
        
        bool PermuteKeyFactory::previous_permutation() {
            // move to previous step
            bool result = make_previous_step(state.begin());
            // update current
            copy_current();
            return result;
        }
        
        bool PermuteKeyFactory::make_previous_step(std::map<int, std::list<int>::iterator>::iterator state_it) {
            if (state.empty())
                return false;
            else {
                if (state_it->second == storage[state_it->first].begin()) {
                    state_it->second = storage[state_it->first].end();
                    state_it->second--;
                    state_it++;
                    if (state_it != state.end())
                        return make_previous_step(state_it);
                    else
                        return false;
                } else {
                    state_it->second--;
                    return true;
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
            PermuteKeyFactory pkf(key);
            while (true) {
                if (m[*pkf.key()] != 0) {
                    os << *pkf.key() << ": " << m[*pkf.key()] << std::endl;
                }
                if (!pkf.next_permutation())
                    break;
            }
            return os;
        }
        
        template std::pair<ProbabilityKey, SolutionSizeType> ProbabilityMatrix::sample<std::mt19937> (std::mt19937&);
        template std::pair<ProbabilityKey, SolutionSizeType> ProbabilityMatrix::sample<std::mt19937> (ProbabilityKey, std::mt19937&);
    }
}