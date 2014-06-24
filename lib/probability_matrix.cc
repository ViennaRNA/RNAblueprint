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

std::size_t MyKeyHash::operator() (const MyKey& k) const {
  // Start with 0 as a hash value   .
  std::size_t hash = 0;

  for (auto elem : k) {
    boost::hash_combine(hash, boost::hash_value(elem.first));
    boost::hash_combine(hash, boost::hash_value(elem.second));
  }

  //std::cerr << "hash is: " << hash << std::endl;

  return hash;
}

ProbabilityMatrix::ProbabilityMatrix (Graph& g) {

  if (debug) {
    std::cerr << "Initializing ProbabilityMatrix..." << std::endl;
  }

  // structure to remember current Ak (attachment vertices)
  std::set<Vertex> currentAk;
  // to store current inner Articulation Points
  std::set<Vertex> currentAi;
  // now start at the outermost ear
  unsigned int k = 1;

  int max_length = 0;
  Graph::children_iterator ear, ear_end;
  for (boost::tie(ear, ear_end) = (g).children(); ear != ear_end; ++ear) {
    int length = boost::num_edges(*ear);
    if (max_length < length)
      max_length = length;
    // also get my (number of ears)
    my++;
  }

  if (debug) {
    std::cerr << "My is: " << my << std::endl;
  }

  // get Pairing matrix for paths, TODO only initialize once for the whole program!
  p = new Pairing(max_length + 1);

  // start at the outermost ear and process inwards
  for (boost::tie(ear, ear_end) = (g).children(); ear != ear_end; ++ear) {

    // before doing anything, update current Ak and Ai
    updateCurrentAkAi(*ear, k, currentAk, currentAi);

    // best is to remember the currentAis at this point
    Aks.push_back(currentAk);

    if (debug) {
      std::cerr << "===========================================" << std::endl
          << "Current k: " << k << std::endl
          << "currentAk:" << std::endl << currentAk << std::endl
          << "currentAi:" << std::endl << currentAi << std::endl;
    }

    // Nk[A6][A10][A1] = sum(AUGC in inner Ap = 9) P[A6][x9][3 pathlength] * P[x9][A10][1] * Nk-1 [x9][A1]

    std::vector<MyKey> key_combinations; // this is what we want to fill in the recursion
    MyKey mykey; // helper to recursively build the posibilities
    std::set<int> cAk; // create a version with ints for recursion
    for (auto ap : currentAk) {
      cAk.insert(boost::get(boost::vertex_color_t(), g.root(), ap));
    }

    // now fill the key_combinations with all kinds of bases
    calculate_combinations(cAk, mykey, key_combinations);

    /*if (debug) {
      std::cerr << "Calculated following combinations of keys to calculate:" << std::endl;
      for (auto thiskey : key_combinations) {
        std::cerr << thiskey << std::endl;
      }
    }*/

    // calculate the probabilities for every key and if not zero, add to matrix
    for (auto thiskey : key_combinations) {
      if (debug) {
        std::cerr << "-------------------------------------------" << std::endl
            << "Calculating probability for key: " << thiskey << std::endl;
      }

      unsigned long long probability = get_probability(thiskey, *ear, currentAk, currentAi, k);

      if (debug) {
        std::cerr << "= " << probability << std::endl;
      }

      if (probability != 0) {
        n[thiskey] = probability;
      }

      // remember max number of sequences/solutions
      if (k == my) {
        nos += probability;
      }
    }


    if (debug) {
      for (auto thiskey : key_combinations) {
        unsigned long long value = get(thiskey);
        if (value != 0) {
          std::cerr << thiskey << " = " << value << std::endl;
        }
      }
    }
    // now going to next ear!
    k++;
  }

  if (debug) {
    std::cerr << "Sub probabilities remembered are: " << std::endl;
    for (auto part : parts) {
      for (auto sub_probability : part) {
        std::cerr << "(" << sub_probability.start << ") --"
            << sub_probability.length << "-- ("
            << sub_probability.end << "); ";
      }
      std::cerr << std::endl;
    }
    std::cerr << "Aks remembered are: " << std::endl;
    for (auto Ak : Aks) {
      std::cerr << Ak << std::endl;
    }
  }
}

ProbabilityMatrix::~ProbabilityMatrix () {
  delete p;
}

void ProbabilityMatrix::updateCurrentAkAi (Graph& g, int k, std::set<Vertex>& currentAk, std::set<Vertex>& currentAi) {
  if (debug) {
    std::cerr << "Updating Ak and Ik for ear " << k << std::endl;
  }
  // Ak and Ai are already stored in graph as a vertex properties
  // write vertex property into currentAk and currentAi
  currentAi.clear();

  BGL_FORALL_VERTICES_T(v, g, Graph) {
    if (g[v].Ak.find(k) != g[v].Ak.end()) {
      currentAk.insert(g.local_to_global(v));
    } else if (g[v].Ai == k) {
      currentAi.insert(g.local_to_global(v));
      // we need to keep Ak from previous glued ears, except those that became Ai this time!
      currentAk.erase(g.local_to_global(v));
    }
  }

  // for the last cycle, Ak is zero, so we need to "define" one Ai as Ak!
  if (currentAk.size() == 0) {
    auto first = currentAi.begin();
    // push the first Ai to Ak
    currentAk.insert(*first);
    // erase it in Ai
    currentAi.erase(first);
  }
}

unsigned long long ProbabilityMatrix::get_probability (MyKey mykey, Graph& g, std::set<Vertex>& Ak, std::set<Vertex>& Ai, unsigned int k) {

  // v--return this--v v--recursion get_sum_of_sum--v v--let's call them sub_probabilities--v   v--prob. last_ear--v
  // Nk[A6][A10][A1]  =   sum(AUGC in inner Ap = 9)   P[A6][x9][3 pathlength] * P[x9][A10][1] * Nk-1 [x9][A1]

  // this is what we calculate here
  unsigned long long max_number_of_sequences = 0;
  // key to get probability from last ear glueing and to store all base combinations for internal articulation points
  MyKey lastkey;
  // container to store all sub-probabilities
  std::vector<SubProbability> sub_probabilities;


  // start to build lastkey
  for (auto v : Ak) {
    bool into_lastear = false;
    if (!g.find_vertex(v).second) {
      // are not in this particular ear
      into_lastear = true;
    } else {
      // check if were Aks before
      Vertex lv = g.global_to_local(v);
      for (unsigned int i = 0; i < k; i++) {
        if (g[lv].Ak.find(i) != g[lv].Ak.end()) {
          into_lastear = true;
        }
      }
    }

    if (into_lastear) {
      // these articulation points are not in this particular ear or were Aks before,
      // therefore we need to look their probability up
      // from last time. so lets generate a key therefore			
      lastkey.insert(std::make_pair(boost::get(boost::vertex_color_t(), g.root(), v), mykey[boost::get(boost::vertex_color_t(), g.root(), v)]));
    }
  }

  // now building this sub_probabilities
  Graph::children_iterator part, part_end;
  for (boost::tie(part, part_end) = (g).children(); part != part_end; ++part) {

    // add to sub_probabilities for this start/end and length
    sub_probabilities.push_back(SubProbability());
    SubProbability &sub_probability = sub_probabilities.back();
    // begin, end and length go into sub_probability
    int i = 0;

    BGL_FORALL_VERTICES_T(v, *part, Graph) {
      if (boost::degree(v, *part) == 1) {
        if (i == 0)
          sub_probability.start = boost::get(boost::vertex_color_t(), *part, v);
        else
          sub_probability.end = boost::get(boost::vertex_color_t(), *part, v);
        i++;
      }
    }
    sub_probability.length = boost::num_edges(*part);
  }
  // remember for later when we have no graph (at backtracing)
  if (parts.size() == k - 1) parts.push_back(sub_probabilities);

  if (debug) {
    for (auto sub_probability : sub_probabilities) {
      std::cerr << "(" << sub_probability.start << ") --"
          << sub_probability.length << "-- ("
          << sub_probability.end << ")\t*\t";
    }
    std::cerr << "N(" << lastkey << ")" << std::endl;
  }

  // calculate sum of sum for all bases colored X (= internal aps)
  // adds base combinatoric to the sub_probabilities and to the lastkey
  if (Ai.size() == 0) {
    max_number_of_sequences += calculate_probability(mykey, lastkey, sub_probabilities);
  } else {
    if (debug) {
      std::cerr << "Make sum of sum: " << k << std::endl;
    }
    make_sum_of_sum(g, Ai, mykey, lastkey, sub_probabilities, k, max_number_of_sequences);
  }

  return max_number_of_sequences;
}

void ProbabilityMatrix::make_sum_of_sum (Graph& g,
    std::set<Vertex>& Ai,
    MyKey& mykey, MyKey& lastkey,
    std::vector<SubProbability>& sub_probabilities,
    unsigned int k,
    unsigned long long& max_number_of_sequences) {

  if (Ai.size() > 0) {
    std::set<Vertex>::iterator it = Ai.begin();
    Vertex v = *it;
    Ai.erase(it);

    for (unsigned int b = 0; b < A_Size; b++) {
      lastkey.insert(std::make_pair(boost::get(boost::vertex_color_t(), g.root(), v), b));
      // recursion starts here
      make_sum_of_sum(g, Ai, mykey, lastkey, sub_probabilities, k, max_number_of_sequences);

      if (Ai.size() == 0) {
        max_number_of_sequences += calculate_probability(mykey, lastkey, sub_probabilities);
      }
      // remove current vertex again to make space for a new base
      lastkey.erase(v);
    }
    // add current vertex again
    Ai.insert(v);
  }
}

int ProbabilityMatrix::get_color_from_key (MyKey& mykey, MyKey& lastkey, int vertex) {
  if (mykey.find(vertex) != mykey.end()) {
    return mykey[vertex];
  } else if (lastkey.find(vertex) != lastkey.end()) {
    return lastkey[vertex];
  } else {
    std::cerr << "Something went wrong! I could not find the color of the vertex in any key (neither last nor current): "
        << vertex << std::endl;
    exit(1);
  }
}

unsigned long long ProbabilityMatrix::calculate_probability (MyKey& mykey, MyKey& lastkey,
    std::vector<SubProbability>& sub_probabilities) {
  // do the actual calculation here!
  unsigned long long multiplied_probabilities = 1;

  // calculate product of sub probabilities
  for (auto sub_probability : sub_probabilities) {
    int startBase = get_color_from_key(mykey, lastkey, sub_probability.start);
    int endBase = get_color_from_key(mykey, lastkey, sub_probability.end);
    // get probability and multiply it
    multiplied_probabilities *= p->get(sub_probability.length, startBase, endBase);

    if (debug) {
      std::cerr << "P(" << enum_to_char(startBase) << ", "
          << enum_to_char(endBase) << ", " << sub_probability.length << ") = "
          << p->get(sub_probability.length, startBase, endBase) << "\t*\t";
    }
  }
  // for k = 0 we should not look up lastkey, as there is no lastkey!
  if (lastkey.size() != 0) {
    // now add probability for last_ear
    unsigned long long last_probability = get(lastkey);
    if (debug) {
      std::cerr << "N{" << lastkey << "} = "
          << last_probability;
    }

    multiplied_probabilities *= last_probability;
  }

  if (debug) {
    std::cerr << "\t=\t" << multiplied_probabilities << std::endl;
  }

  return multiplied_probabilities;
}

void ProbabilityMatrix::calculate_combinations (std::set<int>& Ak, MyKey& mykey, std::vector<MyKey>& key_combinations) {

  if (Ak.size() > 0) {
    std::set<int>::iterator it = Ak.begin();
    int vertex = *it;
    Ak.erase(it);

    for (unsigned int b = 0; b < A_Size; b++) {
      mykey.insert(std::make_pair(vertex, b));
      // recursion starts here
      calculate_combinations(Ak, mykey, key_combinations);

      if (Ak.size() == 0) {
        // remember our generated key
        key_combinations.push_back(mykey);
      }
      // remove current vertex again to make space for a new base
      mykey.erase(vertex);
    }
    // add current vertex again
    Ak.insert(vertex);
  }
}

unsigned long long ProbabilityMatrix::get (MyKey mykey) {

  unsigned long long returnvalue;

  // important for map: if you request with [] an entry will be created for unexisting ones.
  std::unordered_map < MyKey, unsigned long long, MyKeyHash>::const_iterator found = n.find(mykey);

  if (found != n.end()) {
    returnvalue = found->second;
  } else {
    returnvalue = 0;
  }

  return returnvalue;
}

unsigned long long ProbabilityMatrix::get_sum (int k, MyKey mykey, MyKey lastkey, std::vector<MyKey>& key_combinations, std::unordered_map < MyKey, unsigned long long, MyKeyHash>& probabilities) {
  MyKey tempkey; // helper to build all combinations
  std::set<int> cAk; // to get all combinations of keys we need a list of articulation points
  for (auto elem : mykey) {
    if (elem.second == X) {
      cAk.insert(elem.first); // store to get all key combinations
      mykey.erase(elem.first); // delete to in the end get a list of static key-elements
    }
  }
  // now lets calculate all combinations of keys for our points
  if (cAk.size() != 0) {
    calculate_combinations(cAk, tempkey, key_combinations);
  } else {
    key_combinations.push_back(MyKey());
  }

  // need to add the static key-elements in mykey to our keys in the key_combinations
  for (std::vector<MyKey>::iterator it = key_combinations.begin(); it != key_combinations.end(); ++it) {
    (*it).insert(mykey.begin(), mykey.end());
    //std::cerr << *it << " = " << get(*it) << std::endl;
  }

  // now get the sum of all probabililties
  unsigned long long sum = 0;
  for (auto thiskey : key_combinations) {
    // this was sum += get(thiskey) before however we need to take the path probabilities between the Aks into account and
    // multiply them. therefore we need the calculate_probabaility function, which wants this stupid variables
    // at the moment I still need to find a way to create a lastkey where all already colored bases are inside.
    if (k == (int) my) {
      probabilities[thiskey] = get(thiskey);
      sum += get(thiskey);
    } else {
      if (debug) {
        std::cerr << thiskey << std::endl;
      }
      unsigned long long thisprob = calculate_probability(lastkey, thiskey, parts[k]);
      probabilities[thiskey] = thisprob;
      sum += thisprob;
    }
  }
  return sum;
}

std::set<Vertex> ProbabilityMatrix::get_Ak (unsigned int k) {
  if (k <= my) {
    return Aks[k - 1];
  } else {
    std::cerr << "k can't be bigger than my!" << std::endl;
    exit(1);
  }
}

// overload << operator to print mykeys in pair with graph

std::ostream& operator<< (std::ostream& os, MyKey& m) {
  os << "[";
  for (typename MyKey::iterator it = m.begin(); it != m.end(); it++) {
    os << "(" << std::setfill(' ') << std::setw(1) << it->first << "," << enum_to_char(it->second) << ")";
  }
  os << "]";
  return os;
}

// overload << operator to print set of vertices

std::ostream& operator<< (std::ostream& os, std::set<Vertex>& vs) {
  for (auto v : vs) {
    os << v << " ";
  }
  os << std::endl;
  return os;
}