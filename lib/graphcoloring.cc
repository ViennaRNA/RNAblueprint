/* RNAdesign
 * A program for designing RNA molecules.
 *
 * Created on: 13.08.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

// include header
#include "graphcoloring.h"

namespace design {
  namespace detail {
/* 
    template <typename RG>
    void color_blocks (Graph& g, ProbabilityMatrix& pm, RG* rand_ptr) {
      // number of sequences to debug
      if (debug) {
        std::cerr << "Number of sequences for this block: " << pm.number_of_sequences() << std::endl;
      }
      // remember the current key for the next ear.
      ProbabilityKey lastkey;

      // reverse iterate again over all ears to color Aks and all vertices in between
      for (int k = pm.get_my(); k > 0; k--) {
        if (debug) {
          std::cerr << "Start Backtracing at ear " << k << std::endl;
        }
        // get the current Articulation Points
        std::set<Vertex> Ak = pm.get_Ak(k);

        // translate the Ak set into set of ints (no vertex descriptors!)
        ProbabilityKey thiskey;
        for (auto ap : Ak) {
          int vertex = boost::get(boost::vertex_color_t(), g.root(), ap);
          int color = (g.root())[ap].base;
          thiskey.insert(std::make_pair(vertex, color));
        }
        // now do the random coloring of our points
        if (debug) {
          std::cerr << "Try to color this key: " << thiskey << std::endl;
        }
        ProbabilityKey colorkey = color_articulation_points(k, pm, thiskey, lastkey, rand_ptr);
        // remember colorkey for next ear iteration
        lastkey = colorkey;
        if (debug) {
          std::cerr << "Got a colored key: " << colorkey << std::endl;
        }

        // put colors onto graph
        for (auto v : Ak) {
          Vertex lv = g.global_to_local(v);
          if (g[lv].base == N) {
            g[lv].base = colorkey[boost::get(boost::vertex_color_t(), g, lv)];
          } else if (g[lv].base != colorkey[boost::get(boost::vertex_color_t(), g, lv)]) {
            std::cerr << "ERROR: Tried to change a already assigned base of an articulation point." << std::endl;
            exit(1);
          }
          if (debug) {
            std::cerr << "v " << v << ": " << enum_to_char(g[lv].base) << std::endl;
          }
        }
      }

      // now let's color all the vertices on the parts children
      Graph::children_iterator ear, ear_end;
      for (boost::tie(ear, ear_end) = g.children(); ear != ear_end; ++ear) {
        Graph::children_iterator part, part_end;
        for (boost::tie(part, part_end) = (*ear).children(); part != part_end; ++part) {

          print_all_vertex_names(*part, "Coloring a part of the ear");
          color_path_cycle_graph(*part, rand_ptr);
        }
      }
    }
    
    template <typename RG>
    ProbabilityKey color_articulation_points (int k, ProbabilityMatrix& pm, ProbabilityKey& colorkey, ProbabilityKey& lastkey, RG* rand_ptr) {

      ProbabilityKey returnkey;

      // declare random number distribution and get a random number
      std::uniform_real_distribution<float> dist(0, 1);
      // get a random number between 0 and 1.
      float random = dist(*rand_ptr);
      if (debug) {
        std::cerr << "Got a random number: " << random << std::endl;
      }

      // now get all combinations of keys and the sum of all possibilities
      std::vector<ProbabilityKey> key_combinations; // this is what we want to fill now
      ProbabilityMap probabilities;
      unsigned long long sum_of_possibilities = pm.get_sum(k, colorkey, lastkey, key_combinations, probabilities);
      if (debug) {
        std::cerr << "Sum of all possibilities is: " << sum_of_possibilities << std::endl;
      }

      // stochastically take one of the posibilities
      // start at the probability of first possible key and add each other key probability 
      // as long as the random number is bigger.
      unsigned long long sum = 0;
      for (auto thiskey : key_combinations) {
        sum += probabilities[thiskey];
        // if the random number is bigger than our probability, take this base as the first base!
        if (random * sum_of_possibilities < sum) {
          returnkey = thiskey;
          break;
        }
      }
      return returnkey;
    }
    
    template void color_blocks<std::mt19937> (Graph&, ProbabilityMatrix&, std::mt19937*);
    template ProbabilityKey color_articulation_points<std::mt19937> (int, ProbabilityMatrix&, ProbabilityKey&, ProbabilityKey&, std::mt19937*);
*/
  }
}