/* This program reads secundary RNA structures in dot-bracket and
 * builds a graph for a latter ear-decomposition and bipartitness-check
 *
 * Created on: 13.08.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

#ifndef PATHCOLORING_H
#define PATHCOLORING_H

    // include common header with graph definition and global variables
    #include "common.h"

    // include standard library parts
    #include <sstream>

    // include boost components

namespace design {
  namespace detail {
    // typedefs
    typedef matrix< unsigned long long, A_Size, A_Size > rnaMatrix;

    // class definitions
    // Class to get Fibonacci numbers

    class Fibonacci {
    public:
      Fibonacci (unsigned int l);

      unsigned int get (unsigned int n) {
        return numbers[n - 1];
      };
    private:
      std::vector< unsigned int > numbers;
    };

    // Class to get Pairing numbers

    class Pairing {
    public:
      Pairing (unsigned int length);
      unsigned long long get (unsigned int l, unsigned int b1, unsigned int b2);
      unsigned long long get (unsigned int l, unsigned int b1);
      unsigned long long get (unsigned int l);
      //matrix< unsigned int, 4, 4 > getMatrix(unsigned int l) { return p[l]; };
    private:
      std::vector< rnaMatrix > p;
      rnaMatrix multiply (rnaMatrix a, rnaMatrix b);
      unsigned int length;
    };

    // fills the string sequence with random bases, given the first and the last base and the intended length; returns the number of possible solutions
    unsigned long long generate_path_seq (std::deque< int >& sequence, int first, int last, int length);

    // same for cycles (just a wrapper)
    unsigned long long generate_cycle_seq (std::deque< int >& sequence, int first, int length);

    // function that takes a path or cycle (in form of a graph), meassures its length, first and last base assignment, 
    // calls generate_path/cycle_seq, and assignes all the colors to the Graph
    unsigned long long color_path_cycle_graph (Graph& g);

    // this function takes a sequence of bases and assigns them to a graph (must be a path or a cycle!)
    void sequencestring_to_graph (Graph& g, Vertex vertex, std::deque< int >& sequence);
  }
}
#endif
