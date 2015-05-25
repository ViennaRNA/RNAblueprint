/* This file is a boost.test unit test and provides tests for pathcoloring.cc
 *
 *
 * Created on: 23.10.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 * 
 *
 */

// include header
#include "test_common.h"

// include headers containing functions to test
#include "pathcoloring.h"

using namespace design;
using namespace design::detail;

// define heads
namespace PathColoring {

  class TestCase {
  public:
    TestCase (int first, int last, int length, int nos, Sequence sequence);
    int first;
    int last;
    int length;
    int nos;
    Sequence sequence;
  };

  Sequence get_vertex_colors (Graph& g);
  void reset (Graph& g);
}

BOOST_AUTO_TEST_SUITE (PathColoring)

TestCase::TestCase (int first, int last, int length, int nos, Sequence sequence) {
  this->first = first;
  this->last = last;
  this->length = length;
  this->nos = nos;
  this->sequence = sequence;
}

BOOST_AUTO_TEST_CASE (generatePathSeq) {
  // set random generator to a static seed;

  initialize_library(true);

  std::vector < TestCase > testcases;
  testcases.push_back(TestCase(G, G, 2, 2, {G, C, G}));
  testcases.push_back(TestCase(A, A, 0, 1, {A}));
  testcases.push_back(TestCase(Y, Y, 0, 2, {C}));
  testcases.push_back(TestCase(N, G, 3, 5, {U, G, U, G}));
  testcases.push_back(TestCase(U, U, 6, 13, {U, G, U, G, U, A, U}));
  testcases.push_back(TestCase(R, U, 3, 5, {G, U, G, U}));
  testcases.push_back(TestCase(Y, U, 6, 21, {U, G, U, G, C, G, U}));
  testcases.push_back(TestCase(C, N, 4, 5, {C, G, U, G, U}));
  testcases.push_back(TestCase(N, N, 1, 6, {G, U}));
  testcases.push_back(TestCase(N, N, 5, 42, {G, U, G, U, A, U}));
  testcases.push_back(TestCase(-1, A, 6, 5, {A, U, G, U, G, U}));
  testcases.push_back(TestCase(-1, G, 6, 13, {G, U, G, U, G, C}));
  testcases.push_back(TestCase(-1, N, 4, 14, {G, U, G, U}));
  testcases.push_back(TestCase(C, R, 1, 1, {C, G}));
  testcases.push_back(TestCase(R, Y, 1, 3, {G, U}));
  testcases.push_back(TestCase(H, A, 2, 1, {A, U, A}));
  testcases.push_back(TestCase(-1, N, 4, 14, {G, U, G, U}));
  testcases.push_back(TestCase(-1, A, 4, 2, {A, U, G, U}));
  testcases.push_back(TestCase(-1, U, 4, 5, {U, G, U, G}));
  testcases.push_back(TestCase(-1, G, 4, 5, {G, U, G, U}));
  testcases.push_back(TestCase(-1, C, 4, 2, {C, G, U, G}));
  testcases.push_back(TestCase(-1, R, 4, 7, {G, U, G, U}));
  testcases.push_back(TestCase(-1, V, 4, 9, {C, G, U, G}));
  
  for (auto t : testcases) {
    // Variables needed for this function
    int nos;
    Sequence sequence;
    // tell the user what is going on right now
    std::stringstream ss;
    ss << "Checking Sequence (" << enum_to_char(t.first) << ", "
        << enum_to_char(t.last) << ", " << t.length << ")";
    
    std::mt19937 num_gen(1);
    if (t.first >= 0) {
      nos = generate_path_seq(sequence, t.first, t.last, t.length, &num_gen);
    } else {
      nos = generate_cycle_seq(sequence, t.last, t.length, &num_gen);
    }
    // do the comparison of the results
    ss << std::endl << "Got Sequence: " << sequence << " with NOS: " << nos;
    BOOST_TEST_MESSAGE(ss.str());
    BOOST_CHECK(sequence == t.sequence);
    BOOST_CHECK(nos == t.nos);
  }
}

Sequence get_vertex_colors (Graph& g) {
  Sequence sequence;

  BGL_FORALL_VERTICES_T(v, g, Graph) {
    std::cerr << v << ":" << enum_to_char(g[v].base) << std::endl;
    sequence.push_back(g[v].base);
  }
  return sequence;
}

void reset (Graph& g) {
  // reset color

  BGL_FORALL_VERTICES_T(v, g, Graph) {
    g[v].base = N;
  }
}

BOOST_AUTO_TEST_CASE (colorPathGraph) {
  // set random generator to a static seed;
  initialize_library(true);
  // create random generator
  std::mt19937 rand_gen(1);
  // create a graph
  Graph g(10);
  int vertex_name = 0;

  BGL_FORALL_VERTICES_T(v, g, Graph) {
    boost::put(boost::vertex_color_t(), g, v, vertex_name++);
  }
  
  for (unsigned int i = 0; i < boost::num_vertices(g) - 1; i++) {
    boost::add_edge(boost::vertex(i, g), boost::vertex(i + 1, g), g);
  }

  // color this graph!
  BOOST_TEST_MESSAGE("color path-graph");
  rand_gen.seed(1);
  color_path_cycle_graph(g, &rand_gen);
  //print_graph(g, out, "path");
  Sequence sequence{G ,U ,G ,U ,A ,U ,A ,U ,A ,U};
  BOOST_CHECK(sequence == get_vertex_colors(g));

  // color first base and try all over again
  BOOST_TEST_MESSAGE("path_ends_U");
  reset(g);
  g[boost::vertex(boost::num_vertices(g)-1, g)].base = U;
  rand_gen.seed(1);
  color_path_cycle_graph(g, &rand_gen);
  Sequence sequence2{G ,U ,G ,U ,A ,U ,A ,U ,A ,U};
  BOOST_CHECK(sequence2 == get_vertex_colors(g));
  
  // color first base and try all over again
  BOOST_TEST_MESSAGE("path_starts_A");
  reset(g);
  g[boost::vertex(0, g)].base = A;
  rand_gen.seed(1);
  color_path_cycle_graph(g, &rand_gen);
  Sequence sequence1{A ,U ,G ,U ,G ,C ,G ,C ,G ,C};
  BOOST_CHECK(sequence1 == get_vertex_colors(g));

  // color both ends and try all over again
  BOOST_TEST_MESSAGE("path_ends_AC");
  reset(g);
  g[boost::vertex(0, g)].base = A;
  g[boost::vertex(boost::num_vertices(g)-1, g)].base = C;
  rand_gen.seed(1);
  color_path_cycle_graph(g, &rand_gen);
  Sequence sequence3{A ,U ,G ,U ,G ,C ,G ,C ,G ,C};
  BOOST_CHECK(sequence3 == get_vertex_colors(g));

  // make a cycle and color this cycle!
  BOOST_TEST_MESSAGE("cycle");
  reset(g);
  boost::add_edge(boost::vertex(0, g), boost::vertex(boost::num_vertices(g) - 1, g), g);
  rand_gen.seed(1);
  color_path_cycle_graph(g, &rand_gen);
  Sequence sequence4{G ,U ,G ,U ,A ,U ,A ,U ,A ,U};
  BOOST_CHECK(sequence4 == get_vertex_colors(g));

  // set one base and color again this cycle
  BOOST_TEST_MESSAGE("cycle_starts_C");
  reset(g);
  g[boost::vertex(5, g)].base = C;
  rand_gen.seed(1);
  color_path_cycle_graph(g, &rand_gen);
  Sequence sequence5{A ,U ,G ,U ,G ,C ,G ,U ,A ,U};
  BOOST_CHECK(sequence5 == get_vertex_colors(g));
}

BOOST_AUTO_TEST_SUITE_END ()
