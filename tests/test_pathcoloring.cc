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
namespace PathColoring{

    class TestCase
    {
        public :
        TestCase(int first, int last, int length, int nos, Sequence sequence);
        TestCase(Sequence input, int nos, Sequence sequence);
        int first;
        int last;
        int length;
        Sequence input;
        int nos;
        Sequence sequence;
    };

    Sequence get_vertex_colors(Graph & g);
    void reset(Graph & g);
}

BOOST_AUTO_TEST_SUITE(PathColoring)

TestCase::TestCase(int first, int last, int length, int nos, Sequence sequence) {
    this->first = first;
    this->last = last;
    this->length = length;
    this->nos = nos;
    this->sequence = sequence;
}

TestCase::TestCase(Sequence input, int nos, Sequence sequence) {
    this->input = input;
    this->nos = nos;
    this->sequence = sequence;
}

BOOST_AUTO_TEST_CASE(generateSeq) {
    // set random generator to a static seed;
    std::mt19937 rand_gen(1);

    initialize_library(true);

    std::vector < TestCase > testcases;
    testcases.push_back(TestCase(G, G, 2, 2,{G, U, G}));
    testcases.push_back(TestCase(A, A, 0, 1,{A}));
    testcases.push_back(TestCase(N, C, 0, 1,{C}));
    testcases.push_back(TestCase(Y, Y, 0, 2,{C}));
    testcases.push_back(TestCase(N, G, 3, 5,{U, G, U, G}));
    testcases.push_back(TestCase(U, U, 6, 13,{U, G, C, G, U, G, U}));
    testcases.push_back(TestCase(R, U, 3, 5,{G, U, G, U}));
    testcases.push_back(TestCase(Y, U, 6, 21,{C, G, C, G, U, G, U}));
    testcases.push_back(TestCase(C, N, 4, 5,{C, G, U, G, U}));
    testcases.push_back(TestCase(N, N, 1, 6,{U, G}));
    testcases.push_back(TestCase(N, N, 5, 42,{U, A, U, G, U, G}));
    testcases.push_back(TestCase(C, R, 1, 1,{C, G}));
    testcases.push_back(TestCase(R, Y, 1, 3,{G, U}));
    testcases.push_back(TestCase(H, A, 2, 1,{A, U, A}));

    for (auto t : testcases) {
        // tell the user what is going on right now
        std::stringstream ss;
        ss << "-> Checking Sequence (" << enum_to_char(t.first) << ", "
                << enum_to_char(t.last) << ", " << t.length << "):";
        BOOST_TEST_MESSAGE(ss.str());

        rand_gen.seed(1);

        // Build graph
        Graph g(t.length + 1);
        int vertex_name = 0;
        // set the vertex_name

        BGL_FORALL_VERTICES_T(v, g, Graph) {
            boost::put(boost::vertex_color_t(), g, v, vertex_name++);
        }
        // add the edges
        for (unsigned int i = 0; i < boost::num_vertices(g) - 1; i++) {
            boost::add_edge(boost::vertex(i, g), boost::vertex(i + 1, g), g);
        }

        // set the sequence constraints
        g[boost::vertex(0, g)].base = t.first;
        g[boost::vertex(boost::num_vertices(g) - 1, g)].base = t.last;

        //call the function
        boost::multiprecision::mpz_int nos = color_path_graph(g, &rand_gen);
        Sequence sequence = get_vertex_colors(g);
        // do the comparison of the results
        ss.str(std::string());
        ss << std::endl << "Got Sequence: " << sequence << " with NOS: " << nos;
        BOOST_TEST_MESSAGE(ss.str());
        BOOST_CHECK(sequence == t.sequence);
        BOOST_CHECK(nos == t.nos);
    }
}

BOOST_AUTO_TEST_CASE(generateConstraintSeq) {
    // set random generator to a static seed;
    std::mt19937 rand_gen(1);

    initialize_library(true);

    std::vector < TestCase > testcases;
    testcases.push_back(TestCase({G, N, G}, 2,
    {
        G, U, G
    }));
    testcases.push_back(TestCase({N, N, C, N, N}, 4,
    {
        C, G, C, G, C
    }));
    testcases.push_back(TestCase({N, N, G, N, N}, 9,
    {
        A, U, G, U, G
    }));
    testcases.push_back(TestCase({N, N, S, N, N}, 13,
    {
        A, U, G, U, G
    }));

    for (auto t : testcases) {
        // tell the user what is going on right now
        std::stringstream ss;
        ss << "-> Checking Sequence (" << t.input << "):";
        BOOST_TEST_MESSAGE(ss.str());

        rand_gen.seed(1);

        // Build graph
        Graph g(t.input.size());
        int vertex_name = 0;
        // set the vertex_name

        BGL_FORALL_VERTICES_T(v, g, Graph) {
            boost::put(boost::vertex_color_t(), g, v, vertex_name);
            g[v].base = t.input[vertex_name];
            vertex_name++;
        }
        // add the edges
        for (unsigned int i = 0; i < boost::num_vertices(g) - 1; i++) {
            boost::add_edge(boost::vertex(i, g), boost::vertex(i + 1, g), g);
        }

        //call the function
        boost::multiprecision::mpz_int nos = color_path_graph(g, &rand_gen);
        Sequence sequence = get_vertex_colors(g);
        // do the comparison of the results
        ss.str(std::string());
        ss << std::endl << "Got Sequence: " << sequence << " with NOS: " << nos;
        BOOST_TEST_MESSAGE(ss.str());
        BOOST_CHECK(sequence == t.sequence);
        BOOST_CHECK(nos == t.nos);
    }
}

Sequence get_vertex_colors(Graph& g) {
    Sequence sequence;

    BGL_FORALL_VERTICES_T(v, g, Graph) {
        std::cerr << v << ":" << enum_to_char(g[v].base) << std::endl;
        sequence.push_back(g[v].base);
    }
    return sequence;
}

void reset(Graph& g) {
    // reset color

    BGL_FORALL_VERTICES_T(v, g, Graph) {
        g[v].base = N;
    }
}

BOOST_AUTO_TEST_CASE(get_path_pm_test) {
    initialize_library(true);

    std::vector < TestCase > testcases;
    testcases.push_back(TestCase(G, G, 2, 2,{}));
    testcases.push_back(TestCase(A, A, 0, 1,{}));
    testcases.push_back(TestCase(Y, Y, 0, 2,{}));
    testcases.push_back(TestCase(N, G, 3, 5,{}));
    testcases.push_back(TestCase(U, U, 6, 13,{}));
    testcases.push_back(TestCase(R, U, 3, 5,{}));
    testcases.push_back(TestCase(Y, U, 6, 21,{}));
    testcases.push_back(TestCase(C, N, 4, 5,{}));
    testcases.push_back(TestCase(N, N, 1, 6,{}));
    testcases.push_back(TestCase(N, N, 5, 42,{}));
    testcases.push_back(TestCase(C, R, 1, 1,{}));
    testcases.push_back(TestCase(R, Y, 1, 3,{}));
    testcases.push_back(TestCase(H, A, 2, 1,{}));

    for (auto t : testcases) {
        // tell the user what is going on right now
        std::stringstream ss;
        ss << "-> Checking Path (" << enum_to_char(t.first) << ", "
                << enum_to_char(t.last) << ", " << t.length << "):";
        BOOST_TEST_MESSAGE(ss.str());

        // Build graph
        Graph g(t.length + 1);
        int vertex_name = 0;
        // set the vertex_name

        BGL_FORALL_VERTICES_T(v, g, Graph) {
            boost::put(boost::vertex_color_t(), g, v, vertex_name++);
        }
        // add the edges
        for (unsigned int i = 0; i < boost::num_vertices(g) - 1; i++) {
            boost::add_edge(boost::vertex(i, g), boost::vertex(i + 1, g), g);
        }

        // set the sequence constraints
        g[boost::vertex(0, g)].constraint = t.first;
        if (t.first != N) {
            g[boost::vertex(0, g)].special = true;
        }
        g[boost::vertex(boost::num_vertices(g) - 1, g)].constraint = t.last;
        if (t.last != N) {
            g[boost::vertex(boost::num_vertices(g) - 1, g)].special = true;
        }

        //call the function
        ProbabilityMatrix pm = get_path_pm(g);
        // do the comparison of the results
        ss.str(std::string());
        ss << std::endl << "Got PM: " << std::endl << pm << " with MNOS: " << pm.mnos();
        BOOST_TEST_MESSAGE(ss.str());
        BOOST_CHECK(pm.mnos() == t.nos);
    }
}


BOOST_AUTO_TEST_SUITE_END()
