/* This file is a boost.test unit test and provides tests for parsestruct.cc
 *
 *
 * Created on: 02.11.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 * 
 *
 */

// include header
#include "test_common.h"

// include headers containing functions to test
#include "parsestruct.h"

// include std components

// include boost components
#include <boost/graph/iteration_macros.hpp>

// define heads

using namespace design::detail;

BOOST_AUTO_TEST_SUITE(ParseStruct)

/* put this in
0 1 2 3 4 5 6 7 8 9 10111213
( ( ) ( ) ( ) . ) ( ) . ( )
( ) ( ) ( ) ( . ) . ( ) . .
( . . . . . . ) . . . . . .
. . ( . . . ( ) ) . . . . .
. . . . . . . ( . ) . . . .
 */

/* get graph that looks like this 
  5---6---7---9---10---11
  |   |   |
  4   8---0
  |   |   |
  3---2---1       12---13
 */

std::unordered_set<int> getVertexSet(Graph g) {
    std::unordered_set<int> result;

    BGL_FORALL_VERTICES(v, g, Graph) {
        result.insert(boost::get(boost::vertex_color_t(), g, v));
    }
    return result;
}

BOOST_AUTO_TEST_CASE(ParseDotBracket) {

    // parse this dot-bracket
    std::vector< std::string > structures;
    structures.push_back("(()()().)().()");
    structures.push_back("()()()(.).()..");
    structures.push_back("(......)......");
    structures.push_back("..(...()).....");
    structures.push_back(".......(.)....");

    BOOST_TEST_MESSAGE("parse structure to graph");

    Graph g = parse_structures(structures);

    // check if the resulting vertex names are equal to our testcase
    std::unordered_set<int> testcase{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
    BOOST_CHECK(getVertexSet(g) == testcase);

    // check for certain edges
    BOOST_CHECK(boost::edge(int_to_vertex(0, g), int_to_vertex(8, g), g).second);
    BOOST_CHECK(boost::edge(int_to_vertex(2, g), int_to_vertex(3, g), g).second);
    BOOST_CHECK(boost::edge(int_to_vertex(2, g), int_to_vertex(8, g), g).second);
    // check for degree of certain vertices
    BOOST_CHECK(boost::degree(int_to_vertex(2, g), g) == 3);
    BOOST_CHECK(boost::degree(int_to_vertex(0, g), g) == 3);
    BOOST_CHECK(boost::degree(int_to_vertex(13, g), g) == 1);
}

BOOST_AUTO_TEST_CASE(SetConstraints) {
    Graph g(9);
    int vertex_name = 0;

    BGL_FORALL_VERTICES_T(v, g, Graph) {
        boost::put(boost::vertex_color_t(), g, v, vertex_name++);
    }
    
    BOOST_TEST_MESSAGE("set sequence constraints on graph");
    
    std::string constraints = "CNNNYNSNN";
    set_constraints(g, constraints);
    
    BGL_FORALL_VERTICES_T(v, g, Graph) {
        BOOST_CHECK(g[v].constraint == constraints[vertex_to_int(v, g)]);
    }
}

BOOST_AUTO_TEST_SUITE_END()
