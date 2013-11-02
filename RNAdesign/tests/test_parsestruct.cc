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
#include "../parsestruct.h"

// include std components

// include boost components
#include <boost/graph/iteration_macros.hpp>

// define heads

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
	
	Graph result = parse_structures(structures);
	
	// check if the resulting vertex names are equal to our testcase
	std::unordered_set<int> testcase {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
	BOOST_CHECK(getVertexSet(result) == testcase);
	
	// TODO check for certain edges
	// TODO check for degree of certain vertices
}

BOOST_AUTO_TEST_SUITE_END()
