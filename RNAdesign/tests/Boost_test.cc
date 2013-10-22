/* This file tests the generate_path_seq funktion of the file
* pathcoloring.h
*
* Created on: 21.08.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
* compile with g++ -Wall -std=c++11 -lstdc++ -g -o test_pathcoloring ../common.cc ../graphcommon.cc ../printgraph.cc ../pathcoloring.cc test_pathcoloring.cc
*
*/

// Include Boost.Test
#define BOOST_TEST_MODULE Test_Path_Coloring
#include <boost/test/unit_test.hpp>

// include header
#include "../pathcoloring.h"
#include "../common.h"


BOOST_AUTO_TEST_CASE(generate_path_seq) {
	BOOST_CHECK(
		Sequence sequence;
		Sequence check;
		check.push_back(G);
		check.push_back(G);
		generate_path_seq (sequence, G, G, 2);
		sequence == check;
}

