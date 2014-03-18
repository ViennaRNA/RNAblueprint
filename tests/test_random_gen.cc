/* This file is a boost.test unit test and provides tests the random generator
*
*
* Created on: 07.11.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
* 
*
*/

// include header
#include "test_common.h"

// include headers containing functions to test
#include "decompose.h"

BOOST_AUTO_TEST_SUITE(RandomGenerator)

BOOST_AUTO_TEST_CASE(rand_test) {

	BOOST_TEST_MESSAGE("test random number generator");
	
	std::uniform_real_distribution<float> dist(0, 1);
	for (int i = 0; i < 10; i++) {
		rand_gen.seed(1);
		float random = dist(rand_gen);
		BOOST_CHECK(random == 0.417022f);
	}
	
}

BOOST_AUTO_TEST_CASE(rand_test1) {


	BOOST_TEST_MESSAGE("test random number generator again");
	for (int i = 0; i < 10; i++) {
		rand_gen.seed(1);
		std::uniform_real_distribution<float> dist(0, 1);
		float random = dist(rand_gen);
		BOOST_CHECK(random == 0.417022f);
	}
	
}

BOOST_AUTO_TEST_SUITE_END()
