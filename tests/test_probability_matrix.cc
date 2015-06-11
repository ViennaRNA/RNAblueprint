/* This file is a boost.test unit test and provides tests for probability_matrix.cc
 *
 *
 * Created on: 11.06.2015
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 * 
 *
 */

// include header
#include "test_common.h"

// include headers containing functions to test
#include "probability_matrix.h"

// include std components

// include boost components

// define heads

using namespace design;
using namespace design::detail;

BOOST_AUTO_TEST_SUITE (probabilityMatrix)


BOOST_AUTO_TEST_CASE (FirstTest) {
    
    initialize_library(true);
    
    BOOST_TEST_MESSAGE("Creating minimal Matrix and check all content");
    BOOST_CHECK(1 == 1);
}

BOOST_AUTO_TEST_SUITE_END ()
