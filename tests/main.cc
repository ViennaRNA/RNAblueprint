/* This file is a boost.test unit test and should provide tests for many functions of the project
 *
 *
 * @date 23.10.2013
 * @author Stefan Hammer <s.hammer@univie.ac.at>
 * @copyright GPLv3
 *
 * compile with g++ -Wall -std=c++11 -lstdc++ -g -o testRNAblueprint -lboost_unit_test_framework main.cc test_pathcoloring.cc ../pathcoloring.cc ../common.cc ../graphcommon.cc
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "C++ Unit Tests for RNAblueprint"
#include <boost/test/included/unit_test.hpp>

#include "test_common.h"

bool debug = false;
