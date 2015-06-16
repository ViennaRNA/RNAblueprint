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
#include "common.h"

using namespace design::detail;

// include headers containing functions to test

BOOST_AUTO_TEST_SUITE (RandomGenerator)

BOOST_AUTO_TEST_CASE (rand_test) {

  BOOST_TEST_MESSAGE("test random number generator");

  std::mt19937 rand_gen(1);
  std::uniform_real_distribution<float> dist(0, 1);
  for (int i = 0; i < 10; i++) {
    rand_gen.seed(1);
    float random = dist(rand_gen);
    BOOST_CHECK(random == 0.417022f);
  }

}

BOOST_AUTO_TEST_CASE (rand_test1) {


  BOOST_TEST_MESSAGE("test random number generator again");
  
  std::mt19937 rand_gen(1);
  for (int i = 0; i < 10; i++) {
    rand_gen.seed(1);
    std::uniform_real_distribution<float> dist(0, 1);
    float random = dist(rand_gen);
    BOOST_CHECK(random == 0.417022f);
  }

}

BOOST_AUTO_TEST_CASE (rand_library) {
  BOOST_TEST_MESSAGE("test random number generator in library");
  design::initialize_library(false);
  std::vector<std::string> structures = { "." };
  std::mt19937 rand_gen(1);
  design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "A", rand_gen);
  
  std::uniform_real_distribution<float> dist(0, 1);
  float random = dist(*dependency_graph.rand_ptr);
  BOOST_CHECK_CLOSE(random, 0.997185f, 0.0001);
}

BOOST_AUTO_TEST_CASE (rand_library1) {
  BOOST_TEST_MESSAGE("test random number generator in library again");
  design::initialize_library(false);
  std::vector<std::string> structures = { "." };
  std::mt19937 rand_gen(1);
  design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "A", rand_gen);
  
  std::uniform_real_distribution<float> dist(0, 1);
  float random = dist(*dependency_graph.rand_ptr);
  BOOST_CHECK_CLOSE(random, 0.997185f, 0.0001);
}

BOOST_AUTO_TEST_SUITE_END ()
