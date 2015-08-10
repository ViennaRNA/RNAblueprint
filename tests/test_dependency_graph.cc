/* This file is a boost.test unit test and provides tests the internal dependency graph
 *
 * Created on: 03.08.2015
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 * 
 */

// include header
#include "test_common.h"
#include "common.h"
#include "parsestruct.h"

using namespace design::detail;

// include headers containing functions to test

BOOST_AUTO_TEST_SUITE(DependencyGraph)

BOOST_AUTO_TEST_CASE(easy_construction) {

    BOOST_TEST_MESSAGE("test dependency graph construction on an easy example");

    design::initialize_library(true);
    std::vector<std::string> structures = {"(.)"};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NAN", rand_gen);

    BOOST_CHECK(dependency_graph.number_of_sequences() == 6);
}

BOOST_AUTO_TEST_CASE(middle_construction) {

    BOOST_TEST_MESSAGE("test dependency graph construction on an middle example");

    design::initialize_library(true);
    std::vector<std::string> structures = {"()..", ".().", ".(.)"};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNN", rand_gen);

    BOOST_CHECK(dependency_graph.number_of_sequences() == 18);
}

BOOST_AUTO_TEST_CASE(cc_construction) {

    BOOST_TEST_MESSAGE("test dependency graph construction on an cc example");

    design::initialize_library(true);
    std::vector<std::string> structures = {".()..", "..().", "..(.)"};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNN", rand_gen);

    BOOST_CHECK(dependency_graph.number_of_sequences() == 72);
}

BOOST_AUTO_TEST_CASE(cc_construction_constraints) {

    BOOST_TEST_MESSAGE("test dependency graph construction on an cc example with constraints");

    design::initialize_library(true);
    std::vector<std::string> structures = {".()..", "..().", "..(.)"};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNA", rand_gen);

    BOOST_CHECK(dependency_graph.number_of_sequences() == 16);
}

BOOST_AUTO_TEST_CASE(bc_construction) {

    BOOST_TEST_MESSAGE("test dependency graph construction on an bc example");

    design::initialize_library(true);
    std::vector<std::string> structures = {"().().", ".().()", ".(..)."};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNNNN", rand_gen);

    BOOST_CHECK(dependency_graph.number_of_sequences() == 48);
}

BOOST_AUTO_TEST_CASE(bc_sampling) {

    BOOST_TEST_MESSAGE("test sampling on an bc example");

    design::initialize_library(true);
    std::vector<std::string> structures = {"().().", ".().()", ".(..)."};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNNNN", rand_gen);
    BOOST_CHECK(dependency_graph.number_of_sequences() == 48);
    
    dependency_graph.mutate();
    std::cerr << "Sampled sequence: " << dependency_graph.get_sequence_string() << std::endl;
}
/*
BOOST_AUTO_TEST_CASE (cc_sampling) {

  BOOST_TEST_MESSAGE("test sampling on an cc example");

  design::initialize_library(true);
  std::vector<std::string> structures = { "((.().)).", "()()().()", ".()..(..)", ".....()..", "...(..).." };
  std::mt19937 rand_gen(1);
  
  design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNNNNNNN", rand_gen);
  for (int i=0; i<100; i++) {
    dependency_graph.mutate();
  
     std::cerr << "Sampled sequence: " << dependency_graph.get_sequence_string() << std::endl;
  }    
  
  //BOOST_CHECK(random == 0.417022f);
}
 */
BOOST_AUTO_TEST_SUITE_END()
