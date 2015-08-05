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

BOOST_AUTO_TEST_SUITE (GraphCommon)

BOOST_AUTO_TEST_CASE (min_max_degree) {

  BOOST_TEST_MESSAGE("test min max degree function");

  design::initialize_library(true);
  
  Graph graph = parse_structures({ "()()", ".()." });
  
  int max_degree;
  int min_degree;
  std::tie(min_degree, max_degree) = get_min_max_degree(graph);
  BOOST_CHECK(min_degree == 1);
  BOOST_CHECK(max_degree == 2);
}

BOOST_AUTO_TEST_CASE (min_max_degree1) {

  BOOST_TEST_MESSAGE("test min max degree function for single vertex");

  design::initialize_library(true);
  
  Graph graph = parse_structures({ "." });
  
  int max_degree;
  int min_degree;
  std::tie(min_degree, max_degree) = get_min_max_degree(graph);
  BOOST_CHECK(min_degree == 0);
  BOOST_CHECK(max_degree == 0);
}

BOOST_AUTO_TEST_CASE (is_path_boolean) {

  BOOST_TEST_MESSAGE("test min max degree function for single vertex");

  design::initialize_library(true);
  
  Graph graph = parse_structures({ "." });
  
  boost::get_property(graph, boost::graph_name).is_path = true;
  BOOST_CHECK(boost::get_property(graph, boost::graph_name).is_path == true);
}

BOOST_AUTO_TEST_SUITE_END ()
