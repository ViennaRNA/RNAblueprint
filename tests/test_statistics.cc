/* This file is a boost.test unit test and provides tests to check if we sample 
 * with valid statistics
 *
 * Created on: 03.08.2015
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 * 
 */

// include header
#include "test_common.h"
#include "RNAdesign.h"
#include <random>
#include <chrono>

// include headers containing functions to test

BOOST_AUTO_TEST_SUITE(Statistics)

BOOST_AUTO_TEST_CASE(EqualDistribution) {

    BOOST_TEST_MESSAGE("Test if we get a equal base distribution");

    design::initialize_library(false);
    unsigned long seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 rand_gen(seed);

    design::DependencyGraph<std::mt19937>* dependency_graph;
    try{
        dependency_graph = new design::DependencyGraph<std::mt19937>(
        {"...................."}, "", rand_gen);
    }

    catch(std::exception & e) {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    std::unordered_map<char, float> count = {
        {'A', 0},
        {'C', 0},
        {'G', 0},
        {'U', 0}};
    ;

    for (int i = 0; i < 10000; i++) {
        try{
            dependency_graph->set_sequence(); // color the graph and get the sequence
        }

        catch(std::exception & e) {
            std::cerr << e.what() << std::endl;
            exit(EXIT_FAILURE);
        }

        std::string sequence = dependency_graph->get_sequence();
        for (int pos = 0; pos < sequence.length(); pos++) {
            count[sequence[pos]]++;
        }
    }

    for (auto c : count) {
        BOOST_CHECK_CLOSE(c.second / 200000, 0.25000f, 1);
    }
}

BOOST_AUTO_TEST_SUITE_END()
