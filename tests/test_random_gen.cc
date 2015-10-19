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
using namespace design;

// include headers containing functions to test

BOOST_AUTO_TEST_SUITE(RandomGenerator)

BOOST_AUTO_TEST_CASE(rand_test) {

    BOOST_TEST_MESSAGE("test random number generator");

    std::mt19937 rand_gen(1);
    RandomDistType dist(0, SolutionSizeType("713297132980479382748047938274"));
    
    for (int i = 0; i < 10; i++) {
        rand_gen.seed(1);
        SolutionSizeType random = dist(rand_gen);
        BOOST_CHECK(random == SolutionSizeType("245000304800877160769938115407"));
    }

}

BOOST_AUTO_TEST_CASE(rand_test1) {


    BOOST_TEST_MESSAGE("test random number generator again");

    std::mt19937 rand_gen(1);
    for (int i = 0; i < 10; i++) {
        rand_gen.seed(1);
        RandomDistType dist(0, SolutionSizeType("71329713298047938274804793827443536456"));
        SolutionSizeType random = dist(rand_gen);
        BOOST_CHECK(random == SolutionSizeType("9730463849175093333028694892082213"));
    }

}

BOOST_AUTO_TEST_CASE(rand_library) {
    BOOST_TEST_MESSAGE("test random number generator in library");
    design::initialize_library(false);
    std::vector<std::string> structures = {"."};
    std::mt19937 rand_gen(1);
    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "A", rand_gen);

    RandomDistType dist(0, SolutionSizeType("748047938274"));
    SolutionSizeType random = dist(*dependency_graph.rand_ptr);
    std::cerr << random << std::endl;
    BOOST_CHECK(random == SolutionSizeType("703173439372"));
}

BOOST_AUTO_TEST_CASE(rand_library1) {
    BOOST_TEST_MESSAGE("test random number generator in library again");
    design::initialize_library(false);
    std::vector<std::string> structures = {"."};
    std::mt19937 rand_gen(1);
    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "A", rand_gen);

    RandomDistType dist(0, SolutionSizeType("713297132980479382748054325234798237447938"));
    SolutionSizeType random = dist(*dependency_graph.rand_ptr);
    std::cerr << random << std::endl;
    BOOST_CHECK(random == SolutionSizeType("317332826215904285124422076931437098021"));
}

BOOST_AUTO_TEST_SUITE_END()
