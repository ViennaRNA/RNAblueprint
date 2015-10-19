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
#ifdef LIBGMP
    RandomDistType dist(0, SolutionSizeType("713297132980479382748047938274"));
#else
    RandomDistType dist(0, 713297132980479382748047938274.0);
#endif
    
    for (int i = 0; i < 10; i++) {
        rand_gen.seed(1);
        SolutionSizeType random = dist(rand_gen);
#ifdef LIBGMP
        BOOST_CHECK(random == SolutionSizeType("245000304800877160769938115407"));
#else
        BOOST_CHECK_CLOSE(random, 245000304800877160769938115407.0, 0.0001);
#endif
    }

}

BOOST_AUTO_TEST_CASE(rand_test1) {


    BOOST_TEST_MESSAGE("test random number generator again");
    
    for (int i = 0; i < 10; i++) {
        std::mt19937 rand_gen(1);
#ifdef LIBGMP
        RandomDistType dist(0, SolutionSizeType("71329713298047938274804793827443536456"));
#else
        RandomDistType dist(0, 71329713298047938274804793827443536456.0);
#endif
        SolutionSizeType random = dist(rand_gen);
#ifdef LIBGMP
        BOOST_CHECK(random == SolutionSizeType("9730463849175093333028694892082213"));
#else
        BOOST_CHECK_CLOSE(random, 9730463849175093333028694892082213.0, 0.0001);
#endif
    }

}

BOOST_AUTO_TEST_CASE(rand_library) {
    BOOST_TEST_MESSAGE("test random number generator in library");
    design::initialize_library(false);
    std::vector<std::string> structures = {"."};
    std::mt19937 rand_gen(1);
    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "A", rand_gen);
#ifdef LIBGMP
    RandomDistType dist(0, SolutionSizeType("748047938274"));
#else
    RandomDistType dist(0, 748047938274.0);
#endif
    SolutionSizeType random = dist(*dependency_graph.rand_ptr);
    std::cerr << random << std::endl;
#ifdef LIBGMP
    BOOST_CHECK(random == SolutionSizeType("703173439372"));
#else
    BOOST_CHECK_CLOSE(random, 703173439372.0, 0.0001);
#endif
}

BOOST_AUTO_TEST_CASE(rand_library1) {
    BOOST_TEST_MESSAGE("test random number generator in library again");
    design::initialize_library(false);
    std::vector<std::string> structures = {"."};
    std::mt19937 rand_gen(1);
    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "A", rand_gen);
#ifdef LIBGMP
    RandomDistType dist(0, SolutionSizeType("713297132980479382748054325234798237447938"));
#else
    RandomDistType dist(0, 713297132980479382748054325234798237447938.0);
#endif
    SolutionSizeType random = dist(*dependency_graph.rand_ptr);
    std::cerr << random << std::endl;
#ifdef LIBGMP
    BOOST_CHECK(random == SolutionSizeType("317332826215904285124422076931437098021"));
#else
    BOOST_CHECK_CLOSE(random, 317332826215904285124422076931437098021.0, 0.0001);
#endif
}

BOOST_AUTO_TEST_SUITE_END()
