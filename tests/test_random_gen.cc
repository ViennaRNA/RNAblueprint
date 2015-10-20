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
        BOOST_CHECK_CLOSE(random, 7.1128906476233755e+29, 0.0001);
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
    RandomDistType dist(0, 745942039874.8042);
#endif
    SolutionSizeType random = dist(*dependency_graph.rand_ptr);
    std::cerr << random << std::endl;
#ifdef LIBGMP
    BOOST_CHECK(random == SolutionSizeType("703173439372"));
#else
    BOOST_CHECK_CLOSE(random, 743842069983.44971, 0.0001);
#endif
}

#ifdef LIBGMP
BOOST_AUTO_TEST_CASE(int_dist) {
    BOOST_TEST_MESSAGE("test int distribution in library");
    design::initialize_library(false);
    std::vector<std::string> structures = {"."};
    std::mt19937 rand_gen(1);
    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "A", rand_gen);

    RandomDistType dist(0, 10);
    
    std::map<SolutionSizeType, int> count;
    
    for (int i = 0; i < 100000; i++) {
        SolutionSizeType random = dist(*dependency_graph.rand_ptr);
        count[random]++;
    }
    
    for (auto c : count) {
        BOOST_CHECK_CLOSE((double)c.second / 100000, 0.1, 5);
    }
    // 10 must be 0!
    BOOST_CHECK(count[10] == 0);
}
#endif

BOOST_AUTO_TEST_SUITE_END()
