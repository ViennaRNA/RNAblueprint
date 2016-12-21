/* This file is a boost.test unit test and provides tests the random generator
 *
 *
 * @date 07.11.2013
 * @author Stefan Hammer <s.hammer@univie.ac.at>
 * @copyright GPLv3
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
    RandomDistType dist(0, 713297132980479382748047938274.0);
    
    for (int i = 0; i < 10; i++) {
        rand_gen.seed(1);
        SolutionSizeType random = dist(rand_gen);
#ifdef LIBGMP
        BOOST_CHECK_CLOSE(random, 297460595874988724119142727680.0, 0.0001);
#else
        BOOST_CHECK_CLOSE(random, 7.1128906476233755e+29, 0.0001);
#endif
    }

}

BOOST_AUTO_TEST_CASE(rand_copy) {
    //TODO think of a nice way to check if copy construction works
}

BOOST_AUTO_TEST_CASE(int_dist_set) {
    BOOST_TEST_MESSAGE("test int distribution in set_sequence function");
    design::initialize_library(false);
    std::vector<std::string> structures = {"."};
    std::mt19937 rand_gen(1);
    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "N", rand_gen);
    
    std::map<int, int> count;
    // set initial sequence
    dependency_graph.sample();
    
    for (int i = 0; i < 10000; i++) {
        dependency_graph.sample();
        Sequence r = dependency_graph.get_sequence();
        count[r[0]]++;
    }
    
    for (auto c : count) {
        BOOST_CHECK_CLOSE(c.second / 10000.0, 0.25, 3);
    }
    // 10 must be 0!
    BOOST_CHECK(count[4] == 0);
}

BOOST_AUTO_TEST_SUITE_END()
