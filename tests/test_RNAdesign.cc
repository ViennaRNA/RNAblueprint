/* This file is a boost.test unit test and provides tests to check if we sample 
 * with valid statistics
 *
 * Created on: 11.01.2016
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 * 
 */

// include header
#include "test_common.h"
#include "RNAdesign.h"

// include headers containing functions to test

BOOST_AUTO_TEST_SUITE(RNAdesign)

BOOST_AUTO_TEST_CASE(IsBipartite) {

    BOOST_TEST_MESSAGE("Check if we can detect a bipartite graph.");
    
    std::vector<std::string> structures = {"(.)", ".()"};
    bool result = design::graph_is_bipartite(structures);
    BOOST_CHECK(result);
    
    structures.push_back("().");
    result = design::graph_is_bipartite(structures);
    BOOST_CHECK(!result);
}

BOOST_AUTO_TEST_CASE(SequenceStructureCompatible) {

    BOOST_TEST_MESSAGE("Check if we can detect a sequence incompatible to several structures.");
    
    std::vector<std::string> structures = {"(.)"};
    bool result = design::sequence_structure_compatible("NNY", structures);
    BOOST_CHECK(result);
    
    result = design::sequence_structure_compatible("CNW", structures);
    BOOST_CHECK(!result);
}

BOOST_AUTO_TEST_CASE(IncompatiblePositions) {

    BOOST_TEST_MESSAGE("Check if we can get the incompatible sequence positions.");
    
    std::vector<int> result = design::incompatible_sequence_positions("NNY", "(.)");
    std::vector<int> testcase = {};
    BOOST_CHECK(result == testcase);
    
    result = design::incompatible_sequence_positions("CNW", "(.)");
    testcase = {0, 2};
    BOOST_CHECK(result == testcase);
}

BOOST_AUTO_TEST_SUITE_END()
