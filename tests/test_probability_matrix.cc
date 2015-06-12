/* This file is a boost.test unit test and provides tests for probability_matrix.cc
 *
 *
 * Created on: 11.06.2015
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 * 
 *
 */

// include header
#include "test_common.h"

// include headers containing functions to test
#include "probability_matrix.h"

// include std components

// include boost components

// define heads

using namespace design;
using namespace design::detail;

BOOST_AUTO_TEST_SUITE (probabilityMatrix)


BOOST_AUTO_TEST_CASE (PermuteKey1) {
    
    initialize_library(true);
    
    BOOST_TEST_MESSAGE("Test if we can permute a key");
    
    ProbabilityKey pk;
    pk.emplace(0, N);
    pk.emplace(4, S);
    std::stringstream ss;
    ss << "Key to permute:" << std::endl << pk << std::endl;
    BOOST_TEST_MESSAGE(ss.str());
    ss.str(std::string());
    
    std::vector<ProbabilityKey> test = permute_key(pk);
    
    std::vector<ProbabilityKey> control;
    
    for (auto b : base_conversion[N]) {
        for (auto c : base_conversion[S]) {
            ProbabilityKey pk;
            pk.emplace(0, b);
            pk.emplace(4, c);
            control.push_back(pk);
        }
    }
    
    ss << "Control Key:" << std::endl << control << std::endl;
    ss << "Test Key:" << std::endl << test << std::endl;
    BOOST_TEST_MESSAGE(ss.str());
    BOOST_CHECK(test == control);
}

BOOST_AUTO_TEST_CASE (PermuteKey2) {
    
    initialize_library(true);
    
    BOOST_TEST_MESSAGE("Test if we can permute a key with 3");
    
    ProbabilityKey pk;
    pk.emplace(0, N);
    pk.emplace(4, S);
    pk.emplace(12, V);
    std::stringstream ss;
    ss << "Key to permute:" << std::endl << pk << std::endl;
    BOOST_TEST_MESSAGE(ss.str());
    ss.str(std::string());
    
    std::vector<ProbabilityKey> test = permute_key(pk);
    
    std::vector<ProbabilityKey> control;
    
    for (auto b : base_conversion[N]) {
        for (auto c : base_conversion[S]) {
            for (auto d : base_conversion[V]) {
                ProbabilityKey pk;
                pk.emplace(0, b);
                pk.emplace(4, c);
                pk.emplace(12, d);
                control.push_back(pk);
            }
        }
    }
    
    ss << "Control Key:" << std::endl << control << std::endl;
    ss << "Test Key:" << std::endl << test << std::endl;
    BOOST_TEST_MESSAGE(ss.str());
    BOOST_CHECK(test == control);
}

BOOST_AUTO_TEST_CASE (GetPutKey1) {
    
    ProbabilityMatrix pm;
    unsigned long long control = 0;
    
    for (auto b : base_conversion[N]) {
        for (auto c : base_conversion[S]) {
            ProbabilityKey pk;
            pk.emplace(0, b);
            pk.emplace(4, c);
            pm.put(pk, 4);
            control += 4;
        }
    }
    
    ProbabilityKey pk;
    pk.emplace(0, N);
    pk.emplace(4, S);
    
    std::stringstream ss;
    ss << "Key to get:" << std::endl << pk << std::endl;
    BOOST_TEST_MESSAGE(ss.str());
    
    BOOST_CHECK(pm[pk] == control);
}

BOOST_AUTO_TEST_CASE (GetPutKey2) {
    
    ProbabilityMatrix pm;
    unsigned long long control = 0;
    
    for (auto b : base_conversion[N]) {
        for (auto c : base_conversion[S]) {
            ProbabilityKey pk;
            pk.emplace(0, b);
            pk.emplace(4, c);
            pm.put(pk, 4);
            control += 4;
        }
    }
    
    ProbabilityKey pk;
    pk.emplace(0, N);
    pk.emplace(4, N);
    
    std::stringstream ss;
    ss << "Key to get:" << std::endl << pk << std::endl;
    BOOST_TEST_MESSAGE(ss.str());
    
    BOOST_CHECK(pm[pk] == control);
}

BOOST_AUTO_TEST_CASE (GetPutKey3) {
    
    ProbabilityMatrix pm;
    unsigned long long control = 0;
    
    for (auto b : base_conversion[N]) {
        for (auto c : base_conversion[S]) {
            ProbabilityKey pk;
            pk.emplace(0, b);
            pk.emplace(4, c);
            pm.put(pk, 45);
            control += 45;
        }
    }
    
    ProbabilityKey pk;
    pk.emplace(0, N);
    pk.emplace(4, N);
    
    std::stringstream ss;
    ss << "Key to get:" << std::endl << pk << std::endl;
    BOOST_TEST_MESSAGE(ss.str());
    
    BOOST_CHECK(pm[pk] == control);
}

BOOST_AUTO_TEST_CASE (GetNOS1) {
    
    ProbabilityMatrix pm;
    unsigned long long control = 0;
    
    for (auto b : base_conversion[N]) {
        for (auto c : base_conversion[S]) {
            ProbabilityKey pk;
            pk.emplace(0, b);
            pk.emplace(4, c);
            pm.put(pk, 45);
            control += 45;
        }
    }
    
    BOOST_TEST_MESSAGE("Try to get number of sequences for a ProbabilityMatrix");
    
    BOOST_CHECK(pm.mnos() == control);
}

BOOST_AUTO_TEST_CASE (GetNOS2) {
    
    ProbabilityMatrix pm;
    unsigned long long control = 0;
    
    BOOST_TEST_MESSAGE("Try to get number of sequences for empty ProbabilityMatrix");
    
    BOOST_CHECK(pm.mnos() == control);
}
BOOST_AUTO_TEST_SUITE_END ()
