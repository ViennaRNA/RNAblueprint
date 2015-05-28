/* This file is a boost.test unit test and provides tests for parsestruct.cc
 *
 *
 * Created on: 27.05.2015
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 * 
 *
 */

// include header
#include "test_common.h"

// include headers containing functions to test
#include "pairing_matrix.h"

// include std components

// include boost components

// define heads

using namespace design;
using namespace design::detail;

BOOST_AUTO_TEST_SUITE (pairingMatrix)


BOOST_AUTO_TEST_CASE (MinimalMatrix) {
    
    initialize_library(true);
    
    BOOST_TEST_MESSAGE("Creating minimal Matrix and check all content");
    PairingMatrix * p = PairingMatrix::Instance();
    
    BOOST_CHECK(p->get(0, A, A) == 1);
    BOOST_CHECK(p->get(0, U, U) == 1);
    BOOST_CHECK(p->get(0, G, G) == 1);
    BOOST_CHECK(p->get(0, C, C) == 1);
    
    BOOST_CHECK(p->get(1, U, G) == 1);
    BOOST_CHECK(p->get(1, G, U) == 1);
    BOOST_CHECK(p->get(1, U, A) == 1);
    BOOST_CHECK(p->get(1, A, U) == 1);
    BOOST_CHECK(p->get(1, C, G) == 1);
    BOOST_CHECK(p->get(1, G, C) == 1);
    BOOST_CHECK(p->get(1, G, G) == 0);
    BOOST_CHECK(p->get(1, G, A) == 0);
    BOOST_CHECK(p->get(1, A, C) == 0);
    BOOST_CHECK(p->get(1, U, C) == 0);
}

BOOST_AUTO_TEST_CASE (ExtendedMatrix) {
    
    initialize_library(true);
    
    BOOST_TEST_MESSAGE("Creating minimal Matrix and extend to 2");
    PairingMatrix * p = PairingMatrix::Instance();
    
    BOOST_CHECK(p->get(2, A, A) == 1);
    BOOST_CHECK(p->get(2, C, C) == 1);
    BOOST_CHECK(p->get(2, A, G) == 1);
    BOOST_CHECK(p->get(2, A, U) == 0);
}

BOOST_AUTO_TEST_CASE (BigMatrix) {
    
    initialize_library(true);
    
    BOOST_TEST_MESSAGE("Creating big Matrix and ask for things");
    PairingMatrix * p = PairingMatrix::Instance();
    
    BOOST_CHECK(p->get(100, A, A) == 16008811023750101250UL);
    BOOST_CHECK(p->get(56, U, C) == 225851433717UL);
    BOOST_CHECK(p->get(45, U, G) == 1836311903UL);
    BOOST_CHECK(p->get(156, A, A) == 16903146174726102521UL);
}

BOOST_AUTO_TEST_SUITE_END ()