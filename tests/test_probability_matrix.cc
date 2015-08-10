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
    std::set< int > specials = { 0, 4 };
    BOOST_CHECK(pm.getSpecials() == specials);
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
    std::set< int > specials = { 0, 4 };
    BOOST_CHECK(pm.getSpecials() == specials);
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
    std::set< int > specials = { 0, 4 };
    BOOST_CHECK(pm.getSpecials() == specials);
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
    std::set< int > specials = { 0, 4 };
    BOOST_CHECK(pm.getSpecials() == specials);
}

BOOST_AUTO_TEST_CASE (GetNOS2) {
    
    ProbabilityMatrix pm;
    unsigned long long control = 0;
    
    BOOST_TEST_MESSAGE("Try to get number of sequences for empty ProbabilityMatrix");
    
    BOOST_CHECK(pm.mnos() == control);
}

BOOST_AUTO_TEST_CASE (GetSpecials1) {
    
    ProbabilityMatrix pm;
    
    BOOST_TEST_MESSAGE("Try to get specials for empty ProbabilityMatrix");
    
    std::set< int > specials = { };
    BOOST_CHECK(pm.getSpecials() == specials);
}

BOOST_AUTO_TEST_CASE (MultiplyPM1) {
    
    BOOST_TEST_MESSAGE("Multiply ProbabilityMatrix 1");
    
    ProbabilityMatrix x;
    ProbabilityMatrix y;
    ProbabilityMatrix z;
    
    ProbabilityKey xkey;
    for (auto b : base_conversion[N]) {
        xkey[1] = b;
        x.put(xkey, 3*b);
    }
    
    ProbabilityKey ykey;
    for (auto b : base_conversion[N]) {
        for (auto c : base_conversion[N]) {
            ykey[1] = b;
            ykey[2] = c;
            y.put(ykey, 2+b+c);
        }
    }
    
    z = x * y;
    
    ProbabilityKey zkey;
    unsigned long long mnos = 0;
    for (auto b : base_conversion[N]) {
        for (auto c : base_conversion[N]) {
            zkey[1] = b;
            zkey[2] = c;
            BOOST_CHECK(z[zkey] == ((3*b)*(2+b+c)));
            mnos += z[zkey];
        }
    }
    
    std::cerr << "x:" << std::endl << x << std::endl;
    std::cerr << "y:" << std::endl << y << std::endl;
    std::cerr << "z:" << std::endl << z << std::endl;
    
    BOOST_CHECK(z.mnos() == mnos);
}

BOOST_AUTO_TEST_CASE (MultiplyPM2) {
    
    BOOST_TEST_MESSAGE("Multiply ProbabilityMatrix 2");
    
    ProbabilityMatrix x;
    ProbabilityMatrix y;
    ProbabilityMatrix z;
    
    ProbabilityKey xkey;
    for (auto b : base_conversion[N]) {
        xkey[1] = b;
        x.put(xkey, 3*b);
    }
    
    ProbabilityKey ykey;
    for (auto b : base_conversion[N]) {
        ykey[5] = b;
        y.put(ykey, 5+b);
    }
    
    z = x * y;
    
    ProbabilityKey zkey;
    unsigned long long mnos = 0;
    for (auto b : base_conversion[N]) {
        for (auto c : base_conversion[N]) {
            zkey[1] = b;
            zkey[5] = c;
            BOOST_CHECK(z[zkey] == ((3*b)*(5+c)));
            mnos += z[zkey];
        }
    }
    
    std::cerr << "x:" << std::endl << x << std::endl;
    std::cerr << "y:" << std::endl << y << std::endl;
    std::cerr << "z:" << std::endl << z << std::endl;
    
    BOOST_CHECK(z.mnos() == mnos);
}

BOOST_AUTO_TEST_CASE (MultiplyPM3) {
    
    BOOST_TEST_MESSAGE("Multiply ProbabilityMatrix 3");
    
    ProbabilityMatrix x;
    ProbabilityMatrix y;
    ProbabilityMatrix z;
    
    ProbabilityKey xkey;
    for (auto b : base_conversion[N]) {
        xkey[1] = b;
        x.put(xkey, 3*b);
    }
    
    ProbabilityKey ykey;
    for (auto b : base_conversion[N]) {
        ykey[1] = b;
        y.put(ykey, 5+b);
    }
    
    z = x * y;
    
    ProbabilityKey zkey;
    unsigned long long mnos = 0;
    for (auto b : base_conversion[N]) {
        zkey[1] = b;
        BOOST_CHECK(z[zkey] == ((3*b)*(5+b)));
        mnos += z[zkey];
    }
    
    std::cerr << "x:" << std::endl << x << std::endl;
    std::cerr << "y:" << std::endl << y << std::endl;
    std::cerr << "z:" << std::endl << z << std::endl;
    
    BOOST_CHECK(z.mnos() == mnos);
}

BOOST_AUTO_TEST_CASE (MultiplyPM4) {
    
    BOOST_TEST_MESSAGE("Multiply ProbabilityMatrix with empty");
    
    ProbabilityMatrix x;
    ProbabilityMatrix y;
    ProbabilityMatrix z;
    
    ProbabilityKey xkey;
    for (auto b : base_conversion[N]) {
        xkey[1] = b;
        x.put(xkey, 3*b);
    }
    
    ProbabilityKey ykey;
    y.put(ykey, 16);
    
    z = x * y;
    
    ProbabilityKey zkey;
    unsigned long long mnos = 0;
    for (auto b : base_conversion[N]) {
        zkey[1] = b;
        BOOST_CHECK(z[zkey] == 16*(3*b));
        mnos += z[zkey];
    }
    
    std::cerr << "x:" << std::endl << x << std::endl;
    std::cerr << "y:" << std::endl << y << std::endl;
    std::cerr << "z:" << std::endl << z << std::endl;
    
    BOOST_CHECK(z.mnos() == mnos);
}

BOOST_AUTO_TEST_CASE (MultiplyPM5) {
    
    BOOST_TEST_MESSAGE("Multiply two empty ProbabilityMatrix");
    
    ProbabilityMatrix x;
    ProbabilityMatrix y;
    ProbabilityMatrix z;
    
    ProbabilityKey xkey;
    x.put(xkey, 12);
    
    ProbabilityKey ykey;
    y.put(ykey, 24);
    
    z = x * y;
    
    ProbabilityKey zkey;
    unsigned long long mnos = 0;
    BOOST_CHECK(z[zkey] == 12*24);
    mnos += z[zkey];
    
    std::cerr << "x:" << std::endl << x << std::endl;
    std::cerr << "y:" << std::endl << y << std::endl;
    std::cerr << "z:" << std::endl << z << std::endl;
    
    BOOST_CHECK(z.mnos() == mnos);
}

BOOST_AUTO_TEST_CASE (MultiplyPM6) {
    
    BOOST_TEST_MESSAGE("Multiply with complete empty (only constructed) ProbabilityMatrix");
    
    ProbabilityMatrix x;
    ProbabilityMatrix y;
    ProbabilityMatrix z;
    
    ProbabilityKey xkey;
    x.put(xkey, 12);
    
    z = x * y;
    
    BOOST_CHECK(z[xkey] == 12);
    
    std::cerr << "x:" << std::endl << x << std::endl;
    std::cerr << "y:" << std::endl << y << std::endl;
    std::cerr << "z:" << std::endl << z << std::endl;
    
    BOOST_CHECK(z.mnos() == x.mnos());
}

BOOST_AUTO_TEST_CASE (MakeInternal1) {
    
    BOOST_TEST_MESSAGE("make a vertex internal in ProbabilityMatrix");
    
    ProbabilityMatrix x;
    ProbabilityMatrix y;
    
    ProbabilityKey xkey;
    for (auto b : base_conversion[N]) {
        for (auto c : base_conversion[N]) {
            for (auto d : base_conversion[N]) {
                xkey[1] = b;
                xkey[5] = c;
                xkey[7] = d;
                x.put(xkey, d*b+c);
            }
        }
    }
    
    y = make_internal(x, 5);

    std::cerr << "x:" << std::endl << x << std::endl;
    std::cerr << "y:" << std::endl << y << std::endl;
    
    ProbabilityKey testkey;
    testkey[1] = U;
    testkey[7] = C;
    unsigned long long ytest = y[testkey];
    testkey[5] = N;
    unsigned long long xtest = x[testkey];
    BOOST_CHECK(ytest == xtest);
    BOOST_CHECK(x.getSpecials().size() == y.getSpecials().size()+1);
    BOOST_CHECK(x.mnos() == y.mnos());
}

BOOST_AUTO_TEST_CASE (MakeInternal2) {
    
    BOOST_TEST_MESSAGE("make last vertex internal in ProbabilityMatrix");
    
    ProbabilityMatrix x;
    ProbabilityMatrix y;
    
    ProbabilityKey xkey;
    for (auto b : base_conversion[N]) {
        xkey[5] = b;
        x.put(xkey, 3*b+4);
    }
    
    y = make_internal(x, 5);

    std::cerr << "x:" << std::endl << x << std::endl;
    std::cerr << "y:" << std::endl << y << std::endl;
    
    ProbabilityKey testkey;
    unsigned long long ytest = y[testkey];
    testkey[5] = N;
    unsigned long long xtest = x[testkey];
    BOOST_CHECK(ytest == xtest);
    BOOST_CHECK(x.getSpecials().size() == y.getSpecials().size()+1);
    BOOST_CHECK(x.mnos() == y.mnos());
}

BOOST_AUTO_TEST_CASE (RandomlySampleKey1) {
    
    BOOST_TEST_MESSAGE("randomly sample a probability key from a matrix");
    
    std::uniform_real_distribution<float> dist(0, 1);
    std::mt19937 rand_gen(1);
    
    ProbabilityMatrix m;
    
    ProbabilityKey input;
    input[1] = N;
    input[4] = N;
    input[7] = Y;
    
    std::vector<ProbabilityKey> input_keys = permute_key(input);
    for (auto i : input_keys) {
        m.put(i, dist(rand_gen) * 400);
    }
    
    ProbabilityKey constraint;
    constraint[1] = A;
    constraint[4] = V;
    constraint[7] = N;
    
    ProbabilityKey result = m.sample(constraint, &rand_gen);
    BOOST_CHECK(result[1] == A);
    BOOST_CHECK(result[4] == A);
    BOOST_CHECK(result[7] == U);
}

BOOST_AUTO_TEST_CASE (RandomlySampleKey2) {
    
    BOOST_TEST_MESSAGE("randomly sample a probability key from a matrix, no constraints");
    
    std::uniform_real_distribution<float> dist(0, 1);
    std::mt19937 rand_gen(1);
    
    ProbabilityMatrix m;
    
    ProbabilityKey input;
    input[1] = N;
    input[4] = N;
    input[7] = Y;
    
    std::vector<ProbabilityKey> input_keys = permute_key(input);
    for (auto i : input_keys) {
        m.put(i, dist(rand_gen) * 400);
    }
    
    ProbabilityKey result = m.sample(&rand_gen);
    BOOST_CHECK(result[1] == C);
    BOOST_CHECK(result[4] == U);
    BOOST_CHECK(result[7] == U);
}

BOOST_AUTO_TEST_SUITE_END ()
