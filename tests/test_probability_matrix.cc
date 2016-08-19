/* This file is a boost.test unit test and provides tests for probability_matrix.cc
 *
 *
 * @date 11.06.2015
 * @author Stefan Hammer <s.hammer@univie.ac.at>
 * @copyright GPLv3
 *
 * 
 *
 */

// include header
#include "test_common.h"

// include headers containing functions to test
#include "probability_matrix.h"

// include std components
#include <deque>

// include boost components

// define heads

using namespace design;
using namespace design::detail;

BOOST_AUTO_TEST_SUITE(probabilityMatrix)


BOOST_AUTO_TEST_CASE(PermuteKey1) {

    initialize_library(true);

    BOOST_TEST_MESSAGE("Test if we can permute a key");

    ProbabilityKey pk;
    pk.emplace(0, N);
    pk.emplace(4, S);
    std::stringstream ss;
    ss << "Key to permute:" << std::endl << pk << std::endl;
    BOOST_TEST_MESSAGE(ss.str());
    ss.str(std::string());
    
    PermuteKeyFactory pkf(pk);
    std::vector<ProbabilityKey> test;
    while (true) { 
        test.push_back(*pkf.key());
        if (!pkf.next_permutation())
            break;
    }

    std::vector<ProbabilityKey> control;

    for (auto c : base_conversion[S]) {
        for (auto b : base_conversion[N]) {
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

BOOST_AUTO_TEST_CASE(PermuteKey2) {

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
    
    PermuteKeyFactory pkf(pk);
    std::vector<ProbabilityKey> test;
    while (true) { 
        test.push_back(*pkf.key());
        if (!pkf.next_permutation())
            break;
    }

    std::vector<ProbabilityKey> control;

    for (auto d : base_conversion[V]) {
        for (auto c : base_conversion[S]) {
            for (auto b : base_conversion[N]) {
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

BOOST_AUTO_TEST_CASE(PermuteKey3) {

    initialize_library(true);

    BOOST_TEST_MESSAGE("Test if we can permute an empty key");

    ProbabilityKey pk;
    std::stringstream ss;
    ss << "Key to permute:" << std::endl << pk << std::endl;
    BOOST_TEST_MESSAGE(ss.str());
    ss.str(std::string());
    
    PermuteKeyFactory pkf(pk);
    std::vector<ProbabilityKey> test;
    
    while (true) {
        test.push_back(*pkf.key());
        if (!pkf.next_permutation())
            break;
    }

    std::vector<ProbabilityKey> control;
    ProbabilityKey key;
    control.push_back(key);
    ss << "Control Key:" << std::endl << control << std::endl;
    ss << "Test Key:" << std::endl << test << std::endl;
    BOOST_TEST_MESSAGE(ss.str());
    BOOST_CHECK(test == control);
}

BOOST_AUTO_TEST_CASE(PermuteKey4) {

    initialize_library(true);

    BOOST_TEST_MESSAGE("Test if we can permute a key and go back");

    ProbabilityKey pk;
    pk.emplace(0, N);
    pk.emplace(4, S);
    pk.emplace(12, V);
    
    std::stringstream ss;
    ss << "Key to permute:" << std::endl << pk << std::endl;
    BOOST_TEST_MESSAGE(ss.str());
    ss.str(std::string());
    
    std::deque<ProbabilityKey> control;
    std::deque<ProbabilityKey> test;
    
    PermuteKeyFactory pkf(pk); 
    // go forward
    for (int i = 0; i < 14; i++) {
        control.push_back(*pkf.key());
        if (!pkf.next_permutation())
            break;
    }
    
    // and go back again
    while (true) {
        if (!pkf.previous_permutation())
            break;
        test.push_front(*pkf.key());
    }
    
    ss << "Control Key:" << std::endl;
    for (auto elem : control) {
        ss << elem << std::endl;
    }
    ss << "Test Key:" << std::endl;
    for (auto elem : test) {
        ss << elem << std::endl;
    }
    BOOST_TEST_MESSAGE(ss.str());
    BOOST_REQUIRE(test == control);
}

BOOST_AUTO_TEST_CASE(GetPutKey1) {

    ProbabilityMatrix pm;
    SolutionSizeType control = 0;

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
    std::set< int > specials = {0, 4};
    BOOST_CHECK(pm.getSpecials() == specials);
}

BOOST_AUTO_TEST_CASE(GetPutKey2) {

    ProbabilityMatrix pm;
    SolutionSizeType control = 0;

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
    std::set< int > specials = {0, 4};
    BOOST_CHECK(pm.getSpecials() == specials);
}

BOOST_AUTO_TEST_CASE(GetPutKey3) {

    ProbabilityMatrix pm;
    SolutionSizeType control = 0;

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
    std::set< int > specials = {0, 4};
    BOOST_CHECK(pm.getSpecials() == specials);
}

BOOST_AUTO_TEST_CASE(GetNOS1) {

    ProbabilityMatrix pm;
    SolutionSizeType control = 0;

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
    std::set< int > specials = {0, 4};
    BOOST_CHECK(pm.getSpecials() == specials);
}

BOOST_AUTO_TEST_CASE(GetNOS2) {

    ProbabilityMatrix pm;
    SolutionSizeType control = 0;

    BOOST_TEST_MESSAGE("Try to get number of sequences for empty ProbabilityMatrix");

    BOOST_CHECK(pm.mnos() == control);
}

BOOST_AUTO_TEST_CASE(CopyConstructor1) {

    ProbabilityMatrix pm;

    for (auto b : base_conversion[N]) {
        for (auto c : base_conversion[S]) {
            ProbabilityKey pk;
            pk.emplace(0, b);
            pk.emplace(4, c);
            pm.put(pk, 45);
        }
    }

    ProbabilityMatrix check(pm);
    BOOST_TEST_MESSAGE("Try to get number of sequences for a ProbabilityMatrix, copy and compare");

    BOOST_CHECK(pm.mnos() == check.mnos());
    BOOST_CHECK(pm.getSpecials() == check.getSpecials());
}

BOOST_AUTO_TEST_CASE(GetSpecials1) {

    ProbabilityMatrix pm;

    BOOST_TEST_MESSAGE("Try to get specials for empty ProbabilityMatrix");

    std::set< int > specials = {};
    BOOST_CHECK(pm.getSpecials() == specials);
}

BOOST_AUTO_TEST_CASE(MultiplyPM1) {

    BOOST_TEST_MESSAGE("Multiply ProbabilityMatrix 1");

    ProbabilityMatrix x;
    ProbabilityMatrix y;
    ProbabilityMatrix z;

    ProbabilityKey xkey;
    for (auto b : base_conversion[N]) {
        xkey[1] = b;
        x.put(xkey, 3 * b);
    }

    ProbabilityKey ykey;
    for (auto b : base_conversion[N]) {
        for (auto c : base_conversion[N]) {
            ykey[1] = b;
            ykey[2] = c;
            y.put(ykey, 2 + b + c);
        }
    }

    z = x * y;

    ProbabilityKey zkey;
    SolutionSizeType mnos = 0;
    for (auto b : base_conversion[N]) {
        for (auto c : base_conversion[N]) {
            zkey[1] = b;
            zkey[2] = c;
            BOOST_CHECK(z[zkey] == ((3 * b)*(2 + b + c)));
            mnos += z[zkey];
        }
    }

    std::cerr << "x:" << std::endl << x << std::endl;
    std::cerr << "y:" << std::endl << y << std::endl;
    std::cerr << "z:" << std::endl << z << std::endl;

    BOOST_CHECK(z.mnos() == mnos);
}

BOOST_AUTO_TEST_CASE(MultiplyPM2) {

    BOOST_TEST_MESSAGE("Multiply ProbabilityMatrix 2");

    ProbabilityMatrix x;
    ProbabilityMatrix y;
    ProbabilityMatrix z;

    ProbabilityKey xkey;
    for (auto b : base_conversion[N]) {
        xkey[1] = b;
        x.put(xkey, 3 * b);
    }

    ProbabilityKey ykey;
    for (auto b : base_conversion[N]) {
        ykey[5] = b;
        y.put(ykey, 5 + b);
    }

    z = x * y;

    ProbabilityKey zkey;
    SolutionSizeType mnos = 0;
    for (auto b : base_conversion[N]) {
        for (auto c : base_conversion[N]) {
            zkey[1] = b;
            zkey[5] = c;
            BOOST_CHECK(z[zkey] == ((3 * b)*(5 + c)));
            mnos += z[zkey];
        }
    }

    std::cerr << "x:" << std::endl << x << std::endl;
    std::cerr << "y:" << std::endl << y << std::endl;
    std::cerr << "z:" << std::endl << z << std::endl;

    BOOST_CHECK(z.mnos() == mnos);
}

BOOST_AUTO_TEST_CASE(MultiplyPM3) {

    BOOST_TEST_MESSAGE("Multiply ProbabilityMatrix 3");

    ProbabilityMatrix x;
    ProbabilityMatrix y;
    ProbabilityMatrix z;

    ProbabilityKey xkey;
    for (auto b : base_conversion[N]) {
        xkey[1] = b;
        x.put(xkey, 3 * b);
    }

    ProbabilityKey ykey;
    for (auto b : base_conversion[N]) {
        ykey[1] = b;
        y.put(ykey, 5 + b);
    }

    z = x * y;

    ProbabilityKey zkey;
    SolutionSizeType mnos = 0;
    for (auto b : base_conversion[N]) {
        zkey[1] = b;
        BOOST_CHECK(z[zkey] == ((3 * b)*(5 + b)));
        mnos += z[zkey];
    }

    std::cerr << "x:" << std::endl << x << std::endl;
    std::cerr << "y:" << std::endl << y << std::endl;
    std::cerr << "z:" << std::endl << z << std::endl;

    BOOST_CHECK(z.mnos() == mnos);
}

BOOST_AUTO_TEST_CASE(MultiplyPM4) {

    BOOST_TEST_MESSAGE("Multiply ProbabilityMatrix with empty");

    ProbabilityMatrix x;
    ProbabilityMatrix y;
    ProbabilityMatrix z;

    ProbabilityKey xkey;
    for (auto b : base_conversion[N]) {
        xkey[1] = b;
        x.put(xkey, 3 * b);
    }

    ProbabilityKey ykey;
    y.put(ykey, 16);

    z = x * y;

    ProbabilityKey zkey;
    SolutionSizeType mnos = 0;
    for (auto b : base_conversion[N]) {
        zkey[1] = b;
        BOOST_CHECK(z[zkey] == 16 * (3 * b));
        mnos += z[zkey];
    }

    std::cerr << "x:" << std::endl << x << std::endl;
    std::cerr << "y:" << std::endl << y << std::endl;
    std::cerr << "z:" << std::endl << z << std::endl;

    BOOST_CHECK(z.mnos() == mnos);
}

BOOST_AUTO_TEST_CASE(MultiplyPM5) {

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
    SolutionSizeType mnos = 0;
    BOOST_CHECK(z[zkey] == 12 * 24);
    mnos += z[zkey];

    std::cerr << "x:" << std::endl << x << std::endl;
    std::cerr << "y:" << std::endl << y << std::endl;
    std::cerr << "z:" << std::endl << z << std::endl;

    BOOST_CHECK(z.mnos() == mnos);
}

BOOST_AUTO_TEST_CASE(MultiplyPM6) {

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

BOOST_AUTO_TEST_CASE(MakeInternal1) {

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
                x.put(xkey, d * b + c);
            }
        }
    }

    y = make_internal(x, 5);

    std::cerr << "x:" << std::endl << x << std::endl;
    std::cerr << "y:" << std::endl << y << std::endl;

    ProbabilityKey testkey;
    testkey[1] = U;
    testkey[7] = C;
    SolutionSizeType ytest = y[testkey];
    testkey[5] = N;
    SolutionSizeType xtest = x[testkey];
    BOOST_CHECK(ytest == xtest);
    BOOST_CHECK(x.getSpecials().size() == y.getSpecials().size() + 1);
    BOOST_CHECK(x.mnos() == y.mnos());
}

BOOST_AUTO_TEST_CASE(MakeInternal2) {

    BOOST_TEST_MESSAGE("make last vertex internal in ProbabilityMatrix");

    ProbabilityMatrix x;
    ProbabilityMatrix y;

    ProbabilityKey xkey;
    for (auto b : base_conversion[N]) {
        xkey[5] = b;
        x.put(xkey, 3 * b + 4);
    }

    y = make_internal(x, 5);

    std::cerr << "x:" << std::endl << x << std::endl;
    std::cerr << "y:" << std::endl << y << std::endl;

    ProbabilityKey testkey;
    SolutionSizeType ytest = y[testkey];
    testkey[5] = N;
    SolutionSizeType xtest = x[testkey];
    BOOST_CHECK(ytest == xtest);
    BOOST_CHECK(x.getSpecials().size() == y.getSpecials().size() + 1);
    BOOST_CHECK(x.mnos() == y.mnos());
}

BOOST_AUTO_TEST_CASE(RandomlySampleKey1) {

    BOOST_TEST_MESSAGE("randomly sample a probability key from a matrix");

    std::uniform_real_distribution<float> dist(0, 1);
    std::mt19937 rand_gen(1);

    ProbabilityMatrix m;

    ProbabilityKey input;
    input[1] = N;
    input[4] = N;
    input[7] = Y;
    
    PermuteKeyFactory pkf(input);
    while (true) {
        m.put(*pkf.key(), static_cast<SolutionSizeType>(trunc(dist(rand_gen) * 400)));
        if (!pkf.next_permutation())
            break;
    }

    ProbabilityKey constraint;
    constraint[1] = A;
    constraint[4] = V;
    constraint[7] = N;

    std::pair<ProbabilityKey, SolutionSizeType> result = m.sample(constraint, rand_gen);
    BOOST_CHECK(result.first[1] == A);
    BOOST_CHECK(result.first[4] == A);
    BOOST_CHECK(result.first[7] == U);
    BOOST_CHECK(result.second == 630);
}

BOOST_AUTO_TEST_CASE(RandomlySampleKey2) {

    BOOST_TEST_MESSAGE("randomly sample a probability key from a matrix, no constraints");
    
    std::mt19937 rand_gen(1);

    ProbabilityMatrix m;

    ProbabilityKey input;
    input[1] = A;
    input[4] = N;
    input[7] = Y;

    PermuteKeyFactory pkf(input);
    while (true) {
        m.put(*pkf.key(), 16);
        if (!pkf.next_permutation())
            break;
    }
    ProbabilityMatrix c = m;

    std::pair<ProbabilityKey, SolutionSizeType> result = m.sample(rand_gen);
    // remember for every < vertex, a map of < base, count > to get mean value in the end
    std::unordered_map<int, std::unordered_map<int, SolutionSizeType> > result_stats;
    std::unordered_map<int, std::unordered_map<int, double> > stats_check {
        {1, { {A, 1.00}, {C, 0.00}, {G, 0.00}, {U, 0.00} } },
        {4, { {A, 0.25}, {C, 0.25}, {G, 0.25}, {U, 0.25} } },
        {7, { {U, 0.50}, {C, 0.50}, {G, 0.00}, {U, 0.00} } },
    };
    
    for (int i = 0; i < 10000; i++) {
        std::pair<ProbabilityKey, SolutionSizeType> test = m.sample(rand_gen);
        for (auto p : test.first) {
            result_stats[p.first][p.second]++;
        }
    }
    
    for (auto e : stats_check) {
        for (auto b : e.second) {
            BOOST_CHECK_CLOSE((double)result_stats[e.first][b.first] / 10000, b.second, 5);
        }
    }
    BOOST_CHECK(m.mnos() == c.mnos());
    
    
}

BOOST_AUTO_TEST_SUITE_END()
