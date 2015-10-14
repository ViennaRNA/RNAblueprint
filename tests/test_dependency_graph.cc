/* This file is a boost.test unit test and provides tests the internal dependency graph
 *
 * Created on: 03.08.2015
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 * 
 */

// include header
#include "test_common.h"
#include "common.h"
#include "parsestruct.h"

using namespace design::detail;

// include headers containing functions to test

BOOST_AUTO_TEST_SUITE(DependencyGraph)

BOOST_AUTO_TEST_CASE(easy_construction) {

    BOOST_TEST_MESSAGE("test dependency graph construction on an easy example");

    design::initialize_library(true);
    std::vector<std::string> structures = {"(.)"};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NAN", rand_gen);

    BOOST_CHECK(dependency_graph.number_of_sequences() == 6);
}

BOOST_AUTO_TEST_CASE(easy_cc_construction) {

    BOOST_TEST_MESSAGE("test dependency graph construction on an easy cc example");

    design::initialize_library(true);
    std::vector<std::string> structures = {"....."};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNNN", rand_gen);

    BOOST_CHECK(dependency_graph.number_of_sequences() == 1024);
}

BOOST_AUTO_TEST_CASE(middle_construction) {

    BOOST_TEST_MESSAGE("test dependency graph construction on an bc example");

    design::initialize_library(true);
    std::vector<std::string> structures = {"()..", ".().", ".(.)"};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNN", rand_gen);

    BOOST_CHECK(dependency_graph.number_of_sequences() == 18);
}

BOOST_AUTO_TEST_CASE(cc_construction) {

    BOOST_TEST_MESSAGE("test dependency graph construction on an cc example");

    design::initialize_library(true);
    std::vector<std::string> structures = {".()..", "..().", "..(.)"};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNN", rand_gen);

    BOOST_CHECK(dependency_graph.number_of_sequences() == 72);
}

BOOST_AUTO_TEST_CASE(cc_construction_constraints) {

    BOOST_TEST_MESSAGE("test dependency graph construction on an cc example with constraints");

    design::initialize_library(true);
    std::vector<std::string> structures = {".()..", "..().", "..(.)"};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNA", rand_gen);

    BOOST_CHECK(dependency_graph.number_of_sequences() == 16);
}

BOOST_AUTO_TEST_CASE(bc_construction) {

    BOOST_TEST_MESSAGE("test dependency graph construction on an bc example");

    design::initialize_library(true);
    std::vector<std::string> structures = {"().().", ".().()", ".(..)."};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNNNN", rand_gen);

    BOOST_CHECK(dependency_graph.number_of_sequences() == 48);
}

BOOST_AUTO_TEST_CASE(bc_sampling) {

    BOOST_TEST_MESSAGE("test sampling on an bc example");

    design::initialize_library(true);
    std::vector<std::string> structures = {"().().", ".().()", ".(..)."};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNNNN", rand_gen);
    BOOST_CHECK(dependency_graph.number_of_sequences() == 48);
    dependency_graph.set_sequence();
    std::cerr << "Sampled sequence: " << dependency_graph.get_sequence_string() << std::endl;
}

BOOST_AUTO_TEST_CASE(cc_sampling) {

    BOOST_TEST_MESSAGE("test sampling on an cc example");

    design::initialize_library(true);
    std::vector<std::string> structures = {"((.().)).", "()()().()", ".()..(..)", ".....()..", "...(..).."};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNNNNNNN", rand_gen);
    dependency_graph.set_sequence();
    std::cerr << "Sampled sequence: " << dependency_graph.get_sequence_string() << std::endl;
    
    BOOST_CHECK(dependency_graph.number_of_sequences() == 134);
}


BOOST_AUTO_TEST_CASE(output_connected_components) {

    BOOST_TEST_MESSAGE("test the connected_components() function for correct output");

    design::initialize_library(true);
    std::vector<std::string> structures = {"(..)"};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNN", rand_gen);
    std::map< int, std::vector<int> > check_output = {{0, {0, 3}}, {1, {1}}, {2, {2}}};

    BOOST_CHECK(dependency_graph.connected_components() == check_output);
}

BOOST_AUTO_TEST_CASE(set_sequence1) {

    BOOST_TEST_MESSAGE("set a sequence onto the graph");

    design::initialize_library(true);
    std::vector<std::string> structures = {"(..)"};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNN", rand_gen);
    dependency_graph.set_sequence_string("AAAU");
    
    BOOST_CHECK(dependency_graph.get_sequence_string() == "AAAU");
}

BOOST_AUTO_TEST_CASE(set_sequence2) {

    BOOST_TEST_MESSAGE("set a sequence onto the graph with failure");

    design::initialize_library(true);
    std::vector<std::string> structures = {"(..)"};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNN", rand_gen);
    BOOST_REQUIRE_THROW(dependency_graph.set_sequence_string("AAAA"), std::exception);
}

BOOST_AUTO_TEST_CASE(set_sequence3) {

    BOOST_TEST_MESSAGE("set a sequence onto the graph with failure and recover");

    design::initialize_library(true);
    std::vector<std::string> structures = {"(..)"};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNN", rand_gen);
    dependency_graph.set_sequence_string("AGCU");
    BOOST_CHECK(dependency_graph.get_sequence_string() == "AGCU");
    BOOST_REQUIRE_THROW(dependency_graph.set_sequence_string("AAAA"), std::exception);
    BOOST_CHECK(dependency_graph.get_sequence_string() == "AGCU");
}

BOOST_AUTO_TEST_CASE(set_sequence4) {

    BOOST_TEST_MESSAGE("set a sequence onto the graph with failure and recover");

    design::initialize_library(true);
    std::vector<std::string> structures = {"(..)"};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NGNN", rand_gen);
    BOOST_REQUIRE_THROW(dependency_graph.set_sequence_string("AAAU"), std::exception);
}

BOOST_AUTO_TEST_CASE(set_sequence5) {

    BOOST_TEST_MESSAGE("set a sequence onto the path and failure");

    design::initialize_library(true);
    std::vector<std::string> structures = {".."};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "GN", rand_gen);
    BOOST_REQUIRE_THROW(dependency_graph.set_sequence_string("AU"), std::exception);
}

BOOST_AUTO_TEST_CASE(set_sequence6) {

    BOOST_TEST_MESSAGE("set a sequence containing N and die");

    design::initialize_library(true);
    std::vector<std::string> structures = {"().().", ".().()", ".(..)."};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNNNN", rand_gen);
    BOOST_CHECK(dependency_graph.number_of_sequences() == 48);
    BOOST_REQUIRE_THROW(dependency_graph.set_sequence_string("UNNNNG"), std::exception);
}

BOOST_AUTO_TEST_CASE(number_of_sequences_ccID) {

    BOOST_TEST_MESSAGE("check if we can get the nos for ccID");

    design::initialize_library(true);
    std::vector<std::string> structures = {".."};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "WN", rand_gen);
    BOOST_CHECK(dependency_graph.number_of_sequences() == 8);
    std::map< int, std::vector<int> > result = dependency_graph.connected_components();
    std::map< int, std::vector<int> > check = {{0, {0}}, {1, {1}}};
    
    BOOST_CHECK(dependency_graph.number_of_sequences(0) == 2);
    BOOST_CHECK(dependency_graph.number_of_sequences(1) == 4);
}

BOOST_AUTO_TEST_CASE(get_specials1) {

    BOOST_TEST_MESSAGE("get empty list of special vertices");

    design::initialize_library(true);
    std::vector<std::string> structures = {".."};
    std::mt19937 rand_gen(1);
    
    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "WN", rand_gen);
    std::vector< int > check = {};
    BOOST_CHECK(dependency_graph.special_vertices() == check);
}

BOOST_AUTO_TEST_CASE(get_specials2) {

    BOOST_TEST_MESSAGE("get list of special vertices for more complex graph");

    design::initialize_library(true);
    std::vector<std::string> structures = {"().().", ".().()", ".(..)."};
    std::mt19937 rand_gen(1);
    
    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNNNY", rand_gen);
    std::vector< int > check = {1, 4};
    BOOST_CHECK(dependency_graph.special_vertices() == check);
}

BOOST_AUTO_TEST_CASE(mutate_cc_with_ID) {

    BOOST_TEST_MESSAGE("check if we can mutate the cc with ccID");

    design::initialize_library(true);
    std::vector<std::string> structures = {".."};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "WN", rand_gen);
    BOOST_CHECK(dependency_graph.number_of_sequences() == 8);
    dependency_graph.set_sequence();
    std::cerr << dependency_graph.get_sequence_string() << std::endl;
    BOOST_CHECK(dependency_graph.get_sequence_string() == "UG");
    
    boost::multiprecision::mpz_int cnos = dependency_graph.mutate_global(1);
    BOOST_CHECK(cnos == 4);
    BOOST_CHECK(dependency_graph.get_sequence_string() == "UA");
    
    cnos = dependency_graph.mutate_global(0);
    BOOST_CHECK(cnos == 2);
    BOOST_CHECK(dependency_graph.get_sequence_string() == "AA");
    
    cnos = dependency_graph.mutate_global(0);
    BOOST_CHECK(cnos == 2);
    BOOST_CHECK(dependency_graph.get_sequence_string() == "UA");
}

BOOST_AUTO_TEST_CASE(mutate_global1) {

    BOOST_TEST_MESSAGE("check if we can mutate globally");

    design::initialize_library(true);
    std::vector<std::string> structures = {".."};
    std::mt19937 rand_gen(4);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "WN", rand_gen);
    BOOST_CHECK(dependency_graph.number_of_sequences() == 8);
    dependency_graph.set_sequence();
    BOOST_CHECK(dependency_graph.get_sequence_string() == "AG");
    
    boost::multiprecision::mpz_int cnos = dependency_graph.mutate_local_global(1, 0, 0);
    BOOST_CHECK(cnos == 2);
    BOOST_CHECK(dependency_graph.get_sequence_string() == "UG");
    
    cnos = dependency_graph.mutate_local_global(1, 0, 1);
    BOOST_CHECK(cnos == 2);
    BOOST_CHECK(dependency_graph.get_sequence_string() == "UG");
}

BOOST_AUTO_TEST_CASE(mutate_local1) {

    BOOST_TEST_MESSAGE("check if we can mutate locally");

    design::initialize_library(true);
    std::vector<std::string> structures = {".."};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "WN", rand_gen);
    BOOST_CHECK(dependency_graph.number_of_sequences() == 8);
    dependency_graph.set_sequence();
    std::cerr << dependency_graph.get_sequence_string() << std::endl;
    BOOST_CHECK(dependency_graph.get_sequence_string() == "UG");
    
    boost::multiprecision::mpz_int cnos = dependency_graph.mutate_local_global(-1, 0, 0);
    BOOST_CHECK(cnos == 4);
    std::cerr << dependency_graph.get_sequence_string() << std::endl;
    BOOST_CHECK(dependency_graph.get_sequence_string() == "UC");
    
    cnos = dependency_graph.mutate_local_global(-1, 0, 1);
    BOOST_CHECK(cnos == 2);
    std::cerr << dependency_graph.get_sequence_string() << std::endl;
    BOOST_CHECK(dependency_graph.get_sequence_string() == "AC");
}

BOOST_AUTO_TEST_CASE(number_of_sequences_cc) {

    BOOST_TEST_MESSAGE("test dependency graph number_of_sequences(cc_id)");

    design::initialize_library(true);
    std::vector<std::string> structures = {"()....()..()", ".()()..(.)..", ".(.)...()..."};
    
    std::mt19937 rand_gen(1);
    // connected components: cc0 ps 28 #p 5 0,1,2,3,4 cc1 ps 4 #p 1 5 cc2 ps 18 #p 4 6,7,8,9 cc3 ps 6 #p 2 10,11
    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNNNNNNNNNN", rand_gen);
    
    BOOST_CHECK(dependency_graph.number_of_sequences(0) == 28);
    BOOST_CHECK(dependency_graph.number_of_sequences(1) == 4);
    BOOST_CHECK(dependency_graph.number_of_sequences(2) == 18);
    BOOST_CHECK(dependency_graph.number_of_sequences(3) == 6);
    BOOST_CHECK(dependency_graph.number_of_sequences() == 12096);
}

//TODO mutate global local unit tests!
// todo mutate by position unit tests



BOOST_AUTO_TEST_SUITE_END()
