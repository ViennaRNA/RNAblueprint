/* This file is a boost.test unit test and provides tests the internal dependency graph
 *
 * @date 03.08.2015
 * @author Stefan Hammer <s.hammer@univie.ac.at>
 * @copyright GPLv3
 *
 */

// include header
#include "test_common.h"
#include "common.h"
#include "parsestruct.h"

using namespace design::detail;
using namespace design;

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

    auto previous = dependency_graph.get_sequence_string();
    auto cnos = dependency_graph.sample();
    auto actual = dependency_graph.get_sequence_string();
    BOOST_CHECK(actual != previous);
    BOOST_CHECK(cnos == 47);
}

BOOST_AUTO_TEST_CASE(cc_sampling) {

    BOOST_TEST_MESSAGE("test sampling on an cc example");

    design::initialize_library(true);
    std::vector<std::string> structures = {"((.().)).", "()()().()", ".()..(..)", ".....()..", "...(..).."};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNNNNNNN", rand_gen);
    BOOST_CHECK(dependency_graph.number_of_sequences() == 134);
    auto previous = dependency_graph.get_sequence_string();
    auto cnos = dependency_graph.sample();
    auto actual = dependency_graph.get_sequence_string();
    BOOST_CHECK(actual != previous);
    BOOST_CHECK(cnos == 133);

}


BOOST_AUTO_TEST_CASE(output_connected_components) {

    BOOST_TEST_MESSAGE("test the connected_components() function for correct output");

    design::initialize_library(true);
    std::vector<std::string> structures = {"(..)"};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNN", rand_gen);
    // result
    std::map< int, std::vector<int> > check_output = {{0, {0, 3}}, {1, {1}}, {2, {2}}};
    // iterate over connected components
    for (int i = 0; i < dependency_graph.number_of_connected_components(); i++) {
        BOOST_CHECK(dependency_graph.component_vertices(i) == check_output[i]);
    }

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

    std::map< int, std::vector<int> > check = {{0, {0}}, {1, {1}}};
    for (int i = 0; i < dependency_graph.number_of_connected_components(); i++) {
        BOOST_CHECK(dependency_graph.component_vertices(i) == check[i]);
    }

    BOOST_CHECK(dependency_graph.number_of_sequences(0) == 2);
    BOOST_CHECK(dependency_graph.number_of_sequences(1) == 4);
}

BOOST_AUTO_TEST_CASE(get_articulations1) {

    BOOST_TEST_MESSAGE("get empty list of articulation vertices");

    design::initialize_library(true);
    std::vector<std::string> structures = {".."};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "WN", rand_gen);
    std::vector< int > check = {};
    BOOST_CHECK(dependency_graph.articulation_vertices() == check);
}

BOOST_AUTO_TEST_CASE(get_articulations2) {

    BOOST_TEST_MESSAGE("get list of articulation vertices for more complex graph");

    design::initialize_library(true);
    std::vector<std::string> structures = {"().().", ".().()", ".(..)."};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNNNY", rand_gen);
    std::vector< int > check = {1, 4};
    BOOST_CHECK(dependency_graph.articulation_vertices() == check);
}

BOOST_AUTO_TEST_CASE(sample_cc_with_ID) {

    BOOST_TEST_MESSAGE("check if we can sample the cc with ccID");

    design::initialize_library(true);
    std::vector<std::string> structures = {".."};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "WN", rand_gen);
    BOOST_CHECK(dependency_graph.number_of_sequences() == 8);

    for (int i = 0; i < 10; i++) {
        std::string previous = dependency_graph.get_sequence_string();
        SolutionSizeType cnos = dependency_graph.sample_clocal(1);
        std::string actual = dependency_graph.get_sequence_string();

        BOOST_CHECK(cnos == 3);
        BOOST_CHECK(actual[1] != previous[1]);
        BOOST_CHECK(actual[0] == previous[0]);
    }

    for (int i = 0; i < 10; i++) {
        std::string previous = dependency_graph.get_sequence_string();
        SolutionSizeType cnos = dependency_graph.sample_clocal(0);
        std::string actual = dependency_graph.get_sequence_string();

        BOOST_CHECK(cnos == 1);
        BOOST_CHECK(actual[0] != previous[0] && (actual[0] == 'A' || actual[0] == 'U'));
        BOOST_CHECK(actual[1] == previous[1]);
    }
}


BOOST_AUTO_TEST_CASE(sample_cc_with_ID_fix) {

    BOOST_TEST_MESSAGE("check if we can sample the cc with ccID with fixed sequence constraint");

    design::initialize_library(true);
    std::vector<std::string> structures = {".."};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "AN", rand_gen);
    BOOST_CHECK(dependency_graph.number_of_sequences() == 4);

    std::string previous = dependency_graph.get_sequence_string();
    SolutionSizeType cnos = dependency_graph.sample_clocal(0);
    std::string actual = dependency_graph.get_sequence_string();

    BOOST_CHECK(cnos == 0);
    BOOST_CHECK(actual == previous);

}

BOOST_AUTO_TEST_CASE(sample_clocal1) {

    BOOST_TEST_MESSAGE("check if we can sample globally");

    design::initialize_library(true);
    std::vector<std::string> structures = {".."};
    std::mt19937 rand_gen(4);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "WN", rand_gen);
    BOOST_CHECK(dependency_graph.number_of_sequences() == 8);
    BOOST_CHECK(dependency_graph.get_sequence_string()[0] == 'A' || dependency_graph.get_sequence_string()[0] == 'U');

    for (int i = 0; i < 10; i++) {
        std::string previous = dependency_graph.get_sequence_string();
        SolutionSizeType cnos = dependency_graph.sample_local_global(1, 0, 0);
        std::string actual = dependency_graph.get_sequence_string();
        BOOST_CHECK(cnos == 5);
        BOOST_CHECK((actual[0] != previous[0] && actual[1] == previous[1])
                 || (actual[0] == previous[0] && actual[1] != previous[1]));
    }
}

BOOST_AUTO_TEST_CASE(sample_clocal1_fix) {

    BOOST_TEST_MESSAGE("check if we can sample globally with fixed sequence constraint");

    design::initialize_library(true);
    std::vector<std::string> structures = {".."};
    std::mt19937 rand_gen(4);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "AN", rand_gen);
    BOOST_CHECK(dependency_graph.number_of_sequences() == 4);

    std::string previous = dependency_graph.get_sequence_string();
    SolutionSizeType cnos = dependency_graph.sample_local_global(1, 0, 0);
    std::string actual = dependency_graph.get_sequence_string();
    BOOST_CHECK(cnos == 3);
    BOOST_CHECK(actual != previous);
}

BOOST_AUTO_TEST_CASE(sample_clocal1_complete_fix) {

    BOOST_TEST_MESSAGE("check if we can sample globally with completely fixed sequence constraint");

    design::initialize_library(true);
    std::vector<std::string> structures = {".."};
    std::mt19937 rand_gen(4);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "AG", rand_gen);
    BOOST_CHECK(dependency_graph.number_of_sequences() == 1);

    std::string previous = dependency_graph.get_sequence_string();
    SolutionSizeType cnos = dependency_graph.sample_local_global(1, 0, 0);
    std::string actual = dependency_graph.get_sequence_string();
    BOOST_CHECK(cnos == 0);
    BOOST_CHECK(actual == previous);
}

BOOST_AUTO_TEST_CASE(sample_plocal1) {

    BOOST_TEST_MESSAGE("check if we can sample locally");

    design::initialize_library(true);
    std::vector<std::string> structures = {".."};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "WN", rand_gen);
    BOOST_CHECK(dependency_graph.number_of_sequences() == 8);
    BOOST_CHECK(dependency_graph.get_sequence_string()[0] == 'A' || dependency_graph.get_sequence_string()[0] == 'U');

    for (int i = 0; i < 10; i++) {
        std::string previous = dependency_graph.get_sequence_string();
        SolutionSizeType cnos = dependency_graph.sample_local_global(-1, 0, 0);
        std::string actual = dependency_graph.get_sequence_string();
        BOOST_CHECK(cnos == 5);
        BOOST_CHECK((actual[0] != previous[0] && actual[1] == previous[1])
                 || (actual[0] == previous[0] && actual[1] != previous[1]));
    }
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

BOOST_AUTO_TEST_CASE(sample_pos) {

    BOOST_TEST_MESSAGE("sample dependency graph by position");

    design::initialize_library(true);
    std::vector<std::string> structures = { "()....()..()",
                                            ".()()..(.)..",
                                            ".(.)...()..." };
                                           //NNNNNNNGNNNN

    std::mt19937 rand_gen(1);
    // connected components: cc0 ps 28 #p 5 0,1,2,3,4 cc1 ps 4 #p 1 5 cc2 ps 18 #p 4 6,7,8,9 cc3 ps 6 #p 2 10,11
    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNNNNNGNNNN", rand_gen);
    BOOST_CHECK(dependency_graph.number_of_sequences() == 5376);
    // sample CCs which are paths
    for (int i = 0; i<100; i++) {
        // get sequence
        std::string previous = dependency_graph.get_sequence_string();
        BOOST_CHECK(dependency_graph.sample(5) == 3);
        auto actual = dependency_graph.get_sequence_string();
        // check if only position 5 changed
        for (unsigned int j=0; j < previous.length(); j++) {
            if (j != 5 ) {
                BOOST_CHECK(previous[j] == actual[j]);
            } else {
                BOOST_CHECK(previous[j] != actual[j]);
            }
        }
    }

    // sample CCs which are paths
    for (int i = 0; i<10; i++) {
        std::string previous = dependency_graph.get_sequence_string();
        BOOST_CHECK(dependency_graph.sample(10) == 5);
        auto actual = dependency_graph.get_sequence_string();
        // check if only position 10 and 11 changed
        for (int j=0; j < 10; j++) {
            BOOST_CHECK(previous[j] == actual[j]);
        }
        BOOST_CHECK((previous[10] != actual[10]) || (previous[11] != actual[11]));
    }

    // sample more nested subgraphs
    // sample CCs which are paths
    for (int i = 0; i<10; i++) {
        std::string previous = dependency_graph.get_sequence_string();
        BOOST_CHECK(dependency_graph.sample(8) == 1);
        auto actual = dependency_graph.get_sequence_string();
        // check if only position 8 changed
        for (unsigned int j=0; j < previous.length(); j++) {
            if (j != 8 ) {
                BOOST_CHECK(previous[j] == actual[j]);
            } else {
                BOOST_CHECK(previous[j] != actual[j]);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(sample_pos_range) {
    BOOST_TEST_MESSAGE("sample dependency graph by position range");

    design::initialize_library(true);
    std::vector<std::string> structures = { "()....()..()",
                                            ".()()..(.)..",
                                            ".(.)...()..." };
                                           //NNNNNNNGNNNN
                                           //     5
    std::mt19937 rand_gen(1);
    // connected components: cc0 ps 28 #p 5 0,1,2,3,4 cc1 ps 4 #p 1 5 cc2 ps 18 #p 4 6,7,8,9 cc3 ps 6 #p 2 10,11
    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNNNNNGNNNN", rand_gen);

    BOOST_CHECK(dependency_graph.number_of_sequences() == 5376);
    // sample positions which are in one path
    // sample CCs which are paths
    for (int i = 0; i<10; i++) {
        std::string previous = dependency_graph.get_sequence_string();
        BOOST_CHECK(dependency_graph.sample(10, 11) == 5);
        auto actual = dependency_graph.get_sequence_string();
        // check if only position 10 and 11 changed
        for (int j=0; j < 10; j++) {
            BOOST_CHECK(previous[j] == actual[j]);
        }
        BOOST_CHECK((previous[10] != actual[10]) || (previous[11] != actual[11]));
    }
    // sample on different CCs
    for (int i = 0; i<10; i++) {
        std::string previous = dependency_graph.get_sequence_string();
        BOOST_CHECK(dependency_graph.sample(9, 11) == 11);
        auto actual = dependency_graph.get_sequence_string();
        // check if only position 9, 10 and 11 changed
        for (int j=0; j < 9; j++) {
            BOOST_CHECK(previous[j] == actual[j]);
        }
        BOOST_CHECK((previous[9]  != actual[9])  ||
                    (previous[10] != actual[10]) ||
                    (previous[11] != actual[11]));
    }
}

BOOST_AUTO_TEST_CASE(sample_pos_range_fix) {

    BOOST_TEST_MESSAGE("try to sample range with only fixed constraints");

    design::initialize_library(true);
    std::vector<std::string> structures = {"..."};
    std::mt19937 rand_gen(1);

    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "ANN", rand_gen);
    BOOST_CHECK(dependency_graph.number_of_sequences() == 16);

    std::string previous = dependency_graph.get_sequence_string();
    SolutionSizeType cnos = dependency_graph.sample(0, 0);
    std::string actual = dependency_graph.get_sequence_string();
    BOOST_CHECK(cnos == 0);
    BOOST_CHECK(actual == previous);
}

BOOST_AUTO_TEST_CASE(revert_sequence) {
    BOOST_TEST_MESSAGE("test the revert sequence functionality");

    design::initialize_library(true);
    std::vector<std::string> structures = {"....."};
    std::mt19937 rand_gen(1);
    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNNN", rand_gen);

    dependency_graph.set_history_size(10);

    BOOST_CHECK(!dependency_graph.revert_sequence(1));
    std::string sequence = dependency_graph.get_sequence_string();
    BOOST_CHECK(!dependency_graph.revert_sequence(1));

    dependency_graph.sample_clocal(1);
    BOOST_CHECK(dependency_graph.revert_sequence(1));
    BOOST_CHECK(dependency_graph.get_sequence_string() == sequence);

    dependency_graph.sample(0, 3);
    dependency_graph.sample_local_global(1, 0, 0);
    BOOST_CHECK(dependency_graph.revert_sequence(2));
    BOOST_CHECK(dependency_graph.get_sequence_string() == sequence);
    BOOST_CHECK(!dependency_graph.revert_sequence(1));

    dependency_graph.sample();
    dependency_graph.sample();
    BOOST_CHECK(!dependency_graph.revert_sequence(4));

    // now make history big so that we have to forget
    for (int i = 0; i < 10; i++) {
        dependency_graph.sample();
    }
    BOOST_CHECK(!dependency_graph.revert_sequence(12));
    // go all the way back
    while (1) {
        if (!dependency_graph.revert_sequence(1))
        break;
    }
    // check if we did not arrive at initial sequence and indeed forgot something
    BOOST_CHECK(dependency_graph.get_sequence_string() != sequence);

    // now the same with bigger history stack
    dependency_graph.set_history_size(20);
    sequence = dependency_graph.get_sequence_string();
    for (int i = 0; i < 14; i++) {
        dependency_graph.sample();
    }
    while (1) {
        if (!dependency_graph.revert_sequence(1))
        break;
    }
    // check if we did not arrive at initial sequence and indeed forgot something
    BOOST_CHECK(dependency_graph.get_sequence_string() == sequence);
}

BOOST_AUTO_TEST_CASE(set_history_size) {
    BOOST_TEST_MESSAGE("test to set history size function");

    design::initialize_library(true);
    std::vector<std::string> structures = {"....."};
    std::mt19937 rand_gen(1);
    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNNN", rand_gen);

    dependency_graph.set_history_size(10);
    for (int i = 0; i < 20; i++) {
        dependency_graph.sample();
    }
    BOOST_CHECK(dependency_graph.get_history().size() == 10);
    // check if we can make the stack smaller again
    dependency_graph.set_history_size(5);
    BOOST_CHECK(dependency_graph.get_history().size() == 5);
    dependency_graph.sample();
}

BOOST_AUTO_TEST_CASE(set_history_size_zero) {
    BOOST_TEST_MESSAGE("test to set history size to zero");

    design::initialize_library(true);
    std::vector<std::string> structures = {"....."};
    std::mt19937 rand_gen(1);
    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NNNNN", rand_gen);


    BOOST_CHECK_THROW(dependency_graph.set_history_size(0), std::out_of_range);
    BOOST_CHECK(dependency_graph.get_history().size() == 1);
}

BOOST_AUTO_TEST_CASE(get_history) {
    BOOST_TEST_MESSAGE("test the get_history functionality");

    design::initialize_library(true);
    std::vector<std::string> structures = {"..+..."};
    std::mt19937 rand_gen(1);
    design::detail::DependencyGraph<std::mt19937> dependency_graph(structures, "NN+NNN", rand_gen);

    dependency_graph.set_history_size(10);

    for (int i = 0; i < 20; i++) {
        dependency_graph.sample();
    }

    std::vector<std::string> history1 = dependency_graph.get_history();

    std::cout << history1 << std::endl;
    BOOST_CHECK(history1.size() == 10);

    dependency_graph.set_history_size(100);
    for (int i = 0; i < 20; i++) {
        dependency_graph.sample();
    }

    std::vector<std::string> history2 = dependency_graph.get_history();
    std::cout << history2 << std::endl;
    BOOST_CHECK(history2.size() == 30);
    BOOST_CHECK(history1[0] == history2[0]);
    BOOST_CHECK(history1[9] == history2[9]);
}

BOOST_AUTO_TEST_SUITE_END()
