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
#include "graphcommon.h"

using namespace design::detail;

// include headers containing functions to test

BOOST_AUTO_TEST_SUITE(GraphCommon)

BOOST_AUTO_TEST_CASE(min_max_degree) {

    BOOST_TEST_MESSAGE("test min max degree function");

    design::initialize_library(true);

    Graph graph = parse_structures({"()()", ".()."});

    int max_degree;
    int min_degree;
    std::tie(min_degree, max_degree) = get_min_max_degree(graph);
    BOOST_CHECK(min_degree == 1);
    BOOST_CHECK(max_degree == 2);
}

BOOST_AUTO_TEST_CASE(min_max_degree1) {

    BOOST_TEST_MESSAGE("test min max degree function for single vertex");

    design::initialize_library(true);

    Graph graph = parse_structures({"."});

    int max_degree;
    int min_degree;
    std::tie(min_degree, max_degree) = get_min_max_degree(graph);
    BOOST_CHECK(min_degree == 0);
    BOOST_CHECK(max_degree == 0);
}

BOOST_AUTO_TEST_CASE(is_path_boolean) {

    BOOST_TEST_MESSAGE("test min max degree function for single vertex");

    design::initialize_library(true);

    Graph graph = parse_structures({"."});

    boost::get_property(graph, boost::graph_name).is_path = true;
    BOOST_CHECK(boost::get_property(graph, boost::graph_name).is_path == true);
}

BOOST_AUTO_TEST_CASE(vertex_to_int_to_vertex) {
    
    Graph graph = parse_structures({"...."});
    
    for (int i = 0; i < 4; i++) {
        Vertex v = design::detail::int_to_vertex(i, graph);
        BOOST_CHECK( design::detail::vertex_to_int( v, graph) == i );
    }
    // create a subgraph
    Graph& subg = graph.create_subgraph();
    // get vertex descriptor for #1 on graph
    Vertex v = design::detail::int_to_vertex(1, graph);
    // check this vertex
    BOOST_CHECK( design::detail::vertex_to_int( v, graph ) == 1 );
    // add to subgraph and get descriptor on subg
    Vertex u = boost::add_vertex(v, subg);
    // check index on subg
    BOOST_CHECK( design::detail::vertex_to_int( u, subg ) == 1 );
    // check for out of range error
    BOOST_CHECK_THROW( int_to_vertex(6, graph), std::out_of_range );
}

BOOST_AUTO_TEST_SUITE_END()
