/* This file is a boost.test unit test and provides tests the internal dependency graph
 *
 * Created on: 03.08.2015
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 * 
 */

// include boost components
#include <boost/test/unit_test.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

// include header
#include "uninduced_subgraph.hpp"

using namespace boost;

// include headers containing functions to test

BOOST_AUTO_TEST_SUITE(UninducedSubgraph)

BOOST_AUTO_TEST_CASE(simpleGraph) {

    BOOST_TEST_MESSAGE("simple uninduced subgraph");

    typedef uninduced_subgraph< adjacency_list< vecS, vecS, directedS,
        no_property, property< edge_index_t, int > > > Graph;
    typedef Graph::edge_descriptor Edge;
    typedef Graph::vertex_descriptor Vertex;

    const int N = 6;
    Graph G0(N);

    enum { A, B, C, D, E, F};  // for conveniently refering to vertices in G0

    Graph& G1 = G0.create_subgraph(), G2 = G0.create_subgraph();
    enum { A1, B1, C1 }; // for conveniently refering to vertices in G1
    enum { A2, B2, C2, D2 };     // for conveniently refering to vertices in G2
    
    add_vertex(C, G1); // global vertex C becomes local A1 for G1
    add_vertex(E, G1); // global vertex E becomes local B1 for G1
    add_vertex(F, G1); // global vertex F becomes local C1 for G1

    add_vertex(A, G2); // global vertex A becomes local A2 for G2
    add_vertex(B, G2); // global vertex B becomes local B2 for G2
    
    BOOST_CHECK(num_vertices(G0) == 6);
    BOOST_CHECK(num_vertices(G1) == 3);
    BOOST_CHECK(num_vertices(G2) == 2);
    
    // add edges to root graph
    add_edge(A, B, G0);
    add_edge(B, C, G0);
    add_edge(B, D, G0);
    add_edge(E, F, G0);
    
    BOOST_CHECK(num_edges(G0) == 4);
    BOOST_CHECK(num_edges(G1) == 0);
    BOOST_CHECK(num_edges(G2) == 0);
    
    // add edges to G1
    add_edge(A1, B1, G1);
    BOOST_CHECK(num_edges(G0) == 5);
    BOOST_CHECK(num_edges(G1) == 1);
    BOOST_CHECK(num_edges(G2) == 0);
    // num_vertices stays the same
    BOOST_CHECK(num_vertices(G0) == 6);
    BOOST_CHECK(num_vertices(G1) == 3);
    BOOST_CHECK(num_vertices(G2) == 2);
    
    // add edges to G2 with edge descriptor
    Edge g_edge;
    bool present;
    boost::tie(g_edge, present) = edge(E, F, G0);
    BOOST_CHECK(present);
    add_edge(g_edge, G2); // global vertex E becomes local C2 and global vertex F becomes local D2 for G2
    
    BOOST_CHECK(num_edges(G0) == 5);
    BOOST_CHECK(num_edges(G1) == 1);
    BOOST_CHECK(num_edges(G2) == 1);
    
    BOOST_CHECK(num_vertices(G0) == 6);
    BOOST_CHECK(num_vertices(G1) == 3);
    BOOST_CHECK(num_vertices(G2) == 4);
}

BOOST_AUTO_TEST_SUITE_END()
