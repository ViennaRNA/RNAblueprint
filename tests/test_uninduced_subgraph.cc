/* This file is a boost.test unit test and provides tests the internal dependency graph
 *
 * @date 03.08.2015
 * @author Stefan Hammer <s.hammer@univie.ac.at>
 * @copyright GPLv3
 * 
 */

// std lib includes
#include <iostream>

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
    //typedef Graph::edge_descriptor Edge;
    //typedef Graph::vertex_descriptor Vertex;

    const int N = 6;
    Graph G0(N);

    enum { A, B, C, D, E, F};  // for conveniently refering to vertices in G0

    Graph& G1 = G0.create_subgraph();
    Graph& G2 = G1.create_subgraph();
    
    BOOST_CHECK(&G1.parent() == &G0);
    BOOST_CHECK(&G2.parent() == &G1);
    BOOST_CHECK(G0.is_root());
    BOOST_CHECK(!G1.is_root());
    BOOST_CHECK(!G2.is_root());
    
    enum { A1, B1, C1 }; // for conveniently refering to vertices in G1
    enum { A2, B2, C2 };     // for conveniently refering to vertices in G2
    
    add_vertex(C, G1); // global vertex C becomes local A1 for G1
    add_vertex(E, G1); // global vertex E becomes local B1 for G1
    add_vertex(F, G1); // global vertex F becomes local C1 for G1

    add_vertex(C, G2); // global vertex C becomes local A2 for G2
    add_vertex(E, G2); // global vertex E becomes local B2 for G2
    
    BOOST_CHECK(num_vertices(G0) == 6);
    BOOST_CHECK(num_vertices(G1) == 3);
    std::cerr << num_vertices(G1) << std::endl;
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
}

BOOST_AUTO_TEST_CASE(addVertices) {

    BOOST_TEST_MESSAGE("uninduced subgraph add edges");

    typedef uninduced_subgraph< adjacency_list< vecS, vecS, directedS,
        no_property, property< edge_index_t, int > > > Graph;
    //typedef Graph::edge_descriptor Edge;
    typedef Graph::vertex_descriptor Vertex;

    const int N = 3;
    Graph G0(N);
    Graph& G1 = G0.create_subgraph();
    Graph& G2 = G1.create_subgraph();
    
    BOOST_CHECK(&G1.parent() == &G0);
    BOOST_CHECK(&G2.parent() == &G1);
    
    // add vertices to G2
    Vertex n1 = add_vertex(0, G2);
    Vertex n2 = add_vertex(1, G2);
    // check if the global vertex 2 is equal to the returned local vertex
    if (G2.find_vertex(0).second) {
        BOOST_CHECK(G2.find_vertex(0).first == n1);
    } else {
        BOOST_ERROR( "vertex not found!" );
    }
    if (G2.find_vertex(1).second) {
        BOOST_CHECK(G2.find_vertex(1).first == n2);
    } else {
        BOOST_ERROR( "vertex not found!" );
    }
    // and check if this vertex is also present in G1
    if (G1.find_vertex(0).second) {
        BOOST_CHECK(G1.local_to_global(G1.find_vertex(0).first) == 0);
    } else {
        BOOST_ERROR( "vertex not found!" );
    }
    if (G1.find_vertex(0).second) {
        BOOST_CHECK(G1.local_to_global(G1.find_vertex(1).first) == 1);
    } else {
        BOOST_ERROR( "vertex not found!" );
    }
    
    // num_vertices stays the same
    BOOST_CHECK(num_vertices(G0) == 3);
    BOOST_CHECK(num_vertices(G1) == 2);
    BOOST_CHECK(num_vertices(G2) == 2);
    
    // add vertices to G1
    Vertex n3 = add_vertex(2, G1);
    // check if the global vertex 2 is equal to the returned local vertex
    if (G1.find_vertex(2).second) {
        BOOST_CHECK(G1.find_vertex(2).first == n3);
    } else {
        BOOST_ERROR( "vertex not found!" );
    }
    // num_vertices stays the same
    BOOST_CHECK(num_vertices(G0) == 3);
    BOOST_CHECK(num_vertices(G1) == 3);
    BOOST_CHECK(num_vertices(G2) == 2);
    
    // add vertices to G1
    Vertex n4 = add_vertex(G1);
    
    // check if the new local vertex is also in the global graph
    BOOST_CHECK(G0.find_vertex(G1.local_to_global(n4)).second);
    // check if the new local vertex is not in the subgraphs
    BOOST_CHECK(!G2.find_vertex(n4).second);
    
    // num_vertices stays the same
    BOOST_CHECK(num_vertices(G0) == 4);
    BOOST_CHECK(num_vertices(G1) == 4);
    BOOST_CHECK(num_vertices(G2) == 2);
    
    // add vertices to G0
    Vertex n5 = add_vertex(G0);
    
    // check if the new local vertex is not in the subgraphs
    BOOST_CHECK(!G1.find_vertex(n5).second);
    BOOST_CHECK(!G2.find_vertex(n5).second);
    
    // num_vertices stays the same
    BOOST_CHECK(num_vertices(G0) == 5);
    BOOST_CHECK(num_vertices(G1) == 4);
    BOOST_CHECK(num_vertices(G2) == 2);
    
    std::cerr << "All G0 vertices: " << std::endl;
    BGL_FORALL_VERTICES_T(v, G0, Graph) {
        std::cerr << G0.local_to_global(v) << std::endl;
    }
    std::cerr << "All G1 vertices: " << std::endl;
    BGL_FORALL_VERTICES_T(v, G1, Graph) {
        std::cerr << G1.local_to_global(v) << std::endl;
    }
    std::cerr << "All G2 vertices: " << std::endl;
    BGL_FORALL_VERTICES_T(v, G2, Graph) {
        std::cerr << G2.local_to_global(v) << std::endl;
    }
}

BOOST_AUTO_TEST_CASE(addEdge) {

    BOOST_TEST_MESSAGE("uninduced subgraph add edges");

    typedef uninduced_subgraph< adjacency_list< vecS, vecS, directedS,
        no_property, property< edge_index_t, int > > > Graph;
    //typedef Graph::edge_descriptor Edge;
    typedef Graph::vertex_descriptor Vertex;

    const int N = 3;
    Graph G0(N);
    Graph& G1 = G0.create_subgraph();
    Graph& G2 = G1.create_subgraph();
    
    BOOST_CHECK(&G1.parent() == &G0);
    BOOST_CHECK(&G2.parent() == &G1);
    
    // add vertices
    add_vertex(0, G2);
    add_vertex(1, G2);
    BOOST_CHECK(num_vertices(G1) == 2);
    BOOST_CHECK(num_vertices(G2) == 2);
    
    // add edge to G0 which needs propagation
    add_edge(0, 1, G0);
    
    BOOST_CHECK(num_edges(G0) == 1);
    BOOST_CHECK(num_edges(G1) == 0);
    BOOST_CHECK(num_edges(G2) == 0);
    // num_vertices stays the same
    BOOST_CHECK(num_vertices(G0) == 3);
    BOOST_CHECK(num_vertices(G1) == 2);
    BOOST_CHECK(num_vertices(G2) == 2);
    
    // add edge to G0 without propagation
    add_edge(1, 2, G0);
    
    BOOST_CHECK(num_edges(G0) == 2);
    BOOST_CHECK(num_edges(G1) == 0);
    BOOST_CHECK(num_edges(G2) == 0);
    // num_vertices stays the same
    BOOST_CHECK(num_vertices(G0) == 3);
    BOOST_CHECK(num_vertices(G1) == 2);
    BOOST_CHECK(num_vertices(G2) == 2);
    
    // add vertex 2 to G2/G1
    Vertex n = add_vertex(2, G2);
    BOOST_CHECK(G2.local_to_global(n) == 2);
    
    BOOST_CHECK(num_edges(G0) == 2);
    BOOST_CHECK(num_edges(G1) == 0);
    BOOST_CHECK(num_edges(G2) == 0);
    // num_vertices stays the same
    BOOST_CHECK(num_vertices(G0) == 3);
    BOOST_CHECK(num_vertices(G1) == 3);
    BOOST_CHECK(num_vertices(G2) == 3);
    
    std::cerr << "All G0 vertices: " << std::endl;
    BGL_FORALL_VERTICES_T(v, G0, Graph) {
        std::cerr << v << "/" << G0.local_to_global(v) << std::endl;
    }
    std::cerr << "All G1 vertices: " << std::endl;
    BGL_FORALL_VERTICES_T(v, G1, Graph) {
        std::cerr << v << "/" << G1.local_to_global(v) << std::endl;
    }
    std::cerr << "All G2 vertices: " << std::endl;
    BGL_FORALL_VERTICES_T(v, G2, Graph) {
        std::cerr << v << "/" << G2.local_to_global(v) << std::endl;
    }
    std::cerr << "All G0 edges: " << std::endl;
    BGL_FORALL_EDGES_T(e, G0, Graph) {
        std::cerr << source(e, G0) << "->" << target(e, G0) << std::endl;
    }
    std::cerr << "All G1 edges: " << std::endl;
    BGL_FORALL_EDGES_T(e, G1, Graph) {
        std::cerr << source(e, G1) << "->" << target(e, G1) << std::endl;
    }
    std::cerr << "All G2 edges: " << std::endl;
    BGL_FORALL_EDGES_T(e, G2, Graph) {
        std::cerr << source(e, G2) << "->" << target(e, G2) << std::endl;
    }
    
    // add edge to G2 with propagation upwards
    add_edge(0, 2, G2);
    
    BOOST_CHECK(num_edges(G0) == 3);
    BOOST_CHECK(num_edges(G1) == 1);
    BOOST_CHECK(num_edges(G2) == 1);
    // num_vertices stays the same
    BOOST_CHECK(num_vertices(G0) == 3);
    BOOST_CHECK(num_vertices(G1) == 3);
    BOOST_CHECK(num_vertices(G2) == 3);
    
    std::cerr << "All G0 vertices: " << std::endl;
    BGL_FORALL_VERTICES_T(v, G0, Graph) {
        std::cerr << v << "/" << G0.local_to_global(v) << std::endl;
    }
    std::cerr << "All G1 vertices: " << std::endl;
    BGL_FORALL_VERTICES_T(v, G1, Graph) {
        std::cerr << v << "/" << G1.local_to_global(v) << std::endl;
    }
    std::cerr << "All G2 vertices: " << std::endl;
    BGL_FORALL_VERTICES_T(v, G2, Graph) {
        std::cerr << v << "/" << G2.local_to_global(v) << std::endl;
    }
    std::cerr << "All G0 edges: " << std::endl;
    BGL_FORALL_EDGES_T(e, G0, Graph) {
        std::cerr << source(e, G0) << "->" << target(e, G0) << std::endl;
    }
    std::cerr << "All G1 edges: " << std::endl;
    BGL_FORALL_EDGES_T(e, G1, Graph) {
        std::cerr << source(e, G1) << "->" << target(e, G1) << std::endl;
    }
    std::cerr << "All G2 edges: " << std::endl;
    BGL_FORALL_EDGES_T(e, G2, Graph) {
        std::cerr << source(e, G2) << "->" << target(e, G2) << std::endl;
    }
}

BOOST_AUTO_TEST_CASE(addEdgeByDescriptor) {

    BOOST_TEST_MESSAGE("uninduced subgraph add edges by descriptor");

    typedef uninduced_subgraph< adjacency_list< vecS, vecS, directedS,
        no_property, property< edge_index_t, int > > > Graph;
    typedef Graph::edge_descriptor Edge;
    //typedef Graph::vertex_descriptor Vertex;

    const int N = 3;
    Graph G0(N);
    Graph& G1 = G0.create_subgraph();
    Graph& G2 = G1.create_subgraph();
    
    BOOST_CHECK(&G1.parent() == &G0);
    BOOST_CHECK(&G2.parent() == &G1);
    
    // add edge to G0
    add_edge(0, 1, G0);
    
    BOOST_CHECK(num_edges(G0) == 1);
    BOOST_CHECK(num_edges(G1) == 0);
    BOOST_CHECK(num_edges(G2) == 0);
    // num_vertices stays the same
    BOOST_CHECK(num_vertices(G0) == 3);
    BOOST_CHECK(num_vertices(G1) == 0);
    BOOST_CHECK(num_vertices(G2) == 0);
    
    Edge e_global;
    bool exists;
    
    boost::tie(e_global, exists) = edge(0, 1, G0);
    
    BOOST_REQUIRE(exists);
    // add edge to subgraph and propagate up and avoid multi-edge on root
    add_edge(e_global, G2);
    
    BOOST_CHECK(num_edges(G0) == 1);
    BOOST_CHECK(num_edges(G1) == 1);
    BOOST_CHECK(num_edges(G2) == 1);
    // num_vertices stays the same
    BOOST_CHECK(num_vertices(G0) == 3);
    BOOST_CHECK(num_vertices(G1) == 2);
    BOOST_CHECK(num_vertices(G2) == 2);
    
    
    std::cerr << "All G0 vertices: " << std::endl;
    BGL_FORALL_VERTICES_T(v, G0, Graph) {
        std::cerr << v << "/" << G0.local_to_global(v) << std::endl;
    }
    std::cerr << "All G1 vertices: " << std::endl;
    BGL_FORALL_VERTICES_T(v, G1, Graph) {
        std::cerr << v << "/" << G1.local_to_global(v) << std::endl;
    }
    std::cerr << "All G2 vertices: " << std::endl;
    BGL_FORALL_VERTICES_T(v, G2, Graph) {
        std::cerr << v << "/" << G2.local_to_global(v) << std::endl;
    }
    std::cerr << "All G0 edges: " << std::endl;
    BGL_FORALL_EDGES_T(e, G0, Graph) {
        std::cerr << source(e, G0) << "->" << target(e, G0) << std::endl;
    }
    std::cerr << "All G1 edges: " << std::endl;
    BGL_FORALL_EDGES_T(e, G1, Graph) {
        std::cerr << source(e, G1) << "->" << target(e, G1) << std::endl;
    }
    std::cerr << "All G2 edges: " << std::endl;
    BGL_FORALL_EDGES_T(e, G2, Graph) {
        std::cerr << source(e, G2) << "->" << target(e, G2) << std::endl;
    }
}

BOOST_AUTO_TEST_CASE(removeEdge) {

    BOOST_TEST_MESSAGE("uninduced subgraph remove edge");

    typedef uninduced_subgraph< adjacency_list< vecS, vecS, directedS,
        no_property, property< edge_index_t, int > > > Graph;
    typedef Graph::edge_descriptor Edge;
    //typedef Graph::vertex_descriptor Vertex;

    const int N = 3;
    Graph G0(N);
    Graph& G1 = G0.create_subgraph(vertices(G0).first, vertices(G0).second);
    Graph& G2 = G1.create_subgraph(vertices(G1).first, vertices(G1).second);
    
    BOOST_CHECK(&G1.parent() == &G0);
    BOOST_CHECK(&G2.parent() == &G1);
    
    BOOST_CHECK(num_vertices(G0) == 3);
    BOOST_CHECK(num_vertices(G1) == 3);
    BOOST_CHECK(num_vertices(G2) == 3);
    
    // add edge to G0
    add_edge(0, 1, G2);
    add_edge(1, 2, G2);
    
    BOOST_CHECK(num_edges(G0) == 2);
    BOOST_CHECK(num_edges(G1) == 2);
    BOOST_CHECK(num_edges(G2) == 2);
    // num_vertices stays the same
    BOOST_CHECK(num_vertices(G0) == 3);
    BOOST_CHECK(num_vertices(G1) == 3);
    BOOST_CHECK(num_vertices(G2) == 3);
    
    // remove edge from G1 by descriptor
    Edge e_local;
    bool exists;
    boost::tie(e_local, exists) = edge(0, 1, G1);
    
    BOOST_REQUIRE(exists);
    remove_edge(e_local, G1);
    
    BOOST_CHECK(num_edges(G0) == 2);
    BOOST_CHECK(num_edges(G1) == 1);
    BOOST_CHECK(num_edges(G2) == 1);
    
    // remove edge from G1
    remove_edge(1, 2, G1);
    
    BOOST_CHECK(num_edges(G0) == 2);
    BOOST_CHECK(num_edges(G1) == 0);
    BOOST_CHECK(num_edges(G2) == 0);
    
    std::cerr << "All G0 vertices: " << std::endl;
    BGL_FORALL_VERTICES_T(v, G0, Graph) {
        std::cerr << v << "/" << G0.local_to_global(v) << std::endl;
    }
    std::cerr << "All G1 vertices: " << std::endl;
    BGL_FORALL_VERTICES_T(v, G1, Graph) {
        std::cerr << v << "/" << G1.local_to_global(v) << std::endl;
    }
    std::cerr << "All G2 vertices: " << std::endl;
    BGL_FORALL_VERTICES_T(v, G2, Graph) {
        std::cerr << v << "/" << G2.local_to_global(v) << std::endl;
    }
    std::cerr << "All G0 edges: " << std::endl;
    BGL_FORALL_EDGES_T(e, G0, Graph) {
        std::cerr << source(e, G0) << "->" << target(e, G0) << std::endl;
    }
    std::cerr << "All G1 edges: " << std::endl;
    BGL_FORALL_EDGES_T(e, G1, Graph) {
        std::cerr << source(e, G1) << "->" << target(e, G1) << std::endl;
    }
    std::cerr << "All G2 edges: " << std::endl;
    BGL_FORALL_EDGES_T(e, G2, Graph) {
        std::cerr << source(e, G2) << "->" << target(e, G2) << std::endl;
    }
    
}

BOOST_AUTO_TEST_SUITE_END()
