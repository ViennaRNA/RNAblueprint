/* This file is a boost.test unit test and provides tests for decompose.cc
 *
 *
 * Created on: 23.10.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 * 
 *
 */

// include header
#include "test_common.h"

// include headers containing functions to test
#include "decompose.h"
#include "printgraph.h"

// include boost components
#include <boost/graph/iteration_macros.hpp>

// define heads
using namespace design;
using namespace design::detail;

BOOST_AUTO_TEST_SUITE(Decompose)

enum type {
    CC, BC, ED
};

Graph createGraph(int type) {
    Graph g(9);
    int vertex_name = 0;

    BGL_FORALL_VERTICES_T(v, g, Graph) {
        boost::put(boost::vertex_color_t(), g, v, vertex_name++);
    }

    // create a path between 0 and 10
    for (unsigned int i = 0; i < 7; i++) {
        boost::add_edge(boost::vertex(i, g), boost::vertex(i + 1, g), g);
    }

    // create a block
    boost::add_edge(boost::vertex(0, g), boost::vertex(8, g), g);
    boost::add_edge(boost::vertex(0, g), boost::vertex(7, g), g);
    boost::add_edge(boost::vertex(2, g), boost::vertex(8, g), g);
    boost::add_edge(boost::vertex(6, g), boost::vertex(8, g), g);

    if ((type == BC) || (type == CC)) {
        // add a biconnected component
        for (int i = 9; i <= 11; i++) {
            Vertex v = boost::add_vertex(g);
            boost::put(boost::vertex_color_t(), g, v, i);
        }
        boost::add_edge(boost::vertex(7, g), boost::vertex(9, g), g);
        boost::add_edge(boost::vertex(9, g), boost::vertex(10, g), g);
        boost::add_edge(boost::vertex(10, g), boost::vertex(11, g), g);
    }

    if (type == CC) {
        // add a connected component
        Vertex v = boost::add_vertex(g);
        boost::put(boost::vertex_color_t(), g, v, 12);
        v = boost::add_vertex(g);
        boost::put(boost::vertex_color_t(), g, v, 13);
        boost::add_edge(boost::vertex(12, g), boost::vertex(13, g), g);
    }

    /* full graph looks like this now (block plus path as biconnected component and connected component):
      5---6---7---9---10---11
      |   |   |
      4   8---0
      |   |   |
      3---2---1       12---13
     */
    return g;
}

std::unordered_set<int> getVertexSet(Graph& g) {
    std::unordered_set<int> result;

    BGL_FORALL_VERTICES(v, g, Graph) {
        result.insert((int) boost::get(boost::vertex_color_t(), g, v));
    }
    return result;
}

BOOST_AUTO_TEST_CASE(connectedComponents) {
    
    initialize_library(true);
    // create a graph
    Graph g = createGraph(CC);
    BOOST_TEST_MESSAGE("decompose connected components");
    connected_components_to_subgraphs(g);

    int number_of_children = 0;

    Graph::children_iterator child, child_end;
    for (boost::tie(child, child_end) = g.children(); child != child_end; ++child) {
        number_of_children++;
        // for the smaller connected component (12---13)
        if (boost::num_vertices(*child) == 2) {
            // check if both vertices exist and are labeled right
            BOOST_CHECK(boost::get(boost::vertex_color_t(), *child, 0) == 12);
            BOOST_CHECK(boost::get(boost::vertex_color_t(), *child, 1) == 13);
            // check if just one edge exists here
            BOOST_CHECK(boost::num_edges(*child) == 1);
        }
    }
    // check if it is just 2 connected components
    BOOST_CHECK(number_of_children == 2);
}

BOOST_AUTO_TEST_CASE(biconnectedComponents) {
    
    initialize_library(true);
    // get graph
    Graph g = createGraph(BC);
    // print_graph(g, &std::cout, "graph");
    BOOST_TEST_MESSAGE("decompose biconnected components 1");
    biconnected_components_to_subgraphs(g);

    int number_of_children = 0;

    Graph::children_iterator child, child_end;
    for (boost::tie(child, child_end) = g.children(); child != child_end; ++child) {
        number_of_children++;
        // print_graph(*child, &std::cout, "bc");
        // for the smaller connected component (12---13)
        if (boost::num_vertices(*child) == 4) {
            std::unordered_set<int> testcase{7, 9, 10, 11};
            BOOST_CHECK(getVertexSet(*child) == testcase);
            // check if just one edge exists here
            BOOST_CHECK(boost::num_edges(*child) == 3);
        }
    }
    // check if it is just 2 connected components
    BOOST_CHECK(number_of_children == 2);

    // ---------------------------------------------------------------
    // do another one with 2 paths connected on the block
    Graph h = createGraph(BC);
    for (int i = 12; i <= 13; i++) {
        Vertex v = boost::add_vertex(h);
        boost::put(boost::vertex_color_t(), h, v, i);
    }
    boost::add_edge(boost::vertex(0, h), boost::vertex(12, h), h);
    boost::add_edge(boost::vertex(12, h), boost::vertex(13, h), h);

    BOOST_TEST_MESSAGE("decompose biconnected components 2");
    biconnected_components_to_subgraphs(h);

    number_of_children = 0;

    for (boost::tie(child, child_end) = h.children(); child != child_end; ++child) {
        number_of_children++;
        // print_graph(*child, &std::cout, "bc");
        // for the smaller connected component (12---13)
        if (boost::num_vertices(*child) == 4) {
            std::unordered_set<int> testcase{7, 9, 10, 11};
            BOOST_CHECK(getVertexSet(*child) == testcase);
            // check if just one edge exists here
            BOOST_CHECK(boost::num_edges(*child) == 3);
        } else if (boost::num_vertices(*child) == 3) {
            std::unordered_set<int> testcase{0, 12, 13};
            BOOST_CHECK(getVertexSet(*child) == testcase);
            // check if just one edge exists here
            BOOST_CHECK(boost::num_edges(*child) == 2);
        }
    }
    // check if it is just 3 connected components
    BOOST_CHECK(number_of_children == 3);
}

BOOST_AUTO_TEST_CASE(EarDecomposition) {

    // create a graph
    Graph g = createGraph(ED);
    BOOST_TEST_MESSAGE("decompose connected components");
    // set random generator to a static seed;
    initialize_library(true);
    std::mt19937 mt(1);
    ear_decomposition_to_subgraphs(g, mt);

    int number_of_children = 0;

    Graph::children_iterator child, child_end;
    for (boost::tie(child, child_end) = g.children(); child != child_end; ++child) {
        number_of_children++;
        std::unordered_set<int> testcase;
        
        std::pair<Graph::edge_iterator, Graph::edge_iterator> ei = edges(g);
        int ear = (*child)[*ei.first].ear;
        switch (ear) {
            case 1:
                print_graph(*child, &std::cout, "ear-1");
                BOOST_CHECK(boost::num_vertices(*child) == 2);
                // check if both vertices exist and are labeled right
                testcase = {6, 8};
                BOOST_CHECK(getVertexSet(*child) == testcase);
                // check if just two edges exists here
                BOOST_CHECK(boost::num_edges(*child) == 1);
                // check if this is the ear number 0 for all edges

                BGL_FORALL_EDGES(e, *child, Graph) {
                    BOOST_CHECK((*child)[e].ear == 1);
                }
                break;
            case 2:
                print_graph(*child, &std::cout, "ear-2");
                BOOST_CHECK(boost::num_vertices(*child) == 7);
                // check if both vertices exist and are labeled right
                testcase = { 2, 3, 4, 5, 6, 7, 0 };
                BOOST_CHECK(getVertexSet(*child) == testcase);
                // check if just 6 edges exists here
                BOOST_CHECK(boost::num_edges(*child) == 6);
                // check if this is the ear number 2 for all edges

                BGL_FORALL_EDGES(e, *child, Graph) {
                    BOOST_CHECK((*child)[e].ear == 2);
                }
                break;
            case 3:
                print_graph(*child, &std::cout, "ear-3");
                BOOST_CHECK(boost::num_vertices(*child) == 4);
                // check if both vertices exist and are labeled right
                testcase = {0, 1, 2, 8};
                BOOST_CHECK(getVertexSet(*child) == testcase);
                // check if just one edge exists here
                BOOST_CHECK(boost::num_edges(*child) == 4);
                // check if this is the ear number 2 for all edges

                BGL_FORALL_EDGES(e, *child, Graph) {
                    BOOST_CHECK((*child)[e].ear == 3);
                }
                break;
        }
    }
    // check if it is just 3 connected components
    BOOST_CHECK(number_of_children == 3);

    // check if Aks are labeled as special for those two vertices
    std::vector<int> AKs = {0, 2, 6, 8};
    for (auto s : AKs) {
        BOOST_CHECK(g[g.find_vertex(boost::get(boost::vertex_color_t(), g, s)).first].special);
    }
    // check if others are not
    std::vector<int> noAKs = {1, 3, 4, 5, 7};
    for (auto s : noAKs) {
        BOOST_CHECK(!g[g.find_vertex(boost::get(boost::vertex_color_t(), g, s)).first].special);
    }
}

BOOST_AUTO_TEST_CASE(partsBetweenArticulationPoints) {

    // create a graph
    Graph g = createGraph(ED);
    BOOST_TEST_MESSAGE("parts between attachment points to subgraphs");
    // set random generator to a static seed;
    initialize_library(true);
    std::mt19937 mt(1);
    ear_decomposition_to_subgraphs(g, mt);
    print_subgraphs(g, &std::cout, "ears");
    // now the decomposed graph looks different than above, as the random generator went further
    // However, lets call our method on all children:
    int number_of_children = 0;
    Graph::children_iterator ear, ear_end;
    for (boost::tie(ear, ear_end) = g.children(); ear != ear_end; ++ear) {
        parts_between_specials_to_subgraphs(*ear);
        // print_subgraphs(*ear, &std::cout, "ear-parts");
        // there should be sugraphs for all ears now representing the paths and cycles between special nodes (like attatchment points and constraints)
        // or in case of no Ak or Ai, the whole ear is another sugraph
        std::pair<Graph::edge_iterator, Graph::edge_iterator> ei = edges(*ear);
        int k = (*ear)[*ei.first].ear;
        
        Graph::children_iterator child, child_end;
        if (k == 1) {
            boost::tie(child, child_end) = (*ear).children();
            // Ear #1 should not have any children!
            BOOST_CHECK(child == child_end);
        }
        // rest check here:
        for (boost::tie(child, child_end) = (*ear).children(); child != child_end; ++child) {

            // Now check if everything worked out!
            number_of_children++;
            std::unordered_set<int> testcase;
            
            switch (k) {
                case 2:
                    // for ear number 2 there should be two children, one with 3 nodes
                    if (boost::num_vertices(*child) == 3) {
                        // check if both vertices exist and are labeled right
                        testcase = { 6, 7, 0 };
                        BOOST_CHECK(getVertexSet(*child) == testcase);
                        // check if just one edge exists here
                        BOOST_CHECK(boost::num_edges(*child) == 2);
                    // and one with 5 nodes
                    } else if (boost::num_vertices(*child) == 5) {
                        testcase = { 2, 3, 4, 5, 6 };
                        BOOST_CHECK(getVertexSet(*child) == testcase);
                        // check if just one edge exists here
                        BOOST_CHECK(boost::num_edges(*child) == 4);
                    } else {
                        BOOST_FAIL("Found non expected subgraph!");
                    }
                    // check if this is the ear number 2 for all edges
                    BGL_FORALL_EDGES(e, *child, Graph) {
                        BOOST_CHECK((*child)[e].ear == 2);
                    }
                    break;
                case 3:
                    // for ear number 3 (cycle) there should be three children
                    testcase = {0, 1, 2, 8};
                    if (boost::num_vertices(*child) == 2) {
                        std::unordered_set<int> check = getVertexSet(*child);
                        // check if either vertex 2, 8 or 0, 8 are present
                        BOOST_CHECK((check.find(8) != check.end()) && ((check.find(0) != check.end()) != (check.find(2) != check.end())));
                        // check if just one edge exists here
                        BOOST_CHECK(boost::num_edges(*child) == 1);
                    } else if (boost::num_vertices(*child) == 3) {
                        testcase = {0, 1, 2};
                        BOOST_CHECK(getVertexSet(*child) == testcase);
                        // only two edges present here?
                        BOOST_CHECK(boost::num_edges(*child) == 2);
                    }

                    BGL_FORALL_EDGES(e, *child, Graph) {
                        BOOST_CHECK((*child)[e].ear == 3);
                    }
            }
        }
    }

    // check if it is exactly 4 subgraphs in total
    BOOST_CHECK(number_of_children == 5);
}

BOOST_AUTO_TEST_CASE(simpleCircle) {
    Graph g(4);
    int vertex_name = 0;

    BGL_FORALL_VERTICES_T(v, g, Graph) {
        boost::put(boost::vertex_color_t(), g, v, vertex_name++);
    }
    // create a path between 0 and 3 and circle closure
    for (unsigned int i = 0; i < 3; i++) {
        boost::add_edge(boost::vertex(i, g), boost::vertex(i + 1, g), g);
    }
    boost::add_edge(boost::vertex(0, g), boost::vertex(3, g), g);
    
    // decompose into subgraphs
    std::mt19937 mt(1);
    decompose_recursion(g, mt);
    // check
    int number_of_children = 0;
    
    Graph::children_iterator child, child_end;
    for (boost::tie(child, child_end) = g.children(); child != child_end; ++child) {
        number_of_children++;
        // subgraphs should only be paths
        int max_degree;
        int min_degree;
        std::tie(min_degree, max_degree) = get_min_max_degree(*child);
        BOOST_CHECK(max_degree <= 2);
        BOOST_CHECK(min_degree == 1);
    }
    
    BOOST_CHECK(number_of_children == 2);
}

BOOST_AUTO_TEST_SUITE_END()
