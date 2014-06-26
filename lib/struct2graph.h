/* This program reads secundary RNA structures in dot-bracket and
 * builds a graph for a latter ear-decomposition and bipartitness-check
 *
 * Created on: 25.03.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 * Compile with: g++ -std=c++11 -g -lboost_program_options -o struct2graph struct2graph.cc
 */

#ifndef STRUCT2GRAPH_H
#define STRUCT2GRAPH_H

// include standard library parts
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <array>
#include <limits>
#include <chrono>

// include boost components
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/random_spanning_tree.hpp>
#include <boost/graph/iteration_macros.hpp>


namespace design {
  namespace detail {
    // get property with Graph[Vertex].bipartite_color = int;

    struct vertex_property {
      int bipartite_color;
      int color;
    };

    struct edge_property {
      int ear;
      int color;
    };

    namespace boost {

      struct edge_component_t {

        enum {
          num = 555
        };
        typedef edge_property_tag kind;
      };
    }

    //, graph_properties 
    // boost graph template
    typedef boost::adjacency_list_traits< boost::vecS, boost::vecS, boost::undirectedS > Traits;
    typedef boost::subgraph< boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS,
    boost::property< boost::vertex_color_t, int, vertex_property >,
    boost::property< boost::edge_index_t, int, boost::property < boost::edge_component_t, std::size_t, edge_property> > > > Graph;
    typedef Graph::edge_descriptor Edge;
    typedef Graph::vertex_descriptor Vertex;

    // initialise boost command line option parser
    boost::program_options::variables_map init_options (int ac, char* av[]);

    // read the input file into a string
    std::vector<std::string> read_input (std::istream* in);

    // parse the input string into a graph
    Graph parse_graph (std::vector<std::string> structures);

    // print the graph as a gml file to the output
    void print_graph (Graph& g, std::ostream* out, std::string nametag);

    // print all the subgraphs as GML (iterator over subgraphs)
    void print_subgraphs (Graph& g, std::ostream* out, std::string nametag);

    // does the whole graph decomposition, uses methods below
    void decompose_graph (Graph& graph, std::ostream* out);

    // get a vector of all vertices with their component id. finds connected components with DFS
    void connected_components_to_subgraphs (Graph& g);

    // finds connected components with DFS
    void biconnected_components_to_subgraphs (Graph& g);

    // get max degree of a graph
    int get_max_degree (Graph& g);

    // do a Breadth First Search to test for bipartite property
    bool is_bipartite_graph (Graph& g, Vertex startVertex, Edge& ed);

    // do statistics with many spannign trees calculating alpha and beta values using the schieber algorithm
    void do_spanning_tree_stat (Graph& g);

    // do a schieber ear decomposition without any statistics
    void schieber_ear_decomposition (Graph& g);

    // implementation of the schieber ear decomposition
    void ear_decomposition (Graph& g, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex start);

    // get boost random spanning tree for schieber ear decomposition
    void get_random_spanning_tree (Graph& g, std::mt19937& r, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex start);

    // given all the parents of a spanning_tree, find the lca a crossedge (helper for schieber ear decomposition)
    std::pair<Vertex, int> get_lca_distance (Graph& g, std::map<Vertex, Vertex>& parents, Edge e, Vertex r);

    // do a walk in the spanning tree starting at v and ending at the root r -> return a vector with the walk  (helper for schieber ear decomposition)
    std::vector<Vertex> make_tree_walk (std::map<Vertex, Vertex>& parents, Vertex v, Vertex r);

    // calculate spanning tree properties for coloring
    std::pair< unsigned int, unsigned int > calculate_alpha_beta (Graph& g, std::vector<Edge>& crossedges, std::map<int, std::vector<Vertex> >& Ak);

    // write statistic output file
    void print_ab_stat (unsigned int alpha, unsigned int beta, std::map<int, std::vector<Vertex> > Ak, Graph& g, Vertex root, std::vector<Edge>& crossedges);

    // typedefs for ramachandran ear decomposition
    typedef std::pair<Vertex, Vertex> edge_t;
    typedef std::map<edge_t, edge_t> ear_t;
    // struct to remember coloring, time, parents of a vertex

    struct property {
      int color;
      int preorder;
      Vertex parent;
      Vertex low;
      edge_t ear;
    };
    typedef std::map<Vertex, property> ear_propertymap_t;

    // ear decomposition of blocks
    void ramachandran_ear_decomposition (Graph& g);

    // actual dfs for ramachandran ear decomposition (dfs algroithm)
    void ear_dfs (Graph& g, Vertex v, ear_propertymap_t& p, ear_t& ear, unsigned int& counter);
  }
}


#endif
