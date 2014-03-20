/* This program reads secundary RNA structures in dot-bracket and
 * builds a graph for a latter ear-decomposition and bipartitness-check
 *
 * Created on: 25.03.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

// include header
#include "struct2graph.h"

bool debug = false;
bool no_bipartite_check = false;
bool ramachandran = false;
// filenames of graphml files will be starting with this string
std::string outfile = "";
unsigned long seed = 0;
int num_trees = 0;

// overload << operator to print vectors with any content

template <typename T>
std::ostream& operator<< (std::ostream& os, std::vector<T>& vec) {
  int i = 0;
  for (auto elem : vec) {
    os << "(" << i++ << ") " << elem << std::endl;
  }
  return os;
}

// overload << operator to print maps with any content

template <typename U, typename V>
std::ostream& operator<< (std::ostream& os, std::map<U, V>& m) {
  for (typename std::map<U, V>::iterator it = m.begin(); it != m.end(); it++) {
    os << it->first << "," << it->second << std::endl;
  }
  return os;
}

// lexmin implementation for pairs of any type

template <typename T>
std::pair<T, T>& lexmin (std::pair<T, T>& a, std::pair<T, T>& b) {
  if (a.first == b.first) {
    return !(b.second < a.second) ? a : b;
  } else {
    return !(b.first < a.first) ? a : b;
  }
}

//! main program starts here

int main (int ac, char* av[]) {

  // initialize command line options
  boost::program_options::variables_map vm = init_options(ac, av);

  // input handling ( we read from std:in per default and switch to a file if it is given in the --in option
  // std::in will provide a pseudo interface to enter structures directly!
  std::vector<std::string> structures;

  if (vm.count("in")) {
    if (debug) {
      std::cerr << "will read graphml file given in the options." << std::endl;
    }
    std::ifstream* infile = new std::ifstream(vm["in"].as<std::string>(), std::ifstream::in);
    if (infile->is_open()) {
      structures = read_input(infile);
      infile->close();
    } else {
      std::cerr << "Unable to open file";
      return 1;
    }
  } else {
    std::cerr << "Input structures (one per line); @ to quit" << std::endl
        << "....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8" << std::endl;
    structures = read_input(&std::cin); // read infile into array
  }

  std::ostream* out = &std::cout; // out stream

  // variables
  Graph graph = parse_graph(structures); // generate graph from input vector
  *out << "dependency graph:";
  print_graph(graph, out, "root-graph"); // print the graph as GML to a ostream

  decompose_graph(graph, out); // decompose the graph into its connected components, biconnected
  // components and decompose blocks via ear decomposition
  return 0;
}

boost::program_options::variables_map init_options (int ac, char* av[]) {
  // boost option parser
  // http://www.boost.org/doc/libs/1_53_0/doc/html/program_options/tutorial.html
  namespace po = boost::program_options;
  // Group of options that will be allowed only on command line
  po::options_description generic("Generic options");
  generic.add_options()
      ("help,h", "print help message")
      ("debug,v", po::value(&debug)->zero_tokens(), "be debug")
      ;

  // Group of options that will be allowed on command line and in a config file
  po::options_description config("Program options");
  config.add_options()
      ("in,i", po::value<std::string>(), "input file which contains the structures [string]")
      ("out,o", po::value<std::string>(&outfile), "write all (sub)graphs to gml files starting with given name [string]")
      ("seed,s", po::value<unsigned long>(&seed), "random number generator seed [unsigned long]")
      ("ramachandran,r", po::value(&ramachandran)->zero_tokens(), "Use the Ramachandran ear decomposition algorithmus")
      ("stat-trees,t", po::value<int>(&num_trees), "do decomposition statistics: define amount of different spanning trees for every root to calculate [int]")
      ("noBipartiteCheck,b", po::value(&no_bipartite_check)->zero_tokens(), "Don't check if input dependency graph is bipartite")
      ;

  po::positional_options_description p;
  p.add("in", 1).add("out", 2);

  po::options_description cmdline_options;
  cmdline_options.add(generic).add(config);

  po::variables_map vm;
  po::store(po::command_line_parser(ac, av).options(cmdline_options).positional(p).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << cmdline_options << "\n";
    exit(1);
  }
  if (vm.count("out")) {
    if (debug) {
      std::cerr << "graphml files will be written to file." << std::endl;
    }
  }

  return vm;
}

std::vector<std::string> read_input (std::istream* in) {
  // read input file
  std::string line;
  std::vector<std::string> structures;
  while (!in->eof()) {
    getline(*in, line);
    if (line == "@") {
      std::fclose(stdin);
    } else if (line.length() != 0) {
      structures.push_back(line);
    }
  }

  // exit if there is no input
  if (structures.empty()) {
    exit(1);
  }

  if (debug) {
    std::cerr << "Read following structures:" << std::endl;
  }
  // check if structures have equeal length
  unsigned int length = 0;
  for (auto elem : structures) {
    if (debug) {
      std::cerr << elem << std::endl;
    }
    if ((length != elem.length()) && (length != 0)) {
      std::cerr << "Structures have unequal length." << std::endl;
      exit(1);
    }
    length = elem.length();
  }
  return structures;
}

Graph parse_graph (std::vector<std::string> structures) {

  // count the number of positons
  int num_vertices = structures[0].length();
  if (debug) {
    std::cerr << "Generating Graph with " << num_vertices << " vertices." << std::endl;
  }
  Graph g(num_vertices);

  // give the vertices names
  int vertex_name = 0;
  Graph::vertex_iterator v, v_end;
  for (boost::tie(v, v_end) = boost::vertices(g); v != v_end; ++v) {
    boost::put(boost::vertex_color_t(), g, *v, vertex_name++);
  }

  // iterate over structures from input
  for (auto elem : structures) {
    std::vector<int> pair_table(structures[0].length(), 0); // remember position of the open bracket
    unsigned int pos = 0; // position of the character in the structure
    unsigned int open = 0; // remember how many open brackets there are
    // iterate over characters from structure
    while (pos < elem.length()) {
      if (elem[pos] == '(') {
        pair_table[open] = pos;
        if (debug) {
          std::cerr << elem[pos] << ", open count: " << open;
        }
        open++;
      } else if (elem[pos] == ')') {
        open--;
        // check if edge already exists
        bool exists_ab = boost::edge(boost::vertex(pair_table[open], g), boost::vertex(pos, g), g).second;
        bool exists_ba = boost::edge(boost::vertex(pos, g), boost::vertex(pair_table[open], g), g).second;
        if (!exists_ab && !exists_ba) {
          // add edge
          boost::add_edge(boost::vertex(pair_table[open], g), boost::vertex(pos, g), g);
        }
        // reset value
        pair_table[open] = pos;
        if (debug) {
          std::cerr << elem[pos] << ", open count: " << open;
        }
      } else if (elem[pos] != '.') {
        std::cerr << std::endl << "Unknown character in dot bracked representation" << std::endl;
        exit(1);
      }
      // error handling: there can't be more closing brackets than opening ones
      if (open < 0) {
        std::cerr << std::endl << "Unbalanced brackets in make_pair_table" << std::endl;
        exit(1);
      }
      if (debug) {
        std::cerr << " pos count:" << pos << std::endl;
      }
      pos++;
    }
    // error handling: at the end all brackets must be closed again!
    if (open != 0) {
      std::cerr << std::endl << "too few closed brackets in make_pair_table" << std::endl;
      exit(1);
    }
  }

  return g;
}

void print_graph (Graph& g, std::ostream* out, std::string nametag) {

  // print vertex and edge properties from my self defined bundled properties
  boost::dynamic_properties dp;
  dp.property("name", boost::get(boost::vertex_color_t(), g));
  dp.property("bipartite_color", boost::get(&vertex_property::bipartite_color, g));
  dp.property("color", boost::get(&vertex_property::color, g));
  dp.property("ear", boost::get(&edge_property::ear, g));

  if (outfile != "") {
    std::stringstream filename;
    filename << outfile << "-" << nametag << ".graphml";
    std::ofstream graphfile(filename.str());
    if (graphfile.is_open()) {
      boost::write_graphml(graphfile, g, dp, true);
      graphfile << std::endl;
      graphfile.close();
      *out << " done!" << std::endl;
    } else {
      std::cerr << " Unable to create graphml file!" << std::endl;
    }
  } else {
    *out << std::endl;
    boost::write_graphml(*out, g, dp, true);
    *out << std::endl;
    std::cerr << "created graphml!" << std::endl;
  }
}

void print_subgraphs (Graph& g, std::ostream* out, std::string nametag) {
  Graph::children_iterator ci, ci_end;
  int num = 1;
  for (boost::tie(ci, ci_end) = g.children(); ci != ci_end; ++ci) {
    std::stringstream name;
    name << nametag << "-" << num++;
    *out << name.str() << ":";
    print_graph(*ci, out, name.str());
  }
  if (debug) {
    std::cerr << "Printed all subgraphs." << std::endl;
  }
}

void decompose_graph (Graph& graph, std::ostream* out) {

  connected_components_to_subgraphs(graph); // get connected components and make subgraphs

  *out << "subgraphs connected components:" << std::endl;
  // print the just created subgraphs
  print_subgraphs(graph, out, "connected-component");

  // iterate over all subgraphs (connected components)
  Graph::children_iterator ci, ci_end;
  for (boost::tie(ci, ci_end) = graph.children(); ci != ci_end; ++ci) {
    if (!no_bipartite_check) {
      // check if subgraph is bipartite with a simple BFS
      // generate the vertex 0 as vertex_descriptor
      Vertex s = boost::vertex(0, *ci);
      // generate a edge_descriptor for the case that the graph is not bipartite
      Edge ed;
      if (!is_bipartite_graph(*ci, s, ed)) {
        std::cerr << "Graph is not bipartite! Conflict detected on edge " << ed << std::endl;
        exit(1);
      }
    }

    // calculate the max degree of this graph
    int max_degree = get_max_degree(*ci);
    if (debug) {
      std::cerr << "Max degree of subgraph is: " << max_degree << std::endl;
    }

    // split further into biconnected components do ear decomposition
    if (max_degree >= 3) {
      biconnected_components_to_subgraphs(*ci);

      *out << "subgraphs biconnected components:" << std::endl;
      // print the just created subgraphs
      print_subgraphs(*ci, out, "biconnected-component");

      Graph::children_iterator ci_b, ci_b_end;
      for (boost::tie(ci_b, ci_b_end) = (*ci).children(); ci_b != ci_b_end; ++ci_b) {
        // calculate the max degree of this graph
        int max_degree = get_max_degree(*ci_b);
        if (max_degree >= 3) {
          // only for statistics start at all vertices as root for DFS
          if (num_trees > 0) {
            do_spanning_tree_stat(*ci_b);
          } else {

            //TODO starting at 0 does not work atm. maybe underflow of unsigned int/vertex?
            //TODO use do_ear_decompositons (Schieber & Vishkin (1986)) or ear_decomposition (Ramachandran (1992)) one?!
            if (ramachandran) {
              ramachandran_ear_decomposition(*ci_b);
            } else {
              schieber_ear_decomposition(*ci_b);
            }
            *out << "subgraphs ear decomposition:" << std::endl;
            // print the just created subgraphs
            print_subgraphs(*ci_b, out, "decomposed-ear");
          }
        }
      }
    }
  }
}

void connected_components_to_subgraphs (Graph& g) {

  // get list of connected components into the component vector
  // http://www.boost.org/doc/libs/1_53_0/libs/graph/doc/connected_components.html
  // http://www.boost.org/doc/libs/1_53_0/libs/graph/example/connected_components.cpp
  std::vector<int> component(boost::num_vertices(g));
  int num = boost::connected_components(g, &component[0]);

  if (debug) {
    std::cerr << "Number of connected components: " << num << std::endl;
    std::cerr << component << std::endl;
  }

  // iterate over connected component numbers
  for (int i = 0; i != num; ++i) {
    // for each component number generate a new subgraph
    Graph& subg = g.create_subgraph();
    //boost::put(&graph_properties::level, g, "connected_component");
    int vertex = 0;
    // iterate over elements of connected_components
    for (auto elem : component) {
      if (i == elem) {
        // add vertex into current subgraph
        boost::add_vertex(vertex, subg);
      }
      vertex++;
    }
  }
}

void biconnected_components_to_subgraphs (Graph& g) {

  // // get list of biconnected components into the component property map
  // http://www.boost.org/doc/libs/1_53_0/libs/graph/doc/biconnected_components.html
  // http://www.boost.org/doc/libs/1_38_0/libs/graph/example/biconnected_components.cpp
  boost::edge_component_t edge_component;
  boost::property_map < Graph, boost::edge_component_t >::type component = boost::get(edge_component, g);
  unsigned int num = boost::biconnected_components(g, component);
  if (debug) {
    std::cerr << "Number of biconnected components: " << num << std::endl;
  }

  std::vector<Vertex> art_points;
  boost::articulation_points(g, std::back_inserter(art_points));
  if (debug) {
    std::cerr << "Number of articulation points: " << art_points.size() << " ( ";
    for (auto elem : art_points) {
      std::cerr << boost::get(boost::vertex_color_t(), g, elem) << " ";
    }
    std::cerr << ")" << std::endl;
  }

  if (debug) {
    // get graph and iterate over its edges to print connected components table
    typename Graph::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
      std::cerr << *ei << "\t" << "(" << boost::get(boost::vertex_color_t(), g, boost::source(*ei, g)) << ","
          << boost::get(boost::vertex_color_t(), g, boost::target(*ei, g)) << ")"
          << "\tcomponent: " << component[*ei] << std::endl;
    }
  }

  // now need to merge biconnected components that are separated by a articulation point that has a degree == 2 ?!

  // write biconnected components into subgraphs:
  for (unsigned int i = 0; i != num; i++) {
    // for each bicomponent number generate a new subgraph
    Graph& subg = g.create_subgraph();
    //boost::put(&graph_properties::level, g, "biconnected_component");
    // iterate over edges of graph
    //Graph rg = g.root();
    typename Graph::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
      if (i == component[*ei]) {
        // add vertex into current subgraph if not present already
        if (!subg.find_vertex(boost::get(boost::vertex_color_t(), g, boost::target(*ei, g))).second) {
          boost::add_vertex(boost::get(boost::vertex_color_t(), g, boost::target(*ei, g)), subg);
        }
        if (!subg.find_vertex(boost::get(boost::vertex_color_t(), g, boost::source(*ei, g))).second) {
          boost::add_vertex(boost::get(boost::vertex_color_t(), g, boost::source(*ei, g)), subg);
        }
      }
    }
  }

}

int get_max_degree (Graph& g) {
  int max_degree = 0;

  BGL_FORALL_VERTICES_T(v, g, Graph) {
    int current_degree = boost::out_degree(v, g);
    if (current_degree > max_degree) {
      max_degree = current_degree;
    }
  }
  return max_degree;
}

bool is_bipartite_graph (Graph& g, Vertex startVertex, Edge& ed) {
  // This is a Breadth First Search which checks if the graph is bipartit. 
  // If not, returns false and the fills the conflicting edge into the edge_descriptor

  if (debug) {
    std::cerr << "StartVertex is: " << startVertex << std::endl;
    std::cerr << "Number of vertices: " << boost::num_vertices(g) << std::endl;
  }

  // exit value (if bipartite = true, else false)
  bool exit = true;
  // Define A BGL visitor for the BFS algorithm

  class my_bfs_visitor : public boost::default_bfs_visitor {
  public:

    my_bfs_visitor (Edge& ed, bool& exit) : m_ed (ed), m_exit (exit) { }
    Edge& m_ed;
    bool& m_exit;

    enum {
      WHITE, BLACK, GRAY, RED
    };

    void tree_edge (Edge e, Graph g) const {
      if (debug) {
        std::cout << "Detecting Tree edge: " << e << std::endl;
      }
      Vertex u = boost::source(e, g);
      Vertex v = boost::target(e, g);
      if (g[u].bipartite_color == RED) {
        g[v].bipartite_color = BLACK;
      } else {
        g[v].bipartite_color = RED;
      }
    }

    void non_tree_edge (Edge e, Graph g) const {
      if (debug) {
        std::cout << "Detecting Non-Tree edge: " << e << std::endl;
      }
      Vertex u = boost::source(e, g);
      Vertex v = boost::target(e, g);
      if (g[u].bipartite_color == g[v].bipartite_color) {
        if (debug) {
          std::cerr << "u and v have the same color -> not bipartite!" << std::endl;
        }
        m_ed = boost::edge(u, v, g).first;
        // return false if graph is not bipartite
        m_exit = false;
      }
    }
  };

  my_bfs_visitor vis(ed, exit);
  // Do a BGL BFS!
  // http://www.boost.org/doc/libs/1_53_0/libs/graph/doc/breadth_first_search.html
  boost::breadth_first_search(g, startVertex, boost::visitor(vis));
  return exit;
}

void do_spanning_tree_stat (Graph& g) {
  // random generator to make spanning tree sampling
  if (seed == 0) {
    unsigned long clock_seed = std::chrono::system_clock::now().time_since_epoch().count();
    seed = clock_seed;
  }
  std::mt19937 r(seed); // mt19937 is a standard mersenne_twister_engine
  std::cerr << "Using this seed: " << seed << std::endl;

  // start at all vertices of the subgraph as root of the tree

  BGL_FORALL_VERTICES_T(root, g, Graph) {
    // generate num_trees spanning trees for statistics
    for (int i = 1; i != num_trees + 1; i++) {
      // remember predescessor map and all non-tree edges
      std::map<Vertex, Vertex> parents;
      std::vector<Edge> crossedges;

      // get a boost random spanning tree
      get_random_spanning_tree(g, r, parents, crossedges, root);

      // print parents, cross-edges and root vertex
      if (debug) {
        std::cerr << "Root vertex: " << root << std::endl;
        std::cerr << "Spanning tree (vertex, parent) and cross-edges:" << std::endl;
        for (std::map<Vertex, Vertex>::iterator it = parents.begin(); it != parents.end(); ++it) {
          std::cerr << it->first << " => " << it->second << std::endl;
        }
        for (auto elem : crossedges) {
          std::cerr << elem << std::endl;
        }
      }

      // do the actual ear decomposition
      ear_decomposition(g, parents, crossedges, root);

      // calculate the two performance critical variables alpha and beta
      // store attachment vertices per each step of the ear decomposition in Ak
      std::map<int, std::vector<Vertex> > Ak;
      unsigned int alpha;
      unsigned int beta;
      std::tie(alpha, beta) = calculate_alpha_beta(g, crossedges, Ak);

      // write statistics to a file
      print_ab_stat(alpha, beta, Ak, g, root, crossedges);
    }
  }
}

void schieber_ear_decomposition (Graph& g) {
  // remember predescessor map and all non-tree edges
  std::map<Vertex, Vertex> parents;
  std::vector<Edge> crossedges;

  // start Vertex
  Vertex startVertex = boost::vertex((boost::num_vertices(g) - 1), g);

  // random generator to make spanning tree sampling
  if (seed == 0) {
    unsigned long clock_seed = std::chrono::system_clock::now().time_since_epoch().count();
    seed = clock_seed;
  }
  std::mt19937 r(seed); // mt19937 is a standard mersenne_twister_engine
  std::cerr << "Using this seed: " << seed << std::endl;

  // get a boost random spanning tree
  get_random_spanning_tree(g, r, parents, crossedges, startVertex);

  // print parents, cross-edges and root vertex
  if (debug) {
    std::cerr << "Root vertex: " << startVertex << std::endl;
    std::cerr << "Spanning tree (vertex, parent) and cross-edges:" << std::endl;
    for (std::map<Vertex, Vertex>::iterator it = parents.begin(); it != parents.end(); ++it) {
      std::cerr << it->first << " => " << it->second << std::endl;
    }
    for (auto elem : crossedges) {
      std::cerr << elem << std::endl;
    }
  }

  // do the actual ear decomposition
  ear_decomposition(g, parents, crossedges, startVertex);

  // write ears into subgraphs
  int ears;

  BGL_FORALL_EDGES_T(e, g, Graph) {
    if (ears < g[e].ear) {
      ears = g[e].ear;
    }
  }

  for (int i = 0; i != ears + 1; i++) {
    Graph& subg = g.create_subgraph();
    //boost::put(&graph_properties::level, g, "decomposed_ears");
    // iterate over edges of graph
    //Graph rg = g.root();

    BGL_FORALL_EDGES_T(e, g, Graph) {
      if (i == g[e].ear) {
        // add vertex into current subgraph if not present already
        if (!subg.find_vertex(boost::get(boost::vertex_color_t(), g, boost::target(e, g))).second) {
          boost::add_vertex(boost::get(boost::vertex_color_t(), g, boost::target(e, g)), subg);
        }
        if (!subg.find_vertex(boost::get(boost::vertex_color_t(), g, boost::source(e, g))).second) {
          boost::add_vertex(boost::get(boost::vertex_color_t(), g, boost::source(e, g)), subg);
        }
      }
    }
  }
}

void ear_decomposition (Graph& g, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex start) {

  //delca saves (map of distance : (map of edge : lca))
  std::map<int, std::map<Edge, Vertex> > delca;

  // find the lca distances
  for (auto e : crossedges) {
    if (debug) {
      std::cerr << "starting at new chrossedge: " << e << std::endl;
    }
    std::pair<Vertex, int> lcad = get_lca_distance(g, parents, e, start);
    if (debug) {
      std::cerr << "lca " << lcad.first << " has distance " << lcad.second << std::endl;
    }
    delca[lcad.second][e] = lcad.first;
  }

  if (debug) {
    for (std::map<int, std::map<Edge, Vertex> >::iterator it = delca.begin(); it != delca.end(); ++it) {
      std::cerr << it->first << "(";
      for (std::map<Edge, Vertex>::iterator iit = (it->second).begin(); iit != (it->second).end(); ++iit) {
        std::cerr << "[" << iit->first << ", " << iit->second << "]";
      }
      std::cerr << ")" << std::endl;
    }
  }

  // reset ear_integer

  BGL_FORALL_EDGES_T(e, g, Graph) {
    g[e].ear = 0;
  }

  // now start at the biggest distance, at the biggest lexmin crossedge and make walks from the vertices to the lca
  // old values will be overwritten, therefore you get all the ears correctly
  int ear = 0;
  for (std::map<int, std::map<Edge, Vertex> >::reverse_iterator it = delca.rbegin(); it != delca.rend(); ++it) {
    for (std::map<Edge, Vertex>::reverse_iterator iit = (it->second).rbegin(); iit != (it->second).rend(); ++iit) {
      Vertex r = iit->second;
      g[iit->first].ear = ear;
      //std::cerr << "lca is: " << r << "; coloring ear: " << ear << std::endl;
      Vertex i = boost::source(iit->first, g);
      while (i != r) {
        std::map<Vertex, Vertex>::iterator iiit;
        iiit = parents.find(i);
        g[boost::edge(i, iiit->second, g).first].ear = ear;
        //std::cerr << "edge " << boost::edge(i,iiit->second,g).first << " will be ear " << ear << std::endl;
        i = iiit->second;
      }
      i = boost::target(iit->first, g);
      while (i != r) {
        std::map<Vertex, Vertex>::iterator iiit;
        iiit = parents.find(i);
        g[boost::edge(i, iiit->second, g).first].ear = ear;
        //std::cerr << "edge " << boost::edge(i,iiit->second,g).first << " will be ear " << ear << std::endl;
        i = iiit->second;
      }
      ear++;
    }
  }
}

void get_random_spanning_tree (Graph& g, std::mt19937& r, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex start) {

  boost::associative_property_map< std::map<Vertex, Vertex> > pm(parents);
  // call boost random spanning tree here:
  boost::random_spanning_tree(g, r, predecessor_map(pm).root_vertex(start));
  if (debug) {
    std::cerr << "Got a boost random spanning tree..." << std::endl;
  }

  // create the crossedges vector for the later ear-decomposition!
  // clear all values

  BGL_FORALL_EDGES_T(e, g, Graph) {
    g[e].color = 0;
  }
  // mark all tree edges
  for (std::map<Vertex, Vertex>::iterator it = parents.begin(); it != parents.end(); ++it) {
    if (it->first != start) {
      g[boost::edge(it->first, it->second, g).first].color = 1;
    }
  }
  // push all non-tree edges into crossedges

  BGL_FORALL_EDGES_T(e, g, Graph) {
    if (g[e].color == 0) {
      crossedges.push_back(e);
    }
  }
}

std::pair<Vertex, int> get_lca_distance (Graph& g, std::map<Vertex, Vertex>& parents, Edge e, Vertex r) {

  // make walks from the vertices to the root of the tree
  std::vector<Vertex> uwalk = make_tree_walk(parents, boost::target(e, g), r);
  std::vector<Vertex> vwalk = make_tree_walk(parents, boost::source(e, g), r);
  if (debug) {
    for (auto elem : uwalk)
      std::cerr << elem << "->";
    std::cerr << std::endl;
    for (auto elem : vwalk)
      std::cerr << elem << "->";
    std::cerr << std::endl;
  }
  // get the lca from the walks
  Vertex lca = boost::graph_traits<Graph>::null_vertex();
  int distance = -1;
  while (uwalk.back() == vwalk.back()) {
    lca = vwalk.back();
    distance++;
    uwalk.pop_back();
    vwalk.pop_back();
    if ((uwalk.size() == 0) || (vwalk.size() == 0)) {
      break;
    }
  }
  return std::make_pair(lca, distance);
}

std::vector<Vertex> make_tree_walk (std::map<Vertex, Vertex>& parents, Vertex v, Vertex r) {
  std::vector<Vertex> walk;
  Vertex i = v;
  walk.push_back(i);
  while (i != r) {
    std::map<Vertex, Vertex>::iterator it;
    it = parents.find(i);
    walk.push_back(it->second);
    i = it->second;
  }
  return walk;
}

std::pair< unsigned int, unsigned int > calculate_alpha_beta (Graph& g, std::vector<Edge>& crossedges, std::map<int, std::vector<Vertex> >& Ak) {

  // structure to remember Ak (attachment vertices)
  //std::map<int, std::vector<Vertex> > Ak;
  unsigned int alpha = 0;
  unsigned int beta = 0;
  int my = crossedges.size();

  // reset_color

  BGL_FORALL_VERTICES_T(v, g, Graph) {
    g[v].color = 0;
  }

  // iterate over all ear decomposition iterations and inside this loop over all edges of this ear
  // this edges have source and target vertices, iterate over both
  // then iterate over the adjacent out edges of these vertices and get the hightes ear number into maxear
  for (int k = 0; k != my; k++) {

    BGL_FORALL_EDGES_T(e, g, Graph) {
      if (g[e].ear == k) {
        std::vector<Vertex> adjacent_v;
        adjacent_v.push_back(boost::source(e, g));
        adjacent_v.push_back(boost::target(e, g));
        //TODO we calculate maxear twice for many vertices...
        for (auto adja : adjacent_v) {
          int maxear = 0;
          typename Graph::out_edge_iterator ei, ei_end;
          for (boost::tie(ei, ei_end) = boost::out_edges(adja, g); ei != ei_end; ++ei) {
            if (g[*ei].ear > maxear) {
              maxear = g[*ei].ear;
            }
          }

          if (maxear > k) {
            g[adja].color = 1;
          } else {
            g[adja].color = 0;
          }
        }
      }
    }

    //print_graph(g, &std::cout, "test_color");
    // remember Ak for this k
    std::vector< Vertex > thisAk;
    // write colored vertices into thisAk

    BGL_FORALL_VERTICES_T(v, g, Graph) {
      if (g[v].color == 1) {
        thisAk.push_back(v);
      }
    }
    Ak[k] = thisAk;
    if (thisAk.size() > alpha) {
      alpha = thisAk.size();
    }

    if (k > 0) {
      std::vector<Vertex> Akplus1_without_Ak;
      for (auto elem : thisAk) {
        std::vector<Vertex>::iterator iter = find(Ak[k - 1].begin(), Ak[k - 1].end(), elem);
        if (iter == Ak[k - 1].end()) {
          // elem is in thisAk but not in Ak[k-1]
          Akplus1_without_Ak.push_back(elem);
        }
      }
      unsigned int thisbeta = Ak[k - 1].size() + Akplus1_without_Ak.size();
      if (beta < thisbeta) {
        beta = thisbeta;
      }
    }
  }
  return std::make_pair(alpha, beta);
}

void print_ab_stat (unsigned int alpha, unsigned int beta, std::map<int, std::vector<Vertex> > Ak, Graph& g, Vertex root, std::vector<Edge>& crossedges) {

  int my = crossedges.size();

  // write statistic output file
  std::stringstream filename;
  if (outfile != "") {
    filename << outfile << "-statistics.txt";
  } else {
    filename << "statistics.txt";
  }
  std::ofstream statfile(filename.str(), std::ofstream::out | std::ofstream::app);
  if (statfile.is_open()) {
    // print alpha and beta
    statfile << alpha << " " << beta << " ";
    // print Ak sets
    for (int i = 0; i != my; i++) {
      for (auto elem : Ak[i]) {
        if (elem != Ak[i].back()) {
          statfile << boost::get(boost::vertex_color_t(), g, elem) << ",";
        } else {
          statfile << boost::get(boost::vertex_color_t(), g, elem);
        }
      }
      statfile << " ";
    }
    // print root of spanning tree
    statfile << root << " ";
    // also print crossedges to reproduce trees afterwards
    for (auto elem : crossedges) {
      statfile << elem;
    }
    statfile << std::endl;
    statfile.close();
    if (debug) {
      std::cerr << "Statistics written to outfile!" << std::endl;
    }
  } else {
    std::cerr << " Unable to create statistics file!" << std::endl;
  }
}

void ramachandran_ear_decomposition (Graph& g) {
  // blocks need to be decomposed into path. this can be done by Ear Decomposition

  // start Vertex
  Vertex startVertex = boost::vertex((boost::num_vertices(g) - 1), g);

  // map of ear decomposition properties for all vertices as key
  ear_propertymap_t p;
  // map of ear structure.
  ear_t ear;
  // time starts at 0
  unsigned int counter = 0;

  if (debug) {
    std::cout << "StartVertex is: " << startVertex << std::endl;
  }
  // Algorithm from Ramachandran (1992) Parallel Open Ear Decomposition with Applications, page 8/9
  ear_dfs(g, startVertex, p, ear, counter);

  // print out all data-structures at the end
  if (debug) {
    std::cerr << "index\tcolor\tporder\tparent\tlow\tear" << std::endl;
    for (ear_propertymap_t::iterator it = p.begin(); it != p.end(); it++) {
      std::cerr << it->first << "\t" <<
          it->second.color << "\t" <<
          it->second.preorder << "\t" <<
          it->second.parent << "\t" <<
          it->second.low << "\t" <<
          it->second.ear.first << "," <<
          it->second.ear.second << std::endl;
    }
    std::cerr << "index\tear" << std::endl;
    for (ear_t::iterator it = ear.begin(); it != ear.end(); it++) {
      std::cerr << it->first.first << "," <<
          it->first.second << "\t" <<
          it->second.first << "," <<
          it->second.second << std::endl;
    }
    std::cerr << "counter: " << counter << std::endl;
  }

  // create subgraphs from decomposed ears
  std::vector<edge_t> found;
  for (ear_t::iterator it = ear.begin(); it != ear.end(); it++) {
    if (!(std::find(found.begin(), found.end(), it->second) != found.end())) {
      // for each ear create a new subgraph
      Graph& subg = g.create_subgraph();
      //boost::put(&graph_properties::level, g, "decomposed_ears");
      if (debug) {
        std::cerr << "New subgraph for ear (" << it->second.first << "," << it->second.second << ")" << std::endl
            << "Vertices will be included in subgraph: ";
      }
      for (ear_t::iterator ti = it; ti != ear.end(); ti++) {
        if (ti->second == it->second) {
          found.push_back(ti->second);
          // add vertex into current subgraph if not present already
          if (!subg.find_vertex(ti->first.first).second) {
            boost::add_vertex(ti->first.first, subg);
            if (debug) {
              std::cerr << " " << ti->first.first;
            }
          }
          if (!subg.find_vertex(ti->first.second).second) {
            boost::add_vertex(ti->first.second, subg);
            if (debug) {
              std::cerr << " " << ti->first.second;
            }
          }
        }
      }
      if (debug) {
        std::cerr << std::endl;
      }
    }
  }
}

void ear_dfs (Graph& g, Vertex v, ear_propertymap_t& p, ear_t& ear, unsigned int& counter) {

  enum {
    WHITE, BLACK, GRAY
  };
  if (debug) {
    std::cout << "v is: " << v << std::endl;
  }

  // start ear decomposition
  p[v].color = GRAY;
  p[v].preorder = counter;
  counter++;
  p[v].low = boost::num_vertices(g);
  p[v].ear = std::make_pair(boost::num_vertices(g), boost::num_vertices(g));

  // get neighbouring vertices
  typename Graph::out_edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = boost::out_edges(v, g); ei != ei_end; ++ei) {
    if (debug) {
      std::cerr << boost::target(*ei, g) << " is neighbour through edge: " << *ei << std::endl;
    }
    Vertex w = boost::target(*ei, g);
    if (debug) {
      std::cout << "w is: " << w << std::endl;
    }

    if (p[w].color == WHITE) {
      if (debug) {
        std::cout << "w is white" << std::endl;
      }
      p[w].parent = v;
      // start new iteration here
      ear_dfs(g, w, p, ear, counter);
      //TODO: cast low vertex to integer a good idea?
      if ((int) p[w].low >= p[w].preorder) {
        ear[std::make_pair(p[w].parent, w)] = std::make_pair(std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
      } else {
        ear[std::make_pair(p[w].parent, w)] = p[w].ear;
      }
      p[v].low = std::min((int) p[v].low, (int) p[w].low);
      p[v].ear = lexmin(p[v].ear, p[w].ear);
    } else if (p[w].color == GRAY) {
      if (debug) {
        std::cout << "w is gray" << std::endl;
      }
      if (w != p[w].parent) {
        if (debug) {
          std::cout << "found a crossedge: " << v << w << std::endl;
        }
        //TODO: casting vertex in low to integer a bad idea?
        p[v].low = boost::vertex(std::min((int) p[v].low, p[w].preorder), g);
        ear[std::make_pair(w, v)] = std::make_pair(boost::vertex(p[w].preorder, g), boost::vertex(p[v].preorder, g));
        p[v].ear = lexmin(p[v].ear, ear[std::make_pair(w, v)]);
      }
    }
  }
  if (debug) {
    std::cout << "finishing vertex " << v << std::endl;
  }
}


