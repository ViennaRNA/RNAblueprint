// outttakes from struct2graph.cc

// get the spanning tree of our graph
//get_spanning_tree(g, parents, crossedges, startVertex);
//change_spanning_tree(g, r, parents, crossedges, startVertex);


// get spanning tree with DFS
void get_spanning_tree (Graph& g, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex start);

// change spanning tree with sampling a random new one
void change_spanning_tree (Graph& g, std::mt19937& r, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex& root);

void get_spanning_tree (Graph& g, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex start) {

  class my_dfs_visitor : public boost::default_dfs_visitor {
  public:

    my_dfs_visitor (std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges) : p (parents), c (crossedges) { }
    std::map<Vertex, Vertex>& p;
    std::vector<Edge>& c;

    enum {
      WHITE, BLACK, GRAY, RED
    };

    void start_vertex (Vertex s, Graph g) const {
      if (debug) {
        std::cerr << "Start vertex: " << s << std::endl;
      }
    }

    void tree_edge (Edge e, Graph g) const {
      if (debug) {
        std::cerr << "Detecting tree-edge: " << e << std::endl;
      }
      Vertex u = boost::source(e, g);
      Vertex v = boost::target(e, g);
      p[v] = u;
    }

    void forward_or_cross_edge (Edge e, Graph g) const {
      if (debug) {
        std::cerr << "Detecting back-edge: " << e << std::endl;
      }
      //Vertex u = boost::source(e, g);
      //Vertex v = boost::target(e, g);
      c.push_back(e);
    }
  };

  my_dfs_visitor vis(parents, crossedges);

  // Do a BGL DFS!
  // http://www.boost.org/doc/libs/1_53_0/libs/graph/doc/depth_first_search.html
  // Did not work: http://www.boost.org/doc/libs/1_53_0/libs/graph/doc/undirected_dfs.html
  // boost::undirected_dfs(g, boost::visitor(vis), vcolorMap, ecolorMap, rootVertex);
  boost::depth_first_search(g, visitor(vis).root_vertex(start));
}

void change_spanning_tree (Graph& g, std::mt19937& r, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex& start) {

  // take random edge <= tree edges
  // take random old crossedge <= crossedges 
  // calculate lca of old crossedge
  // delete edge in parents
  // add edge to crossedge
  // delete old crossedge in crossedges
  // (add parent of old crossedge and update parents):
  //	walk from old crossedge vertices up to lca (remembering path in vertex)
  //		if reach lca -> do nothing
  //		else if reach one (or is) new crossedge vertex -> invert parents + add new parent for old crossedge
  // return with new parents and crossedges

  // get random crossedgeedge
  std::uniform_int_distribution<int> rand_crossedge(0, crossedges.size() - 1); //(min, max)
  Edge old_crossedge = crossedges[rand_crossedge(r)];
  if (debug) {
    std::cerr << "choosen crossedge: " << old_crossedge << std::endl;
  }
  // calculate lca of old crossedge
  Vertex lca = get_lca_distance(g, parents, old_crossedge, start).first;
  if (debug) {
    std::cerr << "lca from old crossedge: " << lca << std::endl;
  }
  // get current cycle
  std::vector<Vertex> source_walk = make_tree_walk(parents, boost::source(old_crossedge, g), lca);
  std::vector<Vertex> target_walk = make_tree_walk(parents, boost::target(old_crossedge, g), lca);
  // merge cycle paths and delete lca
  target_walk.erase(--target_walk.end());
  source_walk.erase(--source_walk.end());
  for (auto elem : source_walk) {
    target_walk.push_back(elem);
  }
  if (debug) {
    std::cerr << "cycle is:" << std::endl << target_walk << std::endl;
  }
  // get random edge from cycle
  std::uniform_int_distribution<int> rand_tree_edge(0, target_walk.size() - 1); //(min, max)
  int x = rand_tree_edge(r);
  Edge edge = boost::edge(target_walk[x], parents[target_walk[x]], g).first;
  if (debug) {
    std::cerr << "Random Tree Edge is:" << edge << std::endl;
  }
  // erase source target pair of edge in parents
  for (std::map<Vertex, Vertex>::iterator it = parents.begin(); it != parents.end(); ++it) {
    if (((it->second == boost::source(edge, g)) && (it->first == boost::target(edge, g))) ||
        ((it->second == boost::target(edge, g)) && (it->first == boost::source(edge, g)))) {
      parents.erase(it);
    }
  }
  // add edge to crossedge
  crossedges.push_back(edge);
  // delete old crossedge
  std::vector<Edge>::iterator it = std::find(crossedges.begin(), crossedges.end(), old_crossedge);
  crossedges.erase(it);
  // walk from both vertices to add parent of old crossedge and update parents
  std::vector<Vertex> adjacent_v;
  adjacent_v.push_back(boost::source(old_crossedge, g));
  adjacent_v.push_back(boost::target(old_crossedge, g));
  Vertex other = adjacent_v[1];
  for (auto v : adjacent_v) {
    if (debug) {
      std::cerr << "walk from vertex " << v << std::endl;
    }
    if (debug) {
      std::cerr << "other is " << other << std::endl;
    }
    std::vector<Vertex> walk;
    Vertex i = v;
    walk.push_back(i);
    while (i != lca) {
      std::map<Vertex, Vertex>::iterator it;
      it = parents.find(i);
      if (it != parents.end()) {
        walk.push_back(it->second);
      }
      //if reach one (or is) new crossedge vertex -> invert parents + add new parent for old crossedge
      if ((i == boost::source(edge, g)) || (i == boost::target(edge, g))) {
        // invert parents
        Vertex child;
        bool first = true;
        for (auto w : walk) {
          if (!first) {
            parents[w] = child;
          }
          child = w;
          first = false;
        }
        parents[v] = other;
        break;
      }
      i = it->second;
    }
    other = v;
  }
}


/// outtake from decompose.cc (schieber ear decomposition)

// do a schieber ear decomposition without any statistics
void schieber_ear_decomposition (Graph& g);

// implementation of the schieber ear decomposition
void ear_decomposition (Graph& g, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex start);

// get boost random spanning tree for schieber ear decomposition
void get_random_spanning_tree (Graph& g, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex start);

// given all the parents of a spanning_tree, find the lca of a crossedge (helper for schieber ear decomposition)
Vertex get_lca (Graph& g, std::map<Vertex, Vertex>& parents, Edge& e, Vertex& r);

// recursion function for the get_lca_distance function
Vertex lca_recursion (Graph& g, std::map<Vertex, Vertex>& parents, Edge& e, Vertex u, Vertex v, Vertex& r);

// Get the distance from one vertex to another (mostly root) in the spanning tree
int get_distance (Graph& g, std::map<Vertex, Vertex>& parents, Vertex v, Vertex r);

// return the distance of a edge to the root (smallest distance of vertices)
int get_distance (Graph& g, std::map<Vertex, Vertex>& parents, Edge e, Vertex r);

// return a iterator poiting to the parent of a certain vertex in the spannin gtree
Vertex get_parent (std::map<Vertex, Vertex>& parents, Vertex v);

void schieber_ear_decomposition (Graph& g) {
  // remember predescessor map and all non-tree edges
  std::map<Vertex, Vertex> parents;
  std::vector<Edge> crossedges;

  // start Vertex
  Vertex startVertex = boost::vertex((boost::num_vertices(g) - 1), g);

  // get a boost random spanning tree
  get_random_spanning_tree(g, parents, crossedges, startVertex);

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
  int ears = 0;

  BGL_FORALL_EDGES_T(e, g, Graph) {
    if (ears < g[e].ear) {
      ears = g[e].ear;
    }
  }

  for (int i = 0; i != ears + 1; i++) {
    Graph& subg = g.create_subgraph();
    // boost::put(&graph_properties::level, g, "decomposed_ears");
    // iterate over edges of graph
    // Graph rg = g.root();

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

  // detect Articulation Points and push them into the graph as vertex property Ak
  color_articulation_points(g);

  min// calculate the two performance critical variables alpha and beta
  // store attachment vertices per each step of the ear decomposition in Ak
  std::map<int, std::vector<Vertex> > Ak;
  unsigned int alpha;
  unsigned int beta;
  std::tie(alpha, beta) = calculate_alpha_beta(g, crossedges, Ak);
  std::cerr << "Alpha is: " << alpha << std::endl;
  std::cerr << "Beta is: " << beta << std::endl;

}

void ear_decomposition (Graph& g, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex start) {

  //delca saves (map of lca_distance : (map of edge : lca)))
  typedef std::map<int, std::map<Edge, Vertex> > delca_t;
  delca_t delca;

  // find the lca distances
  for (auto e : crossedges) {
    if (debug) {
      std::cerr << "starting at new chrossedge: " << e << std::endl;
    }
    Vertex lca = get_lca(g, parents, e, start);
    int lca_distance = get_distance(g, parents, lca, start);
    if (debug) {
      std::cerr << "lca " << lca << " has distance " << lca_distance << std::endl;
    }
    // need to get the distance from the edge to the root
    // int edge_distance = get_distance(g, parents, e, start);

    delca[lca_distance][e] = lca;
  }

  if (debug) {
    for (auto lca_dist = delca.begin(); lca_dist != delca.end(); ++lca_dist) {
      std::cerr << lca_dist->first << "\t";
      for (auto v_it = (lca_dist->second).begin(); v_it != (lca_dist->second).end(); ++v_it) {
        std::cerr << v_it->first << "\t" << v_it->second << "\n";
      }
    }
  }

  // reset ear_integer

  BGL_FORALL_EDGES_T(e, g, Graph) {
    g[e].ear = 0;
  }

  // now start at the biggest lca_distance, at the biggest edge distane and  at the biggest lexmin crossedge 
  // and make walks from the vertices to the lca old values will be overwritten, therefore you get all the ears correctly
  int ear = 0;
  for (auto lca_dist = delca.begin(); lca_dist != delca.end(); ++lca_dist) {
    for (auto v_it = (lca_dist->second).begin(); v_it != (lca_dist->second).end(); ++v_it) {
      Vertex r = v_it->second;
      g[v_it->first].ear = ear;
      //std::cerr << "lca is: " << r << "; coloring ear: " << ear << std::endl;
      Vertex i = boost::source(v_it->first, g);
      while (i != r) {
        Vertex p = get_parent(parents, i);
        g[boost::edge(i, p, g).first].ear = ear;
        //std::cerr << "edge " << boost::edge(i, p, g).first << " will be ear " << ear << std::endl;
        i = p;
      }
      i = boost::target(v_it->first, g);
      while (i != r) {
        Vertex p = get_parent(parents, i);
        g[boost::edge(i, p, g).first].ear = ear;
        //std::cerr << "edge " << boost::edge(i, p, g).first << " will be ear " << ear << std::endl;
        i = p;
      }
      ear++;
    }
  }
}

void get_random_spanning_tree (Graph& g, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex start) {

  boost::associative_property_map< std::map<Vertex, Vertex> > pm(parents);
  // call boost random spanning tree here:
  boost::random_spanning_tree(g, rand_gen, predecessor_map(pm).root_vertex(start));
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
// TODO this is O(n^2)! 

Vertex get_lca (Graph& g, std::map<Vertex, Vertex>& parents, Edge& e, Vertex& r) {

  // make walks from the vertices to the root of the tree
  Vertex u = boost::source(e, g);
  Vertex v = boost::target(e, g);

  // get the lca from the recursion
  v = lca_recursion(g, parents, e, u, v, r);

  if (debug) {
    std::cerr << "lca is: " << v << std::endl;
  }
  return v;
}

Vertex lca_recursion (Graph& g, std::map<Vertex, Vertex>& parents, Edge& e, Vertex u, Vertex v, Vertex& r) {
  if ((v == r) && (u != r)) {
    v = boost::target(e, g);
    u = get_parent(parents, u);
  }
  // if (debug) { std::cerr << "u is: " << u << "; v is: " << v << std::endl; }
  if (u != v) {
    v = get_parent(parents, v);
    v = lca_recursion(g, parents, e, u, v, r);
  }
  return v;
}

int get_distance (Graph& g, std::map<Vertex, Vertex>& parents, Vertex v, Vertex r) {
  int distance = 0;
  while (v != r) {
    distance++;
    v = get_parent(parents, v);
  }
  return distance;
}

int get_distance (Graph& g, std::map<Vertex, Vertex>& parents, Edge e, Vertex r) {
  // return the distance of a edge to the root (smallest distance of vertices)
  int distance_u = get_distance(g, parents, boost::source(e, g), r);
  int distance_v = get_distance(g, parents, boost::target(e, g), r);
  int return_value;
  distance_u < distance_v ? return_value = distance_u : return_value = distance_v;
  return return_value;
}

Vertex get_parent (std::map<Vertex, Vertex>& parents, Vertex v) {
  // get an vertex iterator and search for the vertex v
  std::map<Vertex, Vertex>::iterator it;
  it = parents.find(v);
  // return null_vertex if we cannot find vertex v in parents.
  Vertex return_vertex;
  it == parents.end() ? return_vertex = boost::graph_traits<Graph>::null_vertex() : return_vertex = it->second;
  return return_vertex;
}

void parts_between_articulation_points_to_subgraphs (Graph& g, int k) {

  // reset edge colors

  BGL_FORALL_EDGES_T(e, g, Graph) {
    g[e].color = 0;
  }

  // find vertex to start our walk
  Vertex start;
  Vertex firstAi; // in case of a cycle take this as start!
  bool is_cycle = true;

  BGL_FORALL_VERTICES_T(v, g, Graph) {
    // reset color
    g[v].color = 0;
    // if degree is one, it is an end
    if (degree_in_ear(v, g, k) == 1) {
      start = v;
      is_cycle = false;
      break;
    } else if (g[v].Ai == k) {
      firstAi = v;
    }
  }
  // in case of the last cycle, we find a Ai to start and taint it 2
  if (is_cycle) {
    start = firstAi;
    g[start].color = 2;
  }

  if (debug) {
    std::cerr << "start is: " << boost::get(boost::vertex_color_t(), g, start) << std::endl;
  }
  // create a new subgraph and a
  // pointer which always points to the newest subgraph added
  Graph *subgptr = &g.create_subgraph();
  // add start vertex to subgraph
  boost::add_vertex(boost::get(boost::vertex_color_t(), g, start), *subgptr);
  // bool to see if we reached our end
  bool end_reached = false;

  while (!end_reached) {

    BGL_FORALL_OUTEDGES_T(start, edge, g, Graph) {
      // if this edge does not belong to our ear, continue to next one
      if (g[edge].ear != k) {
        continue;
      }

      // if this vertex is unvisited, do all the magic
      if (g[edge].color == 0) {
        g[edge].color = 1;
        end_reached = false;
        start = boost::target(edge, g);
        if (debug) {
          std::cerr << "start is: " << boost::get(boost::vertex_color_t(), g, start) << std::endl;
        }
        // add to subgraph
        boost::add_vertex(boost::get(boost::vertex_color_t(), g, start), *subgptr);
        break;
      } else {
        end_reached = true;
      }
    }
    // if the current vertex is a internal articulation point, create new subgraph and add this point again
    // we have to exclude path ends, and cycle ends
    if ((g[start].Ai > 0) && (degree_in_ear(start, g, k) == 2) && (g[start].color != 2)) {
      subgptr = &g.create_subgraph();
      if (debug) {
        std::cerr << "ai on v " << boost::get(boost::vertex_color_t(), g, start) << " is: " << g[start].Ai << std::endl;
      }
      boost::add_vertex(boost::get(boost::vertex_color_t(), g, start), *subgptr);
    }
  }
}

void ramachandran_ear_decomposition (Graph& g) {
  // blocks need to be decomposed into path. this can be done by Ear Decomposition

  // start Vertex
  Vertex startVertex = boost::vertex((boost::num_vertices(g) - 1), g);

  // map of ear structure.
  ear_t ear;
  int ear_nr = -1;

  // do the actual ramachandran ear decomposition
  open_ear_decomposition(g, startVertex, ear);

  // create subgraphs from decomposed ears
  std::vector<edge_t> found;
  for (ear_t::iterator it = ear.begin(); it != ear.end(); it++) {
    // sort the vertex pair by ascending index name
    if (it->second.first > it->second.second) {
      it->second = std::make_pair(it->second.second, it->second.first);
    }

    if (!(std::find(found.begin(), found.end(), it->second) != found.end())) {
      // for each ear create a new subgraph
      Graph& subg = g.create_subgraph();
      ear_nr++;
      //boost::put(&graph_properties::level, g, "decomposed_ears");
      if (debug) {
        std::cerr << "New subgraph for ear (" << g.local_to_global(it->second.first) << ","
            << g.local_to_global(it->second.second) << ")" << std::endl
            << "Vertices will be included in subgraph: ";
      }
      for (ear_t::iterator ti = it; ti != ear.end(); ti++) {
        if (ti->second == it->second) {
          found.push_back(ti->second);
          // add vertex into current subgraph if not present already
          if (!subg.find_vertex(g.local_to_global(ti->first.first)).second) {
            boost::add_vertex(g.local_to_global(ti->first.first), subg);
            if (debug) {
              std::cerr << " " << g.local_to_global(ti->first.first);
            }
          }
          if (!subg.find_vertex(g.local_to_global(ti->first.second)).second) {
            boost::add_vertex(g.local_to_global(ti->first.second), subg);
            if (debug) {
              std::cerr << " " << g.local_to_global(ti->first.second);
            }
          }
          // lable edge properly
          g[boost::edge(ti->first.first, ti->first.second, g).first].ear = ear_nr;
        }
      }
      if (debug) {
        std::cerr << std::endl;
      }
    }
  }
  // detect Articulation Points and push them into the graph as vertex property Ak
  color_articulation_points(g);
}

void color_articulation_points (Graph& g) {
  if (debug) {
    std::cerr << "Write articulation point labels into graph object..." << std::endl;
  }
  // start at the outermost ear and process inwards
  int k = 0;
  Graph::children_iterator ear, ear_end;
  for (boost::tie(ear, ear_end) = g.children(); ear != ear_end; ++ear) {

    BGL_FORALL_VERTICES_T(v, *ear, Graph) {
      // if degree is 1, then it is a end of ear_path and therefore a articulation point of this ear.
      if (degree_in_ear(v, *ear, k) < 2) {
        (*ear)[v].Ak.insert(k);
      } else {
        // if this vertex was an articulation point before and is no end point, add it to inner articulation points.
        if ((*ear)[v].Ak.size() > 0) {
          (*ear)[v].Ai = k;
        }
      }
    }
    // goto next ear
    k++;
  }
}

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

void open_ear_decomposition (Graph& g, Vertex startVertex, ear_t& ear);

void ear_dfs (Graph& g, Vertex v, ear_propertymap_t& p, ear_t& ear, unsigned int& counter);

// lexmin implementation for pairs of any type

template <typename T>
std::pair<T, T>& lexmin (std::pair<T, T>& a, std::pair<T, T>& b) {
  if (a.first == b.first) {
    return !(b.second < a.second) ? a : b;
  } else {
    return !(b.first < a.first) ? a : b;
  }
}

void open_ear_decomposition (Graph& g, Vertex startVertex, ear_t& ear) {

  if (debug) {
    std::cerr << "Starting ramachandran ear decomposition" << std::endl;
  }
  // map of ear decomposition properties for all vertices as key
  ear_propertymap_t p;

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
}

void ear_dfs (Graph& g, Vertex v, ear_propertymap_t& p, ear_t& ear, unsigned int& counter) {

  enum {
    WHITE, BLACK, GRAY
  };
  if (debug) {
    std::cout << "v is: " << g.local_to_global(v) << std::endl;
  }

  // start ear decomposition
  p[v].color = GRAY;
  p[v].preorder = counter;
  counter++;
  p[v].low = boost::num_vertices(g);
  p[v].ear = std::make_pair(boost::num_vertices(g), boost::num_vertices(g));

  // get neighbouring vertices

  BGL_FORALL_ADJ_T(v, w, g, Graph) {
    if (debug) {
      std::cout << "w is: " << g.local_to_global(w) << std::endl;
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
          std::cout << "found a crossedge: (" << g.local_to_global(v) << ","
              << g.local_to_global(w) << ")" << std::endl;
        }
        //TODO: casting vertex in low to integer a bad idea?
        p[v].low = boost::vertex(std::min((int) p[v].low, p[w].preorder), g);
        ear[std::make_pair(w, v)] = std::make_pair(boost::vertex(p[w].preorder, g), boost::vertex(p[v].preorder, g));
        p[v].ear = lexmin(p[v].ear, ear[std::make_pair(w, v)]);
      }
    }
  }
  if (debug) {
    std::cout << "finishing vertex " << g.local_to_global(v) << std::endl;
  }
}



// include boost components
#include <boost/graph/breadth_first_search.hpp>

// do a Breadth First Search to test for bipartite property
bool is_bipartite_graph (Graph& g, Vertex startVertex, Edge& ed);

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