// outttakes from struct2graph.cc

// get the spanning tree of our graph
//get_spanning_tree(g, parents, crossedges, startVertex);
//change_spanning_tree(g, r, parents, crossedges, startVertex);


// get spanning tree with DFS
void get_spanning_tree(Graph& g, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex start);

// change spanning tree with sampling a random new one
void change_spanning_tree(Graph& g, std::mt19937& r, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex& root);

void get_spanning_tree(Graph& g, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex start) {

	class my_dfs_visitor : public boost::default_dfs_visitor {
		public:
		my_dfs_visitor(std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges) : p(parents), c(crossedges) {}
		std::map<Vertex, Vertex>& p;
		std::vector<Edge>& c;
		enum { WHITE, BLACK, GRAY, RED };
		void start_vertex(Vertex s, Graph g) const {
			if (verbose) { std::cerr << "Start vertex: " << s << std::endl; }
		}
		void tree_edge(Edge e, Graph g) const {
			if (verbose) { std::cerr << "Detecting tree-edge: " << e << std::endl; }
			Vertex u = boost::source(e, g);
			Vertex v = boost::target(e, g);
			p[v] = u;
		}
		void forward_or_cross_edge(Edge e, Graph g) const {
			if (verbose) { std::cerr << "Detecting back-edge: " << e << std::endl; }
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

void change_spanning_tree(Graph& g, std::mt19937& r, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex& start) {

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
	std::uniform_int_distribution<int> rand_crossedge(0, crossedges.size()-1);  //(min, max)
	Edge old_crossedge = crossedges[rand_crossedge(r)];
	if (verbose) { std::cerr << "choosen crossedge: " << old_crossedge << std::endl; }
	// calculate lca of old crossedge
	Vertex lca = get_lca_distance(g, parents, old_crossedge, start).first;
	if (verbose) { std::cerr << "lca from old crossedge: " << lca << std::endl; }
	// get current cycle
	std::vector<Vertex> source_walk = make_tree_walk(parents, boost::source(old_crossedge, g), lca);
	std::vector<Vertex> target_walk = make_tree_walk(parents, boost::target(old_crossedge, g), lca);
	// merge cycle paths and delete lca
	target_walk.erase(--target_walk.end());
	source_walk.erase(--source_walk.end());
	for (auto elem : source_walk) {
		target_walk.push_back(elem);
	}
	if (verbose) {	std::cerr << "cycle is:" << std::endl << target_walk << std::endl; }
	// get random edge from cycle
	std::uniform_int_distribution<int> rand_tree_edge(0, target_walk.size()-1);  //(min, max)
	int x = rand_tree_edge(r);
	Edge edge = boost::edge(target_walk[x], parents[target_walk[x]], g).first;
	if (verbose) {	std::cerr << "Random Tree Edge is:" << edge << std::endl; }
	// erase source target pair of edge in parents
	for (std::map<Vertex, Vertex>::iterator it=parents.begin(); it!=parents.end(); ++it) {
		if (((it->second == boost::source(edge,g)) && (it->first == boost::target(edge,g))) ||
		((it->second == boost::target(edge,g)) && (it->first == boost::source(edge,g)))) {
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
		if (verbose) { std::cerr << "walk from vertex " << v << std::endl; }
		if (verbose) { std::cerr << "other is " << other << std::endl; }
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
