/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 25.03.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

// include header
#include "decompose.h"
#include "printgraph.h"
#include "treestatistics.h"
#include "graphcommon.h"

void decompose_graph(Graph& graph, std::ostream* out, int num_trees, bool ramachandran, bool no_bipartite_check) {

	connected_components_to_subgraphs(graph);	// get connected components and make subgraphs

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
		int max_degree = get_min_max_degree(*ci).second;
		if (verbose) { std::cerr << "Max degree of subgraph is: " << max_degree << std::endl; }
		
		// split further into biconnected components do ear decomposition
		if (max_degree >= 3) {
			biconnected_components_to_subgraphs(*ci);
			
			*out << "subgraphs biconnected components:" << std::endl;
			// print the just created subgraphs
			print_subgraphs(*ci, out, "biconnected-component");
			
			Graph::children_iterator ci_b, ci_b_end;
			for (boost::tie(ci_b, ci_b_end) = (*ci).children(); ci_b != ci_b_end; ++ci_b) {
				// calculate the max degree of this graph
				int max_degree = get_min_max_degree(*ci_b).second;
				if (max_degree >= 3) {
					// only for statistics start at all vertices as root for DFS
					if (num_trees > 0) {
						do_spanning_tree_stat(*ci_b, num_trees);
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

void connected_components_to_subgraphs(Graph& g) {
	
	// get list of connected components into the component vector
	// http://www.boost.org/doc/libs/1_53_0/libs/graph/doc/connected_components.html
	// http://www.boost.org/doc/libs/1_53_0/libs/graph/example/connected_components.cpp
	std::vector<int> component(boost::num_vertices(g));
	int num = boost::connected_components(g, &component[0]);
	
	if (verbose) { 
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

void biconnected_components_to_subgraphs(Graph& g) {

	// // get list of biconnected components into the component property map
	// http://www.boost.org/doc/libs/1_53_0/libs/graph/doc/biconnected_components.html
	// http://www.boost.org/doc/libs/1_38_0/libs/graph/example/biconnected_components.cpp
	boost::edge_component_t edge_component;
	boost::property_map < Graph, boost::edge_component_t >::type component = boost::get(edge_component, g);
	unsigned int num = boost::biconnected_components(g, component);
	if (verbose) { std::cerr << "Number of biconnected components: " << num << std::endl; }
	
	std::vector<Vertex> art_points;
	boost::articulation_points(g, std::back_inserter(art_points));
	if (verbose) {	std::cerr << "Number of articulation points: " << art_points.size() << " ( "; 
		for (auto elem : art_points) {
			std::cerr << boost::get(boost::vertex_color_t(), g, elem) << " ";
		}
		std::cerr << ")" << std::endl;	
	}
	
	if (verbose) {
		// get graph and iterate over its edges to print connected components table
		typename Graph::edge_iterator ei, ei_end;
		for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
			std::cerr << *ei << "\t" <<  "(" << boost::get(boost::vertex_color_t(), g, boost::source(*ei, g)) << "," 
			<< boost::get(boost::vertex_color_t(), g, boost::target(*ei, g))<< ")" 
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
				if (!subg.find_vertex(boost::get(boost::vertex_color_t(), g, boost::source(*ei,g))).second) {
					boost::add_vertex(boost::get(boost::vertex_color_t(), g, boost::source(*ei,g)), subg);
				}
			}
		}
	}
	
}

void schieber_ear_decomposition (Graph& g) {
	// remember predescessor map and all non-tree edges
	std::map<Vertex, Vertex> parents;
	std::vector<Edge> crossedges;
	
	// start Vertex
	Vertex startVertex = boost::vertex((boost::num_vertices(g)-1), g);
	
	// get a boost random spanning tree
	get_random_spanning_tree (g, parents, crossedges, startVertex);
	
	// print parents, cross-edges and root vertex
	if (verbose) {
		std::cerr << "Root vertex: " << startVertex << std::endl;
		std::cerr << "Spanning tree (vertex, parent) and cross-edges:" << std::endl;
		for (std::map<Vertex, Vertex>::iterator it=parents.begin(); it!=parents.end(); ++it) {
			std::cerr << it->first << " => " << it->second << std::endl;
		}
		for (auto elem : crossedges) {
			std::cerr << elem << std::endl;
		}
	}
	
	// do the actual ear decomposition
	ear_decomposition (g, parents, crossedges, startVertex);
	
	// calculate the two performance critical variables alpha and beta
	// store attachment vertices per each step of the ear decomposition in Ak
	std::map<int, std::vector<Vertex> > Ak;
	unsigned int alpha;
	unsigned int beta;
	std::tie(alpha, beta) = calculate_alpha_beta(g, crossedges, Ak);
	std::cerr << "Alpha is: " << alpha << std::endl;
	std::cerr << "Beta is: " << beta << std::endl;	
	
	// write ears into subgraphs
	int ears = 0;
	BGL_FORALL_EDGES_T(e, g, Graph) {
		if (ears < g[e].ear) { ears = g[e].ear; }
	}
	
	for (int i = 0; i != ears+1; i++) {
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
				if (!subg.find_vertex(boost::get(boost::vertex_color_t(), g, boost::source(e,g))).second) {
					boost::add_vertex(boost::get(boost::vertex_color_t(), g, boost::source(e,g)), subg);
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
		if (verbose) { std::cerr << "starting at new chrossedge: " << e << std::endl; }
		std::pair<Vertex, int> lcad = get_lca_distance(g, parents, e, start);
		if (verbose) { std::cerr << "lca " << lcad.first << " has distance " << lcad.second << std::endl; }
		delca[lcad.second][e] = lcad.first;
	}
	
	if (verbose) {
		for (std::map<int, std::map<Edge, Vertex> >::iterator it=delca.begin(); it!=delca.end(); ++it) {
			std::cerr << it->first << "(";
			for (std::map<Edge, Vertex>::iterator iit=(it->second).begin(); iit!=(it->second).end(); ++iit) {
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
	for (std::map<int, std::map<Edge, Vertex> >::reverse_iterator it=delca.rbegin(); it!=delca.rend(); ++it) {
		for (std::map<Edge, Vertex>::reverse_iterator iit=(it->second).rbegin(); iit!=(it->second).rend(); ++iit) {
			Vertex r = iit->second;
			g[iit->first].ear = ear;
			//std::cerr << "lca is: " << r << "; coloring ear: " << ear << std::endl;
			Vertex i = boost::source(iit->first,g);
			while (i != r) {
				std::map<Vertex, Vertex>::iterator iiit;
				iiit = parents.find(i);
				g[boost::edge(i,iiit->second,g).first].ear = ear;
				//std::cerr << "edge " << boost::edge(i,iiit->second,g).first << " will be ear " << ear << std::endl;
				i = iiit->second;
			}
			i = boost::target(iit->first,g);
			while (i != r) {
				std::map<Vertex, Vertex>::iterator iiit;
				iiit = parents.find(i);
				g[boost::edge(i,iiit->second,g).first].ear = ear;
				//std::cerr << "edge " << boost::edge(i,iiit->second,g).first << " will be ear " << ear << std::endl;
				i = iiit->second;
			}
			ear++;
		}
	}
}

void get_random_spanning_tree (Graph& g, std::map<Vertex, Vertex>& parents, std::vector<Edge>& crossedges, Vertex start) {
	
	boost::associative_property_map< std::map<Vertex,Vertex> > pm(parents);
	// call boost random spanning tree here:
	boost::random_spanning_tree(g, rand_gen, predecessor_map(pm).root_vertex(start));
	if (verbose) { std::cerr << "Got a boost random spanning tree..." << std::endl; }
	
	// create the crossedges vector for the later ear-decomposition!
	// clear all values
	BGL_FORALL_EDGES_T(e, g, Graph) {
		g[e].color = 0;
	}
	// mark all tree edges
	for (std::map<Vertex, Vertex>::iterator it=parents.begin(); it!=parents.end(); ++it) {
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

std::pair<Vertex, int> get_lca_distance(Graph& g, std::map<Vertex, Vertex>& parents, Edge e, Vertex r) {

	// make walks from the vertices to the root of the tree
	std::vector<Vertex> uwalk = make_tree_walk(parents, boost::target(e, g), r);
	std::vector<Vertex> vwalk = make_tree_walk(parents, boost::source(e, g), r);
	if (verbose) {
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
		if ((uwalk.size() == 0) || (vwalk.size() == 0)) { break; }
	}	
	return std::make_pair(lca,distance);
}

std::vector<Vertex> make_tree_walk(std::map<Vertex, Vertex>& parents, Vertex v, Vertex r) {
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

void ramachandran_ear_decomposition (Graph& g) {
	// blocks need to be decomposed into path. this can be done by Ear Decomposition
	
	// start Vertex
	Vertex startVertex = boost::vertex((boost::num_vertices(g)-1), g);
	
	// map of ear decomposition properties for all vertices as key
	ear_propertymap_t p;
	// map of ear structure.
	ear_t ear;
	// time starts at 0
	unsigned int counter = 0;

	if (verbose) { std::cout << "StartVertex is: " << startVertex << std::endl; }
	// Algorithm from Ramachandran (1992) Parallel Open Ear Decomposition with Applications, page 8/9
	ear_dfs(g, startVertex, p, ear, counter);
	
	// print out all data-structures at the end
	if (verbose) { 
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
			if (verbose) { 	std::cerr << "New subgraph for ear (" << it->second.first << "," << it->second.second << ")" << std::endl 
					<< "Vertices will be included in subgraph: "; }
			for (ear_t::iterator ti = it; ti != ear.end(); ti++) {
				if (ti->second == it->second) {
					found.push_back(ti->second);
					// add vertex into current subgraph if not present already
					if (!subg.find_vertex(ti->first.first).second) {
						boost::add_vertex(ti->first.first, subg);
						if (verbose) { std::cerr << " " << ti->first.first; }
					}
					if (!subg.find_vertex(ti->first.second).second) {
						boost::add_vertex(ti->first.second, subg);
						if (verbose) { std::cerr << " " << ti->first.second; }
					}
				}
			}
			if (verbose) { 	std::cerr << std::endl; }
		}
	}
}

void ear_dfs(Graph& g, Vertex v, ear_propertymap_t& p, ear_t& ear, unsigned int& counter) {
	
	enum { WHITE, BLACK, GRAY };
	if (verbose) { std::cout << "v is: " << v << std::endl; }
	
	// start ear decomposition
	p[v].color = GRAY;
	p[v].preorder = counter;
	counter++;
	p[v].low = boost::num_vertices(g);
	p[v].ear = std::make_pair(boost::num_vertices(g), boost::num_vertices(g));
	
	// get neighbouring vertices
	typename Graph::out_edge_iterator ei, ei_end;
	for (boost::tie(ei, ei_end) = boost::out_edges(v, g);  ei != ei_end; ++ei)
	{
		if (verbose) { std::cerr << boost::target(*ei, g) <<" is neighbour through edge: " << *ei << std::endl; }
		Vertex w = boost::target(*ei, g);
		if (verbose) { std::cout << "w is: " << w << std::endl; }
		
		if (p[w].color == WHITE) {
			if (verbose) { std::cout << "w is white" << std::endl; }
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
			if (verbose) { std::cout << "w is gray" << std::endl; }
			if (w != p[w].parent) {
				if (verbose) { std::cout << "found a crossedge: " << v << w << std::endl; }
				//TODO: casting vertex in low to integer a bad idea?
				p[v].low = boost::vertex(std::min((int) p[v].low, p[w].preorder), g);
				ear[std::make_pair(w, v)] = std::make_pair(boost::vertex(p[w].preorder, g), boost::vertex(p[v].preorder, g));
				p[v].ear = lexmin(p[v].ear, ear[std::make_pair(w, v)]);
			}
		}
	}
	if (verbose) { std::cout << "finishing vertex " << v << std::endl; }
}


