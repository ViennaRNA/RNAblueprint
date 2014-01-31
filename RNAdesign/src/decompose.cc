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

// include boost components
#include <boost/graph/iteration_macros.hpp>

void decompose_graph(Graph& graph, std::ostream* out, int num_trees, bool ramachandran, bool no_bipartite_check) {

	connected_components_to_subgraphs(graph);	// get connected components and make subgraphs

	if (debug) {
		*out << "subgraphs connected components:" << std::endl;
		// print the just created subgraphs
		print_subgraphs(graph, out, "connected-component");
	}
	
	// iterate over all subgraphs (connected components)
	Graph::children_iterator cc, cc_end;
	for (boost::tie(cc, cc_end) = graph.children(); cc != cc_end; ++cc) {
		if (!no_bipartite_check) {
			// check if subgraph is bipartite with a simple BFS
			// generate the vertex 0 as vertex_descriptor
			Vertex s = boost::vertex(0, *cc);
			// generate a edge_descriptor for the case that the graph is not bipartite
			Edge ed;
			if (!is_bipartite_graph(*cc, s, ed)) {
				std::cerr << "Graph is not bipartite! Conflict detected on edge " << ed << std::endl;
				exit(1);
			}
		}
		
		// calculate the max degree of this graph
		int max_degree = get_min_max_degree(*cc).second;
		if (debug) { std::cerr << "Max degree of subgraph is: " << max_degree << std::endl; }
		
		// split further into biconnected components do ear decomposition
		if (max_degree > 2) {
			biconnected_components_to_subgraphs(*cc);
			
			if (debug) {
				*out << "subgraphs biconnected components:" << std::endl;
				// print the just created subgraphs
				print_subgraphs(*cc, out, "biconnected-component");
			}
			
			Graph::children_iterator bc, bc_end;
			for (boost::tie(bc, bc_end) = (*cc).children(); bc != bc_end; ++bc) {
				// calculate the max degree of this graph (biconnected component)
				int max_degree = get_min_max_degree(*bc).second;
				if (max_degree > 2) {
					// only for statistics start at all vertices as root for DFS
					if (num_trees > 0) {
						do_spanning_tree_stat(*bc, num_trees);
					} else {
						//TODO use do_ear_decompositons (Schieber & Vishkin (1986)) or 
						// ear_decomposition (Ramachandran (1992)) one?!
						if (ramachandran) {
							ramachandran_ear_decomposition(*bc);
						} else {
							schieber_ear_decomposition(*bc);
						}
						
						if (debug) {
							*out << "subgraphs ear decomposition:" << std::endl;
							// print the just created subgraphs
							print_subgraphs(*bc, out, "decomposed-ear");
						}
						// now lets push parts between articulation points of an ear to subgraphs
						int k = 0;
						Graph::children_iterator ear, ear_end;
						for (boost::tie(ear, ear_end) = (*bc).children(); ear != ear_end; ++ear) {
							
							parts_between_articulation_points_to_subgraphs (*ear, k++);
							if (debug) {
								*out << "parts between articulation points of ear:" << std::endl;
								// print the just created subgraphs
								print_subgraphs(*ear, out, "articulation-point-parts");
							}
						}
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

void biconnected_components_to_subgraphs(Graph& g) {

	// // get list of biconnected components into the component property map
	// http://www.boost.org/doc/libs/1_53_0/libs/graph/doc/biconnected_components.html
	// http://www.boost.org/doc/libs/1_38_0/libs/graph/example/biconnected_components.cpp
	boost::edge_component_t edge_component;
	boost::property_map < Graph, boost::edge_component_t >::type component = boost::get(edge_component, g);
	unsigned int num = boost::biconnected_components(g, component);
	if (debug) { std::cerr << "Number of biconnected components: " << num << std::endl; }
	
	std::vector<Vertex> art_points;
	boost::articulation_points(g, std::back_inserter(art_points));
	if (debug) {	std::cerr << "Number of articulation points: " << art_points.size() << " ( "; 
		for (auto elem : art_points) {
			std::cerr << boost::get(boost::vertex_color_t(), g, elem) << " ";
		}
		std::cerr << ")" << std::endl;	
	}
	
	if (debug) {
		// iterate over all graph edges to print connected components table
		BGL_FORALL_EDGES_T(e, g, Graph) {
			std::cerr << e << "\t" <<  "(" << boost::get(boost::vertex_color_t(), g, boost::source(e, g)) << "," 
			<< boost::get(boost::vertex_color_t(), g, boost::target(e, g))<< ")" 
			<< "\tcomponent: " << component[e] << std::endl;
		}
	}
	
	// now need to merge biconnected components that are separated by a articulation point that has a degree == 2 !
	for (auto v : art_points) {
		if (boost::degree(v, g) > 2) {
			BGL_FORALL_ADJ_T(v, adj, g, Graph) {
				if ((boost::degree(adj, g) == 2) 
				&& (std::find(art_points.begin(), art_points.end(), adj) != art_points.end())) {
					int nc = -1;
					merge_biconnected_paths(g, v, adj, component, art_points, nc);
				}
			}
		}
	}
	
	// write biconnected components into subgraphs:
	for (unsigned int i = 0; i != num; i++) {
		// only create a subgraph if there is really an edge associated to this component (as we merged many components before)
		bool exists = false;
		BGL_FORALL_EDGES_T(e, g, Graph) {
			if (i == component[e])
				exists = true;
		}
		if (!exists)
			continue;
		
		// for this bicomponent number generate a new subgraph
		Graph& subg = g.create_subgraph();
		//boost::put(&graph_properties::level, g, "biconnected_component");
		// iterate over edges of graph
		//Graph rg = g.root();
		BGL_FORALL_EDGES_T(e, g, Graph) {
			if (i == component[e]) {
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

void merge_biconnected_paths(Graph& g, Vertex p, Vertex v, boost::property_map < Graph, boost::edge_component_t >::type& component, std::vector<Vertex>& art_points, int& nc) {
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
		if (nc == -1)
			nc = component[e];
		else
			component[e] = nc;
	}
	
	BGL_FORALL_ADJ_T(v, adj, g, Graph) {
		if ((adj != p) && (boost::degree(adj, g) == 2) 
		&& (std::find(art_points.begin(), art_points.end(), adj) != art_points.end()))
			merge_biconnected_paths(g, v, adj, component, art_points, nc);
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
	if (debug) {
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
	
	// write ears into subgraphs
	int ears = 0;
	BGL_FORALL_EDGES_T(e, g, Graph) {
		if (ears < g[e].ear) { ears = g[e].ear; }
	}
	
	for (int i = 0; i != ears+1; i++) {
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
				if (!subg.find_vertex(boost::get(boost::vertex_color_t(), g, boost::source(e,g))).second) {
					boost::add_vertex(boost::get(boost::vertex_color_t(), g, boost::source(e,g)), subg);
				}
			}
		}
	}
	
	// detect Articulation Points and push them into the graph as vertex property Ak
	color_articulation_points (g);
	
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
		if (debug) { std::cerr << "starting at new chrossedge: " << e << std::endl; }
		Vertex lca = get_lca(g, parents, e, start);
		int lca_distance = get_distance(g, parents, lca, start);
		if (debug) { std::cerr << "lca " << lca << " has distance " << lca_distance << std::endl; }
		// need to get the distance from the edge to the root
		// int edge_distance = get_distance(g, parents, e, start);
		
		delca[lca_distance][e] = lca;
	}
	
	if (debug) {
		for (auto lca_dist=delca.begin(); lca_dist!=delca.end(); ++lca_dist) {
			std::cerr << lca_dist->first << "\t";
			for (auto v_it=(lca_dist->second).begin(); v_it!=(lca_dist->second).end(); ++v_it) {
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
	for (auto lca_dist=delca.begin(); lca_dist!=delca.end(); ++lca_dist) {
		for (auto v_it=(lca_dist->second).begin(); v_it!=(lca_dist->second).end(); ++v_it) {
			Vertex r = v_it->second;
			g[v_it->first].ear = ear;
			//std::cerr << "lca is: " << r << "; coloring ear: " << ear << std::endl;
			Vertex i = boost::source(v_it->first,g);
			while (i != r) {
				Vertex p = get_parent(parents, i);
				g[boost::edge(i, p, g).first].ear = ear;
				//std::cerr << "edge " << boost::edge(i, p, g).first << " will be ear " << ear << std::endl;
				i = p;
			}
			i = boost::target(v_it->first,g);
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
	
	boost::associative_property_map< std::map<Vertex,Vertex> > pm(parents);
	// call boost random spanning tree here:
	boost::random_spanning_tree(g, rand_gen, predecessor_map(pm).root_vertex(start));
	if (debug) { std::cerr << "Got a boost random spanning tree..." << std::endl; }
	
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
// TODO this is O(n^2)! 
Vertex get_lca(Graph& g, std::map<Vertex, Vertex>& parents, Edge& e, Vertex& r) {

	// make walks from the vertices to the root of the tree
	Vertex u = boost::source(e, g);
	Vertex v = boost::target(e, g);

	// get the lca from the recursion
	v = lca_recursion(g, parents, e, u, v, r);
	
	if (debug) { std::cerr << "lca is: " << v << std::endl; }
	return v;
}

Vertex lca_recursion(Graph& g, std::map<Vertex, Vertex>& parents, Edge& e, Vertex u, Vertex v, Vertex& r) {
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

int get_distance(Graph& g, std::map<Vertex, Vertex>& parents, Vertex v, Vertex r) {
	int distance = 0;
	while (v != r) {
		distance++;
		v = get_parent(parents, v);
	}
	return distance;
}

int get_distance(Graph& g, std::map<Vertex, Vertex>& parents, Edge e, Vertex r) {
	// return the distance of a edge to the root (smallest distance of vertices)
	int distance_u = get_distance(g, parents, boost::source(e, g), r);
	int distance_v = get_distance(g, parents, boost::target(e, g), r);
	int return_value;
	distance_u < distance_v ? return_value = distance_u : return_value = distance_v;
	return return_value;
}

Vertex get_parent(std::map<Vertex, Vertex>& parents, Vertex v) {
	// get an vertex iterator and search for the vertex v
	std::map<Vertex, Vertex>::iterator it;
	it = parents.find(v);
	// return null_vertex if we cannot find vertex v in parents.
	Vertex return_vertex;
	it == parents.end() ? return_vertex = boost::graph_traits<Graph>::null_vertex() : return_vertex = it->second;
	return return_vertex;
}

void color_articulation_points (Graph& g) {
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
				if ((*ear)[v].Ak.size() > 0 ) {
					(*ear)[v].Ai = k;
				}
			}
		}
		// goto next ear
		k++;
	}
}

int degree_in_ear (Vertex& v, Graph& g, int k) {
	int degree = 0;
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
		if (g[e].ear ==	k) {
			degree++;
		}
	}
	return degree;
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
	
	if (debug) { std::cerr << "start is: " << boost::get(boost::vertex_color_t(), g, start) << std::endl; }
	// create a new subgraph and a
	// pointer which always points to the newest subgraph added
	Graph *subgptr = &g.create_subgraph();
	// add start vertex to subgraph
	boost::add_vertex(boost::get(boost::vertex_color_t(), g, start), *subgptr);
	// bool to see if we reached our end
	bool end_reached = false;
	
	while(!end_reached) {
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
				if (debug) { std::cerr << "start is: " << boost::get(boost::vertex_color_t(), g, start) << std::endl; }
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
			if (debug) { std::cerr << "ai on v " << boost::get(boost::vertex_color_t(), g, start) << " is: " << g[start].Ai << std::endl; }
			boost::add_vertex(boost::get(boost::vertex_color_t(), g, start), *subgptr);
		}
	}
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

	if (debug) { std::cout << "StartVertex is: " << startVertex << std::endl; }
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
			if (debug) { 	std::cerr << "New subgraph for ear (" << it->second.first << "," << it->second.second << ")" << std::endl 
					<< "Vertices will be included in subgraph: "; }
			for (ear_t::iterator ti = it; ti != ear.end(); ti++) {
				if (ti->second == it->second) {
					found.push_back(ti->second);
					// add vertex into current subgraph if not present already
					if (!subg.find_vertex(ti->first.first).second) {
						boost::add_vertex(ti->first.first, subg);
						if (debug) { std::cerr << " " << ti->first.first; }
					}
					if (!subg.find_vertex(ti->first.second).second) {
						boost::add_vertex(ti->first.second, subg);
						if (debug) { std::cerr << " " << ti->first.second; }
					}
				}
			}
			if (debug) { 	std::cerr << std::endl; }
		}
	}
	// detect Articulation Points and push them into the graph as vertex property Ak
	color_articulation_points (g);
}

void ear_dfs(Graph& g, Vertex v, ear_propertymap_t& p, ear_t& ear, unsigned int& counter) {
	
	enum { WHITE, BLACK, GRAY };
	if (debug) { std::cout << "v is: " << v << std::endl; }
	
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
		if (debug) { std::cerr << boost::target(*ei, g) <<" is neighbour through edge: " << *ei << std::endl; }
		Vertex w = boost::target(*ei, g);
		if (debug) { std::cout << "w is: " << w << std::endl; }
		
		if (p[w].color == WHITE) {
			if (debug) { std::cout << "w is white" << std::endl; }
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
			if (debug) { std::cout << "w is gray" << std::endl; }
			if (w != p[w].parent) {
				if (debug) { std::cout << "found a crossedge: " << v << w << std::endl; }
				//TODO: casting vertex in low to integer a bad idea?
				p[v].low = boost::vertex(std::min((int) p[v].low, p[w].preorder), g);
				ear[std::make_pair(w, v)] = std::make_pair(boost::vertex(p[w].preorder, g), boost::vertex(p[v].preorder, g));
				p[v].ear = lexmin(p[v].ear, ear[std::make_pair(w, v)]);
			}
		}
	}
	if (debug) { std::cout << "finishing vertex " << v << std::endl; }
}


