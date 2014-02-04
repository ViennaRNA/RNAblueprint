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

void decompose_graph(Graph& graph, std::ostream* out, int num_trees, bool no_bipartite_check) {

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
					//if (num_trees > 0) {
					//	do_spanning_tree_stat(*bc, num_trees);
					//} else {
						// ear_decomposition (Ramachandran (1992))
						ramachandran_ear_decomposition(*bc);
						
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
					//}
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

void ramachandran_ear_decomposition (Graph& g) {
	// blocks need to be decomposed into path. this can be done by Ear Decomposition
	
	// start Vertex
	Vertex startVertex = boost::vertex((boost::num_vertices(g)-1), g);
	
	// map of ear structure.
	ear_t ear;
	int ear_nr = -1;
	
	// do the actual ramachandran ear decomposition
	open_ear_decomposition (g, startVertex, ear);
	
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
			if (debug) { 	std::cerr << "New subgraph for ear (" << it->second.first << "," 
						<< it->second.second << ")" << std::endl 
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
					// lable edge properly
					g[boost::edge(ti->first.first, ti->first.second, g).first].ear = ear_nr;
				}
			}
			if (debug) { 	std::cerr << std::endl; }
		}
	}
	// detect Articulation Points and push them into the graph as vertex property Ak
	color_articulation_points (g);
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
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		// if degree is one, it is an end we can start with
		if (degree_in_ear(v, g, k) == 1) {
			start = v;
			break;
		}
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
		// we have to exclude path ends
		if ((g[start].Ai > 0) && (degree_in_ear(start, g, k) == 2)) {
			subgptr = &g.create_subgraph();
			if (debug) { std::cerr << "ai on v " << boost::get(boost::vertex_color_t(), g, start) << " is: " << g[start].Ai << std::endl; }
			boost::add_vertex(boost::get(boost::vertex_color_t(), g, start), *subgptr);
		}
	}
}

