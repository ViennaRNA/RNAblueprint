/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 25.03.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
* Compile with: g++ -std=c++11 -g -lboost_program_options -o struct2graph struct2graph.cc
*/

// include header
#include "struct2graph.h"

bool verbose = false;
// filenames of graphml files will be starting with this string
std::string outfile = "";

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
std::pair<T, T>& lexmin(std::pair<T, T>& a, std::pair<T, T>& b) {
	if (a.first == b.first) {
		return !(b.second<a.second)?a:b;
	} else {
		return !(b.first<a.first)?a:b;
	}
}

//! main program starts here
int main(int ac, char* av[]) {

	// initialize command line options
	boost::program_options::variables_map vm = init_options(ac, av);
	
	// input handling ( we read from std:in per default and switch to a file if it is given in the --in option
	// std::in will provide a pseudo interface to enter structures directly!
	std::vector<std::string> structures;
	
	if (vm.count("in")) {
		if (verbose) { std::cerr << "will read graphml file given in the options." << std::endl; }
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
		structures = read_input(&std::cin);		// read infile into array
	}
	
	std::ostream* out = &std::cout;			// out stream
	
	// variables
	Graph graph = parse_graph(structures);		// generate graph from input vector
	*out << "dependency graph:";
	print_graph(graph, out, "root-graph");		// print the graph as GML to a ostream
	connected_components_to_subgraphs(graph);	// get connected components and make subgraphs

	*out << "subgraphs connected components:" << std::endl;
	// print the just created subgraphs
	print_subgraphs(graph, out, "connected-component");
	
	// iterate over all subgraphs (connected components)
	Graph::children_iterator ci, ci_end;
	for (boost::tie(ci, ci_end) = graph.children(); ci != ci_end; ++ci) {
		// check if subgraph is bipartite with a simple BFS
		// generate the vertex 0 as vertex_descriptor
		Graph::vertex_descriptor s = boost::vertex(0, *ci);
		// generate a edge_descriptor in case the graph is not bipartite
		Graph::edge_descriptor ed;
		if (!is_bipartite_graph(*ci, s, ed)) {
			std::cerr << "Graph is not bipartite! Conflict detected on edge " << ed << std::endl;
			exit(1);
		}
		
		// calculate the max degree of this graph
		int max_degree = get_max_degree(*ci);
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
				int max_degree = get_max_degree(*ci_b);
				if (max_degree >= 3) {
					// starting at 0 does not work atm. maybe underflow of unsigned int/vertex?
					ear_decomposition(*ci_b, boost::vertex((boost::num_vertices(*ci_b)-1), *ci_b));
					
					*out << "subgraphs ear decomposition:" << std::endl;
					// print the just created subgraphs
					print_subgraphs(*ci_b, out, "decomposed-ear");
				}
			}
		}
	}
	return 0;
}

boost::program_options::variables_map init_options(int ac, char* av[]) {
	// boost option parser
	namespace po = boost::program_options;
	po::options_description desc("Options");
	desc.add_options()
	    ("help", "print help message [boolean]")
	    ("verbose", "be verbose [boolean]")
	    ("in", po::value<std::string>(), "file to open which contains the structures [string]")
	    ("out", po::value<std::string>(), "write all (sub)graphs to gml files starting with given name [string]")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(ac, av, desc), vm);
	po::notify(vm);  

	if (vm.count("help")) {
		std::cout << desc << "\n";
		exit(1);
	}
	if (vm.count("verbose")) {
		verbose = true;
	}
	if (vm.count("out")) {
		if (verbose) { std::cerr << "graphml files will be written to file." << std::endl; }
		outfile = vm["out"].as<std::string>();
	} else {
		if (verbose) { std::cerr << "grapml files go to std-out" << std::endl; }
	}
	
	return vm;
}

std::vector<std::string> read_input(std::istream* in) {
	// read input file
	std::string line;
	std::vector<std::string> structures;
	while (!in->eof()) {
		getline(*in,line);
		if (line == "@") {
			std::fclose(stdin);
		} else if (line.length() != 0) {
			structures.push_back(line);
		}
	}

	if (verbose) { std::cerr << "Read following structures:" << std::endl; }
	// check if structures have equeal length
	unsigned int length = 0;
	for (auto elem : structures) {
		if (verbose) { std::cerr << elem << std::endl; }
		if ((length != elem.length()) && (length != 0)){
			std::cerr << "Structures have unequal length." << std::endl;
			exit(1);
		}
		length = elem.length();
	}
	return structures;
}

Graph parse_graph(std::vector<std::string> structures) {
	// exit if there is no input
	if (structures.empty()) { exit(1); }
	
	// count the number of positons
	int num_vertices = structures[0].length();
	if (verbose) { std::cerr << "Generating Graph with " << num_vertices << " vertices." << std::endl; }
	Graph g(num_vertices);
	
	// give the vertices names
	int vertex_name = 0;
	Graph::vertex_iterator v, v_end;
	for (boost::tie(v,v_end) = boost::vertices(g); v != v_end; ++v) {
		boost::put(boost::vertex_color_t(), g, *v, vertex_name++);
	}

	// iterate over structures from input
	for (auto elem : structures) {
		std::vector<int> pair_table(structures[0].length(), 0);		// remember position of the open bracket
		unsigned int pos = 0;							// position of the character in the structure
		unsigned int open = 0;							// remember how many open brackets there are
		// iterate over characters from structure
		while (pos < elem.length()) {
			if (elem[pos] == '(') {
				pair_table[open] = pos;
				if (verbose) { std::cerr << elem[pos] << ", open count: "<< open; }
				open++;
			} else if (elem[pos] == ')') {
				open--;
				// check if edge already exists
				bool exists_ab = boost::edge(boost::vertex(pair_table[open],g), boost::vertex(pos,g), g).second;
				bool exists_ba = boost::edge(boost::vertex(pos,g), boost::vertex(pair_table[open],g), g).second;
				if (!exists_ab && !exists_ba) {
					// add edge
					boost::add_edge(boost::vertex(pair_table[open],g), boost::vertex(pos,g), g);
				}
				// reset value
				pair_table[open] = pos;
				if (verbose) { std::cerr << elem[pos] << ", open count: "<< open; }
			} else if (elem[pos] != '.') {
				std::cerr << std::endl << "Unknown character in dot bracked representation" << std::endl;
				exit(1);
			}
			// error handling: there can't be more closing brackets than opening ones
			if (open < 0) {
				std::cerr << std::endl << "Unbalanced brackets in make_pair_table" << std::endl;
				exit(1);
			}
			if (verbose) { std::cerr  << " pos count:" << pos << std::endl; }
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

void print_graph(Graph& g, std::ostream* out, std::string nametag) {

	// print vertex and edge properties from my self defined bundled properties
	boost::dynamic_properties dp;
	dp.property("name", boost::get(boost::vertex_color_t(), g));
	dp.property("bipartite_color", boost::get(&vertex_property::bipartite_color, g));
	
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

void print_subgraphs(Graph& g, std::ostream* out, std::string nametag) {
	Graph::children_iterator ci, ci_end;
	int num = 1;
	for (boost::tie(ci, ci_end) = g.children(); ci != ci_end; ++ci) {
		std::stringstream name;
		name << nametag << "-" << num++;
		*out << name.str() << ":";
		print_graph(*ci, out, name.str());
	}
	if (verbose) { std::cerr << "Printed all sugraphs." << std::endl; }
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
	boost::property_map < Graph, boost::edge_component_t >::type component = boost::get(boost::edge_component, g);
	unsigned int num = boost::biconnected_components(g, component);
	if (verbose) { std::cerr << "Number of biconnected components: " << num << std::endl; }
	
	std::vector<Graph::vertex_descriptor> art_points;
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

int get_max_degree(Graph& g) {
	int max_degree = 0;
	typename Graph::vertex_iterator vi, vi_end;
	for (boost::tie(vi, vi_end) = boost::vertices(g);  vi != vi_end; ++vi) {
		int current_degree = boost::out_degree(*vi, g);
		if (current_degree > max_degree) {
			max_degree = current_degree;
		}
	}
	return max_degree;
}

bool is_bipartite_graph(Graph& g, Graph::vertex_descriptor startVertex, Graph::edge_descriptor& ed) {
	// This is a Breadth First Search which checks if the graph is bipartit. 
	// If not, returns false and the fills the conflicting edge into the edge_descriptor
	
	// queue for search stores vertex indexes
	std::vector<Graph::vertex_descriptor> queue;
	// struct to remember coloring
	std::map<Graph::vertex_descriptor, int> color;
	// struct to remember bfs-coloring
	std::map<Graph::vertex_descriptor, int> bfscolor;
	
	enum { WHITE, BLACK, GRAY, RED };
	
	// add start Vertex to queue
	queue.push_back(startVertex);
	color[startVertex] = BLACK;
	if (verbose) { 	std::cerr << "StartVertex is: " << startVertex << std::endl; 
			std::cerr << "Number of vertices: " << boost::num_vertices(g) << std::endl; }

	// do search
	while (!queue.empty()) {
		Graph::vertex_descriptor u = queue.back();
		if (verbose) { std::cerr << "u is: " << u << std::endl; }
		// get neighbouring vertices
		typename Graph::out_edge_iterator ei, ei_end;
		for (boost::tie(ei, ei_end) = boost::out_edges(u, g);  ei != ei_end; ++ei)
		{
			if (verbose) { std::cerr << boost::target(*ei, g) <<" is neighbour through edge: " << *ei << std::endl; }
			Graph::vertex_descriptor v = boost::target(*ei, g);
			if (verbose) { std::cerr << "v is: " << v << std::endl; }
		
			if (bfscolor[v] == WHITE) {
				bfscolor[v] = GRAY;
				queue.push_back(v);
				
				if (color[u] == RED) {
					color[v] = BLACK;
				} else {
					color[v] = RED;
				}
			} else if (color[u] == color[v]) {
				if (verbose) { std::cerr << "u and v have the same color -> not bipartite!" << std::endl; }
				ed = *ei;
				// return true if graph is not bipartite
				return false;
			} else if (color[u] != color[v]) {
				if (verbose) { std::cerr << "u, v have color: " << color[u] << ", " << color[v] << std::endl; }
			}
		}
		bfscolor[u] = BLACK;
		// remove element u from queue
		queue.erase(std::remove(queue.begin(), queue.end(), u), queue.end());
		
		if (verbose) {  std::cerr << "queue is:" << std::endl << queue;
				std::cerr << "color is:" << std::endl << color;
				std::cerr << "bfscolor is:" << std::endl << bfscolor; }
	}
	// return false for bipartite graphs
	return true;
}

void ear_decomposition(Graph& g, Graph::vertex_descriptor startVertex) {
	// blocks need to be decomposed into path. this can be conde by Ear Decomposition
	
	// map of ear decomposition properties for all vertices as key
	ear_propertymap_t p;
	// map of ear structure.
	ear_t ear;
	// time starts at 0
	unsigned int counter = 0;

	if (verbose) { std::cout << "StartVertex is: " << startVertex << std::endl; }
	// Algorithm rom Ramachandran (1992) Parallel Open Ear Decomposition with Applications, page 8/9
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

void ear_dfs(Graph& g, Graph::vertex_descriptor v, ear_propertymap_t& p, ear_t& ear, unsigned int& counter) {
	
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
		Graph::vertex_descriptor w = boost::target(*ei, g);
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
				if (verbose) { std::cout << "found a backedge: " << v << w << std::endl; }
				//TODO: casting vertex in low to integer a bad idea?
				p[v].low = boost::vertex(std::min((int) p[v].low, p[w].preorder), g);
				ear[std::make_pair(w, v)] = std::make_pair(boost::vertex(p[w].preorder, g), boost::vertex(p[v].preorder, g));
				p[v].ear = lexmin(p[v].ear, ear[std::make_pair(w, v)]);
			}
		}
	}
	if (verbose) { std::cout << "finishing vertex " << v << std::endl; }
}

