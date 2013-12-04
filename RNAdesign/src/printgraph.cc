/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 25.03.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

// include header
#include "printgraph.h"

void print_graph (Graph& g, std::ostream* out, std::string nametag) {
	// to convert enums of bases to chars
	std::map<Vertex, char> bases;
	boost::associative_property_map< std::map<Vertex, char> > base_map(bases);
	// to convert list of articulation points to string
	std::map< Vertex, std::string > aps;
	boost::associative_property_map< std::map<Vertex, std::string> > ap_map(aps);
	
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		bases.insert(std::make_pair(v, enum_to_char(g[v].base)));
		
		// convert list of articulation points to string
		std::stringstream ap_stream;
		for (auto elem : g[v].Ak) {
			ap_stream << elem << " ";
		}
		std::string ap_string = ap_stream.str();
		if (ap_string.size() > 1) { ap_string.pop_back(); }					// delete last character
		aps.insert(std::make_pair(v, ap_string));
	}
    
	// print vertex and edge properties from my self-defined bundled properties
	boost::dynamic_properties dp;
	//dp.property("bipartite_color", boost::get(&vertex_property::bipartite_color, g));
	//dp.property("color", boost::get(&vertex_property::color, g));
	dp.property("base", base_map);
	//property("base_enum", boost::get(&vertex_property::base, g));
	dp.property("Ak", ap_map);
	dp.property("Ai", boost::get(&vertex_property::Ai, g));
	dp.property("name", boost::get(boost::vertex_color_t(), g));
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

void print_subgraphs(Graph& g, std::ostream* out, std::string nametag) {
	Graph::children_iterator ci, ci_end;
	int num = 1;
	for (boost::tie(ci, ci_end) = g.children(); ci != ci_end; ++ci) {
		std::stringstream name;
		name << nametag << "-" << num++;
		*out << name.str() << ":";
		print_graph(*ci, out, name.str());
	}
	if (debug) { std::cerr << "Printed all subgraphs." << std::endl; }
}

