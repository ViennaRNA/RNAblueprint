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

void print_graph(Graph& g, std::ostream* out, std::string nametag) {

	// print vertex and edge properties from my self defined bundled properties
	boost::dynamic_properties dp;
	dp.property("bipartite_color", boost::get(&vertex_property::bipartite_color, g));
	dp.property("color", boost::get(&vertex_property::color, g));
	dp.property("base", boost::get(&vertex_property::base, g));
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
	if (verbose) { std::cerr << "Printed all subgraphs." << std::endl; }
}

