/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 25.03.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

#ifndef TREESTATISTICS_H
#define TREESTATISTICS_H

// include common header with graph definition and global variables
#include "common.h"
#include "decompose.h"
#include "printgraph.h"

// include standard library parts
#include <iostream>
#include <sstream>

// include boost components


// do statistics with many spannign trees calculating alpha and beta values using the schieber algorithm
void do_spanning_tree_stat (Graph& g, int num_trees);

// calculate spanning tree properties for coloring
std::pair< unsigned int, unsigned int > calculate_alpha_beta(Graph& g, std::vector<Edge>& crossedges, std::map<int, std::vector<Vertex> >& Ak);

// write statistic output file
void print_ab_stat (unsigned int alpha, unsigned int beta, std::map<int, std::vector<Vertex> > Ak, Graph& g, Vertex root, std::vector<Edge>& crossedges);


#endif
