/* This program reads secundary RNA structures in dot-bracket and
 * builds a graph for a latter ear-decomposition and bipartitness-check
 *
 * Created on: 25.03.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

#ifndef MAIN_H
#define MAIN_H

// include common header with graph definition and global variables
#include "common.h"

// include standard library parts
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>

// include boost components
#include <boost/program_options.hpp>

// initialise boost command line option parser
boost::program_options::variables_map init_options (int ac, char* av[]);

// read the input file into a string
std::vector<std::string> read_input (std::istream* in);

#endif
