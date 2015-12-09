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

// include library DependencyGraph
#include "RNAdesign.h"

// include standard library parts
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <exception>
#include <cstddef>

// include boost components
#include <boost/program_options.hpp>



// initialize boost command line option parser
boost::program_options::variables_map init_options(int ac, char* av[]);

// read the input file into a string
std::tuple<std::vector<std::string>, std::string, std::string > read_input(std::istream* in);

// overload << operator to print vectors with any content

template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T>& vec) {
    int i = 0;
    for (auto elem : vec) {
        os << "(" << i++ << ") " << elem << std::endl;
    }
    return os;
}

#endif
