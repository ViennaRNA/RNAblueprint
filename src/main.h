/* This program reads RNA secondary structures in dot-bracket notation as well as
 * sequence constraints in IUPAC code and fairly samples RNA sequences compatible
 * to both inputs.
 *
 * @date 25.03.2013
 * @author Stefan Hammer <s.hammer@univie.ac.at>
 * @copyright GPLv3
 *
 */

#ifndef MAIN_H
#define MAIN_H

// include library DependencyGraph
#include "RNAblueprint.h"

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
