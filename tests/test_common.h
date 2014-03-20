/* C++ Unit test common header file
 *
 * Created on: 23.10.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

#ifndef TESTCOMMON_H
#define TESTCOMMON_H

// include common header with graph definition and global variables
#include "common.h"

// include standard library parts
#include <sstream>
#include <unordered_set>

// include boost components
#include <boost/test/unit_test.hpp>

// typedefs

// functions
// TODO put getvertexset into a shared file (at the moment there are linker issues)
// std::unordered_set<int> getVertexSet(Graph g);

#endif
