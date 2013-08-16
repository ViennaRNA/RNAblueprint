/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 13.08.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

// include header
#include "graphcoloring.h"

// take highest ear, write all posibilities for open cut points into matrix (together with number of possible sequences for whole ear)

// next ear: find open and fixed cut points
// write into matrix the posibilities of the fixed (internal) cut point(s) from previous matrix and generate the posibilities for new open cut points
// therefore need function that checks the pairing matrix and gives you back the posibility of the last given the first letter of a path.
// eg: length 2, first A -> last has to be U
// length 3, first U -> last letter can be U or A

// in the end add cycle -> all cut points are internal ones.
// now randomly get one combination of cut point bases.

// color the paths in between with pathcoloring.cc
