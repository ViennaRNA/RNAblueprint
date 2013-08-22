/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 25.03.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

// include header
#include "common.h"

std::mt19937 rand_gen;

char enum_to_char(int intletter) {
	char charletter;
	switch ( intletter ) {
		case A:	charletter = 'A';
			break;
		case G:	charletter = 'G';
			break;
		case C:	charletter = 'C';
			break;
		case U:	charletter = 'U';
			break;
	}
	return charletter;
}
