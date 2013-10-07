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
		case X:	charletter = 'X';
			break;
	}
	return charletter;
}

int char_to_enum(char charletter) {
	int intletter;
	switch ( charletter ) {
		case 'A': intletter = A;
			break;
		case 'G': intletter = G;
			break;
		case 'C': intletter = C;
			break;
		case 'U': intletter = U;
			break;
		case 'X': intletter = X;
			break;
	}
	return intletter;
}

std::ostream& operator<< (std::ostream& os, Sequence& sequence) {
	for (auto base : sequence) {
		os << enum_to_char(base);
	}
	return os;
}
