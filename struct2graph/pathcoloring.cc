/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 13.08.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

// include header
#include "pathcoloring.h"

Fibonacci::Fibonacci(unsigned int length)
 : numbers(length) 
{
	// Definition: F1 = 0, F2 = 1, Fn = Fn-1 + Fn-2
	numbers[0] = 0;
	numbers[1] = 1;
	for (unsigned int n = 2; n < length; n++) {
		numbers[n] = numbers[n-1] + numbers[n-2];
	}
}

Pairing::Pairing(unsigned int l)
 : p(l+1)
 , length(l)
{
	// Definition:
	// length 1: p[A][U][1] = 1, p[U][A][1] = 1, p[G][C][1] = 1, p[C][G][1] = 1, p[U][G][1] = 1, p[G][U][1] = 1
	// if length greater than 2, we don't care about the first letter any more -> x
	// length n: 	p[X][A][n] = p[X][U][n-1]
	//		p[X][C][n] = p[X][G][n-1]
	//		p[X][G][n] = p[X][U][n-1] + p[X][C][n-1]
	//		p[X][U][n] = p[X][A][n-1] + p[X][G][n-1]
	
	// fill standard pairing matrix (pathlength = 1)
	p[1][A][U] = 1;
	p[1][U][A] = 1;
	p[1][G][C] = 1;
	p[1][C][G] = 1;
	p[1][G][U] = 1;
	p[1][U][G] = 1;
	
	// fill pathlength up to length (can be done with matrix multiplication of the pairing matrix
	for (unsigned int i = 2; i <= l; i++) {
		p[i] = multiply(p[i-1],p[1]);
	}
	
	if(verbose) {
		std::cerr << "Pairing constructor called and filled" << std::endl;
		for (unsigned int l = 1; l <= l; l++) {
			for (unsigned int i = 0; i <=3; i++) {
				for (unsigned int j = 0; j <=3; j++) {
					std::cerr << get(l, i, j) << ", ";
				}
				std::cerr << std::endl;
			}
			std::cerr << std::endl;
		}
	}
}

rnaMatrix Pairing::multiply(rnaMatrix A, rnaMatrix B) {
	rnaMatrix C;
	int i, j, k;
	int sum;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			sum = 0;
			for (k = 0; k < 4; k++) {
				sum += A[i][k] * B[k][j];
			}
			C[i][j] = sum;
		}
	}
	return C;
}

unsigned int Pairing::get(unsigned int l, unsigned int b1, unsigned int b2) {
	
	if ((l > length) || (b1 > 3) || (b2 > 3))) {
		std::cerr << "Requested a value in pairing matrix which is out of range." << std::endl;
		exit(1);
	}
	return p[l][b1][b2];
}

unsigned int Pairing::get(unsigned int l, unsigned int b1) {
	unsigned int sum = 0;
	for (unsigned int i = 0; i < 4; i++) {
		sum += get(l, b1, i);
	}
	return sum;
}

unsigned int Pairing::get(unsigned int l) {
	unsigned int sum = 0;
	for (unsigned int i = 0; i < 4; i++) {
		sum += get(l, i);
	}
	return sum;
}

unsigned int generate_path_seq (std::string& sequence, int first, int last, unsigned int length) {
	
	// check if length is even number
	if (length % 2 != 0) {
		std::cerr << std::endl << "Length of the path to color is an odd number. This can't be!" << std::endl;
		exit(1);
	}
	
	int number_of_sequences = 0;
	// pairing matrix for every length
	Pairing p(length);
	// delare random number distribution and get a random number
	std::uniform_real_distribution<int> dist(0, 1);
	
	// if begin and end are X, just assign a begining character
	if ((first == X) && (last == X)) {
		number_of_sequences = p.get(length);
		first = floor(dist(rand_gen)*4);
	// if begin or end are X, always start at begining
	} else if ((first == X) && (last != X)){
		first = last;
		last = X;
	}
	
	// start coloring here!
	
	
	/*
	// initialice fibonacci numbers for max length+1
	Fibonacci fibo(length+1);
		
	// calculate total number of possible sequences
	unsigned int number_of_sequences = 2 * (fibo.get(length+1) + fibo.get(length) );
	
	// delare random number distribution and get a random number
	std::uniform_int_distribution<int> dist(0, number_of_sequences);
	
	for (int i = 0; i < 100; i++) {
	
	unsigned int rand = dist(rand_gen);
	//std::cerr << "initial rand is: " << rand << std::endl;
	
	// feed first base to stream
	sequence += first;
	
	// do this as the first one must end with an U or G.
	if (sequence.back() == 'A') {
		sequence += 'U';
	} else if (sequence.back() == 'C') {
		sequence += 'G';
		rand -= fibo.get( length );
	} else {
		rand -= (2 * fibo.get(length));
		if (rand >= fibo.get(length+1)) {
			rand -= fibo.get(length+1);
		}
	}
	
	// extend the sequence
	while (sequence.size() < length) {
	//	std::cerr << "rand is: " << rand << std::endl;
		
		if (rand < fibo.get( length - sequence.size() )) {
			if (*sequence.rbegin() == 'G') { sequence += "CG"; }
			if (*sequence.rbegin() == 'U') { sequence += "AU"; }
		} else {
			rand -= fibo.get( length - sequence.size() );
			(sequence.back() == 'G') ? sequence += 'U' : sequence += 'G';
		}
	//	std::cerr << "Our new sequence: " << sequence << std::endl;
	}
	
	// exchange last character and fit to previous one
	if (sequence.back() != last) {
	//	std::cerr << "Last character does not fit: " << sequence.back() << last << std::endl;
		sequence.back() = last;
		
		if (sequence.back() == 'A') {
			sequence[sequence.length()-2] = 'U';
		} else if (sequence.back() == 'C') {
			sequence[sequence.length()-2] = 'G';
		}
	}
	
	
	std::cerr << sequence << std::endl;
	sequence = "";
	}
	*/
	
	
	
	
	return number_of_sequences;
}

unsigned int generate_cycle_seq (std::string& sequence, int first, unsigned int length) {
	return generate_path_seq (sequence, int first, int first, length);
}
