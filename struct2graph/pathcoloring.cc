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
 , length(l+1)
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
	
	//fill length 0 with probability 1 for same base (important for setting last base)
	p[0][A][A] = 1;
	p[0][U][U] = 1;
	p[0][G][G] = 1;
	p[0][C][C] = 1;
	
	// fill pathlength up to length (can be done with matrix multiplication of the pairing matrix
	for (unsigned int i = 2; i <= l; i++) {
		p[i] = multiply(p[i-1],p[1]);
	}
	
	if(verbose) {
		std::cerr << "Pairing constructor called and filled" << std::endl;
		for (unsigned int k = 0; k <= l; k++) {
			for (unsigned int i = 0; i < 4; i++) {
				for (unsigned int j = 0; j < 4; j++) {
					std::cerr << get(k, i, j) << ", ";
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
	
	// if we request a probability for an unknown (X) character at one or both ends, 
	// return the sum of the probabilities for all characters at this position
	if ((b1 == X) || (b2 == X)) {
		if ((b1 == X) && (b2 == X)) {
			return get(l);
		} else if (b1 == X) {
			return get(l, b2);
		} else if (b2 == X) {
			return get(l, b1);
		}
	}
	
	// check if the requested length is bigger than our initilisation or that a base bigger than 3 is requested
	// -> to avoid segfaults or unknown behaviour!
	if ((l > length) || (b1 > 3) || (b2 > 3)) {
		std::cerr << "Requested a value in pairing matrix which is out of range: p[" << l << "][" << b1 << "][" << b2 << "]" << std::endl;
		exit(1);
	}
	return p[l][b1][b2];
}

unsigned int Pairing::get(unsigned int l, unsigned int b1) {
	// return the sum of all probabilities of the possible characters at position 2
	unsigned int sum = 0;
	for (unsigned int i = 0; i < 4; i++) {
		sum += get(l, b1, i);
	}
	return sum;
}

unsigned int Pairing::get(unsigned int l) {
	// return the sum of all probabilities of the possible characters at position 1 and 2
	unsigned int sum = 0;
	for (unsigned int i = 0; i < 4; i++) {
		sum += get(l, i);
	}
	return sum;
}

unsigned int generate_path_seq (std::string& sequence, int first, int last, int length) {
		
	// pairing matrix for every length
	Pairing p(length+1);
	// set maximum possible muber of sequences for first....last
	int max_number_of_sequences = p.get(length, first, last);
	// delare random number distribution and get a random number
	std::uniform_real_distribution<float> dist(0, 1);
	// number of possible sequences at each possible step
	int number_of_sequences = 0;
	//
	
	// set the first base
	if (first != X) {
		sequence += enum_to_char(first);
		length--;
	}
	
	while (length >= 0) {
		number_of_sequences = p.get(length+1, first, last);
		
		// look in paring matrix for next possible character and remember them
		std::vector< int > posibilities(4);
		for (int i = 0; i < 4; i++) {
			if (p.get(1, first, i) >= 1) {
				posibilities.push_back(i);
			}
		}
		
		// get a rando number between 0 and the max number of seq [0,nos).
		float random = dist(rand_gen)*number_of_sequences;
		
		// stochastically take one of the posibilities
		// start at the probability of first possible character and add each other base probability as long as the random number is bigger.
		int sum = 0;
		for (auto base : posibilities) {
			sum += p.get(length, base, last);
			// if the random number is bigger than our probability, take this base as the first base!
			if (random < sum) {
				sequence += enum_to_char(base);
				first = base;
				length--;
				// dont forget to exit the loop, otherwise will always be first = C;
				break;
			}
		}
	}
	return max_number_of_sequences;
}

unsigned int generate_cycle_seq (std::string& sequence, int first, int length) {

	// max number of sequences to return
	int max_number_of_sequences = 0;
	// check if length is even number
	if (length % 2 != 0) {
		std::cerr << std::endl << "Length of the path to color is an odd number. This can't be!" << std::endl;
		exit(1);
	}
	
	if (first != X) {
		// return a path with same begin and end, but then remove the last character again -> cycle!
		max_number_of_sequences = generate_path_seq (sequence, first, first, length);
		sequence.pop_back();
	} else {
		// initialize fibonacci numbers
		Fibonacci fibo(length+1);
		// max number of sequences is
		max_number_of_sequences = 2*(fibo.get(length + 1)+fibo.get(length - 1)); // is same as 2* Lucas (length)
		// delare random number distribution and get a random number
		std::uniform_real_distribution<float> dist(0, 1); //TODO include 1!!
		float random = dist(rand_gen);
		if (verbose) {
			std::cerr << random << std::endl;
		}
		
		if (random*max_number_of_sequences < fibo.get(length-1)) {
			sequence += 'A';
			generate_path_seq(sequence, U, U, length - 2);
		} else if (random*max_number_of_sequences/2 < fibo.get(length-1)) {
			sequence += 'C';
			generate_path_seq(sequence, G, G, length - 2);
		} else {
			random -= 2*fibo.get(length-1)/max_number_of_sequences;
			if (random*max_number_of_sequences < fibo.get(length + 1)) {
				generate_path_seq(sequence, U, X, length-1);
			} else {
				generate_path_seq(sequence, G, X, length-1);
			}
		}
	}
	return max_number_of_sequences;
}
