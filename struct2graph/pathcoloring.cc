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

unsigned int generate_path_seq (std::string& sequence, char first, char last, unsigned int length) {
	
	// local variables
	std::stringstream seq;			// sequence as stream, for easily appending characters
	
	// check if length is even number
	if (length % 2 != 0) {
		std::cerr << std::endl << "Length of the path to color is an odd number. This can't be!" << std::endl;
		exit(1);
	}
	
	// initialice fibonacci numbers for max length+1
	Fibonacci fibo(length+1);
		
	// calculate total number of possible sequences
	unsigned int number_of_sequences = 2 * (fibo.get(length+1) + fibo.get(length)  );
	
	// delare random number distribution and get a random number
	std::uniform_int_distribution<int> dist(0, number_of_sequences);
	unsigned int rand = dist(rand_gen);
	
	// feed first base to stream
	seq << first;
	
	//do this as the first one must end with an U or G.
	if (seq.str()[-1] == 'A') {
		seq << 'U';
	} else if (seq.str()[-1] == 'C') {
		seq << 'G';
		rand -= fibo.get( length );
	} else {
		rand -= 2 * fibo.get(length);
		if (rand >= fibo.get(length+1)) {
			rand -= fibo.get(length-1);
		}
	}
	
	while (seq.str().size() < length) {
		std::cerr << "length of our sequence: " << seq.str().size() << std::endl;
		std::cerr << "rand is: " << rand << std::endl;
		
		if (rand < fibo.get( length - seq.str().size() )) {
			std::cerr << "rand is: " << rand << " -> lower than fibro: " << fibo.get( length - seq.str().size() ) << std::endl;
			if (seq.str()[-1] == 'G') { seq << "CG"; }
			if (seq.str()[-1] == 'U') { seq << "AU"; }
		} else {
			std::cerr << "is not!" << std::endl;
			rand -= fibo.get( length - seq.str().size() );
			(seq.str()[-1] == 'G') ? seq << 'U' : seq << 'G';
		}
		std::cerr << "Our new sequence: " << seq.str() << std::endl;
	}
	
	sequence = seq.str();
	return number_of_sequences;
}

unsigned int generate_cycle_seq (std::string& sequence, char first, unsigned int length) {
	
	// if it is a cycle, begin and end character are the same
	unsigned int number_of_sequences = generate_path_seq (sequence, first, first, length);
	
	// cut away last base, as cycle ends before
	// sequence substring 0, -1
	
	return number_of_sequences;
}


