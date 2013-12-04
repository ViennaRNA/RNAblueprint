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
#include "graphcommon.h"

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
	
	/*if(debug) {
		std::cerr << "Pairing constructor called and filled" << std::endl;
		for (unsigned int k = 0; k <= l; k++) {
			std::cerr << k << ":" << std::endl;
			for (unsigned int i = 0; i < A_Size; i++) {
				for (unsigned int j = 0; j < A_Size; j++) {
					std::cerr << get(k, i, j) << ", ";
				}
				std::cerr << std::endl;
			}
		}
	}*/
}

rnaMatrix Pairing::multiply(rnaMatrix A, rnaMatrix B) {
	rnaMatrix C;
	int i, j, k;
	long long sum;
	for (i = 0; i < A_Size; i++) {
		for (j = 0; j < A_Size; j++) {
			sum = 0;
			for (k = 0; k < A_Size; k++) {
				sum += A[i][k] * B[k][j];
			}
			C[i][j] = sum;
		}
	}
	return C;
}

unsigned long long Pairing::get(unsigned int l, unsigned int b1, unsigned int b2) {
	
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
	if ((l > length) || (b1 > A_Size-1) || (b2 > A_Size-1)) {
		std::cerr << "Requested a value in pairing matrix which is out of range: p[" << l << "][" << b1 << "][" << b2 << "]" << std::endl;
		exit(1);
	}
	return p[l][b1][b2];
}

unsigned long long Pairing::get(unsigned int l, unsigned int b1) {
	// return the sum of all probabilities of the possible characters at position 2
	unsigned long long sum = 0;
	for (unsigned int i = 0; i < A_Size; i++) {
		sum += get(l, b1, i);
	}
	return sum;
}

unsigned long long Pairing::get(unsigned int l) {
	// return the sum of all probabilities of the possible characters at position 1 and 2
	unsigned long long sum = 0;
	for (unsigned int i = 0; i < A_Size; i++) {
		sum += get(l, i);
	}
	return sum;
}

unsigned long long generate_path_seq (Sequence& sequence, int first, int last, int length) {
		
	// pairing matrix for every length
	Pairing p(length+1);		//TODO initialize only once for the whole program as ist is static content!
	// set maximum possible number of sequences for first....last
	unsigned long long max_number_of_sequences = p.get(length, first, last);
	// delare random number distribution and get a random number
	std::uniform_real_distribution<float> dist(0, 1);
	// number of possible sequences at each possible step
	unsigned long long number_of_sequences = 0;
	// container to remember possible letters after our current letter
	std::vector< int > posibilities;
	
	// set the first base
	if (first != X) {
		sequence.push_back(first);
		length--;
	}
	
	if (debug) {
		std::cerr << "Max Number of Sequences is: " << max_number_of_sequences << std::endl;
		std::cerr << "Length is: " << length << std::endl;
		std::cerr << "Sequence is: " << sequence << std::endl;
		std::cerr << "First is: " << enum_to_char(first) << std::endl;
		std::cerr << "Last is: " << enum_to_char(last) << std::endl;
	}
	
	while (length >= 0) {
		if (first == X) {
			number_of_sequences = p.get(length, first, last);
		} else {
			number_of_sequences = p.get(length+1, first, last);
		}
		// look in paring matrix for next possible character and remember them
		posibilities.clear();
		for (int i = 0; i < A_Size; i++) {
			if (p.get(1, first, i) >= 1) {
				posibilities.push_back(i);
			}
		}
		
		// get a random number between 0 and 1.
		float random = dist(rand_gen);
		
		// stochastically take one of the posibilities
		// start at the probability of first possible character and add each other base probability as long long as the random number is bigger.
		unsigned long long sum = 0;
		for (auto base : posibilities) {
			sum += p.get(length, base, last);
			// if the random number is bigger than our probability, take this base as the first base!
			if (random*number_of_sequences < sum) {
				sequence.push_back(base);
				// our new begin is the chosen base.
				first = base;
				length--;
				// dont forget to exit the loop, otherwise will always be first = C;
				break;
			}
		}
		if (sum == 0) {
			std::cerr << std::endl << "The requested sequence cannot be colored! Conflict at: " 
				<< enum_to_char(first) << ", " << enum_to_char(last) << ", " << length << std::endl;
			exit(1);
		}
		
		if (debug) {
			std::cerr << "Number of Sequences is: " << number_of_sequences << std::endl;
			std::cerr << "Random is: " << random*number_of_sequences << std::endl;
			std::cerr << "Possibilities is: " << posibilities << std::endl;
			std::cerr << "Sequence is: " << sequence << std::endl;
			std::cerr << "--->" << std::endl;
			std::cerr << "First is: " << enum_to_char(first) << std::endl;
			std::cerr << "Last is: " << enum_to_char(last) << std::endl;
			std::cerr << "Length is: " << length << std::endl;
		}
	}
	return max_number_of_sequences;
}

unsigned long long generate_cycle_seq (Sequence& sequence, int first, int length) {

	// max number of sequences to return
	long long max_number_of_sequences = 0;
	// check if length is even number
	if (length % 2 != 0) {
		std::cerr << std::endl << "Length of the cycle to color is an odd number. This can't be!" << std::endl;
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
		std::uniform_real_distribution<float> dist(0, 1);
		float random = dist(rand_gen);
		
		if (debug) { std::cerr << random << std::endl; }
		
		// if random number is smaller than fibo(n-1)/2Lucas(n) -> add an A and color the rest with U,U,n-2
		if (random*max_number_of_sequences < fibo.get(length-1)) {
			sequence.push_back(A);
			generate_path_seq(sequence, U, U, length - 2);
		// if random number is smaller than fibo(n-1)/Lucas(n) -> add an C and color the rest with G,G,n-2
		} else if (random*max_number_of_sequences/2 < fibo.get(length-1)) {
			sequence.push_back(C);
			generate_path_seq(sequence, G, G, length - 2);
		// else -> change random number, choose either G or U and color the rest with G/U,X,n-1
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

unsigned long long color_path_cycle_graph (Graph& g) {
	unsigned long long max_number_of_sequences = 0;
	
	// check if given graph is indeed a path with max_degree = 2 and two ends with degree = 1;
	int max_degree;
	int min_degree;
	std::tie(min_degree, max_degree) = get_min_max_degree(g);
	
	if (max_degree > 2) {
		std::cerr << std::endl << "This graph is no cycle or path (max degree > 2). I can't color this!" << std::endl;
		exit(1);
	}
	
	// find out the length of the path
	unsigned int length = boost::num_edges(g);
	
	// find out the degree of the first and the last base and the base of all non 'X' colored
	std::vector< Vertex > ends;
	std::vector< Vertex > colored_bases;
	Sequence sequence;
	
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		// remember ends of the path
		if (boost::out_degree(v, g) == 1) {
			ends.push_back(v);
		} else if (g[v].base != X) {
			// remember non X bases which are no ends
			colored_bases.push_back(v);
		}
		
		// reset vertex.color tag to 0 -> must be done for sequencestring_to_graph
		g[v].color = 0;
	}
	
	if (ends.size() == 2) {
		// it is a path!
		// check if any of the non-end vertices have assigend colors
		if (colored_bases.size() > 0) {
			std::cerr << std::endl << "This path already is partly colored in between. I can't color this!" << std::endl;
			exit(1);
		}
		// call generate_path_seq and color the vertices accordingly
		max_number_of_sequences = generate_path_seq (sequence, g[ends[0]].base, g[ends[1]].base, length);
		// assign this sequence of bases to the graph
		sequencestring_to_graph(g, ends[0], sequence);
			
	} else if (max_degree == 2 && min_degree == 2 && ends.size() == 0) {
		// it is a cycle (all vertices degree 2)
		if (colored_bases.size() == 1) {
			// start to color at exact this vertex
			max_number_of_sequences = generate_cycle_seq (sequence, g[colored_bases[0]].base, length);
			// assign this sequence of bases to the graph
			sequencestring_to_graph(g, colored_bases[0], sequence);
		} else if (colored_bases.size() == 0) {
			// start to color at any vertex with X
			max_number_of_sequences = generate_cycle_seq (sequence, g[boost::vertex(0, g)].base, length);
			// assign this sequence of bases to the graph
			sequencestring_to_graph(g, boost::vertex(0, g), sequence);
		} else {
			std::cerr << std::endl << "This cycle already is partly colored in between. I can't color this!" << std::endl;
			exit(1);
		}
	} else if (length == 0 && boost::num_vertices(g) == 1) {
		// its a single vertex!
		max_number_of_sequences = generate_path_seq (sequence, g[boost::vertex(0, g)].base, X, length);
		sequencestring_to_graph(g, boost::vertex(0, g), sequence);
	} else {
		// this is no path  - more than two "ends"
		std::cerr << std::endl << "This graph is no cycle or path. I can't color this!" << std::endl;
		exit(1);
	}
	
	return max_number_of_sequences;
}

void sequencestring_to_graph(Graph& g, Vertex vertex, Sequence& sequence) {
	// check if sequence is empty
	if (sequence.size() == 0) {
		std::cerr << std::endl << "Tried to color a path/cycle with a base sequence of zero length. This is impossible!" << std::endl;
		exit(1);
	}

	// check if we are going to overwrite an already existing assignment with a different one.
	if (g[vertex].base != X) {
		if (g[vertex].base != sequence.front()) {
			std::cerr << "Tried to color following vertex with a base, but it is already colored with another one! " 
				<< boost::get(boost::vertex_color_t(), g, vertex) << "/" << enum_to_char(g[vertex].base) << ", new color: " << enum_to_char(sequence.front()) << std::endl;
			exit(1);
		}
	}

	// assign the current vertex with the first character
	g[vertex].base = sequence.front();
	// delete the first character
	sequence.pop_front();
	// mark as done (color = 1)
	g[vertex].color = 1;
	
	// iterate over all adjacent vertices and find the uncolored one.
	// recursively call this function on this vertex again
	BGL_FORALL_ADJ_T(vertex, adj, g, Graph) {
		if (g[adj].color == 0) {
			sequencestring_to_graph(g, adj, sequence);
		}
	}
}


