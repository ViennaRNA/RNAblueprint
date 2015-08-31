/* RNAdesign
 * A program for designing RNA molecules.
 *
 * Created on: 27.05.2015
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

// include header
#include "pairing_matrix.h"

namespace design {
    namespace detail {

        Fibonacci::Fibonacci(unsigned int length)
        : numbers(length) {
            // Definition: F1 = 0, F2 = 1, Fn = Fn-1 + Fn-2
            numbers[0] = 0;
            numbers[1] = 1;
            for (unsigned int n = 2; n < length; n++) {
                numbers[n] = numbers[n - 1] + numbers[n - 2];
            }
        }

        PairingMatrix* PairingMatrix::Instance() {
            if (!_instance) {
                if (debug) { std::cerr << "Initialize new pairing matrix" << std::endl; }
                _instance = new PairingMatrix();
            }
            return _instance;
        }

        PairingMatrix * PairingMatrix::_instance = nullptr;

        PairingMatrix::~PairingMatrix() {
            delete _instance;
            _instance = nullptr;
        }

        PairingMatrix::PairingMatrix()
        : p(2), length(1) {
            // Definition:
            // length 1: p[A][U][1] = 1, p[U][A][1] = 1, p[G][C][1] = 1, p[C][G][1] = 1, p[U][G][1] = 1, p[G][U][1] = 1
            // if length greater than 2, we don't care about the first letter any more -> x
            // length n: 	p[N][A][n] = p[N][U][n-1]
            //		p[N][C][n] = p[N][G][n-1]
            //		p[N][G][n] = p[N][U][n-1] + p[N][C][n-1]
            //		p[N][U][n] = p[N][A][n-1] + p[N][G][n-1]

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
            
            /*if (debug) {
                std::cerr << "Pairing constructor called and filled" << std::endl;
                for (unsigned int k = 0; k <= l; k++) {
                    std::cerr << std::endl << k << ":" << std::endl;
                    std::cerr << "  " << enum_to_char(0) << " " << enum_to_char(1) << " " << enum_to_char(2) << " " << enum_to_char(3) << " " << std::endl;
                    for (unsigned int i = 0; i < A_Size; i++) {
                        std::cerr << enum_to_char(i) << " ";
                        for (unsigned int j = 0; j < A_Size; j++) {
                            std::cerr << get(k, i, j) << " ";
                        }
                        std::cerr << std::endl;
                    }
                }
            }*/
        }
        
        void PairingMatrix::extend(int newlength) {
            if (debug) { std::cerr << "Extending pairing matrix from length " << length << " to " << newlength << std::endl; }
            // fill pathlength up to length (can be done with matrix multiplication of the pairing matrix
            for (unsigned int i = length+1; i <= newlength; i++) {
                p.push_back( multiply(p[i - 1], p[1]) );
            }
            length = newlength;
        }

        rnaMatrix PairingMatrix::multiply(rnaMatrix A, rnaMatrix B) {
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

        double PairingMatrix::get(unsigned int l, unsigned int b1, unsigned int b2) {

            // if we request a probability for an unknown (N) character at one or both ends, 
            // return the sum of the probabilities for all characters at this position

            //std::cerr << "b1 is " << enum_to_char(b1) << b1 << ", b2 is " << enum_to_char(b2) << b2 << std::endl;

            if ((b1 >= A_Size) || (b2 >= A_Size)) {
                if ((b1 >= A_Size) && (b2 >= A_Size)) {
                    //std::cerr << "b1>=Abet; b2>=Abet" << std::endl;
                    double sum = 0;
                    for (auto i : base_conversion[b2]) {
                        sum += get(l, b1, i);
                    }
                    return sum;

                } else if (b1 >= A_Size) {
                    //std::cerr << "b1>=Abet" << std::endl;
                    double sum = 0;
                    for (auto i : base_conversion[b1]) {
                        sum += get(l, i, b2);
                    }
                    return sum;

                } else if (b2 >= A_Size) {
                    //std::cerr << "b2>=Abet" << std::endl;
                    return get(l, b2, b1);
                }
            } else {
                //std::cerr << "b1<Abet; b2<Abet" << std::endl;
                if (l > length) {
                    // we need to extend the pairing matrix!
                    extend(l);
                } else if ((b1 >= A_Size) || (b2 >= A_Size)) {
                    // check if the requested length is bigger than our initialization or that a base bigger than 3 is requested
                    // -> to avoid segfaults or unknown behavior!
                    std::stringstream ss;
                    ss << "Requested a value in pairing matrix which is out of range: p[" << l << "][" << b1 << "][" << b2 << "]";
                    throw( new std::out_of_range( ss.str() ));
                }
                return p[l][b1][b2];
            }
        }
    }
}
