/* RNAdesign
 * A program for designing RNA molecules.
 *
 * Created on: 25.03.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

// include header
#include "common.h"

namespace design {
    namespace detail {

        bool debug = false;
        bool * debug_ptr = &debug;

        //A, C, G, U, R, Y, K, M, S, W, B, D, H, V, N

        char enum_to_char(int intletter) {
            char charletter;
            switch (intletter) {
                case A: charletter = 'A';
                    break;
                case C: charletter = 'C';
                    break;
                case G: charletter = 'G';
                    break;
                case U: charletter = 'U';
                    break;
                case R: charletter = 'R';
                    break;
                case Y: charletter = 'Y';
                    break;
                case K: charletter = 'K';
                    break;
                case M: charletter = 'M';
                    break;
                case S: charletter = 'S';
                    break;
                case W: charletter = 'W';
                    break;
                case B: charletter = 'B';
                    break;
                case D: charletter = 'D';
                    break;
                case H: charletter = 'H';
                    break;
                case V: charletter = 'V';
                    break;
                case N: charletter = 'N';
                    break;
                default: charletter = '-';
                    break;
            }
            return charletter;
        }

        int char_to_enum(char charletter) {
            int intletter;
            switch (charletter) {
                case 'A': intletter = A;
                    break;
                case 'C': intletter = C;
                    break;
                case 'G': intletter = G;
                    break;
                case 'U': intletter = U;
                    break;
                case 'R': intletter = R;
                    break;
                case 'Y': intletter = Y;
                    break;
                case 'K': intletter = K;
                    break;
                case 'M': intletter = M;
                    break;
                case 'S': intletter = S;
                    break;
                case 'W': intletter = W;
                    break;
                case 'B': intletter = B;
                    break;
                case 'D': intletter = D;
                    break;
                case 'H': intletter = H;
                    break;
                case 'V': intletter = V;
                    break;
                case 'N': intletter = N;
                    break;
                default: 
                    throw std::out_of_range("This is not a valid IUPAC notation!");
                    break;
            }
            return intletter;
        }

        std::ostream& operator<<(std::ostream& os, Sequence& sequence) {
            for (auto base : sequence) {
                os << enum_to_char(base);
            }
            return os;
        }
    }
}