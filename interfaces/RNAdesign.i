/* RNAdesign.i */
%module RNAdesign
%{

/* Includes the header in the wrapper code */
#include "../lib/RNAdesign.h"

%}

%include "std_string.i"
%include "std_vector.i"

namespace std {
    %template() vector<std::string>;
}

namespace design {
    void initialize_library(bool debug);
    
    template<typename R>
    class DependencyGraph {
        public:
            DependencyGraph(std::vector<std::string> structures, std::string constraints, unsigned long seed);
            DependencyGraph(std::vector<std::string> structures, std::string constraints);
            DependencyGraph(std::vector<std::string> structures);
            ~DependencyGraph();
            std::string get_sequence();
            void mutate();
            void mutate(int position);
            unsigned long long number_of_sequences();
    };
    
    %template(DependencyGraphMT) DependencyGraph<std::mt19937>;
}
