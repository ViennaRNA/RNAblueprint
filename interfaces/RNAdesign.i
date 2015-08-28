/* RNAdesign.i */
%module RNAdesign
%{

/* Includes the header in the wrapper code */
#include "../lib/RNAdesign.h"

%}

%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"

namespace std {
    %template() std::vector<std::string>;
    %template() std::map< int, std::vector<int> >;
}

%include "exception.i"
%exception{

    try {
        $action
    }

    catch(const std::exception & e) {
        SWIG_exception_fail(SWIG_RuntimeError, e.what());
    }
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
            void set_sequence(std::string sequence);
            void set_sequence();
            unsigned long long mutate_local(int min_num_pos, int max_num_pos);
            unsigned long long mutate_global(int min_num_pos, int max_num_pos);
            unsigned long long mutate(int position);
            unsigned long long mutate(int start, int end);
            unsigned long long number_of_sequences();
            unsigned long long number_of_sequences(int connected_component_ID);
            std::map< int, std::vector<int> > connected_components();
            std::vector< int > special_vertices();
    };
    
    %template(DependencyGraphMT) DependencyGraph<std::mt19937>;
}
