/* RNAdesign.i */

/*
 * include exception handling from C++ to scripting language
 */
%include "exception.i"
%exception{

  try {
    $action
      }

  catch(const std::exception & e) {
    SWIG_exception_fail(SWIG_RuntimeError, e.what());
  }
 }

/*
 * include the RNAdesign header it self
 */

%module RNAdesign
%{

  /* Includes the header in the wrapper code */
#include "../lib/RNAdesign.h"

  %}

/*
 * stuff necessary to get the RNAdesign DependencyGraph stuff working
 */

%include "std_string.i"
%include "std_vector.i"

namespace std {
   %template(StringVector) std::vector<std::string>;
   %template(IntVector) std::vector<int>;
 }


namespace design {
  void initialize_library(bool debug);
  
  template<typename R>
    class DependencyGraph {
  public:
    DependencyGraph(std::vector<std::string> structures, std::string constraints, unsigned long seed);
    DependencyGraph(std::vector<std::string> structures, std::string constraints);
    DependencyGraph(std::vector<std::string> structures);
    DependencyGraph(const DependencyGraph& copy);
    ~DependencyGraph();
    std::string get_graphml();
    std::string get_graphml(int connected_component_ID);
    std::string get_sequence();
    void set_sequence(std::string sequence);
    void set_sequence();
    bool revert_sequence();
    bool revert_sequence(unsigned int jump);
    double mutate_local(int min_num_pos, int max_num_pos);
    double mutate_local();
    double mutate_global(int min_num_pos, int max_num_pos);
    double mutate_global(int connected_component_ID);
    double mutate_global();
    double mutate(int position);
    double mutate(int start, int end);
    double number_of_sequences();
    double number_of_sequences(int connected_component_ID);
    int number_of_connected_components();
    std::vector< int > component_vertices(int connected_component_ID);
    std::vector< int > special_vertices();
  };
    
  %template(DependencyGraphMT) DependencyGraph<std::mt19937>;
}


/*
 * START WITH THE STUFF NECESSARY TO HANDLE vector<int> AND map<int, vector<int>>
 */

/*
 * extend the standard vector to get the vector itself
 *
 * TODO It should return the vector it self, actually I make a copy of
 * it and it is used realy strange in Perl @vector =
 * @{$ref->{$key}->getvec($ref->{$key})}
 
%extend std::vector<int>{
  std::vector<int> getvec(const std::vector<int> & v){
    std::vector<int> w(v);
    return w;
  }
 }
*/
/*
 * remap the standard perl functions of the map<int, vector<int>> to be used as hash
 
%rename(FETCH)  std::map<int, std::vector<int> >::get;
%rename(STORE)  std::map<int, std::vector<int> >::set;
%rename(EXISTS) std::map<int, std::vector<int> >::has_key;
%rename(DELETE) std::map<int, std::vector<int> >::del;
%rename(SCALAR) std::map<int, std::vector<int> >::size;
%rename(CLEAR)  std::map<int, std::vector<int> >::clear;
*/

/*
 * implement iteration support for the map<int, vector<int>> construct
 
%{
#include <map>
#include <string>
#include <sstream>
  // For iteration support, will leak if iteration stops before the end ever.
  static std::map<void*, std::map<int, std::vector<int> >::const_iterator> iterstate;
  const char *current(std::map<int, std::vector<int> >& map) {
    std::map<void*, std::map<int, std::vector<int> >::const_iterator>::iterator it = iterstate.find(&map);
    
    if (it != iterstate.end() && map.end() == it->second) {
      // clean up entry in the global map
      iterstate.erase(it);
      it = iterstate.end();
    }
    if (it == iterstate.end()){
      return NULL;
    }
    else{
      //used to convert the int we get from the it->second->first to a string which is actually need for iteration
      std::string string = static_cast<std::ostringstream*>( &(std::ostringstream() << it->second->first) )->str();
      return string.c_str();
    }
  }
 
  %}
*/

/*
 * implement the convenient Perl keys function

%extend std::map<int, std::vector<int> > {
  std::map<int, std::vector<int> > *TIEHASH() {
    return $self;
  }
  const char *FIRSTKEY() {
    iterstate[$self] = $self->begin();
    return current(*$self);
  }
  const char *NEXTKEY(const std::string&) {
    ++iterstate[$self];
    return current(*$self);
  }
 }
 */