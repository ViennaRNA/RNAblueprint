#ifndef MOD_LIB_GRAPHMORPHISM_FILTEREDWRAPPER_H
#define	MOD_LIB_GRAPHMORPHISM_FILTEREDWRAPPER_H

#include <mod/Error.h>
#include <mod/lib/PairToRangeAdaptor.h>

#include <boost/graph/filtered_graph.hpp>

#include <utility>

namespace mod {
namespace lib {
namespace GraphMorphism {

template<typename Graph>
struct FilteredWrapper {
	using base_traits = boost::graph_traits<Graph>;
	// Graph
	using vertex_descriptor = typename base_traits::vertex_descriptor;
	using edge_descriptor = typename base_traits::edge_descriptor;
	using directed_category = typename base_traits::directed_category;
	using edge_parallel_category = typename base_traits::edge_parallel_category;
	using traversal_category = typename base_traits::traversal_category;

	// IncidenceGraph
	typedef typename boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;
	typedef typename boost::graph_traits<Graph>::degree_size_type degree_size_type;

	// BidirectionalGraph
	typedef typename boost::graph_traits<Graph>::in_edge_iterator in_edge_iterator;

	// AdjacencyGraph
	// adjacency_iterator

	// VertexListGraph
	// vertex_iterator
	// vertices_size_type

	// EdgeListGraph
	typedef typename boost::graph_traits<Graph>::edge_iterator edge_iterator;
	typedef typename boost::graph_traits<Graph>::edges_size_type edges_size_type;
	typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iterator;

	//graph_traits<filtered_graph>::edge_iterator 
	//
	//The type for the iterators returned by edges(), which is:
	//filter_iterator<EdgePredicate, graph_traits<Graph>::edge_iterator>
	//The iterator is a model of MultiPassInputIterator.
	//graph_traits<filtered_graph>::adjacency_iterator 
	//
	//The type for the iterators returned by adjacent_vertices(). The adjacency_iterator models the same iterator concept as out_edge_iterator.

	typedef typename boost::graph_traits<Graph>::vertices_size_type vertices_size_type;

	//
	//The type used for dealing with the number of vertices in the graph.
	//graph_traits<filtered_graph>::edge_size_type 
	//
	//The type used for dealing with the number of edges in the graph.
	//graph_traits<filtered_graph>::degree_size_type 
public:

	FilteredWrapper(const Graph & g)
	: g(g), map(num_vertices(g), std::numeric_limits<vertices_size_type>::max()) {
		vertices_size_type count = 0;
		//		std::cout << "FilterWrapper" << std::endl;

		for(typename boost::graph_traits<Graph>::vertex_descriptor v : asRange(vertices(g))) {
			vertices_size_type vId = get(boost::vertex_index_t(), g, v);
			map[vId] = count;
			//			std::cout << "Map: " << vId << " => " << count << std::endl;
			//			std::cout << "ReverseMap: " << reverseMap.size() << " => " << vId << std::endl;
			reverseMap.push_back(vId);
			count++;
		}
	}
public:
	const Graph &g;
	std::vector<vertices_size_type> map;
	std::vector<vertices_size_type> reverseMap;
};

template<typename Graph>
FilteredWrapper<Graph> makeFilteredWrapper(const Graph &g) {
	return FilteredWrapper<Graph > (g);
}

template<typename Graph>
struct FilteredWrapperIndexMap {
	typedef typename boost::graph_traits<Graph>::vertices_size_type VSizeType;

	FilteredWrapperIndexMap() : g(nullptr) { }

	explicit FilteredWrapperIndexMap(const FilteredWrapper<Graph> &g) : g(&g) { }

	VSizeType operator[](typename boost::graph_traits<Graph>::vertex_descriptor v) const {
		assert(g);
		VSizeType vId = get(boost::vertex_index_t(), g->g, v);
		if(g->map[vId] == std::numeric_limits<VSizeType>::max()) {
			//			std::cout << "FilteredWrapperIndexMap; request for " << v << "(" << vId << ") is not defined" << std::endl;
			//			std::cout << "Map is:" << std::endl << "from:";
			//			for(unsigned int i = 0; i < g->map.size(); i++) std::cout << "\t" << i;
			//			std::cout << std::endl << "to:  ";
			//			for(const auto &p : g->map) {
			//				std::cout << "\t";
			//				if(p == g->g.null_vertex()) std::cout << "-(-)";
			//				else std::cout << p << "(" << get(boost::vertex_index_t(), g->g, p) << ")";
			//			}
			//			std::cout << std::endl;
		}
		assert(g->map[vId] != std::numeric_limits<VSizeType>::max());
		return g->map[vId];
	}
private:
	const FilteredWrapper<Graph> *g;
};

} // namespace GraphMorphism
} // namespace lib
} // namespace mod
namespace boost {

template<typename InnerGraph>
struct graph_traits<mod::lib::GraphMorphism::FilteredWrapper<InnerGraph> > {
	typedef mod::lib::GraphMorphism::FilteredWrapper<InnerGraph> Graph;
	// Graph
	typedef typename Graph::vertex_descriptor vertex_descriptor;
	typedef typename Graph::edge_descriptor edge_descriptor;
	typedef typename Graph::directed_category directed_category;
	typedef typename Graph::edge_parallel_category edge_parallel_category;
	typedef typename Graph::traversal_category traversal_category;

	static vertex_descriptor null_vertex() {
		return graph_traits<InnerGraph>::null_vertex();
	}

	// IncidenceGraph
	typedef typename Graph::out_edge_iterator out_edge_iterator;
	typedef typename Graph::degree_size_type degree_size_type;

	// BidirectionalGraph
	typedef typename Graph::in_edge_iterator in_edge_iterator;

	// AdjacencyGraph
	// adjacency_iterator

	// VertexListGraph
	// vertex_iterator
	// vertices_size_type

	// EdgeListGraph
	typedef typename Graph::edge_iterator edge_iterator;
	typedef typename Graph::edges_size_type edges_size_type;
	typedef typename Graph::vertex_iterator vertex_iterator;

	//graph_traits<filtered_graph>::edge_iterator 
	//
	//The type for the iterators returned by edges(), which is:
	//filter_iterator<EdgePredicate, graph_traits<Graph>::edge_iterator>
	//The iterator is a model of MultiPassInputIterator.
	//graph_traits<filtered_graph>::adjacency_iterator 
	//
	//The type for the iterators returned by adjacent_vertices(). The adjacency_iterator models the same iterator concept as out_edge_iterator.

	typedef typename Graph::vertices_size_type vertices_size_type;

	//
	//The type used for dealing with the number of vertices in the graph.
	//graph_traits<filtered_graph>::edge_size_type 
	//
	//The type used for dealing with the number of edges in the graph.
	//graph_traits<filtered_graph>::degree_size_type 
};

} // namespace boost
namespace mod {
namespace lib {
namespace GraphMorphism {

// IncidenceGraph

template<typename Graph>
std::pair<typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::out_edge_iterator,
typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::out_edge_iterator>
out_edges(typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::vertex_descriptor v,
		const mod::lib::GraphMorphism::FilteredWrapper<Graph> &g) {
	return out_edges(v, g.g);
}

template<typename Graph>
typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::vertex_descriptor
source(typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::edge_descriptor e, const mod::lib::GraphMorphism::FilteredWrapper<Graph> &g) {
	return source(e, g.g);
}

template<typename Graph>
typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::vertex_descriptor
target(typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::edge_descriptor e, const mod::lib::GraphMorphism::FilteredWrapper<Graph> &g) {
	return target(e, g.g);
}

template<typename Graph>
typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::degree_size_type
out_degree(typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::vertex_descriptor v,
		const mod::lib::GraphMorphism::FilteredWrapper<Graph> &g) {
	return out_degree(v, g.g);
}

// BidirectionalGraph

template<typename Graph>
std::pair<typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::in_edge_iterator,
typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::in_edge_iterator>
in_edges(typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::vertex_descriptor v, const mod::lib::GraphMorphism::FilteredWrapper<Graph> &g) {
	return in_edges(v, g.g);
}

template<typename Graph>
typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::degree_size_type
in_degree(typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::vertex_descriptor v,
		const mod::lib::GraphMorphism::FilteredWrapper<Graph> &g) {
	return in_degree(v, g.g);
}
//degree(e, g)	degree_size_type
// AdjacencyGraph 
//adjacent_vertices(v, g)	std::pair<adjacency_iterator, adjacency_iterator>

// VertexListGraph

template<typename Graph>
std::pair<typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::vertex_iterator,
typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::vertex_iterator>
vertices(const mod::lib::GraphMorphism::FilteredWrapper<Graph> &g) {
	return vertices(g.g);
}

template<typename Graph>
typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::vertices_size_type
num_vertices(const mod::lib::GraphMorphism::FilteredWrapper<Graph> &g) {
	return std::distance(vertices(g.g).first, vertices(g.g).second);
}

// EdgeListGraph

template<typename Graph>
std::pair<typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::edge_iterator,
typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::edge_iterator>
edges(const mod::lib::GraphMorphism::FilteredWrapper<Graph> &g) {
	return edges(g.g);
}

template<typename Graph>
typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::edges_size_type
num_edges(const mod::lib::GraphMorphism::FilteredWrapper<Graph> &g) {
	return std::distance(edges(g.g).first, edges(g.g).second);
}

// AdjacencyMatrix 

template<typename Graph>
std::pair<typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::edge_descriptor, bool>
edge(typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::vertex_descriptor u,
		typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::vertex_descriptor v,
		const mod::lib::GraphMorphism::FilteredWrapper<Graph> &g) {
	return edge(u, v, g.g);
}

} // namespace GraphMorphism
} // namespace lib
} // namespace mod
namespace boost {
// PropertyGraph

template<typename Graph>
struct property_traits<mod::lib::GraphMorphism::FilteredWrapperIndexMap<Graph> > {
	typedef typename graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::vertices_size_type value_type;
	typedef typename graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::vertices_size_type reference;
	typedef typename graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::vertex_descriptor key_type;
	typedef readable_property_map_tag category;
};

template <typename Graph>
struct vertex_property_type<mod::lib::GraphMorphism::FilteredWrapper<Graph> > : vertex_property_type<Graph> {
};

template <typename Graph>
struct edge_property_type<mod::lib::GraphMorphism::FilteredWrapper<Graph> > : edge_property_type<Graph> {
};

template <typename Graph>
struct graph_property_type<mod::lib::GraphMorphism::FilteredWrapper<Graph> > : graph_property_type<Graph> {
};

template<typename Graph, typename Property>
struct property_map<mod::lib::GraphMorphism::FilteredWrapper<Graph>, Property> {
	typedef typename property_map<Graph, Property>::type type;
	typedef typename property_map<Graph, Property>::const_type const_type;
};

template<typename Graph>
struct property_map<mod::lib::GraphMorphism::FilteredWrapper<Graph>, vertex_index_t> {
	typedef mod::lib::GraphMorphism::FilteredWrapperIndexMap<Graph> type;
	typedef type const_type;
};

} // namespace boost
namespace mod {
namespace lib {
namespace GraphMorphism {

template<typename Graph>
typename boost::property_traits<mod::lib::GraphMorphism::FilteredWrapperIndexMap<Graph> >::value_type
get(const mod::lib::GraphMorphism::FilteredWrapperIndexMap<Graph> &map,
		typename boost::property_traits<mod::lib::GraphMorphism::FilteredWrapperIndexMap<Graph> >::key_type k) {
	return map[k];
}

template<typename Graph, typename PropertyTag>
typename boost::property_map<mod::lib::GraphMorphism::FilteredWrapper<Graph>, PropertyTag>::const_type
get(PropertyTag t, const mod::lib::GraphMorphism::FilteredWrapper<Graph> &g) {
	return get(t, g.g);
}

template<typename Graph>
mod::lib::GraphMorphism::FilteredWrapperIndexMap<Graph> get(boost::vertex_index_t, const mod::lib::GraphMorphism::FilteredWrapper<Graph> &g) {
	return mod::lib::GraphMorphism::FilteredWrapperIndexMap<Graph > (g);
}

template<typename Graph, typename PropertyTag, typename VertexOrEdge>
typename boost::property_traits<typename boost::property_map<mod::lib::GraphMorphism::FilteredWrapper<Graph>, PropertyTag>::const_type>::reference
get(PropertyTag t, const mod::lib::GraphMorphism::FilteredWrapper<Graph> &g, VertexOrEdge ve) {
	return get(get(t, g), ve);
}

// Other

template<typename Graph>
inline typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::vertex_descriptor
vertex(typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::vertices_size_type n,
		const mod::lib::GraphMorphism::FilteredWrapper<Graph> &g) {

	for(typename boost::graph_traits<mod::lib::GraphMorphism::FilteredWrapper<Graph> >::vertex_descriptor v : asRange(vertices(g))) {
		if(get(boost::vertex_index_t(), g, v) == n) return v;
	}
	MOD_ABORT;
}

} // namespace GraphMorphism
} // namespace lib
} // namespace mod

#endif	/* MOD_LIB_GRAPHMORPHISM_FILTEREDWRAPPER_H */