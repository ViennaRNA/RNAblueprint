/* 
 * File:   uninduced_subgraph.hpp
 * Author: Stefan Hammer <jango@tbi.univie.ac.at>
 *
 * Created on October 6, 2015, 2:33 PM
 */

#ifndef UNINDUCED_SUBGRAPH_HPP
#define	UNINDUCED_SUBGRAPH_HPP

#include <boost/graph/subgraph.hpp>

namespace boost {

// Invariants of an induced subgraph:
//   - If vertex u is in subgraph g, then u must be in g.parent().
//   - If edge e is in subgraph g, then e must be in g.parent().

// In contrast to the subgraph class, this is NOT TRUE:
//   - If edge e=(u,v) is in the root graph, then edge e
//     is also in any subgraph that contains both vertex u and v.

// The Graph template parameter must have a vertex_index and edge_index
// internal property. It is assumed that the vertex indices are assigned
// automatically by the graph during a call to add_vertex(). It is not
// assumed that the edge vertices are assigned automatically, they are
// explicitly assigned here.
    
    
template <typename Graph>
class uninduced_subgraph : public subgraph<Graph> {

public:
    
    uninduced_subgraph(const graph_property_type& p) : subgraph<Graph>(p) {  }

    uninduced_subgraph(vertices_size_type n, const graph_property_type& p = graph_property_type()) : subgraph<Graph>(n, p) {  }
    
    // copy constructor
    uninduced_subgraph(const uninduced_subgraph& x)
        : m_parent(x.m_parent), m_edge_counter(x.m_edge_counter)
        , m_global_vertex(x.m_global_vertex), m_global_edge(x.m_global_edge)
    {
        if(x.is_root())
        {
         m_graph = x.m_graph;
        }
        // Do a deep copy (recursive).
        // Only the root graph is copied, the subgraphs contain
        // only references to the global vertices they own.
        typename uninduced_subgraph<Graph>::children_iterator i,i_end;
        boost::tie(i,i_end) = x.children();
        for(; i != i_end; ++i)
        {
            uninduced_subgraph<Graph> child = this->create_subgraph();
            child = *i;
            vertex_iterator vi,vi_end;   
            boost::tie(vi,vi_end) = vertices(*i);
            for (;vi!=vi_end;++vi)  
            {
                add_vertex(*vi,child);
            }
            edge_iterator ei,ei_end;
            boost::tie(ei,ei_end) = edges(*i);
            for (;ei!=ei_end;++ei)  
            {
                add_edge(*ei,child);
            }
       }
    }
    
    // Create a subgraph with the specified edge set.
    template <typename EdgeIterator>
    uninduced_subgraph<Graph>& create_subgraph(EdgeIterator first, EdgeIterator last) {
        m_children.push_back(new uninduced_subgraph<Graph>());
        m_children.back()->m_parent = this;
        for(; first != last; ++first) {
            add_edge(*first, *m_children.back());
        }
        return *m_children.back();
    }
    
    // Create a subgraph with the specified edge set.
    template <typename VertexIterator, typename EdgeIterator>
    uninduced_subgraph<Graph>& create_subgraph(VertexIterator v_first, VertexIterator v_last, EdgeIterator e_first, EdgeIterator e_last) {
        m_children.push_back(new uninduced_subgraph<Graph>());
        m_children.back()->m_parent = this;
        for(; v_first != v_last; ++v_first) {
            add_vertex(*v_first, *m_children.back());
        }
        for(; e_first != e_last; ++e_first) {
            add_edge(*e_first, *m_children.back());
        }
        return *m_children.back();
    }

};

//===========================================================================
// Functions special to the Subgraph Class

template <typename G>
typename uninduced_subgraph<G>::vertex_descriptor
add_vertex(typename uninduced_subgraph<G>::vertex_descriptor u_global,
           uninduced_subgraph<G>& g)
{
    BOOST_ASSERT(!g.is_root());
    typename uninduced_subgraph<G>::vertex_descriptor u_local;

    u_local = add_vertex(g.m_graph);
    g.m_global_vertex.push_back(u_global);
    g.m_local_vertex[u_global] = u_local;
    
    return u_local;
}

//===========================================================================
// Functions required by the MutableGraph concept

namespace detail {

    template <typename Vertex, typename Graph>
    std::pair<typename uninduced_subgraph<Graph>::edge_descriptor, bool>
    add_edge_recur_up(Vertex u_global, Vertex v_global,
                      const typename Graph::edge_property_type& ep,
                      uninduced_subgraph<Graph>& g, uninduced_subgraph<Graph>* orig)
    {
        if(g.is_root()) {
            typename uninduced_subgraph<Graph>::edge_descriptor e_global;
            bool inserted;
            boost::tie(e_global, inserted) = add_edge(u_global, v_global, ep, g.m_graph);
            put(edge_index, g.m_graph, e_global, g.m_edge_counter++);
            g.m_global_edge.push_back(e_global);
            //children_add_edge(u_global, v_global, e_global, g.m_children, orig);
            return std::make_pair(e_global, inserted);
        } else {
            typename uninduced_subgraph<Graph>::edge_descriptor e_global;
            bool inserted;
            boost::tie(e_global, inserted) = 
                    add_edge_recur_up(u_global, v_global, ep, *g.m_parent, orig);
            
            // insert edge into the current uninduced_subgraph on the way down
            g.local_add_edge(g.global_to_local(u_global), g.global_to_local(v_global), e_global);
            return std::make_pair(e_global, inserted);
        }
    }

} // namespace detail

// Add an edge to the uninduced_subgraph g, specified by the local vertex descriptors u
// and v. In addition, the edge will be added to any (all) parent uninduced_subgraphs.

template <typename G>
std::pair<typename uninduced_subgraph<G>::edge_descriptor, bool>
add_edge(typename uninduced_subgraph<G>::vertex_descriptor u,
         typename uninduced_subgraph<G>::vertex_descriptor v,
         const typename G::edge_property_type& ep,
         uninduced_subgraph<G>& g)
{
    if (g.is_root()) {
        // u and v are really global
        return detail::add_edge_recur_up(u, v, ep, g, &g);
    } else {        
        typename uninduced_subgraph<G>::edge_descriptor e_local, e_global;
        bool inserted;
        boost::tie(e_global, inserted) =
            detail::add_edge_recur_up(g.local_to_global(u),
                                      g.local_to_global(v),
                                      ep, g, &g);
        e_local = g.global_to_local(e_global);
        return std::make_pair(e_local, inserted);
    }
}

template <typename G>
std::pair<typename uninduced_subgraph<G>::edge_descriptor, bool>
add_edge(typename uninduced_subgraph<G>::vertex_descriptor u,
         typename uninduced_subgraph<G>::vertex_descriptor v,
         uninduced_subgraph<G>& g)
{ return add_edge(u, v, typename G::edge_property_type(), g); }


template <typename G>
std::pair<typename uninduced_subgraph<G>::edge_descriptor, bool>
add_edge(typename uninduced_subgraph<G>::edge_descriptor e_global,
           uninduced_subgraph<G>& g)
{
    BOOST_ASSERT(!g.is_root());
    typename uninduced_subgraph<G>::vertex_descriptor u_local, v_local;
    typename uninduced_subgraph<G>::edge_descriptor e_local;
    // add vertices to this subgraph first
    u_local = add_vertex(source(e_global, g.root()), g);
    v_local = add_vertex(target(e_global, g.root()), g);
    // create this local edge
    e_local = add_edge(u_local, v_local, g);

    return e_local;
}

template <typename G>
void
remove_edge(typename uninduced_subgraph<G>::vertex_descriptor u,
            typename uninduced_subgraph<G>::vertex_descriptor v,
            uninduced_subgraph<G>& g)
{
    if(g.is_root()) {
        detail::remove_edge_recur_down(u, v, g);
    } else {
        detail::remove_edge_recur_down(g.local_to_global(u),
                                     g.local_to_global(v), g);
    }
}

template <typename G>
void
remove_edge(typename uninduced_subgraph<G>::edge_descriptor e, uninduced_subgraph<G>& g)
{
    remove_edge(source(e, g), target(e, g), g);
}

}

#endif	/* UNINDUCED_SUBGRAPH_HPP */
