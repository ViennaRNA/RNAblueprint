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
    typedef graph_traits<Graph> Traits;
    typedef std::list<uninduced_subgraph<Graph>*> ChildrenList;
public:
    // Graph requirements
    typedef typename Traits::vertex_descriptor         vertex_descriptor;
    typedef typename Traits::edge_descriptor           edge_descriptor;
    typedef typename Traits::directed_category         directed_category;
    typedef typename Traits::edge_parallel_category    edge_parallel_category;
    typedef typename Traits::traversal_category        traversal_category;

    // IncidenceGraph requirements
    typedef typename Traits::out_edge_iterator         out_edge_iterator;
    typedef typename Traits::degree_size_type          degree_size_type;

    // AdjacencyGraph requirements
    typedef typename Traits::adjacency_iterator        adjacency_iterator;

    // VertexListGraph requirements
    typedef typename Traits::vertex_iterator           vertex_iterator;
    typedef typename Traits::vertices_size_type        vertices_size_type;

    // EdgeListGraph requirements
    typedef typename Traits::edge_iterator             edge_iterator;
    typedef typename Traits::edges_size_type           edges_size_type;

    typedef typename Traits::in_edge_iterator          in_edge_iterator;

    typedef typename edge_property_type<Graph>::type   edge_property_type;
    typedef typename vertex_property_type<Graph>::type vertex_property_type;
    typedef subgraph_tag                               graph_tag;
    typedef Graph                                      graph_type;
    typedef typename graph_property_type<Graph>::type  graph_property_type;
    
    uninduced_subgraph() : subgraph<Graph>() {  }
    
    uninduced_subgraph(const graph_property_type& p) : subgraph<Graph>(p) {  }

    uninduced_subgraph(vertices_size_type n, const graph_property_type& p = graph_property_type()) : subgraph<Graph>(n, p) {  }
    
    // copy constructor
    uninduced_subgraph(const uninduced_subgraph& x) {
        subgraph<Graph>::m_parent = x.m_parent;
        subgraph<Graph>::m_edge_counter = x.m_edge_counter;
        subgraph<Graph>::m_global_vertex = x.m_global_vertex;
        subgraph<Graph>::m_global_edge = x.m_global_edge;
        
        if(x.is_root()) {
            subgraph<Graph>::m_graph = x.m_graph;
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
                add_vertex(*vi, child);
            }
            edge_iterator ei,ei_end;
            boost::tie(ei,ei_end) = edges(*i);
            for (;ei!=ei_end;++ei)  
            {
                add_edge(*ei, child);
            }
       }
    }
    
    // Create a subgraph
    uninduced_subgraph<Graph>& create_subgraph() {
        m_children.push_back(new uninduced_subgraph<Graph>());
        m_children.back()->m_parent = this;
        return *m_children.back();
    }

    // Create a subgraph with the specified vertex set.
    template <typename VertexIterator>
    uninduced_subgraph<Graph>& create_subgraph(VertexIterator first, VertexIterator last) {
        m_children.push_back(new uninduced_subgraph<Graph>());
        m_children.back()->m_parent = this;
        for(; first != last; ++first) {
            add_vertex(*first, *m_children.back());
        }
        return *m_children.back();
    }
    
    // Create a subgraph with the specified edge set.
    template <typename EdgeIterator>
    uninduced_subgraph<Graph>& create_subgraph(std::pair<EdgeIterator, EdgeIterator> its) {
        m_children.push_back(new uninduced_subgraph<Graph>());
        m_children.back()->m_parent = this;
        for(; its.first != its.second; ++its.first) {
            add_edge(*its.first, *m_children.back());
        }
        return *m_children.back();
    }
    
    // Create a subgraph with the specified edge set.
    template <typename VertexIterator, typename EdgeIterator>
    uninduced_subgraph<Graph>& create_subgraph(std::pair<VertexIterator, VertexIterator> v_its, std::pair<EdgeIterator, EdgeIterator> e_its) {
        m_children.push_back(new uninduced_subgraph<Graph>());
        m_children.back()->m_parent = this;
        for(; v_its.first != v_its.second; ++v_its.first) {
            add_vertex(*v_its.first, *m_children.back());
        }
        for(; e_its.first != e_its.second; ++e_its.first) {
            add_edge(*e_its.first, *m_children.back());
        }
        return *m_children.back();
    }
    
    // Return the children subgraphs of this graph/subgraph.
    // Use a list of pointers because the VC++ std::list doesn't like
    // storing incomplete type.
    typedef indirect_iterator<
        typename ChildrenList::const_iterator
      , uninduced_subgraph<Graph>
      , std::bidirectional_iterator_tag
    >
    children_iterator;

    typedef indirect_iterator<
        typename ChildrenList::const_iterator
      , uninduced_subgraph<Graph> const
      , std::bidirectional_iterator_tag
    >
    const_children_iterator;

    std::pair<const_children_iterator, const_children_iterator> children() const {
      return std::make_pair(const_children_iterator(m_children.begin()),
                            const_children_iterator(m_children.end()));
    }

    std::pair<children_iterator, children_iterator> children() {
      return std::make_pair(children_iterator(m_children.begin()),
                            children_iterator(m_children.end()));
    }
    
    // Return the parent graph.
    uninduced_subgraph& parent() { return static_cast<uninduced_subgraph&>(*subgraph<Graph>::m_parent); }
    const uninduced_subgraph& parent() const { static_cast<uninduced_subgraph&>(*subgraph<Graph>::m_parent); }

public: // Needs new declaration
    ChildrenList m_children;
    
};

//===========================================================================
// Functions special to the Subgraph Class

namespace detail {
    
    template <typename G>
    typename uninduced_subgraph<G>::vertex_descriptor
    add_vertex_recur_up(typename uninduced_subgraph<G>::vertex_descriptor u_global,
                uninduced_subgraph<G>& g)
    {
        if (!g.is_root()) {
            if (!g.find_vertex(u_global).second) {
                typename uninduced_subgraph<G>::vertex_descriptor u_local;
                
                detail::add_vertex_recur_up(u_global, g.parent());
                
                u_local = add_vertex(g.m_graph);
                g.m_global_vertex.push_back(u_global);
                g.m_local_vertex[u_global] = u_local;
                
                return u_local;
            } else {
                return g.find_vertex(u_global).first;
            }
        } else {
            return u_global;
        }
    }
    
} // namespace detail

template <typename G>
typename uninduced_subgraph<G>::vertex_descriptor
add_vertex(typename uninduced_subgraph<G>::vertex_descriptor u_global,
           uninduced_subgraph<G>& g)
{
    BOOST_ASSERT(!g.is_root());
    typename uninduced_subgraph<G>::vertex_descriptor u_local;

    u_local = detail::add_vertex_recur_up(u_global, g);
    
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
            return std::make_pair(e_global, inserted);
        } else {
            typename uninduced_subgraph<Graph>::edge_descriptor e_global;
            bool inserted;
            boost::tie(e_global, inserted) = 
                    add_edge_recur_up(u_global, v_global, ep, static_cast<uninduced_subgraph<Graph>&>(*g.m_parent), orig);
            
            // insert edge into the current uninduced_subgraph on the way down
            std::cerr << "added local edge: " << u_global << "->" << v_global << std::endl;
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
typename uninduced_subgraph<G>::edge_descriptor
add_edge(typename uninduced_subgraph<G>::edge_descriptor e_global,
           uninduced_subgraph<G>& g)
{
    if(g.is_root()) {
        return e_global;
    } else {
        typename uninduced_subgraph<G>::edge_descriptor e_local;
        typename uninduced_subgraph<G>::vertex_descriptor u_local, v_local;
        bool exists;
        boost::tie(e_local, exists) = g.find_edge(e_global);

        if (!exists) {
            // recursion up!
            e_local = add_edge(e_global, static_cast<uninduced_subgraph<G>&>(*g.m_parent));
            // add vertices to this subgraph first
            u_local = add_vertex(source(e_global, g.root()), g);
            v_local = add_vertex(target(e_global, g.root()), g);
            // insert edge into the current uninduced_subgraph on the way down
            e_local = g.local_add_edge(u_local, v_local, e_global);
        }
        return e_local;
    }
}
    
namespace detail {
    template <typename Vertex, typename Graph>
    void remove_edge_recur_down(Vertex u_global, Vertex v_global,
                                uninduced_subgraph<Graph>& g)
    {
        Vertex u_local, v_local;
        u_local = g.m_local_vertex[u_global];
        v_local = g.m_local_vertex[v_global];
        remove_edge(u_local, v_local, g.m_graph);
        children_remove_edge(u_global, v_global, g.m_children);
    }
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
    typename subgraph<G>::edge_descriptor e_global = g.local_to_global(e);
    remove_edge(e, g.m_graph); // kick edge from this graph
    detail::children_remove_edge<G>(e_global, g.m_children);
    
}

}

#endif	/* UNINDUCED_SUBGRAPH_HPP */
