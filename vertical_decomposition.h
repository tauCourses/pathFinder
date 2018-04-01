#ifndef VERTICAL_DECOMPOSITION_H
#define VERTICAL_DECOMPOSITION_H

#include <CGAL/Arr_vertical_decomposition_2.h>

// ------------------------------------------------------------------------
// Add a vertical segment from the given vertex to some other arrangement
// feature.
template <typename Arrangement_2, typename Kernel>
typename Arrangement_2::Halfedge_const_handle
add_vertical_segment(Arrangement_2& arr, typename Arrangement_2::Vertex_handle v,
                     CGAL::Object obj, Kernel& ker)
{
  typedef typename Arrangement_2::X_monotone_curve_2        Segment_2;
  
  Segment_2                                      seg;
  typename Arrangement_2::Vertex_const_handle    vh;
  typename Arrangement_2::Halfedge_const_handle  hh;
  typename Arrangement_2::Face_const_handle      fh;
  typename Arrangement_2::Vertex_handle          v2;

  if (CGAL::assign(vh, obj)) {
    // The given feature is a vertex.
    seg = Segment_2(v->point(), vh->point()); 
    v2 = arr.non_const_handle(vh);
  }
  else if (CGAL::assign(hh, obj)) {
    // The given feature is a halfedge. We ignore fictitious halfedges.
    if (hh->is_fictitious())
      return typename Arrangement_2::Halfedge_const_handle();

    // Check whether v lies in the interior of the x-range of the edge (in
    // which case this edge should be split).
    const typename Kernel::Compare_x_2 cmp_x = ker.compare_x_2_object();
    if (cmp_x(v->point(), hh->target()->point()) == CGAL::EQUAL) {
      // In case the target of the edge already has the same x-coordinate as
      // the vertex v, just connect these two vertices.
      seg = Segment_2(v->point(), hh->target()->point());
      v2 = arr.non_const_handle(hh->target());
    }
    else {
      // Compute the vertical projection of v onto the segment associated
      // with the halfedge. Split the edge and connect v with the split point.
      typedef typename Kernel::Line_2 Line;
      Line supp_line(hh->source()->point(), hh->target()->point());
      Line vert_line(v->point(), Point_2(v->point().x(), v->point().y() + 1));
      typename Arrangement_2::Point_2  point;
      CGAL::assign(point, ker.intersect_2_object()(supp_line, vert_line));
      seg = Segment_2(v->point(), point);
      arr.split_edge(arr.non_const_handle(hh),
                     Segment_2(hh->source()->point(), point),
                     Segment_2(point, hh->target()->point()));
      v2 = arr.non_const_handle(hh->target());
    }
  }
  // Ignore faces and empty objects.
  else return typename Arrangement_2::Halfedge_const_handle();

  // Add the vertical segment to the arrangement using its two end vertices.
  return arr.insert_at_vertices(seg, v, v2);
}

// ------------------------------------------------------------------------
// Construct the vertical decomposition of the given arrangement.
template <typename Arrangement_2, typename Kernel>
void vertical_decomposition(Arrangement_2& arr, Kernel& ker)
{
  typedef std::pair<typename Arrangement_2::Vertex_const_handle,
                    std::pair<CGAL::Object, CGAL::Object> > Vd_entry;

  // For each vertex in the arrangment, locate the feature that lies
  // directly below it and the feature that lies directly above it.
  std::list<Vd_entry>            vd_list;
  CGAL::decompose(arr, std::back_inserter(vd_list));

  // Go over the vertices (given in ascending lexicographical xy-order),
  // and add segements to the feautres below and above it.
  const typename Kernel::Equal_2          equal = ker.equal_2_object();
  typename std::list<Vd_entry>::iterator  it, prev = vd_list.end();
  for (it = vd_list.begin(); it != vd_list.end(); ++it) {
    // If the feature above the previous vertex is not the current vertex,
    // add a vertical segment to the feature below the vertex.
    typename Arrangement_2::Vertex_const_handle v;
    if ((prev == vd_list.end()) ||
        !CGAL::assign(v, prev->second.second) ||
        !equal(v->point(), it->first->point()))
      add_vertical_segment(arr, arr.non_const_handle(it->first),
                           it->second.first, ker);
    // Add a vertical segment to the feature above the vertex.
    add_vertical_segment(arr, arr.non_const_handle(it->first),
                         it->second.second, ker);
    prev = it;
  }
}

#endif
