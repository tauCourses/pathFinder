//
// Created by t-idkess on 18-Mar-18.
//

#ifndef INC_2_3_CGAL_DEFINES_H
#define INC_2_3_CGAL_DEFINES_H

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Point_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/centroid.h>
#include <CGAL/Arr_batched_point_location.h>


typedef typename CGAL::Gmpq Number_type;
typedef typename CGAL::Cartesian<Number_type> Kernel;
typedef typename Kernel::FT FT;
typedef typename Kernel::Point_2 Point_2;
typedef typename Kernel::Segment_2 Segment_2;
typedef typename Kernel::Vector_2 Vector_2;
typedef typename CGAL::Polygon_2<Kernel> Polygon_2;
typedef typename CGAL::Polygon_with_holes_2<Kernel>         Polygon_with_holes_2;
typedef typename CGAL::Polygon_with_holes_2<Kernel>::Hole_const_iterator 
                                                            HoleCIter;
typedef typename CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>> 
                                                            Arrangement_2;
typedef typename Polygon_2::Vertex_iterator                 VrtxIter;
typedef typename Polygon_2::Vertex_const_iterator           VrtxCIter;
typedef typename Polygon_2::Edge_const_iterator             EdgeIter;
typedef CGAL::Arr_naive_point_location<Arrangement_2>       NaivePL;

typedef CGAL::Arr_point_location_result<Arrangement_2>  PointLocationResult;
typedef std::pair<Point_2, PointLocationResult::Type>   QueryResult;
typedef Arrangement_2::Vertex_const_handle              ArrVrtxCHandle;
typedef Arrangement_2::Halfedge_const_handle            ArrHedgeCHandle;
typedef Arrangement_2::Face_const_handle                ArrFaceCHandle;
typedef Arrangement_2::Face_handle                      ArrFaceHandle;
typedef Arrangement_2::Halfedge                         HEdge;
typedef Arrangement_2::Halfedge_const_handle            ArrHedgeCHandle;
#endif //INC_2_3_CGAL_DEFINES_H
