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


typedef typename CGAL::Gmpq Number_type;
typedef typename CGAL::Cartesian<Number_type> Kernel;
typedef typename Kernel::FT FT;
typedef typename Kernel::Point_2 Point_2;
typedef typename Kernel::Segment_2 Segment_2;
typedef typename Kernel::Vector_2 Vector_2;
typedef typename CGAL::Polygon_2<Kernel> Polygon_2;
typedef typename CGAL::Polygon_with_holes_2<Kernel> Polygon_with_holes_2;
typedef CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>> Arrangement_2;
typedef typename Polygon_2::Vertex_iterator VertexIterator;

#endif //INC_2_3_CGAL_DEFINES_H
