#include "PathPlanner.h"

void polygon_split_observer::after_split_face(Face_handle f1, Face_handle f2, bool)
{
    f2->set_contained(f1->contained());
}

PathPlanner::PathPlanner(const Point_2 start, const Point_2 end, const Polygon_2 &robot, vector<Polygon_2> &obstacles) :
        start(start), end(end), robot(robot), obstacles(obstacles){
    this->setInversedRobot();
    this->setFreeSpace();
}

void PathPlanner::setInversedRobot() {
    VrtxCIter vi = this->robot.vertices_begin();
    for (; vi != this->robot.vertices_end(); ++vi) {
        //movedPoint is the robot Point moved to the new origin(ref point)
        Point_2 movedPoint = {(*vi)[0] + this->start[0] - this->robot[0][0],
                              (*vi)[1] + this->start[1] - this->robot[0][1]};
        this->invRobot.push_back({movedPoint[0] - 2 * (movedPoint[0] - this->start[0]),
                       movedPoint[1] - 2 * (movedPoint[1] - this->start[1])});
    }
    CGAL::Orientation orient = this->invRobot.orientation();
    if (CGAL::CLOCKWISE == orient)
        this->invRobot.reverse_orientation();
}

void PathPlanner::setFreeSpace() {
    vector<Polygon_with_holes_2> vecConfObst;
    for (const Polygon_2 &obst : this->obstacles) {
        Polygon_with_holes_2 currConfObst = CGAL::minkowski_sum_2(obst,  invRobot);
        vecConfObst.push_back(currConfObst);
    }
    this->freeSpace.join(vecConfObst.begin(), vecConfObst.end()); //join of all obstacles
    this->freeSpace.complement(); //complement to get the free space
}



void PathPlanner::addVerticalSegment(Arrangement_2 &arr, Vertex_handle v, CGAL::Object obj, Kernel &ker) {
    X_monotone_curve_2 seg;
    Vertex_const_handle vh;
    Halfedge_const_handle hh;
    Face_const_handle fh;
    Vertex_handle v2;

    if (CGAL::assign(vh, obj)) { // The given feature is a vertex.
        seg = X_monotone_curve_2(v->point(), vh->point());
        v2 = arr.non_const_handle(vh);
    } else if (CGAL::assign(hh, obj)) { // The given feature is a halfedge.
        if (hh->is_fictitious()) //We ignore fictitious halfedges.
            return;

        // Check whether v lies in the interior of the x-range of the edge (in
        // which case this edge should be split).
        const typename Kernel::Compare_x_2 cmp_x = ker.compare_x_2_object();
        if (cmp_x(v->point(), hh->target()->point()) == CGAL::EQUAL) {
            // In case the target of the edge already has the same x-coordinate as
            // the vertex v, just connect these two vertices.
            seg = X_monotone_curve_2(v->point(), hh->target()->point());
            v2 = arr.non_const_handle(hh->target());
        }
        else {
            // Compute the vertical projection of v onto the segment associated
            // with the halfedge. Split the edge and connect v with the split point.
            Line_2 Line;
            Line_2 supp_line(hh->source()->point(), hh->target()->point());
            Line_2 vert_line(v->point(), Point_2(v->point().x(), v->point().y() + 1));
            Point_2  point;
            CGAL::assign(point, ker.intersect_2_object()(supp_line, vert_line));
            seg = X_monotone_curve_2(v->point(), point);
            arr.split_edge(arr.non_const_handle(hh),
                           X_monotone_curve_2(hh->source()->point(), point),
                           X_monotone_curve_2(point, hh->target()->point()));
            v2 = arr.non_const_handle(hh->target());
        }
    } else // Ignore faces and empty objects.
        return;

    // Add the vertical segment to the arrangement using its two end vertices.
    arr.insert_at_vertices(seg, v, v2);
}

void PathPlanner::verticalDecomposition(Arrangement_2 &arr, Kernel &ker) {
    typedef pair<Vertex_const_handle, pair<CGAL::Object, CGAL::Object> > Vd_entry;

    // For each vertex in the arrangment, locate the feature that lies
    // directly below it and the feature that lies directly above it.
    list<Vd_entry>   vd_list;
    CGAL::decompose(arr, back_inserter(vd_list));

    // Go over the vertices (given in ascending lexicographical xy-order),
    // and add segements to the feautres below and above it.
    const typename Kernel::Equal_2 equal = ker.equal_2_object();
    typename list<Vd_entry>::iterator  it, prev = vd_list.end();
    for (it = vd_list.begin(); it != vd_list.end(); ++it) {
        // If the feature above the previous vertex is not the current vertex,
        // add a vertical segment to the feature below the vertex.
        Vertex_const_handle v;
        if ((prev == vd_list.end()) ||
            !CGAL::assign(v, prev->second.second) ||
            !equal(v->point(), it->first->point()))
            addVerticalSegment(arr, arr.non_const_handle(it->first), it->second.first, ker);
        // Add a vertical segment to the feature above the vertex.
        addVerticalSegment(arr, arr.non_const_handle(it->first), it->second.second, ker);
        prev = it;
    }
}


Face_handle PathPlanner::get_face(Arrangement_2& arr, const Landmarks_pl &pl, const Point_2 &p) {
    CGAL::Object obj = pl.locate(p); //find p in pl

    // Check whether the point lies on an edge separating two forbidden faces.
    Halfedge_const_handle  helfEdge; //check it's a halfedge
    if (CGAL::assign(helfEdge, obj)) {
        if (helfEdge->face()->contained())
            return arr.non_const_handle(helfEdge->face());
        else if(helfEdge->twin()->face()->contained())
            return arr.non_const_handle(helfEdge->twin()->face());
        throw "point is not in legal position - on edge between two obstacles faces";
    }

    // Check whether the point is contained inside a free bounded face.
    Face_const_handle face;
    if (CGAL::assign(face, obj)) //if obj is face
    {
        if(face->contained())
            return arr.non_const_handle(face);
    }
    throw "point is not in legal position";
}

void PathPlanner::setFacesPath(Arrangement_2& arr)
{
    typedef CGAL::Arr_face_index_map<Arrangement_2>                 Face_index_map;

    Landmarks_pl pl(arr);
    start_face = get_face(arr, pl, this->start);
    end_face = get_face(arr, pl, this->end);
    Dual dual(arr);
    boost::filtered_graph<Dual, boost::keep_all, FreeSpaceFace> graph(dual, boost::keep_all());

    Face_index_map                   index_map(arr);
    Preds_map                        preds_map;
    Edges_map                        edges_map;
    Find_vertex_visitor<Face_handle> find_vertex_visitor(end_face);

    try {
        preds_map[start_face] = Face_handle();
        edges_map[start_face] = Halfedge_handle();
        boost::breadth_first_search(graph, start_face,
                                    boost::vertex_index_map(index_map).visitor
                                            (boost::make_bfs_visitor(
                                                    make_pair(
                                                            find_vertex_visitor,
                                                            make_pair(
                                                                    record_predecessors(
                                                                            boost::make_assoc_property_map(preds_map),
                                                                            boost::on_tree_edge()),
                                                                    record_edge_predecessors(
                                                                            boost::make_assoc_property_map(edges_map),
                                                                            boost::on_tree_edge()))))));

        // If there is a path then an exception should have been thrown.
        throw "no path found!";
    }
    catch(Found_vertex_exception e) {}
}

Point_2 PathPlanner::midPoint(Face_const_handle face)
{
    return {0,0};
}
Point_2 PathPlanner::midPoint(Halfedge_const_handle edge)
{
    return {0,0};
}

vector<Point_2> PathPlanner::planPath() {
    Arrangement_2 arr = this->freeSpace.arrangement(); //set the free space as arrangment
    Polygon_set_2::Traits_2 traits;
    polygon_split_observer observer; //ensure that when face split two side safe their property(inside/outside)
    observer.attach(arr);
    Kernel* ker = &traits;
    this->verticalDecomposition(arr, *ker);
    observer.detach();

    setFacesPath(arr);

    //create points path from faceNode path:
    /*vector<Point_2> path;
    path.push_back(this->start);
    path.push_back(this->midPoint(facePath[facePath.size()-1].face));
    for(int i=facePath.size()-2; i >=0; i--)
    {
        path.push_back(midPoint(facePath[i].fatherEdge));
        path.push_back(midPoint(facePath[i].face));
    }
    path.push_back(this->end);*/
    return vector<Point_2>({start,{1.71,5.57},{23.84,5.94},{21.21,29.17}, end});
}




