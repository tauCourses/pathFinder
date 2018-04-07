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
        this->invRobot.push_back(Point_2(- ((*vi)[0] - this->robot[0][0]), - ((*vi)[1] - this->robot[0][1])));
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

    Vertex_const_handle vertex;
    if (CGAL::assign(vertex, obj)) {
        Face_iterator it = arr.faces_begin();
        for(;it!=arr.faces_end();it++)
        {
            if(it == arr.unbounded_face())
                continue;

            ccb_haledge_circulator first = it->outer_ccb();
            ccb_haledge_circulator circ = first;
            do {
                Halfedge_const_handle temp = circ;
                if(temp->source()->point() == vertex->point())
                {
                    if(it->contained())
                        return arr.non_const_handle(it);
                    else
                        break;
                }
            } while (++circ != first);
        }
        throw "point is not in a legal position - on a vertex between obstacles faces";
    }

    Halfedge_const_handle  helfEdge; //check it's a halfedge
    if (CGAL::assign(helfEdge, obj)) {
        if (helfEdge->face()->contained())
            return arr.non_const_handle(helfEdge->face());
        else if(helfEdge->twin()->face()->contained())
            return arr.non_const_handle(helfEdge->twin()->face());
        throw "point is not in a legal position - on an edge between two obstacles faces";
    }

    // Check whether the point is contained inside a free bounded face.
    Face_const_handle face;
    if (CGAL::assign(face, obj)) //if obj is face
    {
        if(face->contained())
            return arr.non_const_handle(face);
    }
    throw "point is not in a legal position - inside an obstacle face";
}

void PathPlanner::addFacesToQueue(Arrangement_2 &arr, Face_handle face) {
    ccb_haledge_circulator first = face->outer_ccb();
    ccb_haledge_circulator circ = first;
    do {
        Halfedge_const_handle temp = circ;
        Face_handle twinFace = arr.non_const_handle(temp->twin()->face());
        if(!twinFace->contained())
            continue;
        auto search = lastFaceMap.find(twinFace);
        if(search != lastFaceMap.end()) //if face already exist - we already scan it in BFS no need to do anything
            continue;

        lastFaceMap[twinFace] = face;
        lastEdgeMap[twinFace] = arr.non_const_handle(temp);
        queue.push_back(twinFace);
    } while (++circ != first);
}

void PathPlanner::setFacesPath(Arrangement_2& arr) // run BFS from start_face to end_face
{
    Landmarks_pl pl(arr);
    start_face = get_face(arr, pl, this->start);
    end_face = get_face(arr, pl, this->end);

    lastFaceMap[start_face] = Face_handle();
    lastEdgeMap[start_face] = Halfedge_handle();

    this->queue.push_back(start_face);
    while(!queue.empty())
    {
        Face_handle face = queue.front();
        queue.pop_front();
        if(face == end_face)
            return;

        this->addFacesToQueue(arr, face);
    }
    throw "no path found!";
}

vector<Point_2> PathPlanner::reversedPath(Arrangement_2& arr, Kernel& ker){ //create path from BFS results
    vector<Point_2> path;
    Kernel::Construct_midpoint_2  midp = ker.construct_midpoint_2_object();
    path.push_back(this->end);

    Face_handle f = this->end_face;
    do {
        path.push_back(point_in_vertical_trapezoid(f, arr, ker));
        Halfedge_handle he = lastEdgeMap[f];
        if (he != Halfedge_handle())
            path.push_back(midp(he->source()->point(), he->target()->point()));
        f = this->lastFaceMap[f];
    } while (f != Face_handle());
    path.push_back(this->start);

    return path;
}

vector<Point_2> PathPlanner::planPath() {
    Arrangement_2 arr = this->freeSpace.arrangement(); //set the free space as arrangment
    this->addFrame(arr); //add frame around the arrangement

    Polygon_set_2::Traits_2 traits;
    polygon_split_observer observer; //ensure that when face split two side safe their property(inside/outside)
    observer.attach(arr);
    Kernel* ker = &traits;
    this->verticalDecomposition(arr, *ker);
    observer.detach();

    setFacesPath(arr);
    vector<Point_2> path, reversedPath = this->reversedPath(arr, *ker);
    for(int i = static_cast<int>(reversedPath.size())-1; i >= 0; i--)
        path.push_back(Point_2(reversedPath[i]));

    return path;
}

Point_2 PathPlanner::point_in_vertical_trapezoid(Face_const_handle f, const Arrangement_2& arr, const Kernel& ker)
{
    const typename Arrangement_2::Traits_2::Is_vertical_2 is_vertical =
            arr.traits()->is_vertical_2_object();
    const typename Kernel::Construct_midpoint_2 midpoint =
            ker.construct_midpoint_2_object();

    Halfedge_const_handle he1, he2;
    CGAL::Arr_halfedge_direction direction;
    bool found = false;
    ccb_haledge_circulator first = f->outer_ccb();
    ccb_haledge_circulator circ = first;
    do {
        if (!is_vertical(circ->curve())) {
            // The current edge is not vertical: assign it as either he1 or he2.
            if (!found) {
                he1 = circ;
                direction = he1->direction();
                found = true;
                continue;
            }
            if (circ->direction() != direction) {
                he2 = circ;
                break;
            }
        }
    } while (++circ != first);

    // Take the midpoint of the midpoints of he1 and he2.
    return midpoint(midpoint(he1->source()->point(), he1->target()->point()),
                    midpoint(he2->source()->point(), he2->target()->point()));
}

void PathPlanner::addFrame(Arrangement_2 &arr) {
    FT mostLeft = this->start.x() < this->end.x() ? this->start.x() : this->end.x();
    FT mostRight = this->start.x() < this->end.x() ? this->end.x() : this->start.x();
    FT mostUp = this->start.y() < this->end.y() ? this->end.y() : this->start.y();
    FT mostDown = this->start.y() < this->end.y() ? this->start.y() : this->end.y();

    for (Arr_VrtxCIter vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
    {
        if(vit->point()[0] < mostLeft)
            mostLeft = vit->point().x();
        if(vit->point()[0] > mostRight)
            mostRight = vit->point().x();
        if(vit->point()[1] < mostDown)
            mostDown = vit->point().y();
        if(vit->point()[1] > mostUp)
            mostUp = vit->point().y();
    }
    Point_2     upperLeft(mostLeft - 1, mostUp + 1),
                upperRight(mostRight +1 , mostUp + 1),
                lowerRight(mostRight + 2, mostDown -1),
                lowerLeft(mostLeft - 2, mostDown -1);

    Segment_2   upperBound(upperLeft, upperRight),
                rightBound(upperRight, lowerRight),
                lowerBound(lowerRight, lowerLeft),
                leftBound(lowerLeft, upperLeft);

    Halfedge_handle tempEdge = arr.insert_in_face_interior(upperBound, arr.unbounded_face());
    Vertex_handle startVertex = tempEdge->source();
    tempEdge = arr.insert_from_left_vertex(rightBound, tempEdge->target());
    tempEdge = arr.insert_from_right_vertex(lowerBound, tempEdge->target());
    tempEdge = arr.insert_at_vertices(leftBound, tempEdge->target(), startVertex);

    tempEdge->twin()->face()->set_contained(true);
    arr.unbounded_face()->set_contained(false);
}