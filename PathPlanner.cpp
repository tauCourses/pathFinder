#include "PathPlanner.h"


FaceNode::FaceNode(FT distance, Face_handle face, Face_handle prevFace, Halfedge_handle prevEdge) :
        distance(distance), face(face), prevFace(prevFace), prevEdge(prevEdge) {}

bool CmpfaceNodePtrs::operator()(const FaceNode *lhs, const FaceNode *rhs) const {
    if(lhs->distance != rhs->distance)
        return lhs->distance < rhs->distance;

    if(lhs->prevEdge == Halfedge_handle())
        return rhs->prevEdge != Halfedge_handle();
    else if(rhs->prevEdge == Halfedge_handle())
        return lhs->prevEdge != Halfedge_handle();
    if (lhs->prevEdge->source()->point() != rhs->prevEdge->source()->point())
        return lhs->prevEdge->source()->point() < rhs->prevEdge->source()->point();
    return lhs->prevEdge->target()->point() < rhs->prevEdge->target()->point();
}

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

Point_2 PathPlanner::midPoint(Point_2 a, Point_2 b)
{
    return {(a.x() + b.x()) /2, (a.y()+b.y())/2 };
}

FT PathPlanner::edgeDistance(Halfedge_handle a, Halfedge_handle b)
{
    FT distance;
    Point_2 a_point, b_point = midPoint(b->source()->point(), b->target()->point());
    if(a == Halfedge_handle())
        a_point = Point_2(this->start);
    else
        a_point = midPoint(a->source()->point(), a->target()->point());

    distance = (a_point.x() - b_point.x()) * (a_point.x() - b_point.x()) +
               (a_point.y() - b_point.y()) * (a_point.y() - b_point.y());

    return sqrt(CGAL::to_double(distance));;
}

void PathPlanner::addFacesToQueue(Arrangement_2 &arr, FaceNode* faceNode) {
    ccb_haledge_circulator first = faceNode->face->outer_ccb();
    ccb_haledge_circulator circ = first;
    do {
        Halfedge_handle tempEdge = arr.non_const_handle(circ);
        Face_handle twinFace = arr.non_const_handle(tempEdge->twin()->face());

        if(!twinFace->contained())
            continue;
        auto search = facesMap.find(twinFace);
        if(search != facesMap.end()) //if face already exist try to improve
        {
            FaceNode* temp = &(facesMap[twinFace]);
            if(temp->processed)
                continue;
            FT tempDistance = faceNode->distance + edgeDistance(faceNode->prevEdge, tempEdge);
            if(tempDistance < temp->distance)
            {
                queue.erase(temp); //remove from set beacuse it's in wrong position
                temp = &(facesMap[twinFace]);
                temp->distance = tempDistance;
                temp->prevFace = faceNode->face;
                temp->prevEdge = tempEdge;
                queue.insert(temp); //insert in the right position
            }
        } else {
            FT tempDistance = faceNode->distance + edgeDistance(faceNode->prevEdge, tempEdge);
            facesMap[twinFace] = FaceNode(tempDistance, Face_handle(twinFace), Face_handle(faceNode->face), Halfedge_handle(tempEdge));
            this->queue.insert(&(facesMap[twinFace]));
        }
    } while (++circ != first);
}

void PathPlanner::setFacesPath(Arrangement_2& arr) // run BFS from start_face to end_face
{
    Landmarks_pl pl(arr);
    start_face = get_face(arr, pl, this->start);
    end_face = get_face(arr, pl, this->end);

    facesMap[start_face] = FaceNode(0,start_face, Face_handle(),  Halfedge_handle());
    this->queue.insert(&(facesMap[start_face]));
    while(!queue.empty())
    {
        FaceNode* faceNode = *(queue.begin());
        if(faceNode->face == end_face)
            return;
        faceNode->processed = true;
        this->addFacesToQueue(arr, faceNode);
        queue.erase(faceNode);
    }
    throw "no path found!";
}

vector<Point_2> PathPlanner::reversedPath(Arrangement_2& arr, Kernel& ker){ //create path from BFS results
    vector<Point_2> path;
    Kernel::Construct_midpoint_2  midp = ker.construct_midpoint_2_object();
    path.push_back(this->end);

    Face_handle f = this->end_face;
    do {
        FaceNode* temp = &(facesMap[f]);
        Halfedge_handle he = temp->prevEdge;
        if (he != Halfedge_handle())
            path.push_back(midp(he->source()->point(), he->target()->point()));
        f = temp->prevFace;
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

void PathPlanner::addFrame(Arrangement_2 &arr) {
    FT mostLeft = this->start.x() < this->end.x() ? this->start.x() : this->end.x();
    FT mostRight = this->start.x() < this->end.x() ? this->end.x() : this->start.x();
    FT mostUp = this->start.y() < this->end.y() ? this->end.y() : this->start.y();
    FT mostDown = this->start.y() < this->end.y() ? this->start.y() : this->end.y();

    for (Arr_VrtxCIter vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
    {
        if(vit->point().x() < mostLeft)
            mostLeft = vit->point().x();
        if(vit->point().x() > mostRight)
            mostRight = vit->point().x();
        if(vit->point().y() < mostDown)
            mostDown = vit->point().y();
        if(vit->point().y() > mostUp)
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
