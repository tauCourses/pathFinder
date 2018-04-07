#ifndef PATHFINDER_PATHPLANNER_H
#define PATHFINDER_PATHPLANNER_H

#include <vector>
#include <list>
#include <map>
#include <set>
#include <math.h>

#include "CGAL_defines.h"

using namespace std;

class FaceNode
{
public:
    explicit FaceNode(FT distance, Face_handle face, Face_handle prevFace, Halfedge_handle prevEdge);
    FaceNode() = default;

    FT distance;
    Face_handle face;
    Face_handle prevFace;
    Halfedge_handle prevEdge;
    bool processed = false;
};

struct CmpfaceNodePtrs
{
    bool operator()(const FaceNode* lhs, const FaceNode* rhs) const
    {
        return lhs->distance < rhs->distance;
    }
};

class polygon_split_observer : public CGAL::Arr_observer<Arrangement_2>
{
    void after_split_face(Face_handle f1, Face_handle f2, bool) override;
};

class PathPlanner {
private:

    const Point_2 &start;
    const Point_2 &end;
    const Polygon_2 &robot;
    vector<Polygon_2> &obstacles;

    Polygon_2 invRobot;
    Polygon_set_2 freeSpace;

    //start and end faces:
    Face_handle start_face;
    Face_handle end_face;

    //bfs maps:
    set<FaceNode*, CmpfaceNodePtrs> queue; //use set because need to delete efficiently
    map<Face_handle, FaceNode> facesMap;

    void setInversedRobot();
    void setFreeSpace();
    void addFrame(Arrangement_2 &arr);

    void verticalDecomposition(Arrangement_2 &arr, Kernel &ker);
    void addVerticalSegment(Arrangement_2 &arr, Vertex_handle v, CGAL::Object obj, Kernel &ker);

    Face_handle get_face(Arrangement_2& arr, const Landmarks_pl &pl, const Point_2 &p);
    void setFacesPath(Arrangement_2& arr);
    void addFacesToQueue(Arrangement_2 &arr, FaceNode* faceNode);
    FT edgeDistance(Halfedge_handle a, Halfedge_handle b);
    Point_2 midPoint(Point_2 a, Point_2 b);

    vector<Point_2> reversedPath(Arrangement_2& arr, Kernel& ker);

    void printFace(Face_handle face);
    void printArr(Arrangement_2 &arr);

public:
    PathPlanner(const Point_2 start, const Point_2 end, const Polygon_2 &robot, vector<Polygon_2> &obstacles);
    vector<Point_2> planPath();

};


#endif //PATHFINDER_PATHPLANNER_H
