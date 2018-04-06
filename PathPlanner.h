#ifndef PATHFINDER_PATHPLANNER_H
#define PATHFINDER_PATHPLANNER_H

#include <vector>
#include <list>
#include <map>

#include "CGAL_defines.h"

/*#include <boost/property_map/property_map.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/filtered_graph.hpp>*/


using namespace std;

class Less_than_handle {
public:
    template <typename Type>
    bool operator()(Type s1, Type s2) const { return (&(*s1) < &(*s2)); }
};

typedef map<Face_handle, Face_handle, Less_than_handle>     Preds_map;
typedef map<Face_handle, Halfedge_handle, Less_than_handle> Edges_map;


class polygon_split_observer : public CGAL::Arr_observer<Arrangement_2>
{
    void after_split_face(Face_handle f1, Face_handle f2, bool);
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
    list<Face_handle> queue;
    Preds_map preds_map;
    Edges_map edges_map;

    void setInversedRobot();
    void setFreeSpace();

    void verticalDecomposition(Arrangement_2 &arr, Kernel &ker);
    void addVerticalSegment(Arrangement_2 &arr, Vertex_handle v, CGAL::Object obj, Kernel &ker);

    Face_handle get_face(Arrangement_2& arr, const Landmarks_pl &pl, const Point_2 &p);
    void setFacesPath(Arrangement_2& arr);

    Point_2 point_in_vertical_trapezoid(Face_const_handle f, const Arrangement_2& arr, const Kernel& ker);
    vector<Point_2> reversedPath(Arrangement_2& arr, Kernel& ker);

    void printFace(Face_handle face);

    void addFaces(Arrangement_2& arr, Face_handle face);

public:
    PathPlanner(const Point_2 start, const Point_2 end, const Polygon_2 &robot, vector<Polygon_2> &obstacles);

    vector<Point_2> planPath();

    void addFrame(Arrangement_2 &arr);

    void printArr(Arrangement_2 &arr);


};


#endif //PATHFINDER_PATHPLANNER_H
