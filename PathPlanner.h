#ifndef PATHFINDER_PATHPLANNER_H
#define PATHFINDER_PATHPLANNER_H

#include <vector>
#include "CGAL_defines.h"

#include <boost/property_map/property_map.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/filtered_graph.hpp>


using namespace std;

class Found_vertex_exception : public std::exception {};

template <typename Vertex> class Find_vertex_visitor {
private:
    Vertex m_goal;
public:
    typedef boost::on_finish_vertex            event_filter;

    Find_vertex_visitor(Vertex v) : m_goal(v) {}

    template <class Graph> void operator()(Vertex v, const Graph& g)
    { if (v == m_goal) throw Found_vertex_exception(); }
};


class Less_than_handle {
    template <typename Type>
    bool operator()(Type s1, Type s2) const { return (&(*s1) < &(*s2)); }
};

class FreeSpaceFace {
    bool operator()(Face_handle f) const
    { f->contained(); }
    bool operator()(Face_const_handle f) const
    { f->contained(); }
};

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
    typedef std::map<Face_handle, Face_handle, Less_than_handle>     Preds_map;
    typedef std::map<Face_handle, Halfedge_handle, Less_than_handle> Edges_map;

    void setInversedRobot();
    void setFreeSpace();

    void verticalDecomposition(Arrangement_2 &arr, Kernel &ker);
    void addVerticalSegment(Arrangement_2 &arr, Vertex_handle v, CGAL::Object obj, Kernel &ker);

    Face_handle get_face(Arrangement_2& arr, const Landmarks_pl &pl, const Point_2 &p);
    void setFacesPath(Arrangement_2& arr);

    Point_2 midPoint(Face_const_handle face);
    Point_2 midPoint(Halfedge_const_handle edge);

public:
    PathPlanner(const Point_2 start, const Point_2 end, const Polygon_2 &robot, vector<Polygon_2> &obstacles);

    vector<Point_2> planPath();

};


#endif //PATHFINDER_PATHPLANNER_H
