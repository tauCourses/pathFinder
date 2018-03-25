#include "Path.h"

Path::Path(const vector<Point_2> &path) : _path(path) {}

Polygon_2 Path::hidamariSketchRobot(const Segment_2 &segment, const Polygon_2 &robot) const {
    Polygon_2 robot1(robot);
    if (!robot1.is_clockwise_oriented()) robot1.reverse_orientation();
    Polygon_2 res;
    Vector_2 DirectionVector = (segment.target() - segment.source()).perpendicular(CGAL::COUNTERCLOCKWISE);
    FT maxValue, minValue;
    Polygon_2::iterator maxPoint, minPoint;
    bool first = true;
    auto currentVertex = robot1.vertices_begin();
    for (; currentVertex != robot1.vertices_end(); ++currentVertex) {
        FT value = (*currentVertex - CGAL::ORIGIN) * DirectionVector;
        if (first || minValue > value) {
            minValue = value;
            minPoint = currentVertex;
        }
        if (first || maxValue < value) {
            maxValue = value;
            maxPoint = currentVertex;
        }
        first = false;
    }
    auto firstVertexVector = (*robot1.vertices_begin() - CGAL::ORIGIN) - (segment.source() - CGAL::ORIGIN);
    currentVertex = robot1.vertices_begin();
    while (currentVertex != maxPoint) ++currentVertex;
    res.push_back(*currentVertex - firstVertexVector);
    Vector_2 currentShift = (segment.target() - segment.source());
    res.push_back(*currentVertex - firstVertexVector + currentShift);
    ++currentVertex;
    while (currentVertex != maxPoint) {
        if (currentVertex == robot1.vertices_end()) {
            currentVertex = robot1.vertices_begin();
            continue;
        }
        if (currentVertex == minPoint) {
            res.push_back(*currentVertex + currentShift - firstVertexVector);
            currentShift = {0,0};
        }
        res.push_back(*currentVertex + currentShift - firstVertexVector);
        ++currentVertex;
    }
    return res;
}

vector<Polygon_2> Path::strechRobotToSegment(const Segment_2 &segment, const Polygon_2 &robot) {
    list<CGAL::Partition_traits_2<Kernel>::Polygon_2> partition_polys;
    Polygon_2 robot1(robot);
    if (robot1.is_clockwise_oriented()) robot1.reverse_orientation();
    CGAL::approx_convex_partition_2(robot1.vertices_begin(),
                                    robot1.vertices_end(),
                                    std::back_inserter(partition_polys));
    vector<Polygon_2> realPartitions;
    for (auto &p:partition_polys) {
        Polygon_2 temp;
        for (auto it = p.vertices_begin(); it != p.vertices_end(); ++it) temp.push_back({it->x(), it->y()});
        realPartitions.push_back(hidamariSketchRobot(segment, temp));
    }
    return realPartitions;
}

bool Path::verifyLine(const Segment_2 &segment, const Polygon_2 &robot,
                      vector<Polygon_2> &obstacles) {
    vector<Polygon_2> wideRobot = strechRobotToSegment(segment, robot);
    bool res = true;
    for (size_t i = 0; i < obstacles.size(); ++i) {
        if (obstacles.at(i).is_clockwise_oriented()) obstacles.at(i).reverse_orientation();
        bool res1 = false;
        for (auto &p :wideRobot) {
            if (p.is_clockwise_oriented()) p.reverse_orientation();
            if (CGAL::do_intersect(p, obstacles.at(i))) {
                res1 = true;
                break;
            }
        }
        if (!obstacles.at(i).is_clockwise_oriented())obstacles.at(i).reverse_orientation();
        if (res1) {
            _badObstacles.push_back((int) i);
            res = false;
        }
    }
    return res;
}

bool Path::verify(const Point_2 &start, const Point_2 &end, const Polygon_2 &robot,
                  vector<Polygon_2> &obstacles) {
    if (_path.empty()) {
        cout << "FAILURE: result path is empty!" << endl;
        return false;
    }
    if (start != _path.front() || end != _path.back()) {
        if (start != _path.front()) cout << "FAILURE: result path start is not equal to the start point!" << endl;
        if (start != _path.front()) cout << "FAILURE: result path end is not equal to the end point!" << endl;
        _path.clear();
        return false;
    }
    for (size_t i = 0; i < _path.size() - 1; ++i)
        if (!verifyLine({_path.at(i), _path.at(i + 1)}, robot, obstacles)) {
            _badPath.push_back(_path.at(i));
            _badPath.push_back(_path.at(i + 1));
            while (_path.size() > i + 1) _path.pop_back();
            return false;
        }
    return true;
}

