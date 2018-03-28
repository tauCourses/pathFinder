#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <boost/timer.hpp>
#include "CGAL_defines.h"

#include "Path.h"

using namespace std;


Point_2 loadPoint_2(std::ifstream &is) {
    Kernel::FT x, y;
    is >> x >> y;
    Point_2 point(x, y);
    return point;
}

Polygon_2 loadPolygon(ifstream &is) {
    size_t polygon_size = 0;
    is >> polygon_size;
    Polygon_2 ret;
    while (polygon_size--)
        ret.push_back(loadPoint_2(is));
    return ret;
}

vector<Polygon_2> loadPolygons(ifstream &is) {
    size_t number_of_polygons = 0;
    is >> number_of_polygons;
    vector<Polygon_2> ret;
    while (number_of_polygons--)
        ret.push_back(loadPolygon(is));
    return ret;
}

Polygon_2 getInverseRobot(const Polygon_2 &robot, const Point_2 &refPoint)
{
    Polygon_2 ret;
    for (VertexIterator vi = robot.vertices_begin(); vi != robot.vertices_end(); ++vi) {
        //movedPoint is the robot Point moved to the new origin(ref point)
        Point_2 movedPoint = {(*vi)[0] + refPoint[0] - robot[0][0] , (*vi)[1] + refPoint[1] - robot[0][1]};
        ret.push_back({movedPoint[0] - 2 * (movedPoint[0] - refPoint[0]), movedPoint[1] - 2 * (movedPoint[1] - refPoint[1])});
    }

    return ret;
}
vector<Point_2> findPath(const Point_2 &start,
                         const Point_2 &end,
                         const Polygon_2 &robot,
                         vector<Polygon_2> &obstacles) {
    Polygon_2 inverseRobot = getInverseRobot(robot, start);

    vector<Polygon_with_holes_2> obstaclesC;
    for(Polygon_2& obstacle: obstacles)
        obstaclesC.push_back(CGAL::minkowski_sum_2(inverseRobot, obstacle));

    cout << "first obstacleC " << obstaclesC[0] << endl;
    Arrangement_2 arr; //is it the right data structure?
    insert(arr, obstaclesC[0]); // not working!
    //create trapzohid
    //createGraph
    //findPathInGraph
    //return pathToPoints


    cout << "inverseRobot " << inverseRobot << endl;
    cout << "start - " << start << endl;
    cout << "first robot point " << robot[0] << endl;
    cout << "end " << end << endl;
    return vector<Point_2>({start,{1.71,5.57},{23.84,5.94},{21.21,29.17}, end});
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        cerr << "[USAGE]: inputRobot inputObstacles outputFile" << endl;
        return 1;
    }

    ifstream inputRobotFile(argv[1]), inputObstaclesFile(argv[2]);
    if (!inputRobotFile.is_open() || !inputObstaclesFile.is_open()) {
        if (!inputRobotFile.is_open()) cerr << "ERROR: Couldn't open file: " << argv[1] << endl;
        if (!inputObstaclesFile.is_open()) cerr << "ERROR: Couldn't open file: " << argv[2] << endl;
        return -1;
    }

    auto startPoint = loadPoint_2(inputRobotFile);
    auto endPoint = loadPoint_2(inputRobotFile);
    if (startPoint == endPoint) {
        cout << "Can't have the same startPoint endPoint" << endl;
        return -1;
    };
    auto robot = loadPolygon(inputRobotFile);
    inputRobotFile.close();

    auto obstacles = loadPolygons(inputObstaclesFile);
    inputObstaclesFile.close();
    boost::timer timer;
    auto result = Path(findPath(startPoint, endPoint, robot, obstacles));
    auto secs = timer.elapsed();
    cout << "Path created:      " << secs << " secs" << endl;
    cout << "Path validation:   " << ((result.verify(startPoint, endPoint, robot, obstacles)) ? "Success!" : "Failure")
         << endl;

    ofstream outputFile;
    outputFile.open(argv[3]);
    if (!outputFile.is_open()) {
        cerr << "ERROR: Couldn't open file: " << argv[3] << endl;
        return -1;
    }
    outputFile << result << endl;
    outputFile.close();
    return 0;
}