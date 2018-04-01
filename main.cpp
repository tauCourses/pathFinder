#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <boost/timer.hpp>
#include "CGAL_defines.h"

#include "Path.h"
#include "vertical_decomposition.h"

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
    CGAL::Orientation orient = ret.orientation();
    if( CGAL::CLOCKWISE == orient )
      ret.reverse_orientation();
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

//----------------------------------------------------------------------------
void
print_polygon( const Polygon_2& pgn )
{
  cout << pgn.size() << " ";
  VrtxCIter vit = pgn.vertices_begin();
  for( ; vit != pgn.vertices_end(); ++vit )
  {
    double x = CGAL::to_double((*vit)[0]);
    double y = CGAL::to_double((*vit)[1]);
    cout << x << " " << y << " ";
  }
  cout << endl;
}

//----------------------------------------------------------------------------
void
print_polygon( const Polygon_with_holes_2& pwh )
{
  if( !pwh.is_unbounded())
  {
    cout << "{ Outer boundary " << endl;
    print_polygon (pwh.outer_boundary());
  }
  else
    cout << "{ Unbounded polygon." << endl;
  unsigned int k = 1;
  std::cout << pwh.number_of_holes() << " holes:" << std::endl;
  for (HoleCIter hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit, ++k)
  {
    cout << " Hole #" << k << " = " << endl;
    print_polygon (*hit);
  }
  cout << "}" << endl;
}

//----------------------------------------------------------------------------
Polygon_2 getInversedRobot( const Polygon_2& robot,
                            const Point_2&   refPoint )
{
  Polygon_2 ret;
  VrtxCIter vi = robot.vertices_begin();
  for( ; vi != robot.vertices_end(); ++vi)
  {
    //movedPoint is the robot Point moved to the new origin(ref point)
    Point_2 movedPoint = { (*vi)[0] + refPoint[0] - robot[0][0],
                           (*vi)[1] + refPoint[1] - robot[0][1] };
    ret.push_back( { movedPoint[0] - 2 * (movedPoint[0] - refPoint[0]),
                     movedPoint[1] - 2 * (movedPoint[1] - refPoint[1])});
  }
  CGAL::Orientation orient = ret.orientation();
  if( CGAL::CLOCKWISE == orient )
    ret.reverse_orientation();
  return ret;
}

//----------------------------------------------------------------------------
bool insertPolygonIntoArng( Arrangement_2&              arr,
                            const Polygon_with_holes_2& pwh )
{
  if( pwh.is_unbounded() )
    return false;

  const Polygon_2& plg = pwh.outer_boundary();
  CGAL::insert(arr, plg.edges_begin(), plg.edges_end() );
  for( HoleCIter hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit )
  {
    const Polygon_2& plg = *hit;
    CGAL::insert(arr, plg.edges_begin(), plg.edges_end() );
  }
  return true;
}

//----------------------------------------------------------------------------
pair<const ArrFaceCHandle*, const ArrFaceCHandle*>
getPathFaces( const Point_2&        start,
              const Point_2&        end,
              const Arrangement_2&  arr,
              bool&                 bSame )
{
  list<Point_2> pts;
  pts.push_back(start);
  pts.push_back(end);
  std::list<QueryResult>  results;
  locate(arr, pts.begin(), pts.end(), std::back_inserter(results));
  std::list<QueryResult>::const_iterator it;
  bool bProperFace[2];
  int i = 0;
  const ArrFaceCHandle* pResFace = nullptr;
  for( it = results.begin(); it != results.end(); ++it, ++i)
  {
    //cout << "The point (" << it->first << ") is located ";
    if( pResFace = boost::get<ArrFaceCHandle>(&(it->second)) )
      bProperFace[i] = true;
    else
      bProperFace[i] = false;
  }
  it = results.begin();
  ++it;
  bool bSameObject = it->second == results.begin()->second;
  //if( bSameObject && bProperFace[1] )
  //  cout << "inside "
  //       << (((*f)->is_unbounded()) ? "the unbounded" : "a bounded")
  //       << " face." << endl;
  pair<const ArrFaceCHandle*, const ArrFaceCHandle*> res( boost::get<ArrFaceCHandle>(&(results.begin()->second)) ,
                                                          boost::get<ArrFaceCHandle>(&(it->second))  );

  bSame = (bSameObject && bProperFace[1]);
  return res;
}

//----------------------------------------------------------------------------
// Iterating through a DCEL face
// https://doc.cgal.org/latest/Arrangement_on_surface_2/index.html
// Example 2.5
void print_ccb (Arrangement_2::Ccb_halfedge_const_circulator circ)
{
  Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
  std::cout << "(" << curr->source()->point() << ")";
  do {
    const HEdge& he = *curr;//->handle();
    std::cout << " [" << he.curve() << "] "
              << "(" << he.target()->point() << ")";
  } while (++curr != circ);
  std::cout << std::endl;
}


bool print_arr_face(const ArrFaceCHandle&  f)
{
  // Print the outer boundary.
  if (f->is_unbounded())
    std::cout << "Unbounded face. " << std::endl;
  else {
    std::cout << "Outer boundary: ";
    print_ccb (f->outer_ccb());
  }
  // Print the boundary of each of the holes.
  Arrangement_2::Hole_const_iterator hi;
  int index = 1;
  for (hi = f->holes_begin(); hi != f->holes_end(); ++hi, ++index) {
    std::cout << " Hole #" << index << ": ";
    print_ccb (*hi);
  }
}
//----------------------------------------------------------------------------
vector<Point_2> findPath(const Point_2&      start,
                         const Point_2&      end,
                         const Polygon_2&    robot,
                         vector<Polygon_2>&  obstacles)
{
  Polygon_2 invRobot = getInversedRobot( robot, start );
  vector<Polygon_with_holes_2> vecConfObst;
  for(const Polygon_2& obst : obstacles )
  {
    Polygon_with_holes_2 currConfObst = CGAL::minkowski_sum_2( obst,
                                                               invRobot );
    vecConfObst.push_back( currConfObst );
    //print_polygon(currConfObst.outer_boundary());
  }

  CGAL::Arr_segment_traits_2<Kernel> traits;
  Arrangement_2 arr(&traits);
  for( const Polygon_with_holes_2& pwh : vecConfObst)
    insertPolygonIntoArng( arr, pwh );

  bool bSameFace = false;
  pair<const ArrFaceCHandle*, const ArrFaceCHandle*> arrFaces = 
                                    getPathFaces(start, end, arr, bSameFace );
  if( !bSameFace )
    return vector<Point_2>();

  cout << "Query points are in the same Arrangement face. "
       << "A solution exists" << endl;
  //
  //... to be continued...
  //
  Kernel* kernel = &traits;
  vertical_decomposition(arr, *kernel);

  pair<const ArrFaceCHandle*, const ArrFaceCHandle*> vertArrFaces = 
                                     getPathFaces(start, end, arr, bSameFace);
  //print_arr_face(*(vertArrFaces.first));
  //print_arr_face(*(vertArrFaces.second));
  //cout << *(vertArrFaces.first) << endl;
  cout << "==========================" << endl;
  cout << arr << endl;
  return vector<Point_2>({start,{1.71,5.57},{23.84,5.94},{21.21,29.17}, end});
}
//----------------------------------------------------------------------------

int main(int argc, char *argv[]) {
    if (argc != 4) {
        cerr << "[USAGE]: inputRobot inputObstacles outputFile" << endl;
        return 1;
    }

    ifstream inputRobotFile(argv[1]), inputObstaclesFile(argv[2]);
    if (!inputRobotFile.is_open() || !inputObstaclesFile.is_open())
    {
        if (!inputRobotFile.is_open())
            cerr << "ERROR: Couldn't open file: " << argv[1] << endl;
        if (!inputObstaclesFile.is_open())
            cerr << "ERROR: Couldn't open file: " << argv[2] << endl;
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
