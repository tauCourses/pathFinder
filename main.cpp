#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <boost/timer.hpp>
#include "CGAL_defines.h"

#include "Path.h"
#include "vertical_decomposition.h"

using namespace std;

void print_ccb (Arrangement_2::Ccb_halfedge_const_circulator circ);

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
list<ArrFaceCHandle>
getPathFaces( const Point_2&        start,
              const Point_2&        end,
              const Arrangement_2&  arr )
{
  list<Point_2> pts;
  pts.push_back(start);
  pts.push_back(end);
  std::list<QueryResult>  queryRes;
  locate(arr, pts.begin(), pts.end(), std::back_inserter(queryRes));
  list<QueryResult>::const_iterator it;
  int i = 0;
  list<ArrFaceCHandle>  res;
  for( it = queryRes.begin(); it != queryRes.end(); ++it, ++i)
  {
    //cout << "The point (" << it->first << ") is located ";
    if( const ArrFaceCHandle* pResFace = boost::get<ArrFaceCHandle>(&(it->second)) )
    {
      res.push_back( *pResFace );
      cout << "Face found" << endl;
    }
    else if( const ArrHedgeCHandle* e = boost::get<ArrHedgeCHandle>(&(it->second)))
    {
      ArrFaceCHandle hLeftFace  = (*e)->face();
      ArrFaceCHandle hRightFace = (*e)->twin()->face();
      res.push_back( hLeftFace );
      res.push_back( hRightFace );
      const Point_2 src = (*e)->source()->point();
      const Point_2 trg = (*e)->target()->point();
      double x0 = CGAL::to_double(src[0]);
      double y0 = CGAL::to_double(src[1]);
      double x1 = CGAL::to_double(trg[0]);
      double y1 = CGAL::to_double(trg[1]);
      cout << "on an edge: [" << x0 << ", " << y0 
           << "] -> [" << x1 << ", " << y1 << "]" << std::endl;
    }
    else if (const ArrVrtxCHandle* v = 
                                 boost::get<ArrVrtxCHandle>(&(it->second)))
    { 
      cout << "on "
                << (((*v)->is_isolated()) ? "an isolated" : "a")
                << " vertex: " << (*v)->point() << std::endl;       
    }
  }
  return res;
}

//----------------------------------------------------------------------------
void updateMinMax( const Point_2& pt, 
                   double& minx, double& miny,
                   double& maxx, double& maxy  )
{
  double x = CGAL::to_double(pt[0]);
  double y = CGAL::to_double(pt[1]);
  minx = min(x, minx);
  maxx = max(x, maxx);
  miny = min(y, miny);
  maxy = max(x, maxy);
}

void updateBoundaryRanges( ArrCCBHedgeCCirc circ, 
                           double& minx, double& miny,
                           double& maxx, double& maxy )
{
  ArrCCBHedgeCCirc curr = circ;
  do {
    const HEdge& he = *curr;//->handle();
    const Point_2& pt = he.source()->point();
    updateMinMax( pt, minx, miny, maxx, maxy );
  } while (++curr != circ);
}
//----------------------------------------------------------------------------
void addOuterRectBoundary( Arrangement_2&  arr, 
                           const Polygon_2& robot,
                           const Point_2&  start, 
                           const Point_2&  end )
{
  ArrFaceCHandle ff = arr.fictitious_face();
  Bbox_2 robotBbox = robot.bbox();
  double minx = robotBbox.xmin();
  double miny = robotBbox.ymin();
  double maxx = robotBbox.xmax();
  double maxy = robotBbox.ymax();
  double rx = fabs(maxx - minx);
  double ry = fabs(maxy - miny);

  updateMinMax( start, minx, miny, maxx, maxy );
  updateMinMax( end, minx, miny, maxx, maxy );

  Arrangement_2::Hole_const_iterator  hi;
  for( hi = ff->holes_begin(); hi != ff->holes_end(); ++hi ) 
    updateBoundaryRanges( *hi, minx, miny, maxx, maxy );
  minx -= rx;
  maxx += rx;
  miny -= ry;
  maxy += ry; 
  Polygon_2 outb;
  for( int i = 0; i < 4; ++i )
    outb.push_back( Point_2( (i < 2) ? minx:maxx, 
                             (i == 1 || i == 2) ? maxy:miny ) );
  CGAL::insert(arr, outb.edges_begin(), outb.edges_end());
}
//----------------------------------------------------------------------------
// Iterating through a DCEL face boundary
void print_ccb( ArrCCBHedgeCCirc circ )
{
  ArrCCBHedgeCCirc curr = circ;
  int k = 1;
  do{
    ++k;
  }while(++curr != circ);
  cout << k << " ";
  const Point_2& pt = curr->source()->point();
  cout << CGAL::to_double(pt[0]) << " " 
       << CGAL::to_double(pt[1]) << " ";
  do {
    const HEdge& he = *curr;//->handle();
    const Point_2& pt = he.target()->point();
    cout << CGAL::to_double(pt[0]) << " " 
         << CGAL::to_double(pt[1]) << " ";
  } while (++curr != circ);
  std::cout << std::endl;
}


bool print_arr_face(const ArrFaceCHandle&  f)
{
  // Print the outer boundary.
  if (f->is_unbounded())
    cout << "Unbounded face. " << std::endl;
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
list<ArrFaceCHandle>::const_iterator 
getTwoInstanceFace( const list<ArrFaceCHandle>& faces )
{
  list<ArrFaceCHandle>::const_iterator i = faces.begin();
  for( ; i != faces.end(); ++i )
  {
    list<ArrFaceCHandle>::const_iterator j = i;
    for( ++j; j != faces.end(); ++j )
      if( *i == *j )
        return i;
  }
  return faces.end();
}

//----------------------------------------------------------------------------
Point_2 getCenter( Arrangement_2::Face CurrFace )
{
  ArrCCBHedgeCCirc circ = CurrFace.outer_ccb();
  cout << "Here2" << endl;

  ArrCCBHedgeCCirc curr = circ;
  vector<Point_2> pts;
  do {
    const HEdge he = *curr;//->handle();
    cout << "Here 3" << endl;
    Arrangement_2::Vertex v = *(he.target());
    Point_2 pt = v.point();
    cout << "Here 4" << endl;
    pts.push_back( pt );
  } while (++curr != circ);
  return CGAL::centroid( pts.begin(), pts.end() );  
}
//----------------------------------------------------------------------------
void createGraph( Arrangement_2&   arrConfSpace, 
                  Arrangement_2&   arrTrpzSplit,
                  ArrFaceCHandle&  hPathFace )
{
  list<Point_2> TrpzCntrs;
  Arrangement_2::Face_const_iterator iTrpzFaces = arrTrpzSplit.faces_begin(); 
  cout << "arrTrpzSplit.faces() = " << arrTrpzSplit.number_of_faces() << endl;
  for( ; iTrpzFaces != arrTrpzSplit.faces_end(); ++iTrpzFaces )
  {
    Arrangement_2::Face f = *iTrpzFaces;
    TrpzCntrs.push_back( getCenter( f ) );
  }
  cout << "TrpzCntrs = " << TrpzCntrs.size() << endl;
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
  Arrangement_2 arrConfSpace(&traits);
  for( const Polygon_with_holes_2& pwh : vecConfObst)
    insertPolygonIntoArng( arrConfSpace, pwh );

  list<ArrFaceCHandle> arrFaces = getPathFaces( start, end, arrConfSpace );
  list<ArrFaceCHandle>::const_iterator iPathFace = getTwoInstanceFace(arrFaces);
  bool bPathExists = iPathFace != arrFaces.end();
  if( !bPathExists )
    return vector<Point_2>();

  ArrFaceCHandle hPathFace = *iPathFace;
  bool bPathFaceIsUnbounded = hPathFace->is_unbounded();
  if( bPathFaceIsUnbounded )
      addOuterRectBoundary( arrConfSpace, robot, start, end );

  //print_arr_face(*(arrFaces.first));
  //print_arr_face(*(arrFaces.second));

  cout << "Query points are in the same Arrangement face. "
       << "A solution exists" << endl;
  //
  //... to be continued...
  //
  Kernel* kernel = &traits;
  Arrangement_2 arrTrpzSplit(arrConfSpace);
  vertical_decomposition(arrTrpzSplit, *kernel);

  list<ArrFaceCHandle> vertArrFaces = getPathFaces(start, end, arrTrpzSplit);
  print_arr_face(*vertArrFaces.begin());
  print_arr_face(*(++vertArrFaces.begin()));
  print_arr_face(*(++(++vertArrFaces.begin())));

  createGraph(arrConfSpace, arrTrpzSplit, hPathFace );
  //print_arr_face(*(vertArrFaces.second));
  //cout << *(vertArrFaces.first) << endl;
  cout << "==========================" << endl;
  //cout << arr << endl;
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
