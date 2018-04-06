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
extern int dijkstra_compute_path(char*, int, int, int, vector<int>& );

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


void print_arr_face(const ArrFaceCHandle&  f)
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

void print_arrangement( const Arrangement_2& arr )
{
  Arrangement_2::Face_const_iterator iFace = arr.faces_begin();
  for( ; iFace != arr.faces_end(); ++iFace )
  {
    if ((*iFace).is_unbounded())
      continue;
    ArrFaceCHandle f = iFace;
    print_arr_face( f );
  }
}

void print_segment( const Segment_2& s )
{
  const Point_2& p1 = s.source();
  const Point_2& p2 = s.target();
  double x1 = CGAL::to_double( p1[0] );
  double y1 = CGAL::to_double( p1[1] );
  double x2 = CGAL::to_double( p2[0] );
  double y2 = CGAL::to_double( p2[1] );
  cout << "[" << x1 << ", " << y1 << " - " << x2 << ", " << y2 << "]" << endl;
}

void print_point( const Point_2& p)
{
  double x1 = CGAL::to_double( p[0] );
  double y1 = CGAL::to_double( p[1] );
  cout << "(" << x1 << ", " << y1 << ")" << endl;
}

//----------------------------------------------------------------------------
Polygon_2 getInversedRobot( const Polygon_2& robot,
                            const Point_2&   refPoint )
{
  Polygon_2 ret;
  VrtxCIter vi = robot.vertices_begin();
/*  
  for( ; vi != robot.vertices_end(); ++vi)
  {
    //movedPoint is the robot Point moved to the new origin(ref point)
    Point_2 movedPoint = { (*vi)[0] + refPoint[0] - robot[0][0],
                           (*vi)[1] + refPoint[1] - robot[0][1] };
    ret.push_back( { movedPoint[0] - 2 * (movedPoint[0] - refPoint[0]),
                     movedPoint[1] - 2 * (movedPoint[1] - refPoint[1])});
  }
*/
  VrtxCIter iVrtx = robot.vertices_begin();
  Kernel::FT xOff = robot[0][0];// - refPoint[0];
  Kernel::FT yOff = robot[0][1];// - refPoint[1];  
  for(; iVrtx != robot.vertices_end(); ++ iVrtx )
  {
    ret.push_back({ -(*iVrtx)[0] + xOff, -(*iVrtx)[1] + yOff });
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
/*
bool insertCcbIntoArng( Arrangement_2& arr,
                        ArrCCBHedgeCCirc circ,
                        list<Segment_2>& ForbiddenEdges )
{
  vector<Segment_2> vecSegm;
  ArrCCBHedgeCCirc curr = circ;
  do {
    const HEdge& he = *curr;
    Point_2 pt1 = he.source()->point();
    Point_2 pt2 = he.target()->point();
    Segment_2 s(pt1, pt2);
    if( pt1[0] == pt2[0] )
        ForbiddenEdges.push_back(s);
    vecSegm.push_back( s );
  } while (++curr != circ);
  CGAL::insert( arr, vecSegm.begin(), vecSegm.end());
  return true;
}

bool insertFaceIntoArng( Arrangement_2&  arr,
                         ArrFaceCHandle  hFace,
                         list<Segment_2>& ForbiddenEdges )
{
  if( hFace->is_unbounded() )
    return false;
  insertCcbIntoArng( arr, hFace->outer_ccb(), ForbiddenEdges );
  
  Arrangement_2::Hole_const_iterator hi;
  for( hi = hFace->holes_begin(); hi != hFace->holes_end(); ++hi) 
    insertCcbIntoArng( arr, *hi, ForbiddenEdges );
  return true; 
}
*/

void getForbiddenEdges( ArrFaceCHandle  hFace,
                        list<Segment_2>& ForbiddenEdges )
{
  Arrangement_2::Hole_const_iterator hi;
  cout << "=== Forbidden Segements === "<< endl;
  for( hi = hFace->holes_begin(); hi != hFace->holes_end(); ++hi) 
  {
    ArrCCBHedgeCCirc circ = *hi;
    ArrCCBHedgeCCirc curr = circ;
    do {
      const HEdge& he = *curr;
      Point_2 pt1 = he.source()->point();
      Point_2 pt2 = he.target()->point();
      Segment_2 s(pt1, pt2);
      if( pt1[0] == pt2[0] )
      {
        ForbiddenEdges.push_back(s);
        print_segment( s );
      }
    } while (++curr != circ);
  }
}
//----------------------------------------------------------------------------
list<ArrFaceCHandle>
getPathFaces( const Point_2&        start,
              const Point_2&        end,
              const Arrangement_2&  arr,
              int*                  pPivotIdx = NULL )
{
  list<Point_2> pts;
  pts.push_back(start);
  pts.push_back(end);
  std::list<QueryResult>  queryRes;
  locate(arr, pts.begin(), pts.end(), std::back_inserter(queryRes));
  list<QueryResult>::const_iterator it;
  int i, nPivotIdx = 0;
  list<ArrFaceCHandle>  res;
  for( it = queryRes.begin(); it != queryRes.end(); ++it, ++i)
  {
    bool bPushBack = it->first == end;
    if( const ArrFaceCHandle* pResFace = boost::get<ArrFaceCHandle>(&(it->second)) )
    {
      ArrFaceCHandle hFace = *pResFace;
      if( bPushBack )
      {
        nPivotIdx = res.size();
        res.push_back( hFace );
      }
      else
      {
        res.insert( res.begin(), hFace );
        nPivotIdx = 1;
      }
    }
    else if( const ArrHedgeCHandle* e = boost::get<ArrHedgeCHandle>(&(it->second)))
    {
      ArrFaceCHandle hLeftFace  = (*e)->face();
      ArrFaceCHandle hRightFace = (*e)->twin()->face();
      if( bPushBack )
      {
        nPivotIdx = res.size();
        res.push_back( hLeftFace );
        res.push_back( hRightFace );
      }
      else
      {
        res.insert( res.begin(), hLeftFace );
        res.insert( res.begin(), hRightFace );
        nPivotIdx = 2;
      }
    }
    else if (const ArrVrtxCHandle* v = 
                                 boost::get<ArrVrtxCHandle>(&(it->second)))
    { 
      cout << "on "
                << (((*v)->is_isolated()) ? "an isolated" : "a")
                << " vertex: " << (*v)->point() << std::endl;
      cout << "Not implemented yet" << endl;
    }
  }
  if( pPivotIdx != NULL )
    *pPivotIdx = nPivotIdx;
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
  maxy = max(y, maxy);
}

void updateBoundaryRanges( ArrCCBHedgeCCirc circ, 
                           double& minx, double& miny,
                           double& maxx, double& maxy )
{
  ArrCCBHedgeCCirc curr = circ;
  do {
    const HEdge& he = *curr;
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
    outb.push_back( Point_2( (i == 0 || i == 3) ? minx:maxx, 
                             (i == 0 || i == 1) ? miny:maxy ) );
  CGAL::insert(arr, outb.edges_begin(), outb.edges_end());
  // Add two additional rectangulars on the left and right.
  Polygon_2 LeftOut;
  LeftOut.push_back( Point_2( minx - rx, miny ) );
  LeftOut.push_back( Point_2( minx, miny ) );
  LeftOut.push_back( Point_2( minx, maxy ) );
  LeftOut.push_back( Point_2( minx - rx, maxy ) );
  CGAL::insert(arr, LeftOut.edges_begin(), LeftOut.edges_end());
  Polygon_2 RightOut;
  RightOut.push_back( Point_2( maxx, miny ) );
  RightOut.push_back( Point_2( maxx + rx, miny ) );
  RightOut.push_back( Point_2( maxx + rx, maxy ) );
  RightOut.push_back( Point_2( maxx, maxy ) );
  CGAL::insert(arr, RightOut.edges_begin(), RightOut.edges_end());
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
/*
Point_2 getCenter( const ArrFaceCHandle& hCurrFace )
{
  ArrCCBHedgeCCirc circ = hCurrFace->outer_ccb();
  ArrCCBHedgeCCirc curr = circ;
  double x = 0.0;
  double y = 0.0;
  int n = 0;
  vector<Point_2> pts;
  do {
    const HEdge& he = *curr;
    Point_2 pt = he.source()->point();
    x += CGAL::to_double( pt[0] );
    y += CGAL::to_double( pt[1] );
    ++n;
  } while (++curr != circ);
  return Point_2( x/n, y/n);  
}
*/
//----------------------------------------------------------------------------
// Assumption: [p1, p2] is vertical, i.e. p1[0] == p2[0], 
// and all the forbidden edges are.
bool isForbidden( const list<Segment_2>& ForbiddenEdges, 
                  const Point_2&    p1, 
                  const Point_2&    p2 )
{
  list<Segment_2>::const_iterator iSeg = ForbiddenEdges.begin();
  for( ; iSeg != ForbiddenEdges.end(); ++iSeg )
  {
    const Segment_2& seg = *iSeg;
    const Point_2& s = seg.source();
    const Point_2& t = seg.target();
    double yu = max(CGAL::to_double(s[1]), CGAL::to_double(t[1]) );
    double yd = min(CGAL::to_double(s[1]), CGAL::to_double(t[1]) );
    if( p1[0] == s[0] )
    {
      // Vertical collision.
      double p1y = CGAL::to_double( p1[1] );
      double p2y = CGAL::to_double( p2[1] );
      bool b1 = yd < p1y && p1y < yu;
      bool b2 = yd < p2y && p2y < yu; 
      if( b1 || b2 )
      {
        cout << "=== Forbidden case: "<< endl;
        print_point( p1 );
        print_point( p2 );
        print_segment( *iSeg );
        cout << " ==== " << endl;
        return true;
      }
    }
  }
  return false;
}
//----------------------------------------------------------------------------
// Assumption: Faces are bounded.

bool twoFacesHaveCommonEdge( const ArrFaceCHandle& hFace1,
                             const ArrFaceCHandle& hFace2, 
                             double* px, double* py,
                             const list<Segment_2>& ForbiddenEdges )
{
  if( hFace1->is_unbounded() || hFace2->is_unbounded() )
      return false;

  ArrCCBHedgeCCirc circ1 = hFace1->outer_ccb();
  ArrCCBHedgeCCirc circ2 = hFace2->outer_ccb();
  ArrCCBHedgeCCirc curr1 = circ1;
  do 
  {
    ArrHedgeCHandle hHe1 = curr1;
    ArrCCBHedgeCCirc curr2 = circ2;
    do
    {
      ArrHedgeCHandle hHe2 = curr2;
      Point_2 p1s = hHe1->source()->point();
      Point_2 p2t = hHe2->target()->point();
      Point_2 p1t = hHe1->target()->point();
      Point_2 p2s = hHe2->source()->point();
 
      if( p1s == p2t && p1t == p2s && p1s[0] == p1t[0] )
      {
        if( isForbidden( ForbiddenEdges, p1s, p1t ) )
          continue;
        //It has to be vertical
        double x1 = CGAL::to_double( p1s[0] );
        double y1 = CGAL::to_double( p1s[1] );
        double x2 = CGAL::to_double( p1t[0] );
        double y2 = CGAL::to_double( p1t[1] );
        *px = (x1 + x2)/ 2.0;
        *py = (y1 + y2)/ 2.0;
        return true;
      }
    }while( ++curr2 != circ2 );
  }while( ++curr1 != circ1 );
  return false;
}
//----------------------------------------------------------------------------
void prepareDataInvokeDijkstra( const list<ArrFaceCHandle>& TrpzInPathFace,
                                const list<int>&            StartIdxs,
                                const list<int>&            EndIdxs,
                                const list<Segment_2>&      ForbiddenEdges,
                                vector<Point_2>&            ResPath )
{
  int nMtrxSize = TrpzInPathFace.size();
  list<ArrFaceCHandle>::const_iterator fi = TrpzInPathFace.begin();
  int i = 0, j = 0;
  char* pNeigMtrx = new char[nMtrxSize * nMtrxSize];
  Point_2* pNeigCntrs = new Point_2[nMtrxSize * nMtrxSize];
  for( ; fi != TrpzInPathFace.end(); ++fi, ++i )
  {
    list<ArrFaceCHandle>::const_iterator si = fi;
    for( ++si, j = i+1; si != TrpzInPathFace.end(); ++si, ++j )
    {
      double cx, cy;
      Point_2 CurrCntr;
      bool bHaveCommon = twoFacesHaveCommonEdge( *fi, *si, &cx, &cy, 
                                                 ForbiddenEdges );
      pNeigMtrx[i*nMtrxSize + j] = bHaveCommon;
      pNeigMtrx[j*nMtrxSize + i] = bHaveCommon;
       if( bHaveCommon )
      {
        Point_2 CurrCntr(cx, cy);
        cout << "Intersection point  ";
        print_point(CurrCntr);
        pNeigCntrs[i*nMtrxSize + j] = CurrCntr;
        pNeigCntrs[j*nMtrxSize + i] = CurrCntr;
      }
    }
  }
  list<int>::const_iterator iStartIdx = StartIdxs.begin();
  vector<int> RecordPath;
  for( ; iStartIdx != StartIdxs.end(); ++iStartIdx )
  {
    list<int>::const_iterator iEndIdx = EndIdxs.begin();
    for(; iEndIdx != EndIdxs.end(); ++iEndIdx )
    {
     cout << "Start Trpz Idx = " << *iStartIdx << endl;
      cout << "End Trpz Idx = " << *iEndIdx << endl;
      vector<int> CurrPath;
      int nCurrLen = 0;
      if( pNeigMtrx[ *iStartIdx * nMtrxSize + *iEndIdx ] )
      {
        CurrPath.push_back( *iStartIdx );
        CurrPath.push_back( *iEndIdx );
      }
      else if( *iStartIdx != *iEndIdx )
      {
        nCurrLen = dijkstra_compute_path( pNeigMtrx, nMtrxSize, 
                                          *iStartIdx, *iEndIdx, CurrPath );
        if( nCurrLen > 0 )
        {
          cout << "Path Found" << endl;
          vector<int>::const_iterator iPrevIdx = CurrPath.begin();
          vector<int>::const_iterator iCurrIdx = CurrPath.begin();
          for( ++iCurrIdx; iCurrIdx != CurrPath.end(); ++iCurrIdx, ++iPrevIdx )
          {
            Point_2 CurrPt = pNeigCntrs[*iPrevIdx * nMtrxSize + *iCurrIdx];
            print_point( CurrPt );
          }
        }

      }
      if( 0 == RecordPath.size() || nCurrLen < RecordPath.size() )
      {
        RecordPath.clear();
        RecordPath.insert( RecordPath.begin(), CurrPath.begin(), CurrPath.end() );
      }
    }
  }
  vector<int>::const_iterator iPrevIdx = RecordPath.begin();
  vector<int>::const_iterator iCurrIdx = RecordPath.begin();
  for( ++iCurrIdx; iCurrIdx != RecordPath.end(); ++iCurrIdx, ++iPrevIdx )
  {
    Point_2 CurrPt = pNeigCntrs[*iPrevIdx * nMtrxSize + *iCurrIdx];
    ResPath.push_back(CurrPt);
  }
  delete[] pNeigCntrs;
  delete[] pNeigMtrx;
}
//----------------------------------------------------------------------------
void applyGraphSearch( const Arrangement_2&        arrTrpzSplit,
                       const list<Segment_2>&      ForbiddenEdges,
                       const list<ArrFaceCHandle>& StartEndFaces,
                       int                         nPivotIdx,
                       vector<Point_2>&            ResPath )
{
  list<int>            StartIdxs;
  list<int>            EndIdxs;
  list<ArrFaceCHandle> AllTrpz;
  Arrangement_2::Face_const_iterator iTrpzFaces = arrTrpzSplit.faces_begin();
  for( ; iTrpzFaces != arrTrpzSplit.faces_end(); ++iTrpzFaces )
  {
    if ((*iTrpzFaces).is_unbounded())
      continue;
    ArrFaceCHandle f = iTrpzFaces;
    AllTrpz.push_back( f );
  }
  int idx = 0;
  int k = 0;
  list<ArrFaceCHandle>::const_iterator iSEF = StartEndFaces.begin();
  list<ArrFaceCHandle>::const_iterator iAllTrpz = AllTrpz.begin();
  for( ; iAllTrpz != AllTrpz.end(); ++iAllTrpz, ++idx )
  {
    if( k >= nPivotIdx )
        break;
    if( *iAllTrpz == *iSEF )
    {
      StartIdxs.push_back( idx );
      ++k;
      ++iSEF;
    }
  }
  iAllTrpz = AllTrpz.begin();
  for( idx = 0; iAllTrpz != AllTrpz.end(); ++iAllTrpz, ++idx )
  {
    if( iSEF ==  StartEndFaces.end() )
        break;
    if( *iAllTrpz == *iSEF )
    {
      EndIdxs.push_back( idx );
      ++k;
      ++iSEF;
    }
  }
  prepareDataInvokeDijkstra( AllTrpz, StartIdxs, EndIdxs, ForbiddenEdges, ResPath );
}
//----------------------------------------------------------------------------
vector<Point_2> findPath(const Point_2&      start,
                         const Point_2&      end,
                         const Polygon_2&    robot,
                         vector<Polygon_2>&  obstacles)
{
  if( 0 == obstacles.size() )
    return vector<Point_2>({start,end});

  cout << "Going from ";
  print_point(start);
  cout << " to ";
  print_point(end);
  Polygon_2 invRobot = getInversedRobot( robot, start );
  print_polygon( invRobot );
  cout << "^^^ Inv Robot ^^^" << endl; 
  vector<Polygon_with_holes_2> vecConfObst;
  for(const Polygon_2& obst : obstacles )
  {
    Polygon_with_holes_2 currConfObst = CGAL::minkowski_sum_2( obst,
                                                               invRobot );
    vecConfObst.push_back( currConfObst );
  }

  Polygon_set_2 FreeSpace;
  FreeSpace.join(vecConfObst.begin(), vecConfObst.end());
  FreeSpace.complement();
  Arrangement_2 arrConfSpace = FreeSpace.arrangement();

  Polygon_set_2::Traits_2 traits;
  /*
  CGAL::Arr_segment_traits_2<Kernel> traits;
  Arrangement_2 arrConfSpace(&traits);
  for( const Polygon_with_holes_2& pwh : vecConfObst)
    insertPolygonIntoArng( arrConfSpace, pwh );
  */
  print_arrangement( arrConfSpace );
  list<ArrFaceCHandle> arrFaces = getPathFaces( start, end, arrConfSpace );
  list<ArrFaceCHandle>::const_iterator iPathFace = getTwoInstanceFace(arrFaces);
  bool bPathExists = iPathFace != arrFaces.end();
  cout << "bPathExists = " << (bPathExists ? "true":"false") << endl;

  if( !bPathExists )
    return vector<Point_2>();

  ArrFaceCHandle hPathFace = *iPathFace;
  bool bPathFaceIsUnbounded = hPathFace->is_unbounded();
  if( bPathFaceIsUnbounded )
      addOuterRectBoundary( arrConfSpace, robot, start, end );

  arrFaces.clear();
  arrFaces = getPathFaces( start, end, arrConfSpace );
  iPathFace = getTwoInstanceFace(arrFaces);
  hPathFace = *iPathFace;

  print_arr_face(hPathFace);
  cout << "=== Path Face ends here == " << endl;

  print_arrangement( arrConfSpace );
  cout << "=== arrConfSpace end here === " << endl;

  list<Segment_2> ForbiddenEdges;
  getForbiddenEdges( hPathFace, ForbiddenEdges );


/*
  arrFaces.clear();
  arrFaces = getPathFaces( start, end, arrConfSpace );
  iPathFace = getTwoInstanceFace(arrFaces);
  hPathFace = *iPathFace;

  Arrangement_2 arrConfSpaceClean(&traits);
  list<Segment_2> ForbiddenEdges;
  insertFaceIntoArng( arrConfSpaceClean, hPathFace, ForbiddenEdges );


  print_arrangement( arrConfSpaceClean );
  cout << "=== arrConfSpaceClean end here === " << endl;
*/

  Kernel* kernel = &traits;
  Arrangement_2 arrTrpzSplit(arrConfSpace ); //Clean);
  vertical_decomposition(arrTrpzSplit, *kernel);

  print_arrangement( arrTrpzSplit );
  cout << "=== arrTrpzSplit end here === " << endl;



  int nPivotIdx = 0;
  list<ArrFaceCHandle> StartEndFaces = getPathFaces( start, end, 
                                                     arrTrpzSplit,
                                                     &nPivotIdx );

  cout << "=== Start Face === "<< endl;
  print_arr_face( *StartEndFaces.begin() );
  cout << "=== End Face === "<< endl;
  print_arr_face( *(++StartEndFaces.begin() ) );

  vector<Point_2> ResPath;
  applyGraphSearch( arrTrpzSplit, ForbiddenEdges, 
                    StartEndFaces, nPivotIdx, ResPath );
  ResPath.insert( ResPath.begin(), start );
  ResPath.push_back(end);
  return ResPath;
  //return vector<Point_2>({start,{1.71,5.57},{23.84,5.94},{21.21,29.17}, end});
}
//----------------------------------------------------------------------------

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
    cout << "Path validation:   " << ((result.verify(startPoint, endPoint, robot, obstacles)) ? "Success!" : "FAILURE")
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
