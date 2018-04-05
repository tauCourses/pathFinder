#include <boost/config.hpp>
 
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/properties.hpp>
 
#include <boost/property_map/property_map.hpp>
 
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>
 
int dijkstra_compute_path( char* pNeigMtrx, int nMtrxSize, 
                           int nSrcIdx, int nDstIdx, std::vector<int>& ResPath )
{
  typedef float Weight;
  typedef boost::property<boost::edge_weight_t, Weight> WeightProperty;
  typedef boost::property<boost::vertex_name_t, std::string> NameProperty;
 
  typedef boost::adjacency_list < boost::listS, boost::vecS, boost::directedS,
    NameProperty, WeightProperty > Graph;
 
  typedef boost::graph_traits < Graph >::vertex_descriptor Vertex;
 
  typedef boost::property_map < Graph, boost::vertex_index_t >::type IndexMap;
  typedef boost::property_map < Graph, boost::vertex_name_t >::type NameMap;
 
  typedef boost::iterator_property_map < Vertex*, IndexMap, Vertex, Vertex& > PredecessorMap;
  typedef boost::iterator_property_map < Weight*, IndexMap, Weight, Weight& > DistanceMap;
 
 
  // Create a graph
  Graph g;
  std::vector<Vertex> AllVrt; 
  int i,j;
  for( i = 0; i < nMtrxSize; ++i )
  {
    char buf[4] = {'\0', '\0', '\0', '\0' };
    std::ostringstream str_buf(buf);
    str_buf << i; 
    Vertex v = boost::add_vertex(str_buf.str(), g);
    AllVrt.push_back(v);
  }
 
  Weight weight1 = 1;

  for( i = 0; i < nMtrxSize; ++i )
  {
    for( j = 0; j < nMtrxSize; ++j )
    {
      if( pNeigMtrx[i*nMtrxSize+j] )
        boost::add_edge(AllVrt[i], AllVrt[j], weight1, g);
    }
  }
  
  // Create things for Dijkstra
  std::vector<Vertex> predecessors(boost::num_vertices(g)); // To store parents
  std::vector<Weight> distances(boost::num_vertices(g)); // To store distances
 
  IndexMap indexMap = boost::get(boost::vertex_index, g);
  PredecessorMap predecessorMap(&predecessors[0], indexMap);
  DistanceMap distanceMap(&distances[0], indexMap);
 
  // Compute shortest paths from v0 to all vertices, and store the output in predecessors and distances
  // boost::dijkstra_shortest_paths(g, v0, boost::predecessor_map(predecessorMap).distance_map(distanceMap));
  // This is exactly the same as the above line - it is the idea of "named parameters" - you can pass the
  // prdecessor map and the distance map in any order.
  boost::dijkstra_shortest_paths(g, AllVrt[nSrcIdx], 
                                 boost::distance_map(distanceMap).predecessor_map(predecessorMap));
 
  // Output results
  NameMap nameMap = boost::get(boost::vertex_name, g);
 
  typedef std::vector<Graph::edge_descriptor> PathType;
 
  PathType path;
 
  Vertex v = AllVrt[nDstIdx]; // We want to start at the destination and work our way back to the source
  for(Vertex u = predecessorMap[v]; // Start by setting 'u' to the destintaion node's predecessor
      u != v; // Keep tracking the path until we get to the source
      v = u, u = predecessorMap[v]) // Set the current vertex to the current predecessor, and the predecessor to one level up
  {
    std::pair<Graph::edge_descriptor, bool> edgePair = boost::edge(u, v, g);
    Graph::edge_descriptor edge = edgePair.first;
 
    path.push_back( edge );
  }
 
  std::istringstream strSrc(nameMap[boost::source(*(path.rbegin()), g)]);
  int nStartIdx;
  strSrc >> nStartIdx;
  ResPath.push_back(nStartIdx);
  for(PathType::reverse_iterator pathIterator = path.rbegin(); pathIterator != path.rend(); ++pathIterator)
  {
    std::istringstream strTrg(nameMap[boost::target(*pathIterator, g)]);
    int nCurTrgIdx;
    strTrg >> nCurTrgIdx;
    ResPath.push_back( nCurTrgIdx );
  }
  return ResPath.size();
}
 
