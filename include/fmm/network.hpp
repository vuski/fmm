/**
 * Fast map matching.
 *
 * Network class
 *
 * @author: Can Yang
 * @version: 2017.11.11
 * 
 * @revision : vuski
 * @version: 2021.11.29
 */

#define _CRT_SECURE_NO_WARNINGS

#ifndef FMM_NETWORK_HPP
#define FMM_NETWORK_HPP

#include <iostream>
#include <math.h> // Calulating probability
#include <iomanip>
#include <algorithm> // Partial sort copy
#include <unordered_set> // Partial sort copy
#include <unordered_map>
#include <map>
#include <string>
#include <stdexcept>

#include "fmm/type.hpp"
#include "fmm/util.hpp"
#include "fmm/geom_algorithm.hpp"
#include "fmm/transition_graph.hpp"
#include "fmm/gps_reader.hpp"

#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/stringbuffer.h"

#include "THST/RTree.h"

//using namespace FMM;

namespace FMM {

/**
 * Road network class
 */
struct NodeQ {
    unsigned int index;
    double value;
};

struct cmp {
    bool operator()(NodeQ& a, NodeQ& b) {
        return a.value > b.value;
    }
};


struct GraphCSR {
    int numNode;
    int numEdge;
    std::vector<unsigned int> rowOffset;
    std::vector<unsigned int> colIndex;
    std::vector<double> value;
};


 /**
*  Road edge property.
*/
struct EdgeProperty
{
    EdgeIndex index; /**< Index of the edge */
    double cost; /**< cost of the edge */
};



/**
    * Predecessor Map. It stores for each node, the previous node
    * visited, which is part of the shortest path routing result.
    */
typedef std::unordered_map<NodeIndex, NodeIndex> PredecessorMap;


/**
    * Successor Map. It stores for each node, the next node
    * visited, which is part of the bidirectional shortest path routing result.
    */
typedef std::unordered_map<NodeIndex, NodeIndex> SuccessorMap;

/**
    * Distance map. It stores for each node, the distance visited from a source
    * node, which is part of the shortest path routing result.
    */

typedef std::unordered_map<int, double> DistanceMap;

struct EdgeRtreeObj {
    spatial::Box2<double> bbox;
    Edge* eg;
};

struct Indexable{
    const double* min(const EdgeRtreeObj& value) const { return value.bbox.min; }
    const double* max(const EdgeRtreeObj& value) const { return value.bbox.max; }
};


  
class Network {
public:

    static size_t calcDuration_bg;
    static size_t calcDuration_thst;
    bool costDefault;
    static size_t calcCount;
    static size_t calcDuration;

  /**
   *  Constructor of Network
   *
   *  @param filename: the path to a network file in ESRI shapefile format
   *  @param id_name: the name of the id field
   *  @param source_name: the name of the source field
   *  @param target_name: the name of the target field
   *  @param mode: mode name, only applies to OSM network
   *
   */
  Network(const std::string &filename,
          const std::string &id_name = "id",
          const std::string &source_name = "source",
          const std::string &target_name = "target"
        ) {

      if (FMM::check_file_extension(filename, "shp"))
      {
          //read_ogr_file(filename,id_name,source_name,target_name);
      }
      else if (FMM::check_file_extension(filename, "geojson"))
      {
          read_geojson_file(filename, id_name, source_name, target_name);
      }
      else
      {
          std::string message = "Network file not supported ";
          SPDLOG_CRITICAL(message);
          throw std::runtime_error(message);
      }
  };

  Network(const std::string& filename,
      const std::string& id_name,
      const std::string& source_name,
      const std::string& target_name,
      const std::string& cost_name_,
      int srid_) {
      cost_name = cost_name_;
      SPDLOG_INFO("cost_name_ : {} -> {}", cost_name_, cost_name);
      srid = srid_;
      if (FMM::check_file_extension(filename, "geojson"))
      {
          read_geojson_file(filename, id_name, source_name, target_name);
      }

	
  };

  Network() {

  };

  void init(const std::string& filename,
	  const std::string& id_name,
	  const std::string& source_name,
	  const std::string& target_name,
	  const std::string& cost_name_,
	  int srid_)
  {
	  cost_name = cost_name_;
	  SPDLOG_INFO("cost_name_ : {} -> {}", cost_name_, cost_name);
	  srid = srid_;
	  if (FMM::check_file_extension(filename, "geojson"))
	  {
		  read_geojson_file(filename, id_name, source_name, target_name);
	  }
  }

  /**
   * Get number of nodes in the network
   * @return number of nodes
   */
  int get_node_count() const {
      return node_id_vec.size();
  }
  /**
   * Get number of edges in the network
   * @return number of edges
   */
  int get_edge_count() const {
      return edges.size();
  }

  const Edge& get_edgeByID(EdgeID id) const {
      return edges[get_edge_index(id)];
  };

  const Edge& get_edgeByIndex(EdgeIndex index) const {
      return edges[index];
  }

  /**
   * Get edges in the network
   * @return a constant reference to the edges
   */
  const std::vector<Edge> &get_edges() const {
      return edges;
  }
  /**
   * Get edge ID from index
   * @param index index of edge
   * @return edge ID
   */
  EdgeID get_edge_id(EdgeIndex index) const {
      return index < edges.size() ? edges[index].id : -1;
  }

  /**
   * Get edge index from ID
   * @param id edge id
   * @return edge index
   */
  inline EdgeIndex get_edge_index(EdgeID id) const {
      return edge_map.at(id);
  }


  inline EdgeIndex get_edge_index(NodeIndex source, NodeIndex target) const
  {
      auto& tg = nodeSTtoEgIndexMap.find(source)->second;
      return tg.find(target)->second;
  };
  /**
   * Get node ID from index
   * @param index index of node
   * @return node ID
   */
  NodeID get_node_id(NodeIndex index) const {
      return index < num_vertices ? node_id_vec[index] : -1;
  }
  /**
   * Get node index from id
   * @param id node id
   * @return node index
   */
  NodeIndex get_node_index(NodeID id) const {
      return node_map.at(id);
  }
  /**
   * Get node geometry from index
   * @param index node index
   * @return point of a node
   */
  FMM::Point get_node_geom_from_idx(NodeIndex index) const {
      return vertex_points[index];
  }


  /**
   *  Search for k nearest neighboring (KNN) candidates of a
   *  trajectory within a search radius
   *
   *  @param trajectory: input trajectory
   *  @param k: the number of candidates
   *  @param radius: the search radius
   *  @return a 2D vector of Candidates containing
   *  the candidates selected for each point in a trajectory
   *
   */
  FMM::Traj_Candidates search_tr_cs_knn(
    FMM::Trajectory &trajectory, std::size_t k, double radius) const {

      return search_tr_cs_knn(trajectory.geom, k, radius);
  }


  void searchRtreeRecursive(Point_Candidates& pcs, 
      const double px, const double py,
      double radius, const int& i) const
  {
      spatial::Box2<double> searchBox = {
              { px - radius, py - radius } ,
              { px + radius, py + radius } };
      std::vector<EdgeRtreeObj> results;
      rtree2.query(spatial::intersects<2>(searchBox.min, searchBox.max), std::back_inserter(results));

      int Nitems = results.size();
      for (unsigned int j = 0; j < Nitems; ++j) {
          // Check for detailed intersection
          // The two edges are all in OGR_linestring
          Edge* edge = results[j].eg;
          double offset;
          double dist;
          double closest_x, closest_y;
          double total_length;

          //현재 점과 선분과의 거리, 선분 시작점에서 선따라 얼마나 떨어진 점인지, 그 점은 무엇인지 구한다.
          linear_referencing(
              px, py, edge->geom,
              dist, offset, closest_x, closest_y,
              total_length); //계산을 위해 새로 추가함

          if (dist <= radius) { //bbox뿐만 아니라 거리가 반경 안에 들어와야 한다.
            // index, offset, dist, edge, pseudo id, point
              Candidate c = { 0, //인덱스는 일단 0으로 넣고 요 밑에서 업데이트
                             offset, //물리적 오프셋 입력
                             (edge->cost) * (offset / total_length), //여기에서 물리적 거리를 cost로 변경
                             dist,
                             edge,
                             Point(closest_x, closest_y) };
              pcs.push_back(c);
          }
      } //한 점에 대한 모든 후보군을 구한다.

      //중간에라도 하나도 없으면 뱉어낸다.
      //SPDLOG_INFO("Candidate count point :{}  {} (filter to k)", i, pcs.size());
      if (pcs.empty()) {
          //SPDLOG_WARN("Candidate not found for point {} {} {}", i, px, py);
          searchRtreeRecursive(pcs, px, py, radius * 1.5, i);;
      }
   
  }

  /**
   * Search for k nearest neighboring (KNN) candidates of a
   * linestring within a search radius
   *
   * @param geom
   * @param k number of candidates
   * @param radius search radius
   * @return a 2D vector of Candidates containing
   * the candidates selected for each point in a linestring
   */
  FMM::Traj_Candidates search_tr_cs_knn(const FMM::LineString &geom,
                                            std::size_t k,
                                            double radius)const {
      int NumberPoints = geom.get_num_points();
      Traj_Candidates tr_cs(NumberPoints);
      unsigned int current_candidate_index = num_vertices;
      //num_vertices는 원래 네트워크의 vertex 번호

      for (int i = 0; i < NumberPoints; ++i) { //입력한 궤적의 점 개수만큼 돈다.
        // SPDLOG_DEBUG("Search candidates for point index {}",i);
        // Construct a bounding boost_box
          double px = geom.get_x(i);
          double py = geom.get_y(i);
          Point_Candidates pcs;

          //auto t2 = std::chrono::high_resolution_clock::now();
          
          searchRtreeRecursive(pcs, px, py, radius, i); //못찾으면 찾을때까지. 새로운 거리를 리턴
          //auto t3 = std::chrono::high_resolution_clock::now();

          //calcDuration_thst += std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();

          //rtree2.Search(minxy, maxxy, MySearchCallback);

         




          // KNN part K-nearest neighbors
          // 후보군보다 개수가 많으면 소트해서 가장 가까운 것들 순으로 k개를 취한다.
          if (pcs.size() <= k) {
              tr_cs[i] = pcs;
          }
          else {
              tr_cs[i] = Point_Candidates(k);
              std::partial_sort_copy(
                  pcs.begin(), pcs.end(),
                  tr_cs[i].begin(), tr_cs[i].end(),
                  candidate_compare);
          }
          for (int m = 0; m < tr_cs[i].size(); ++m) {
              tr_cs[i][m].index = current_candidate_index + m; //인덱스는 현재 edge 총 숫자에 후보군 숫자를 더한다.
          }
          current_candidate_index += tr_cs[i].size();
          //SPDLOG_INFO("current_candidate_index {}\t{}", i,current_candidate_index);
      }
      return tr_cs;
  }

  /**
   * Get edge geometry
   * @param edge_id edge id
   * @return Geometry of edge
   */
  const FMM::LineString &get_edge_geom(EdgeID edge_id) const {
      return edges[get_edge_index(edge_id)].geom;
  }
  /**
   * Extract the geometry of a complete path, whose two end segment will be
   * clipped according to the input trajectory
   * @param traj input trajectory
   * @param complete_path complete path
   */
  FMM::LineString complete_path_to_geometry(
    const FMM::LineString &traj,
    const C_Path &complete_path) const {
      // if (complete_path->empty()) return nullptr;
      LineString line;
      if (complete_path.empty()) return line;
      int Npts = traj.get_num_points();
      int NCsegs = complete_path.size();
      if (NCsegs == 1) {
          double dist;
          double firstoffset;
          double lastoffset;
          const LineString& firstseg = get_edge_geom(complete_path[0]);
          linear_referencing_physically(traj.get_x(0), traj.get_y(0), firstseg,
              dist, firstoffset);
          linear_referencing_physically(traj.get_x(Npts - 1), traj.get_y(Npts - 1),
              firstseg, dist, lastoffset);
          LineString firstlineseg = cutoffseg_unique(firstseg, firstoffset,
              lastoffset);
          append_segs_to_line(&line, firstlineseg, 0);
      }
      else {
          const LineString& firstseg = get_edge_geom(complete_path[0]);
          const LineString& lastseg = get_edge_geom(complete_path[NCsegs - 1]);
          double dist;
          double firstoffset;
          double lastoffset;
          linear_referencing_physically(traj.get_x(0), traj.get_y(0), firstseg,
              dist, firstoffset);
          linear_referencing_physically(traj.get_x(Npts - 1), traj.get_y(Npts - 1),
              lastseg, dist, lastoffset);
          LineString firstlineseg = cutoffseg(firstseg, firstoffset, 0); //마지막 변수가 0이면 처음을 끊고
          LineString lastlineseg = cutoffseg(lastseg, lastoffset, 1); //마지막 변수가 1이면 뒤를 끊는다.
          append_segs_to_line(&line, firstlineseg, 0); //마지막 변수가 0이면 첫 점부터 끝점까지
          if (NCsegs > 2) {
              for (int i = 1; i < NCsegs - 1; ++i) {
                  const LineString& middleseg = get_edge_geom(complete_path[i]);
                  append_segs_to_line(&line, middleseg, 1); //마지막 변수가 1이면 둘째 점부터 끝점까지 (앞과 연속으로 이어지므로)
              }
          }
          append_segs_to_line(&line, lastlineseg, 1); //마지막 변수가 1이면 둘째 점부터 끝점까지 (앞과 연속으로 이어지므로)
      }
      return line;
  }

  /**
   * Get all node geometry
   * @return a vector of points
   */
  const std::vector<FMM::Point> &get_vertex_points() const {
      return vertex_points;
  }
  /**
   * Get node geometry
   * @param u node index
   * @return geometry of node
   */
  const FMM::Point &get_vertex_point(NodeIndex u) const {
      return vertex_points[u];
  }
  /**
   * Extract the geometry of a route in the network
   * @param path a route stored with edge ID
   * @return the geometry of the route
   */
  FMM::LineString route2geometry(const std::vector<EdgeID> &path) const {
      LineString line;
      if (path.empty()) return line;
      // if (complete_path->empty()) return nullptr;
      int NCsegs = path.size();
      for (int i = 0; i < NCsegs; ++i) {
          EdgeIndex e = get_edge_index(path[i]);
          const LineString& seg = edges[e].geom;
          if (i == 0) {
              append_segs_to_line(&line, seg, 0);
          }
          else {
              append_segs_to_line(&line, seg, 1);
          }
      }
      //SPDLOG_DEBUG("Path geometry is {}",line.exportToWkt());
      return line;
  }
  /**
   * Extract the geometry of a route in the network
   * @param path a route stored with edge Index
   * @return the geometry of the route
   */
  FMM::LineString route2geometry(
    const std::vector<EdgeIndex> &path) const {
      LineString line;
      if (path.empty()) return line;
      // if (complete_path->empty()) return nullptr;
      int NCsegs = path.size();
      for (int i = 0; i < NCsegs; ++i) {
          const LineString& seg = edges[path[i]].geom;
          if (i == 0) {
              append_segs_to_line(&line, seg, 0);
          }
          else {
              append_segs_to_line(&line, seg, 1);
          }
      }
      //SPDLOG_DEBUG("Path geometry is {}",line.exportToWkt());
      return line;
  }
  /**
   * Compare two candidate according to their GPS error
   * @param a candidate a
   * @param b candidate b
   * @return true if a.dist<b.dist
   */
  static bool candidate_compare(const Candidate &a, const Candidate &b) {
      if (a.dist != b.dist) {
          return a.dist < b.dist;
      }
      else {
          return a.edge->index < b.edge->index;
      }
  }

  unsigned int get_num_vertices() const {
      return num_vertices;
  }

  //const Network& get_network() {
  //    return this;
  //}

  void initNetworkGraph()      
  {
      const std::vector<Edge>& edges = get_edges();
      SPDLOG_INFO("Construct graph from network edges start");
      num_vertices = get_node_count();
      //CSR 형식으로 변환해둔다.
      convertNetworkToCSRGraph(edges, get_node_count());
      SPDLOG_INFO("CSR 그래프로 변환 완료");
  }


  /**
  * Dijkstra Shortest path query from source to target
  * @param source
  * @param target
  * @return a vector of edge index representing the path from source to target
  */

  size_t single_source_upperbound_dijkstra_new(
      NodeIndex s,
      double delta,
      PredecessorMap& predecessor,
      DistanceMap& distance) const {
      using std::priority_queue;
      priority_queue<NodeQ, std::vector<NodeQ>, cmp> Q;

      size_t searchCount = 0;

      NodeIndex scanBegin, scanEnd;
      unsigned int i;

      //Heap Q;
      // Initialization
      Q.push({ s, 0 });
      predecessor[s] = s;
      distance[s] = 0;

      double endTime = 0;
      // Dijkstra search
      while (!Q.empty()) {

          NodeQ node = Q.top();
          Q.pop(); //바로 빼버려도 상관 없음

          if (distance.find(node.index) != distance.end()) //이미 있으면
          {
              //현재 top으로 올라온 노드의 기록 값이, 중간에 기록된 값보다 크다면 지우고 건너뛴다.
              //같은 노드가 중간 계산으로 삭제되지 않고 중복 추가 되기 때문에 이런 일이 발생한다.
              //굳이 지우지 않아도 되지만, 중복 프로세스를 방지할 수 있다.
              if (distance[node.index] < node.value) continue;
          }

          NodeIndex currentNode = node.index;
          double currentValue = node.value;

          //if (currentValue > delta) break; //delta로 설정해놓은 옵션에 도달하면 멈춘다.
          scanBegin = graphStr.rowOffset[currentNode];
          scanEnd = graphStr.rowOffset[currentNode + 1];


          for (i = scanBegin; i < scanEnd; i++)
          {
              //searchCount++;
              NodeIndex nextNode = graphStr.colIndex[i];
              endTime = currentValue + graphStr.value[i];

              if (endTime > delta) continue; //제한보다 크면 건너뛴다.

              //auto iter = dmap->find(v);
              if (distance.find(nextNode) == distance.end())
              { //없으면 새로 넣는다.
                  predecessor[nextNode] = currentNode;  //
                  distance[nextNode] = endTime; //
                  Q.push({ nextNode, endTime });
              }
              else { //있으면
                  if (distance[nextNode] > endTime) //지금 계산한 것이 전에 도달했던 계산시간보다 적으면, 즉, 최소시간으로 업데이트
                  {
                      predecessor[nextNode] = currentNode;  //업데이트
                      distance[nextNode] = endTime; //업데이트
                      Q.push({ nextNode, endTime }); //겹쳐도 그냥 추가한다.
                  }
              }
          }
      }

      return searchCount;
  };



  /**
   * Calculate heuristic distance from p1 to p2,which is used in Astar routing.
   * @param p1
   * @param p2
   * @return the Euclidean distance from p1 to p2
   */
  double calc_heuristic_dist(
      const Point& p1, const Point& p2) const 
  {
      double hdist = sqrt((p2.x - p1.x) * (p2.x - p1.x) +
          (p2.y - p1.y) * (p2.y - p1.y));
      //SPDLOG_INFO("calc_heuristic_dist network.costDefault :{}", costDefault);
      return costDefault ? hdist : hdist / 5000;// 1666; 분속 meter/min 기준
      //time으로 cost를 둘 때는, 
      // 네트워크에 존재하는 가장 빠른 속도보다 같거나 빠르게 둔다.(분속 : 1분당 이동거리 m 기준)
      //현재 시속 300km = 분속 5000m 로 둠.
  }





  /**
 * AStar Shortest path query from source to target
 * @param source
 * @param target
 * @return a vector of edge index representing the path from source to target
 */
  std::vector<EdgeIndex> shortest_path_astar(NodeIndex source,
      NodeIndex target) const 
  {

      auto startT = std::chrono::high_resolution_clock::now();

      SPDLOG_TRACE("Shortest path astar starts");
      if (source == target) return {};

      bool found = false;

      using std::priority_queue;
      priority_queue<NodeQ, std::vector<NodeQ>, cmp> Q;

      NodeIndex scanBegin, scanEnd;
      unsigned int i;

      //const std::vector<Point>& vertex_points = get_vertex_points();

      PredecessorMap predecessor;
      DistanceMap distance;
      // Initialization
      double h = calc_heuristic_dist(vertex_points[source], vertex_points[target]);

      Q.push({ source, h });
      predecessor[source] = source;
      distance[source] = 0;

      double endTime = 0;

      // Dijkstra search
      while (!Q.empty()) {

          NodeQ node = Q.top();
          Q.pop();


          NodeIndex currentNode = node.index;
          double currentValue = node.value;
          if (currentNode == target) {
              found = true;
              break;
          }


          if (distance.find(currentNode) != distance.end()) //이미 있으면
          {
              //현재 top으로 올라온 노드의 기록 값이, 중간에 기록된 값보다 크다면 지우고 건너뛴다.
              //같은 노드가 중간 계산으로 삭제되지 않고 중복 추가 되기 때문에 이런 일이 발생한다.
              //굳이 지우지 않아도 되지만, 중복 프로세스를 방지할 수 있다.
              h = calc_heuristic_dist(vertex_points[currentNode], vertex_points[target]);
              if (distance[currentNode] + h < currentValue) continue;
          }

          scanBegin = graphStr.rowOffset[currentNode];
          scanEnd = graphStr.rowOffset[currentNode + 1];

          for (i = scanBegin; i < scanEnd; i++)
          {
              //searchCount++;
              NodeIndex nextNode = graphStr.colIndex[i];
              endTime = currentValue + graphStr.value[i];

              h = calc_heuristic_dist(vertex_points[nextNode], vertex_points[target]);

              if (distance.find(nextNode) == distance.end())
              { //없으면 새로 넣는다.
                  predecessor[nextNode] = currentNode;  //
                  distance[nextNode] = endTime; //
                  Q.push({ nextNode, endTime + h });
              }
              else { //있으면
                  if (distance[nextNode] > endTime) //지금 계산한 것이 전에 도달했던 계산시간보다 적으면, 즉, 최소시간으로 업데이트
                  {
                      predecessor[nextNode] = currentNode;  //업데이트
                      distance[nextNode] = endTime; //업데이트
                      Q.push({ nextNode, endTime + h }); //겹쳐도 그냥 추가한다.
                  }
              }
          }
      }

      auto endT = std::chrono::high_resolution_clock::now();
      //calcCount++;
       // calcDuration += std::chrono::duration_cast<std::chrono::microseconds>(endT - startT).count();
      // Backtrack from target to source
      if (found) return back_track(source, target, predecessor, distance);
      else return {};
  }



  bool shortest_path_astar_dist(
      NodeIndex source, NodeIndex target, double& dist) const {
      auto startT = std::chrono::high_resolution_clock::now();

      SPDLOG_TRACE("Shortest path astar starts");
      if (source == target) return false;

      bool found = false;

      using std::priority_queue;
      priority_queue<NodeQ, std::vector<NodeQ>, cmp> Q;

      NodeIndex scanBegin, scanEnd;
      unsigned int i;

      //const std::vector<Point>& vertex_points = get_vertex_points();

      //PredecessorMap predecessor;
      DistanceMap distance;
      // Initialization
      double h = calc_heuristic_dist(vertex_points[source], vertex_points[target]);

      Q.push({ source, h });
      distance[source] = 0;

      double endTime = 0;

      // Dijkstra search
      while (!Q.empty()) {

          NodeQ node = Q.top();
          Q.pop();

          NodeIndex currentNode = node.index;
          double currentValue = node.value;
          if (currentNode == target) {
              found = true;
              break;
          }

          if (distance.find(currentNode) != distance.end()) //이미 있으면
          {
              //현재 top으로 올라온 노드의 기록 값이, 중간에 기록된 값보다 크다면 지우고 건너뛴다.
              //같은 노드가 중간 계산으로 삭제되지 않고 중복 추가 되기 때문에 이런 일이 발생한다.
              //굳이 지우지 않아도 되지만, 중복 프로세스를 방지할 수 있다.
              h = calc_heuristic_dist(vertex_points[currentNode], vertex_points[target]);
              if (distance[currentNode] + h < currentValue) continue;
          }



          scanBegin = graphStr.rowOffset[currentNode];
          scanEnd = graphStr.rowOffset[currentNode + 1];

          for (i = scanBegin; i < scanEnd; i++)
          {
              //searchCount++;
              NodeIndex nextNode = graphStr.colIndex[i];
              endTime = currentValue + graphStr.value[i];

              h = calc_heuristic_dist(vertex_points[nextNode], vertex_points[target]);

              if (distance.find(nextNode) == distance.end())
              { //없으면 새로 넣는다.
                 // predecessor[nextNode] = currentNode;  //
                  distance[nextNode] = endTime; //
                  Q.push({ nextNode, endTime + h });
              }
              else { //있으면
                  if (distance[nextNode] > endTime) //지금 계산한 것이 전에 도달했던 계산시간보다 적으면, 즉, 최소시간으로 업데이트
                  {
                      //predecessor[nextNode] = currentNode;  //업데이트
                      distance[nextNode] = endTime; //업데이트
                      Q.push({ nextNode, endTime + h }); //겹쳐도 그냥 추가한다.
                  }
              }
          }
      }

      auto endT = std::chrono::high_resolution_clock::now();
      calcCount++;
      calcDuration += std::chrono::duration_cast<std::chrono::microseconds>(endT - startT).count();

      if (found) {
          //std::cout << "성공. 거리 : " << dist << "\t소요시간 : " << clock() - timeT << std::endl;
          dist = distance[target];
          return true;
      }
      else {
          //std::cout << "실패" << std::endl;
          dist = std::numeric_limits<double>::infinity();
          return false;
          
      }

  }



  /**
   * Backtrack the routing result to find a path from source to target
   * @param source
   * @param target
   * @param pmap predecessor map
   * @param dmap distance map
   * @return a vector of edge index representing the path from source to target
   */
  std::vector<EdgeIndex> back_track(NodeIndex source,
      NodeIndex target,
      PredecessorMap& pmap,
      DistanceMap& dmap) const {
      // SPDLOG_TRACE("Backtrack starts");
      if (dmap.find(target) == dmap.end()) {
          return {};
      }
      else {
          std::vector<EdgeIndex> path;
          NodeIndex v = target;
          NodeIndex u = pmap[v];
          while (v != source) {
              //double cost = dmap.at(v) - dmap.at(u);
              // SPDLOG_TRACE("u {} d {} v {} d {} cost {}",
              //              get_node_id(u), dmap.at(u),
              //              get_node_id(v), dmap.at(v), cost);
              path.push_back(get_edge_index(u, v));
              v = u;
              u = pmap[v];
          }
          std::reverse(path.begin(), path.end());
          return path;
      }
  }


private:


  void read_geojson_file(const std::string& filename,
        const std::string& id_name,
        const std::string& source_name,
        const std::string& target_name) {

      using std::cout;
      using std::endl;

      costDefault = cost_name.compare("") == 0 ? true : false;
      SPDLOG_INFO("cost는 기본으로? : {}, {}", costDefault, cost_name);
      //srid = 5179;
      rapidjson::Document doc;

      {
          using namespace rapidjson;
          char* buffer;
          {
              clock_t t2 = clock();
              int count;
              FILE* fp = NULL;
              fopen_s(&fp, filename.c_str(), "r");

              //read size of file
              fseek(fp, 0, SEEK_END);
              int64_t size = ftell(fp);

              buffer = (char*)malloc(size + 1);
              memset(buffer, 0, size + 1);

              fseek(fp, 0, SEEK_SET);
              count = fread(buffer, size, 1, fp);
              //처리	
              fclose(fp);
              clock_t t3 = clock();
              std::cout << " 읽기시간 : " << t3 - t2 << std::endl;
          }

          std::cout << "문자열 길이:" << strlen(buffer) << std::endl;
          doc.Parse(buffer);

          std::cout << "parsed" << std::endl;
          free(buffer);
      }


      using namespace rapidjson;
      const Value& features = doc["features"];
      cout << "geojson 멤버 수:" << features.Size() << endl;
      unsigned int featureSize = features.Size();

      const char* idStr = id_name.c_str();
      const char* sourceStr = source_name.c_str();
      const char* targetStr = target_name.c_str();
      const char* costStr = cost_name.c_str();
      SPDLOG_INFO("필드 이름: {}, {}, {}, {}", idStr, sourceStr, targetStr, costStr);

      std::unordered_map<size_t, Edge> keymap;

      for (unsigned int i = 0; i < featureSize; i++)//도형 하나하나 루프
      {
          const Value& feature = features[i];
          const Value& properties = feature["properties"];
          std::string linetype = feature["geometry"]["type"].GetString();

          const Value& linestring = linetype == "MultiLineString" ?
              feature["geometry"]["coordinates"].GetArray()[0] :
              feature["geometry"]["coordinates"];//폴리곤들의 배열    


          //좌표값을 읽는다.
          LineString geom;
          unsigned int coordSize = linestring.Size();

          if (coordSize > 1) { //점이 하나로 이루어져 있을 때는 버린다.                
              for (unsigned int k = 0; k < coordSize; k++) {
                  //cout << linestring[k].GetArray()[0].GetDouble() << "\t"
                  //    << linestring[k].GetArray()[1].GetDouble() << endl;
                  geom.add_point(linestring[k].GetArray()[0].GetDouble(),
                      linestring[k].GetArray()[1].GetDouble());
              } //k                 
          }

          //속성을 읽는다.
          EdgeID id = properties[idStr].GetInt();
          NodeID source = properties[sourceStr].GetInt();
          NodeID target = properties[targetStr].GetInt();
          double cost = (costDefault ?
              geom.get_length() :
              properties[costStr].GetDouble()
              );

          //double timeDist = properties["timeFT"].GetDouble();

          //고유번호 따기. 기존에 존재하면 받아온다.
          NodeIndex s_idx, t_idx;
          if (node_map.find(source) == node_map.end()) {
              s_idx = node_id_vec.size();
              node_id_vec.push_back(source);
              node_map.insert({ source, s_idx });
              vertex_points.push_back(geom.get_point(0));
          }
          else {
              s_idx = node_map[source];
          }

          if (node_map.find(target) == node_map.end()) {
              t_idx = node_id_vec.size();
              node_id_vec.push_back(target);
              node_map.insert({ target, t_idx });
              int npoints = geom.get_num_points();
              vertex_points.push_back(geom.get_point(npoints - 1));
          }
          else {
              t_idx = node_map[target];
          }


          Edge eg = {
              0, //edge의 고유 키값. 일단 0으로
              id,
              s_idx,
              t_idx,
              cost,
              geom
          };

          //현재 넣은 네트워크 기준으로 키값 만듬
          size_t nodeKey = s_idx * 10000000ULL + t_idx;
          if (keymap.find(nodeKey) == keymap.end())
          { //없으면
              keymap[nodeKey] = eg;
          }
          else {
              SPDLOG_INFO("source-target node 중복 네트워크 \"source\": {}, \"target\" :{},",
                  source, target);
              if (keymap[nodeKey].cost > eg.cost) { //이러면 안되지만, 만약 키가 존재할 경우, cost가 작은 경우에만 업데이트 해준다.
                  keymap[nodeKey] = eg;
              }
          }

          //struct Edge
          //{
          //    EdgeIndex index; /**< Index of an edge, which is continuous [0,N-1] */
          //    EdgeID id; /**< Edge ID, can be discontinuous integers */
          //    NodeIndex source; /**< source node index */
          //    NodeIndex target; /**< target node index */
          //    double length; /**< length of the edge polyline */
          //    FMM::CORE::LineString geom; /**< the edge geometry */
          //};





      } //i 도형 하나하나

      //중복을 체크했으므로 이제 집어넣는다.
      EdgeIndex index = 0;
      for (auto eg : keymap)
      {
          Edge e = eg.second;
          e.index = index;

          edges.push_back(e);
          //SPDLOG_INFO("index : {} \t id : {} \t s_idx : {} \t t_idx {} \t length : {} \t",
          //    index, id, s_idx, t_idx, geom.get_length());
          edge_map.insert({ e.id, index });

          //소스와 타겟으로 엣지 인덱스 구하기 위해 넣는다.
          nodeSTtoEgIndexMap[e.source][e.target] = index;
          index++;
      }


      num_vertices = node_id_vec.size();
      SPDLOG_INFO("Number of edges {} nodes {}", edges.size(), num_vertices);
      build_rtree_index();
      SPDLOG_INFO("Read network done.");
  }    // Network constructor


  /**
   * Concatenate a linestring segs to a linestring line, used in the
   * function complete_path_to_geometry
   *
   * @param line: linestring which will be updated
   * @param segs: segs that will be appended to line
   * @param offset: the number of points skipped in segs.
   */
  static void append_segs_to_line(FMM::LineString *line,
                                  const FMM::LineString &segs,
                                  int offset = 0) {
      int Npoints = segs.get_num_points();
      for (int i = 0; i < Npoints; ++i) {
          if (i >= offset) {
              line->add_point(segs.get_x(i), segs.get_y(i));
          }
      }
  }

  /**
   * Build rtree for the network
   */
  
  void build_rtree_index() {

      // Build an rtree for candidate search
      SPDLOG_DEBUG("Create THST rtree");
      // create some Items
      for (std::size_t i = 0; i < edges.size(); ++i) {

          Edge* edge = &edges[i];
          double x1, y1, x2, y2;
          boundingbox_geometry(edge->geom, x1, y1, x2, y2);

          spatial::Box2<double> bbox = { {x1,y1},{x2,y2} };
          EdgeRtreeObj ero = { bbox, edge };

          rtree2.insert(ero);

      }
      SPDLOG_DEBUG("Create boost rtree done");
  }

  void convertNetworkToCSRGraph(const std::vector<Edge>& edges,
      unsigned int nodeSize) {
      //초기화
      graphStr.numNode = nodeSize;
      graphStr.numEdge = edges.size();
      graphStr.rowOffset.resize(graphStr.numNode + 1);
      graphStr.colIndex.resize(graphStr.numEdge);
      graphStr.value.resize(graphStr.numEdge);


      using std::map;
      std::unordered_set<size_t> nodeSet;
      size_t maxNode = 0;

      //결과 map을 graphStr에 저장하기
      map<NodeIndex, map<NodeIndex, double>> data; //오름차순으로 기본 정렬 됨
      for (const Edge& edge : edges) {
          data[edge.source][edge.target] = edge.cost;


          size_t nodeKey = edge.source * 1000000ULL + edge.target;
          if (nodeSet.find(nodeKey) != nodeSet.end())
          {
              SPDLOG_INFO("중복. 이러면 안됨 \"source\": {}, \"target\" :{},",
                  edge.source, edge.target);
          }
          else {
              nodeSet.insert(nodeKey);
          }
          if (edge.source > maxNode) maxNode = edge.source;
          if (edge.target > maxNode) maxNode = edge.target;
      }
      SPDLOG_INFO("source-target combination : {}, maxNode : {}",
          nodeSet.size(),
          maxNode);


      map<NodeIndex, map<NodeIndex, double>>::iterator iter1;
      map<NodeIndex, double>::iterator iter2;
      unsigned int index = 0; //전체 데이터 번호. 끝까지 가면 edge+1이 된다.
      unsigned int rowIndex = 0; //실제 진행에 따른 row 번호
      for (iter1 = data.begin(); iter1 != data.end(); iter1++) {

          while (rowIndex <= (iter1->first)) { //비어있는 row를 고려한다.
              graphStr.rowOffset[rowIndex] = index;
              rowIndex++;
          }
          map<NodeIndex, double>& row = iter1->second;
          for (iter2 = row.begin(); iter2 != row.end(); iter2++) {
              graphStr.colIndex[index] = iter2->first;
              graphStr.value[index] = iter2->second;
              index++;
          }
      }
      graphStr.rowOffset[rowIndex] = index; //마지막 추가

      SPDLOG_INFO("csr graph nodes : {}\t edges : {}\t {}",
          graphStr.numNode, graphStr.value.size(), index);
  }


  int srid;   // Spatial reference id
  spatial::RTree<double, EdgeRtreeObj, 2, 64, 32, Indexable> rtree2;

  std::vector<Edge> edges;   // all edges in the network
  std::unordered_map<NodeIndex, std::unordered_map<NodeIndex, EdgeIndex>> nodeSTtoEgIndexMap;

  NodeIDVec node_id_vec;
  unsigned int num_vertices;
  NodeIndexMap node_map;
  EdgeIndexMap edge_map;
  std::vector<FMM::Point> vertex_points;
  std::string cost_name = ""; //기본값으로 ""이면 길이를 계산, 아니면 cost를 별도설정한다.

  //그래프 속성
  GraphCSR graphStr;
  static constexpr double DOUBLE_MIN = 1.e-6;

}; // Network

} // FMM
#endif /* FMM_NETWORK_HPP */