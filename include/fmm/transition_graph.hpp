/**
 * Fast map matching.
 *
 * Definition of a TransitionGraph, which is a wrapper of trajectory
 * candidates, raw trajectory and UBODT.
 * This class is designed for optimal path inference where
 * Viterbi algorithm is implemented.
 *
 * @author: Can Yang
 * @version: 2020.01.31
 * 
 * @revision : vuski
 * @version: 2021.11.29
 */

#ifndef FMM_TRANSITION_GRAPH_HPP
#define FMM_TRANSITION_GRAPH_HPP

#include <float.h>


#include "fmm/type.hpp"


using namespace FMM;


namespace FMM
{


    /**
 * Class related with map matching
 */

 /**
  * %Candidate edge matched to a GPS point
  */
struct Candidate
{
    //고유 처리 인덱스
    FMM::NodeIndex index; /**< The index is defined for a specific
                            candidate the index starting from N where N is the
                            numeber of vertices in the graph */
    double offsetLength; //offset distance from the start of polyline to p' */
    double offsetCost; //최단거리를 구할 때 사용하므로 cost offset으로 변경/**< offset distance from the start of polyline to p' */
    double dist; /**< distance from original point p to map matched point p' */
    Edge* edge;  /**< candidate edge */
    FMM::Point point; /**< boost point */ //선 위의 점 dist에 참고
};

typedef std::vector<Candidate> Point_Candidates; /**< Point candidates */
typedef std::vector<Point_Candidates> Traj_Candidates;
/**< trajectory  candidates */

typedef std::vector<const Candidate*> OptCandidatePath;
/**< Optimal candidates*/

typedef std::vector<FMM::EdgeID> O_Path; /**< Optimal path, edge id matched to
each point in the trajectory */

typedef std::vector<FMM::EdgeID> C_Path; /**< Complete path, ids of
a sequence of topologically connected edges.*/

/**
    * A candidate matched to a point
    */
struct MatchedCandidate {
    Candidate c; /**< Candidate matched to a point */
    double ep;  /**< emission probability */
    double tp;  /**< transition probability to previous matched candidate */
    double sp_dist; /**< shortest path distance to previous matched candidate */
};

/**
    * A vector of candidates, used for representing the matching result of
    * a trajectory.
    */
typedef std::vector<MatchedCandidate> MatchedCandidatePath;

/**
    * Map matched result representation
    */
struct MatchResult {
    int id; /**< id of the trajectory to be matched */
    MatchedCandidatePath opt_candidate_path; /**< A vector of candidate matched
    to each point of a trajectory. It is stored in order to export more
    detailed map matching information. */
    O_Path opath; /**< the optimal path,
                                containing id of edges matched to each
                                point in a trajectory */
    C_Path cpath; /**< the complete path, containing ids of a sequence of
                        topologically connected edges traversed by the
                        trajectory.  */
    std::vector<int> indices; /**< index of opath edge in cpath */
    FMM::LineString mgeom; /**< the geometry of the matched path */
};



/**
 * A node in the transition graph
 */
struct TGNode {
  const Candidate *c; /**< Candidate */
  TGNode *prev; /**< previous optimal candidate */
  double ep; /**< emission probability */
  double tp; /**< transition probability from previous optimal candidate */
  double cumu_prob; /**< current node's accumulative probability */
  double sp_dist; /**< shorest path distance from previous optimal
                       candidate to current node */
};

/**
 * A layer of nodes in the transition graph.
 */
typedef std::vector<TGNode> TGLayer;
/**
 * The optimal path of nodes in the transition graph
 */
typedef std::vector<const TGNode*> TGOpath;

/**
 * Transition graph class in HMM.
 *
 * The class stores the underlying transition graph of a HMM, which stores
 * the probabilities of candidates matched to a trajectories.
 */
class TransitionGraph
{
public:

    std::vector<bool> connectedToNext; //다음 단계의 점들과 이어졌는지 여부 update_tg의 결과로 나온다.
  /**
   * Transition graph constructor.
   *
   * A transition graph is created to store the probability of
   * emission and transition in HMM model where optimal path will be
   * infered.
   *
   * @param tc        Trajectory candidates
   * @param gps_error GPS error
   */
  TransitionGraph(const Traj_Candidates &tc, double gps_error) {
      //궤적 각각의 점 하나당 네트워크 후보군 집합
      for (const std::vector<Candidate>& cs : tc)
          //for (auto cs = tc.begin(); cs != tc.end(); ++cs) 
      {
          std::vector<TGNode> layer; //TGLayer는 vector<TGNode>

          //한점에 대한 후보군들 루프
          for (const Candidate& iter : cs) //auto iter = cs->begin(); iter != cs->end(); ++iter) 
          {
              double ep = calc_ep(iter.dist, gps_error);
              layer.push_back(TGNode{
                  &iter, // const Candidate* c; /**< Candidate */
                  nullptr,//    TGNode* prev; /**< previous optimal candidate */
                  ep,//    double ep; /**< emission probability */ 실제 떨어진 거리로 확률을 계산 한다. 
                  0,//    double tp; /**< transition probability from previous optimal candidate */
                -std::numeric_limits<double>::infinity(), //    double cumu_prob; /**< current node's accumulative probability */
                  0 //    double sp_dist; /**< shorest path distance from previous optimal
                  });
          }
          layers.push_back(layer);
      }
      if (!tc.empty()) {
          reset_layer(&(layers[0]));
      }
  }

  /**
   * Calculate transition probability
   * @param  sp_dist Shortest path distance between two candidates
   * @param  eu_dist Euclidean distance between two candidates
   * @return transition probability in HMM
   */
  static double calc_tp(double sp_dist,double eu_dist) {
      //두 수치가 근접한 정도로 전이확률을 구한다.
      //return std::min(sp_dist, eu_dist) / std::max(sp_dist, eu_dist);
      //이렇게 하면, 실제 시간과 매칭시키게 되므로, 얼핏 생각하면 맞는 것 갖지만,
      //소요시간이 길어질 때, 빙빙 도는 경로를 택한다. 즉, 막힐 때 도는 경로를 택한다.
      //아래 식으로 해야, 막혀서 분자가 커져서 1 이상이 나와도 그냥 그대로 간다.

      return std::min(1.0, eu_dist / sp_dist);

      //return eu_dist>=sp_dist ? //궤적길이가 매칭된 길이보다 큰가?
      //    1.0 :
      //    eu_dist/sp_dist; //궤적길이가 매칭된 길이보다 짧을 때
  }

  /**
   * Calculate emission probability
   * @param  dist  The actual gps error from observed point to matched point
   * @param  error The GPS sensor's accuracy
   * @return  emission probability in HMM
   */
  static double calc_ep(double dist,double error) {
      //double a = std::max(dist-500.0, 0.0) / error; //오차 500m까지는 무시하고 진행.
      double a = dist / error;
      return exp(-0.5 * a * a); //표준편차 공식에서 지수값만. 앞의 상수는 안쓴다
  }

  /**
   * Reset all the proability data stored in a layer of the transition graph
   * @param layer A layer in the transition graph
   */
  void reset_layer(TGLayer *layer) {
      for (auto iter = layer->begin(); iter != layer->end(); ++iter) {
          iter->cumu_prob = log(iter->ep);
          iter->prev = nullptr;
      }
  }

  /**
   * Find the optimal candidate in a layer with
   * the highest accumulative probability.
   * @param  layer [description]
   * @return  pointer to the optimal node in the transition graph
   */
  const TGNode *find_optimal_candidate(const TGLayer &layer) {
      const TGNode* opt_c = nullptr;
      double final_prob = -std::numeric_limits<double>::infinity();
      for (auto c = layer.begin(); c != layer.end(); ++c) {
          if (final_prob < c->cumu_prob) {
              final_prob = c->cumu_prob;
              opt_c = &(*c);
          }
      }
      return opt_c;
  }
  /**
   * Backtrack the transition graph to find an optimal path
   * @return An optimal path connecting the first layer with last layer and
   * has the highest accumulative probability value.
   */
  TGOpath backtrack(int trid) {
      SPDLOG_TRACE("Backtrack on transition graph");

      std::vector<int> endvec; //중간에 끊겼든, 모두 이어졌든 끊긴 지점을 모두 구한다.
      for (int i = connectedToNext.size()-1 ; i >=1 ; i--)
      {
          //지금 점은 다음이 끊어져 있고, 하나 전 점은 지금과 이어져 있으면.
          if (!connectedToNext[i] && connectedToNext[i-1]) 
          { 
              endvec.push_back(i);
              //if (trid == 20) SPDLOG_INFO("{} / {}", i, connectedToNext[i]);
          }          
      }
      if (endvec.empty()) endvec.push_back(0); //하나도 없으면 기본으로 하나
      //for (int tt : endvec)   SPDLOG_INFO("끊긴지점 : {}", tt);


      TGOpath opath;
      
      for (int endPoint : endvec) //끝지점마다 돈다.
      {


          TGNode* track_cand = nullptr;
          double final_prob = -std::numeric_limits<double>::infinity();

          //끊긴 지점의 노드를 참조한다.
          std::vector<TGNode>& last_layer = layers[endPoint];


          //마지막 노드에서 가장 확률이 높은 노드를 구한다. 
          //이미 cumulated 되어 있으므로, 그냥 고르면 된다.
          //즉 여기서 최종 결정
          for (auto c = last_layer.begin(); c != last_layer.end(); ++c) {
              if (final_prob < c->cumu_prob) {
                  final_prob = c->cumu_prob;
                  track_cand = &(*c);
              }
          }

          //이제 뒤로 되짚어가면서 각각의  gps 점들에 매칭된 edge를 구한다.

          int i = endPoint; // layers.size();
          //SPDLOG_INFO("PROB : {}", final_prob);
          if (final_prob > -std::numeric_limits<double>::infinity())
          {
              opath.push_back(track_cand);
              --i;
              SPDLOG_TRACE("Optimal candidate {} edge id {} sp {} tp {} cp {}",
                  i, track_cand->c->edge->id, track_cand->sp_dist, track_cand->tp,
                  track_cand->cumu_prob);
              // Iterate from tail to head to assign path

              while ((track_cand = track_cand->prev) != nullptr) {
                  opath.push_back(track_cand);
                  --i;
                  SPDLOG_TRACE("Optimal candidate {} edge id {} sp {} tp {} cp {}",
                      i, track_cand->c->edge->id, track_cand->sp_dist, track_cand->tp,
                      track_cand->cumu_prob);
              }
              //맨 뒤의 점에서 출발하여, 뒤에서 짚어나가면서 추가하고
              //마지막에 reverse해서 정 방향을 만들어 return한다.
             
          }
      }
      std::reverse(opath.begin(), opath.end());
      SPDLOG_TRACE("Backtrack on transition graph done");
      return opath;
  }
  /**
   * Get a reference to the inner layers of the transition graph.
   */
  std::vector<TGLayer> &get_layers() {
      return layers;
  }


  void print_optimal_info() {
      int N = layers.size();
      if (N < 1) return;
      const TGNode* global_opt_node = nullptr;
      for (int i = N - 1; i >= 0; --i) {
          const TGNode* local_opt_node = find_optimal_candidate(layers[i]);
          if (global_opt_node != nullptr) {
              global_opt_node = global_opt_node->prev;
          }
          else {
              global_opt_node = local_opt_node;
          }
          SPDLOG_TRACE("Point {} global opt {} local opt {}",
              i, (global_opt_node == nullptr) ? -1 : global_opt_node->c->edge->id,
              (local_opt_node == nullptr) ? -1 : local_opt_node->c->edge->id);
      }
  };
private:
  // candidates of a trajectory
  std::vector<TGLayer> layers; //전체 size는 gps 궤적 점의 개수만큼이고,
  
  //TGLayer각각에는 각 점들의 후보군이 담겨 있다.
};


}


/**
 * Utility functions for writing data to std stream
 */
namespace std {

    /**
     * Print a vector of values
     * @tparam T vector element type
     * @param os stream to write result
     * @param vec input vector
     * @return the stream with comma separated values of the vector written
     */
    template<typename T>
    std::ostream& operator<<(std::ostream& os,
        const std::vector<T>& vec) {
        if (!vec.empty()) {
            std::copy(vec.begin(), vec.end() - 1,
                std::ostream_iterator<T>(os, ","));
            os << vec.back();
        }
        return os;
    }

    /**
     * Write trajectory candidate to a stream
     * @param os stream to write
     * @param tr_cs trajectory candidate
     * @return the stream with trajectory candidate information written
     */
    std::ostream& operator<<(std::ostream& os,
        const FMM::Traj_Candidates& tr_cs) {
        os << "\nCandidate "
            << std::fixed << std::setw(4) << "step" << ";"
            << std::fixed << std::setw(6) << "index" << ";"
            << std::fixed << std::setw(8) << "offset" << ";"
            << std::fixed << std::setw(8) << "distance" << ";"
            << std::fixed << std::setw(8) << "edge_id" << '\n';
        for (auto tr_cs_iter = tr_cs.begin();
            tr_cs_iter != tr_cs.end(); ++tr_cs_iter) {
            for (auto p_cs_iter = tr_cs_iter->begin();
                p_cs_iter != tr_cs_iter->end();
                ++p_cs_iter) {
                os << "Candidate "
                    << std::fixed << std::setw(4) << std::distance(tr_cs.begin(),
                        tr_cs_iter) << ";"
                    << std::fixed << std::setw(6) << p_cs_iter->index << ";"
                    << std::fixed << std::setw(8) << p_cs_iter->offsetCost << ";"
                    << std::fixed << std::setw(8) << p_cs_iter->dist << ";"
                    << std::fixed << std::setw(8) << p_cs_iter->edge->id << '\n';
            }
        }
        return os;
    }


    /**
     * Write optimal candidate path into a stream
     * @param os stream to write
     * @param opath optimal candidate path
     * @return the stream with candidate path information written
     */
    std::ostream& operator<<(std::ostream& os,
        const FMM::OptCandidatePath& opath) {
        for (int i = 0; i < opath.size(); ++i) {
            // std::cout <<"Write edge "<< i <<" edge "<< opath[i]->edge->id <<"\n";
            os << opath[i]->edge->id;
            if (i != opath.size() - 1)
                os << ",";
        }
        return os;
    }

    /**
     * Write a point into a stream
     * @param os stream to write
     * @param geom point
     * @return the stream with wkt point written.
     */
    std::ostream& operator<<(std::ostream& os,
        const FMM::Point& geom) {
        os << std::setprecision(12) << "POINT(" << geom.x << " " << geom.y << ")";
    };

} // namespace std


#endif /* FMM_TRANSITION_GRAPH_HPP */