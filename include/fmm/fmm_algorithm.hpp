/**
 * Fast map matching.
 *
 * fmm algorithm implementation
 *
 * @author: Can Yang
 * @version: 2020.01.31
 * 
 * @revision : vuski
 * @version: 2021.11.29
 */

#ifndef FMM_FMM_ALGORITHM_H_
#define FMM_FMM_ALGORITHM_H_

#include "fmm/util.hpp"
#include "fmm/type.hpp"
#include "fmm/geom_algorithm.hpp"
#include "fmm/transition_graph.hpp"

#include "fmm/network.hpp"
#include "fmm/ubodt.hpp"
#include "fmm/gps_reader.hpp"
#include "fmm/mm_writer.hpp"

#include <string>


#include "cxxopts/cxxopts.hpp"

using namespace FMM;


namespace FMM {

    /**
     * Configuration class for fmm algorithm
     */
    struct FastMapMatchConfig {
        /**
         * Constructor of FastMapMatch configuration
         * @param k_arg the number of candidates
         * @param r_arg the search radius, in map unit, which is the same as
         * GPS data and network data.
         * @param gps_error the gps error, in map unit
         *
         */
        FastMapMatchConfig(int k_arg, double r_arg,
            double gps_error,
            double reverse_tolerance) :
            k(k_arg), radius(r_arg), gps_error(gps_error),
            reverse_tolerance(reverse_tolerance) {
        };

        int k; /**< Number of candidates */
        double radius; /**< Search radius*/
        double gps_error; /**< GPS error */
        double reverse_tolerance;
        /**
         * Check if the configuration is valid or not
         * @return true if valid
         */
         //reverse tolerance는 0~1사이가 되어야 하는 것 같음
        bool validate() const {
            if (gps_error <= 0 || radius <= 0 || k <= 0 || reverse_tolerance < 0
                || reverse_tolerance>1) {
                SPDLOG_CRITICAL(
                    "Invalid mm parameter k {} r {} gps error {} reverse_tolerance {}",
                    k, radius, gps_error, reverse_tolerance);
                return false;
            }
            return true;
        }
        /**
         * Print information about this configuration
         */
        void print() const {
            SPDLOG_INFO("FMMAlgorithmConfig");
            SPDLOG_INFO("k {} radius {} gps_error {} reverse_tolerance {}",
                k, radius, gps_error, reverse_tolerance);
        };

        /**
         * Load configuration from argument data
         * @param arg_data argument data
         * @return a %FastMapMatchConfig object
         */
        static FastMapMatchConfig load_from_arg(
            const cxxopts::ParseResult& arg_data) {
            int k = arg_data["candidates"].as<int>();
            double radius = arg_data["radius"].as<double>();
            double gps_error = arg_data["error"].as<double>();
            double reverse_tolerance = arg_data["reverse_tolerance"].as<double>();
            return FastMapMatchConfig{ k, radius, gps_error, reverse_tolerance };
        };
        /**
         * Register arguments to an option object
         */
        static void register_arg(cxxopts::Options& options) {
            options.add_options()
                ("k,candidates", "Number of candidates",
                    cxxopts::value<int>()->default_value("8"))
                ("r,radius", "Search radius",
                    cxxopts::value<double>()->default_value("300.0"))
                ("reverse_tolerance", "Ratio of reverse movement allowed",
                    cxxopts::value<double>()->default_value("0.0"))
                ("e,error", "GPS error",
                    cxxopts::value<double>()->default_value("50.0"));
        }

        /**
         * Register help information to a string stream
         */
        static void register_help(std::ostringstream& oss) {
            oss << "-k/--candidates (optional) <int>: Number of candidates (8)\n";
            oss << "-r/--radius (optional) <double>: search "
                "radius (network data unit) (300)\n";
            oss << "-e/--error (optional) <double>: GPS error "
                "(network data unit) (50)\n";
            oss << "--reverse_tolerance (optional) <double>: proportion "
                "of reverse movement allowed on an edge\n";
        };

    };

    /**
     * Fast map matching algorithm/model.
     *
     *
     */
    class FastMapMatch {
    public:

        //소요시간 측정용
        static size_t loopCount;
        static size_t calcDuration;
        static size_t updateTGDuration;
        static size_t ubodtCount;
        /**
         * Constructor of Fast map matching model
         * @param network road network
         * @param graph road network graph
         * @param ubodt Upperbounded origin destination table
         */
        FastMapMatch(const Network& network,
            std::shared_ptr<UBODT> ubodt)
            : network_(network), ubodt_(ubodt) {
        };


        FastMapMatch()
        {

        };


        void init(Network network,
            std::shared_ptr<UBODT> ubodt)
        {
            network_ = network;
            ubodt_ = ubodt;
        }

        Network getNetwork()
        {
            return network_;
        }

        std::shared_ptr<UBODT> getUBODT()
        {
            return ubodt_;
        }
        /**
         * Match a trajectory to the road network
         * @param  traj   input trajector data
         * @param  config configuration of map matching algorithm
         * @return map matching result
         */
        MatchResult match_traj(const Trajectory& traj,
            const FastMapMatchConfig& config) {
            //SPDLOG_INFO("reverse_tolerance : {}", config.reverse_tolerance);
            SPDLOG_DEBUG("Count of points in trajectory {}", traj.geom.get_num_points());
            SPDLOG_DEBUG("Search candidates");

            //config에서 설정한 반경으로 찾는데, 현재 cost가 거리든 시간이든 상관 없다. 그냥 거리로 계산한다.
            auto startTime = std::chrono::high_resolution_clock::now();
            Traj_Candidates tc = network_.search_tr_cs_knn(
                traj.geom, config.k, config.radius);
            //찾은 선분들의 총 숫자는 전체  edge 숫자에 각 궤적의 점들마다 radius 안의 가까운 점들이 얼마나 있는지 그 숫자를 모두 더한 수
            auto endTime = std::chrono::high_resolution_clock::now();
            calcDuration += std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();

            SPDLOG_DEBUG("Trajectory candidate {}", tc);
            if (tc.empty()) return MatchResult{};
            SPDLOG_DEBUG("Generate transition graph");

            //gps error를 고려한, 앞의 찾기 결과의 변형. \
            //뒤의 계산을 위해 새로 TGLayer 구조체를 만들고 값들을 비워둔다.
            TransitionGraph tg(tc, config.gps_error);

            SPDLOG_DEBUG("Update cost in transition graph");
            // The network will be used internally to update transition graph
            update_tg(&tg, traj, config.reverse_tolerance);

            SPDLOG_DEBUG("Optimal path inference");
            //결정은 마지막 점에 할당된 확률들을 비교하면서 함수 시작하자마자 한번에 찾는다.
            //그 다음에는 역방향 추적해가면서 하나씩
            TGOpath tg_opath = tg.backtrack(traj.id);
            //o path는 gps 점의 개수만큼만 존재한다.


            //여기부터는 결과를 정리하는 부분
            SPDLOG_DEBUG("Optimal path size {}", tg_opath.size());
            MatchedCandidatePath matched_candidate_path(tg_opath.size());
            std::transform(tg_opath.begin(), tg_opath.end(),
                matched_candidate_path.begin(),
                [](const TGNode* a) {
                    return MatchedCandidate{
                      *(a->c), a->ep, a->tp, a->sp_dist
                    };
                });
            O_Path opath(tg_opath.size());
            std::transform(tg_opath.begin(), tg_opath.end(),
                opath.begin(),
                [](const TGNode* a) {
                    return a->c->edge->id;
                });
            std::vector<int> indices;
            const std::vector<FMM::Edge>& edges = network_.get_edges();
            C_Path cpath = ubodt_->construct_complete_path(traj.id, tg_opath, edges,
                &indices,
                config.reverse_tolerance
            );
            if (cpath.empty())  //만약 없으면, ubodt에 존재하지 않아서 그러함.
            {
                cpath = construct_complete_path_astar(traj.id, tg_opath, edges,
                    &indices,
                    config.reverse_tolerance);
            }


            //SPDLOG_DEBUG("Opath is {}", opath);
            //SPDLOG_DEBUG("Indices is {}", indices);
            //SPDLOG_DEBUG("Complete path is {}", cpath);
            LineString mgeom = network_.complete_path_to_geometry(
                traj.geom, cpath);
            return MatchResult{
              traj.id, matched_candidate_path, opath, cpath, indices, mgeom };
        }

        /**
         * Match a wkt linestring to the road network.
         * @param wkt WKT representation of a trajectory
         * @param config Map matching configuration
         * @return Map matching result in POD format used in Python API
         */
         //PYTHON::PyMatchResult match_wkt(
         //    const std::string &wkt,const FastMapMatchConfig &config);


         ///**
         // * Match GPS data stored in a file
         // * @param  gps_config    [description]
         // * @param  result_config [description]
         // * @param  config        [description]
         // * @return a string storing information about running time and
         // * statistics.
         // */
         //std::string match_gps_file(
         //  const FMM::CONFIG::GPSConfig &gps_config,
         //  const FMM::CONFIG::ResultConfig &result_config,
         //  const FastMapMatchConfig &config,
         //  bool use_omp = true
         //);

    protected:
        /**
         * Get shortest path distance between two candidates
         * @param  ca from candidate
         * @param  cb to candidate
         * @return  shortest path value
         */
        double get_sp_dist(const Candidate* ca,
            const Candidate* cb,
            double reverse_tolerance) {
            double sp_dist = 0;
            if (ca->edge->id == cb->edge->id && ca->offsetCost <= cb->offsetCost) {
                sp_dist = cb->offsetCost - ca->offsetCost;
            }
            else if (ca->edge->id == cb->edge->id &&
                ca->offsetCost - cb->offsetCost < ca->edge->cost * reverse_tolerance)
            { //역방향으로 진행했지만, 해당 선분의 역방향 진행 허용치율 안쪽에 들어와 있으면, 0
                sp_dist = 0;
            }
            else if (ca->edge->target == cb->edge->source)
            { //연속된 선분에 매칭되어 있으면  -----.---o----.-----
              // Transition on the same OD nodes
                sp_dist = ca->edge->cost - ca->offsetCost + cb->offsetCost;
            }
            else
            { //인접한 것도 아니고 떨어져 있으면, 기존에 구한 최단거리등을 참고로 해서 선분 따라 거리를 구한다.
                Record* r = ubodt_->look_up(ca->edge->target, cb->edge->source);

                // No sp path exist from O to D.
                if (r == nullptr) //UBODT에 없으면
                {
                    //일단 nullptr이면 다시 최단거리를 찾는다.
                    //그래도 없으면 정말 안 닿는 곳이므로 infinity를 리턴한다.
                    double cost;
                    //bool isSPexist = graph_.shortest_path_dijkstra_dist(ca->edge->target, cb->edge->source, cost);
                    bool isSPexist = network_.shortest_path_astar_dist(ca->edge->target, cb->edge->source, cost);
                    if (isSPexist) {
                        sp_dist = cost + ca->edge->cost - ca->offsetCost + cb->offsetCost;
                    }
                    else {
                        //SPDLOG_WARN("최단거리가 존재하지 않음!!");
                        sp_dist = std::numeric_limits<double>::infinity();
                    }
                }
                else {
                    // calculate original SP distance
                    sp_dist = r->cost + ca->edge->cost - ca->offsetCost + cb->offsetCost;
                    ubodtCount++;
                }

            }
            return sp_dist;
        }
        /**
         * Update probabilities in a transition graph
         * @param tg transition graph
         * @param traj raw trajectory
         * @param config map match configuration
         */
        void update_tg(TransitionGraph* tg,
            const Trajectory& traj,
            double reverse_tolerance = 0) {
            auto startTime = std::chrono::high_resolution_clock::now();
            SPDLOG_DEBUG("Update transition graph");

            //tg는 traj를 이용해서 만든 것이지만, 밑에 계속 traj가 필요한 것 같다.
            //전체 점에 대한 후보군 집합을 꺼낸다. private이므로 getter이용
            std::vector<TGLayer>& layers = tg->get_layers();

            //궤적의 선분별 거리를 계산해둔다. 
            std::vector<double> trj_costs = network_.costDefault ?
                cal_eu_dist(traj.geom) : //기본 옵션이면 궤적 선분 길이로
                cal_timestamp_cost(traj.geom, traj.timestamps); //timestamp가 있으면, 시간 차로

            int N = layers.size(); //이 개수는 궤적의 점 개수와 일치
            tg->connectedToNext.resize(N, false); //연결안됨으로 초기화

            //궤적 점 하나하나를 돌면서 점수를 업데이트 한다.
            //처음부터 돈다.
            //https://ratsgo.github.io/data%20structure&algorithm/2017/11/14/viterbi/
            //위의 포스팅을 참고하면 이해가 쉽다.
            //각각의 단계에서 후행노드 기준으로,
            //하나의 후행노드는 가장 확률이 높은 선행 노드를 결정하게 된다.
            //후행노드는 모두 값을 갖게 되며
            //선행노드는 중복이나 미할당이 있을 수 있다.
            //1-2단계까지 나가면, 0-1단계에서 살아남은 1단계의 점들 몇개는
            //2단계에서 지목받지 못함으로써 더 이상 확률이 존재하지 않을 수도 있다.
            for (int i = 0; i < N - 1; ++i) {
                //SPDLOG_INFO("------------------------------------------------Update layer {} ", i);
                bool connected = false;
                update_layer(
                    i,         //궤적 점 번호
                    layers[i],       //i 번째 궤적점의 후보군
                    layers[i + 1],   //i+1 번째 궤적점의 후보군
                    trj_costs[i],        //궤적 선분의 길이
                    reverse_tolerance,  //거꾸로 갈 때 패널티
                    connected);        //그래서 연결되었는가? 초기에는 false
                //매칭 점수를 구하는 과정
                if (connected) {
                    tg->connectedToNext[i] = true; //연결관계를 기록한다.
                }
                else {
                    //tg->connectedToNext[i] = false; --> 초기값이 false
                    //SPDLOG_WARN("Traj {} unmatched as point {} and {} not connected", traj.id, i, i + 1);

                    //tg->print_optimal_info();
                    //break;
                    tg->reset_layer(&(layers[i + 1])); //끊겼으므로 리셋한다.
                }
            }
            SPDLOG_DEBUG("Update transition graph done");
            auto endTime = std::chrono::high_resolution_clock::now();
            size_t durTotal = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
            updateTGDuration += durTotal;

        }
        /**
         * Update probabilities between two layers a and b in the transition graph
         * @param level   the index of layer a
         * @param la_ptr  layer a
         * @param lb_ptr  layer b next to a
         * @param eu_dist Euclidean distance between two observed point
         * @param connected the variable is set to false if the layer is not connected
         * with the next layer
         */
        void update_layer(int level, TGLayer& la_ptr, TGLayer& lb_ptr,
            double trj_cost, double reverse_tolerance,
            bool& connected) {
            // SPDLOG_TRACE("Update layer");
            //std::vector<TGNode>& lb = lb_ptr;

              //이중 루프를 통해, gps 궤적하나와 그 다음에 대해서
              //후보군 쌍들을 각각 하나씩 차례로 훑어가면서
              //가장 큰 확률로 결정짓는다.
            bool layer_connected = false;
            for (TGNode& iter_a : la_ptr) //= la_ptr->begin(); iter_a != la_ptr->end(); ++iter_a) 
            {
                //NodeIndex source = iter_a.c->index;

                for (TGNode& iter_b : lb_ptr) //for (auto iter_b = lb_ptr->begin(); iter_b != lb_ptr->end(); ++iter_b) 
                {
                    //물리적 거리가 아닌 cost 기준으로 모두 수정
                    double sp_cost = get_sp_dist(
                        iter_a.c, iter_b.c,
                        reverse_tolerance);

                    if (sp_cost == std::numeric_limits<double>::infinity()) continue; //무한이면 건너뛴다.


                    //전이 확률도 시간 개념으로 바꿔줌
                    //전이 확률 = 직선거리 / 최단거리. 1을 넘을 수는 없다.
                    //SPDLOG_INFO("network_.costDefault :{}", network_.costDefault);

                    //이제 timestamp를 이용해서 eu_cost를 계산해야 한다.
                    ///그러기 위해서는 먼저 trajectory에 timestamp를 넣어야 한다.


                    //물리적 거리로 할지, 소요시간으로 할지는 앞에서 이미 계산해서 이 함수로 넘어왔다.
                    double tp = TransitionGraph::calc_tp(sp_cost, trj_cost);
                    //if (tp < 0.8) continue;
                    //double tp =
                    //if (tp < 0.8 || tp * iter_b.ep < 0.8) {
                    //    //SPDLOG_INFO("prob : {}", tp * iter_b.ep);
                    //    continue; //지나치게 돌아가면 해당 궤적은 건너 뛴다.
                    //}

                    double temp = iter_a.cumu_prob
                        + log(tp)         //이것저것 고려한 trasition prob. 계산. 다음 단계로의 전이 확률
                        + log(iter_b.ep); //떨어진 거리 정도에 따라 계산 //방사 확률

                    SPDLOG_TRACE("L {} f {} t {} sp {} dist {} tp {} ep {} fcp {} tcp {}",
                        level, iter_a->c->edge->id, iter_b->c->edge->id,
                        sp_dist, eu_dist, tp, iter_b->ep, iter_a->cumu_prob,
                        temp);

                    //선행 노드 값에 후행 노드의 이번 경우를 고려한 결과가
                    //다른 후행 노드 값보다 크면, 업데이트
                    //SPDLOG_INFO("a.cumu_prob : {} \t b.cumu_prob : {}", iter_a.cumu_prob, iter_b.cumu_prob);
                    if (temp >= iter_b.cumu_prob)
                    {
                        if (temp > -std::numeric_limits<double>::infinity()) {
                            layer_connected = true;
                        }
                        iter_b.cumu_prob = temp;
                        iter_b.prev = &(iter_a); //앞의 노드를 이번 계산 결과 고려해서 바꿔치기
                        iter_b.tp = tp;
                        iter_b.sp_dist = sp_cost;
                    }
                    loopCount++; //루프 수 구하기
                } //for (TGNode& iter_b : lb_ptr)  
            }//for (TGNode& iter_a : la_ptr)
          //if (connected!=nullptr){
            connected = layer_connected;

            //}
            // SPDLOG_TRACE("Update layer done");
        }


        C_Path construct_complete_path_astar(
            int traj_id, const TGOpath& path,
            const std::vector<FMM::Edge>& edges,
            std::vector<int>* indices,
            double reverse_tolerance
        ) const {
            C_Path cpath;
            if (!indices->empty()) indices->clear();
            if (path.empty()) return cpath;
            int N = path.size();
            cpath.push_back(path[0]->c->edge->id);
            int current_idx = 0;
            indices->push_back(current_idx);
            SPDLOG_TRACE("Insert index {}", current_idx);
            for (int i = 0; i < N - 1; ++i) {
                const Candidate* a = path[i]->c;
                const Candidate* b = path[i + 1]->c;
                SPDLOG_DEBUG("Check point {} a {} b {}", i, a->edge->id, b->edge->id);
                if ((a->edge->id != b->edge->id) || (a->offsetCost - b->offsetCost >
                    a->edge->cost * reverse_tolerance)) {
                    // segs stores edge index
                    bool connected = true;
                    std::vector<EdgeIndex> segs = network_.shortest_path_astar(a->edge->target, b->edge->source);
                    //std::vector<EdgeIndex> segs = graph_.shortest_path_dijkstra(a->edge->target, b->edge->source);
                    // No transition exist in UBODT

                    if (segs.empty() && a->edge->target != b->edge->source) {
                        SPDLOG_DEBUG("Edges not found connecting a b");
                        SPDLOG_DEBUG("reverse movement {} tolerance {}",
                            a->offset - b->offset, a->edge->length * reverse_tolerance);
                        SPDLOG_WARN("Traj {} unmatched as edge {} L {} offset {}"
                            " and edge {} L {} offset {} disconnected",
                            traj_id, a->edge->id, a->edge->cost, a->offsetCost,
                            b->edge->id, b->edge->cost, b->offsetCost);
                        connected = false; //아무일 안함
                        //indices->clear();
                        //return C_Path();
                    }

                    if (segs.empty()) {
                        SPDLOG_DEBUG("Edges ab are adjacent");
                    }
                    else {
                        SPDLOG_DEBUG("Edges connecting ab are {}", segs);
                    }

                    for (int e : segs) {
                        cpath.push_back(edges[e].id);
                        ++current_idx;
                    }

                    cpath.push_back(b->edge->id);
                    ++current_idx;
                    indices->push_back(current_idx);
                    SPDLOG_TRACE("Insert index {}", current_idx);
                }
                else {
                    indices->push_back(current_idx);
                    SPDLOG_TRACE("Insert index {}", current_idx);
                }
            }
            return cpath;
        }



    private:
        Network network_;
        std::shared_ptr<UBODT> ubodt_;
    };

}

#endif //FMM_FMM_ALGORITHM_H_