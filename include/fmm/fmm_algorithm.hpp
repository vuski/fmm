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
         //reverse tolerance�� 0~1���̰� �Ǿ�� �ϴ� �� ����
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

        //�ҿ�ð� ������
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

            //config���� ������ �ݰ����� ã�µ�, ���� cost�� �Ÿ��� �ð��̵� ��� ����. �׳� �Ÿ��� ����Ѵ�.
            auto startTime = std::chrono::high_resolution_clock::now();
            Traj_Candidates tc = network_.search_tr_cs_knn(
                traj.geom, config.k, config.radius);
            //ã�� ���е��� �� ���ڴ� ��ü  edge ���ڿ� �� ������ ���鸶�� radius ���� ����� ������ �󸶳� �ִ��� �� ���ڸ� ��� ���� ��
            auto endTime = std::chrono::high_resolution_clock::now();
            calcDuration += std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();

            SPDLOG_DEBUG("Trajectory candidate {}", tc);
            if (tc.empty()) return MatchResult{};
            SPDLOG_DEBUG("Generate transition graph");

            //gps error�� �����, ���� ã�� ����� ����. \
            //���� ����� ���� ���� TGLayer ����ü�� ����� ������ ����д�.
            TransitionGraph tg(tc, config.gps_error);

            SPDLOG_DEBUG("Update cost in transition graph");
            // The network will be used internally to update transition graph
            update_tg(&tg, traj, config.reverse_tolerance);

            SPDLOG_DEBUG("Optimal path inference");
            //������ ������ ���� �Ҵ�� Ȯ������ ���ϸ鼭 �Լ� �������ڸ��� �ѹ��� ã�´�.
            //�� �������� ������ �����ذ��鼭 �ϳ���
            TGOpath tg_opath = tg.backtrack(traj.id);
            //o path�� gps ���� ������ŭ�� �����Ѵ�.


            //������ʹ� ����� �����ϴ� �κ�
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
            if (cpath.empty())  //���� ������, ubodt�� �������� �ʾƼ� �׷���.
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
            { //���������� ����������, �ش� ������ ������ ���� ���ġ�� ���ʿ� ���� ������, 0
                sp_dist = 0;
            }
            else if (ca->edge->target == cb->edge->source)
            { //���ӵ� ���п� ��Ī�Ǿ� ������  -----.---o----.-----
              // Transition on the same OD nodes
                sp_dist = ca->edge->cost - ca->offsetCost + cb->offsetCost;
            }
            else
            { //������ �͵� �ƴϰ� ������ ������, ������ ���� �ִܰŸ����� ����� �ؼ� ���� ���� �Ÿ��� ���Ѵ�.
                Record* r = ubodt_->look_up(ca->edge->target, cb->edge->source);

                // No sp path exist from O to D.
                if (r == nullptr) //UBODT�� ������
                {
                    //�ϴ� nullptr�̸� �ٽ� �ִܰŸ��� ã�´�.
                    //�׷��� ������ ���� �� ��� ���̹Ƿ� infinity�� �����Ѵ�.
                    double cost;
                    //bool isSPexist = graph_.shortest_path_dijkstra_dist(ca->edge->target, cb->edge->source, cost);
                    bool isSPexist = network_.shortest_path_astar_dist(ca->edge->target, cb->edge->source, cost);
                    if (isSPexist) {
                        sp_dist = cost + ca->edge->cost - ca->offsetCost + cb->offsetCost;
                    }
                    else {
                        //SPDLOG_WARN("�ִܰŸ��� �������� ����!!");
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

            //tg�� traj�� �̿��ؼ� ���� ��������, �ؿ� ��� traj�� �ʿ��� �� ����.
            //��ü ���� ���� �ĺ��� ������ ������. private�̹Ƿ� getter�̿�
            std::vector<TGLayer>& layers = tg->get_layers();

            //������ ���к� �Ÿ��� ����صд�. 
            std::vector<double> trj_costs = network_.costDefault ?
                cal_eu_dist(traj.geom) : //�⺻ �ɼ��̸� ���� ���� ���̷�
                cal_timestamp_cost(traj.geom, traj.timestamps); //timestamp�� ������, �ð� ����

            int N = layers.size(); //�� ������ ������ �� ������ ��ġ
            tg->connectedToNext.resize(N, false); //����ȵ����� �ʱ�ȭ

            //���� �� �ϳ��ϳ��� ���鼭 ������ ������Ʈ �Ѵ�.
            //ó������ ����.
            //https://ratsgo.github.io/data%20structure&algorithm/2017/11/14/viterbi/
            //���� �������� �����ϸ� ���ذ� ����.
            //������ �ܰ迡�� ������ ��������,
            //�ϳ��� ������� ���� Ȯ���� ���� ���� ��带 �����ϰ� �ȴ�.
            //������� ��� ���� ���� �Ǹ�
            //������� �ߺ��̳� ���Ҵ��� ���� �� �ִ�.
            //1-2�ܰ���� ������, 0-1�ܰ迡�� ��Ƴ��� 1�ܰ��� ���� ���
            //2�ܰ迡�� ������� �������ν� �� �̻� Ȯ���� �������� ���� ���� �ִ�.
            for (int i = 0; i < N - 1; ++i) {
                //SPDLOG_INFO("------------------------------------------------Update layer {} ", i);
                bool connected = false;
                update_layer(
                    i,         //���� �� ��ȣ
                    layers[i],       //i ��° �������� �ĺ���
                    layers[i + 1],   //i+1 ��° �������� �ĺ���
                    trj_costs[i],        //���� ������ ����
                    reverse_tolerance,  //�Ųٷ� �� �� �г�Ƽ
                    connected);        //�׷��� ����Ǿ��°�? �ʱ⿡�� false
                //��Ī ������ ���ϴ� ����
                if (connected) {
                    tg->connectedToNext[i] = true; //������踦 ����Ѵ�.
                }
                else {
                    //tg->connectedToNext[i] = false; --> �ʱⰪ�� false
                    //SPDLOG_WARN("Traj {} unmatched as point {} and {} not connected", traj.id, i, i + 1);

                    //tg->print_optimal_info();
                    //break;
                    tg->reset_layer(&(layers[i + 1])); //�������Ƿ� �����Ѵ�.
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

              //���� ������ ����, gps �����ϳ��� �� ������ ���ؼ�
              //�ĺ��� �ֵ��� ���� �ϳ��� ���ʷ� �Ⱦ�鼭
              //���� ū Ȯ���� �������´�.
            bool layer_connected = false;
            for (TGNode& iter_a : la_ptr) //= la_ptr->begin(); iter_a != la_ptr->end(); ++iter_a) 
            {
                //NodeIndex source = iter_a.c->index;

                for (TGNode& iter_b : lb_ptr) //for (auto iter_b = lb_ptr->begin(); iter_b != lb_ptr->end(); ++iter_b) 
                {
                    //������ �Ÿ��� �ƴ� cost �������� ��� ����
                    double sp_cost = get_sp_dist(
                        iter_a.c, iter_b.c,
                        reverse_tolerance);

                    if (sp_cost == std::numeric_limits<double>::infinity()) continue; //�����̸� �ǳʶڴ�.


                    //���� Ȯ���� �ð� �������� �ٲ���
                    //���� Ȯ�� = �����Ÿ� / �ִܰŸ�. 1�� ���� ���� ����.
                    //SPDLOG_INFO("network_.costDefault :{}", network_.costDefault);

                    //���� timestamp�� �̿��ؼ� eu_cost�� ����ؾ� �Ѵ�.
                    ///�׷��� ���ؼ��� ���� trajectory�� timestamp�� �־�� �Ѵ�.


                    //������ �Ÿ��� ����, �ҿ�ð����� ������ �տ��� �̹� ����ؼ� �� �Լ��� �Ѿ�Դ�.
                    double tp = TransitionGraph::calc_tp(sp_cost, trj_cost);
                    //if (tp < 0.8) continue;
                    //double tp =
                    //if (tp < 0.8 || tp * iter_b.ep < 0.8) {
                    //    //SPDLOG_INFO("prob : {}", tp * iter_b.ep);
                    //    continue; //����ġ�� ���ư��� �ش� ������ �ǳ� �ڴ�.
                    //}

                    double temp = iter_a.cumu_prob
                        + log(tp)         //�̰����� ����� trasition prob. ���. ���� �ܰ���� ���� Ȯ��
                        + log(iter_b.ep); //������ �Ÿ� ������ ���� ��� //��� Ȯ��

                    SPDLOG_TRACE("L {} f {} t {} sp {} dist {} tp {} ep {} fcp {} tcp {}",
                        level, iter_a->c->edge->id, iter_b->c->edge->id,
                        sp_dist, eu_dist, tp, iter_b->ep, iter_a->cumu_prob,
                        temp);

                    //���� ��� ���� ���� ����� �̹� ��츦 ����� �����
                    //�ٸ� ���� ��� ������ ũ��, ������Ʈ
                    //SPDLOG_INFO("a.cumu_prob : {} \t b.cumu_prob : {}", iter_a.cumu_prob, iter_b.cumu_prob);
                    if (temp >= iter_b.cumu_prob)
                    {
                        if (temp > -std::numeric_limits<double>::infinity()) {
                            layer_connected = true;
                        }
                        iter_b.cumu_prob = temp;
                        iter_b.prev = &(iter_a); //���� ��带 �̹� ��� ��� ����ؼ� �ٲ�ġ��
                        iter_b.tp = tp;
                        iter_b.sp_dist = sp_cost;
                    }
                    loopCount++; //���� �� ���ϱ�
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
                        connected = false; //�ƹ��� ����
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