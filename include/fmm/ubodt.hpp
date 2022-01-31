/**
 * Fast map matching.
 *
 * Upperbounded origin destination table
 *
 * @author: Can Yang
 * @version: 2020.01.31
 * 
  * @revision : vuski
 * @version: 2021.11.29
 */

#ifndef FMM_UBODT_H_
#define FMM_UBODT_H_

#define _CRT_SECURE_NO_WARNINGS

#include <fstream>
#include <stdexcept>

#include "fmm/util.hpp"
#include "fmm/type.hpp"
#include "fmm/network.hpp"
#include "fmm/transition_graph.hpp"


#include <omp.h>

using namespace FMM;


namespace FMM {


    /**
     * %Record type of the upper bounded origin destination table
     */
    struct Record {
        NodeIndex source; /**< source node*/
        NodeIndex target; /**< target node*/
        NodeIndex first_n; /**< next node visited from source to target */
        NodeIndex prev_n; /**< last node visited before target */
        EdgeIndex next_e; /**< next edge visited from source to target */
        double cost; /**< distance from source to target */
        Record* next; /**< the next record stored in hashtable */
    };



    /**
     * Upperbounded origin destination table
     */
    class UBODT {
    public:
        UBODT(const UBODT&) = delete;
        UBODT& operator=(const UBODT&) = delete;
        /**
         * Constructor of UBODT from bucket number and multiplier
         * @param buckets_arg    Bucket number
         * @param multiplier_arg A multiplier used for querying, recommended to be
         * the number of nodes in the graph.
         */
        UBODT(int buckets_arg, int multiplier_arg) :
            buckets(buckets_arg), multiplier(multiplier_arg) {
            SPDLOG_TRACE("Intialization UBODT with buckets {} multiplier {}",
                buckets, multiplier);
            hashtable = (Record**)malloc(sizeof(Record*) * buckets);
            for (int i = 0; i < buckets; i++) {
                hashtable[i] = nullptr;
            }
            SPDLOG_TRACE("Intialization UBODT finished");
        }

        ~UBODT() {
            /* Clean hashtable */
            SPDLOG_TRACE("Clean UBODT");
            int i;
            for (i = 0; i < buckets; ++i) {
                Record* head = hashtable[i];
                Record* curr;
                while ((curr = head) != nullptr) {
                    head = head->next;
                    free(curr);
                }
            }
            // Destory hash table pointer
            free(hashtable);
            SPDLOG_TRACE("Clean UBODT finished");
        }
        /**
         * Look up the row according to a source node and a target node
         * @param  source source node
         * @param  target target node
         * @return  A row in the ubodt if the od pair is found, otherwise nullptr
         * is returned.
         */
        Record* look_up(NodeIndex source, NodeIndex target) const {
            unsigned int h = cal_bucket_index(source, target);
            Record* r = hashtable[h];
            while (r != nullptr) {
                if (r->source == source && r->target == target) {
                    return r;
                }
                else {
                    r = r->next; //하나의 버킷에 담은 그 다음 내용으로
                }
            }
            return r;
        }

        /**
         * Look up a shortest path (SP) containing edges from source to target.
         * In case that SP is not found, empty is returned.
         * @param  source source node
         * @param  target target node
         * @return  a shortest path connecting source to target
         */
        std::vector<EdgeIndex> look_sp_path(NodeIndex source,
            NodeIndex target) const {
            std::vector<EdgeIndex> edges;
            if (source == target) { return edges; }
            Record* r = look_up(source, target);
            // No transition exist from source to target
            if (r == nullptr) { return edges; }
            while (r->first_n != target) {
                edges.push_back(r->next_e);
                r = look_up(r->first_n, target);
            }
            edges.push_back(r->next_e);
            return edges;
        }


        /**
         * Construct the complete path (a vector of edge ID) from an optimal path
         * (a vector of optimal nodes in the transition graph)
         *
         * @param path an optimal path
         * @param edges a vector of edges
         * @param indices the index of each optimal edge in the complete path
         * @return a complete path (topologically connected).
         * If there is a large gap in the optimal
         * path implying complete path cannot be found in UBDOT,
         * an empty path is returned
         */
        C_Path construct_complete_path(
            int traj_id, const TGOpath& path,
            const std::vector<Edge>& edges,
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
                    std::vector<EdgeIndex> segs = look_sp_path(a->edge->target, b->edge->source);

                    // No transition exist in UBODT
                    if (segs.empty() && a->edge->target != b->edge->source) {
                        SPDLOG_DEBUG("Edges not found connecting a b");
                        SPDLOG_DEBUG("reverse movement {} tolerance {}",
                            a->offset - b->offset, a->edge->cost * reverse_tolerance);
                        SPDLOG_DEBUG("Traj {} unmatched as edge {} L {} offset {}"
                            " and edge {} L {} offset {} disconnected",
                            traj_id, a->edge->id, a->edge->cost, a->offset,
                            b->edge->id, b->edge->cost, b->offset);

                        indices->clear();
                        return C_Path();
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

        /**
         * Get the upperbound of the UBODT
         * @return upperbound value
         */
        double get_delta() const {
            return delta;
        }
        /**
         * Find the bucket index for an OD pair
         * @param  source origin/source node
         * @param  target destination/target node
         * @return  bucket index
         */
        unsigned int cal_bucket_index(NodeIndex source,
            NodeIndex target) const {
            return (source * multiplier + target) % buckets;
        }


        /**
         *  Insert a record into the hash table
         * @param r a record to be inserted
         */
        void insert(Record* r) {
            //int h = (r->source*multiplier+r->target)%buckets ;

            int h = cal_bucket_index(r->source, r->target);
            //첫번째 경우에는 초기값으로 nullptr이 들어 있으므로, r->next에는 null포인터가 들어간다.
            //두 번째 경우에는 기존에 저장된 값들이 r-> next로 밀려난다.
            r->next = hashtable[h];

            //그리고 아래에 비로소 r에 저장된 값들이 들어간다.
            hashtable[h] = r;
            if (r->cost > delta) delta = r->cost;
            ++num_rows;
        }

        inline long long get_num_rows() {
            return num_rows;
        };

        static std::shared_ptr<UBODT> read_ubodt_binary_new(const std::string& filename, Network& network)
        {


            SPDLOG_INFO("Reading UBODT file (binary format) from {}", filename);

            std::vector<Record> data;
            {
                clock_t t2 = clock();
                FILE* fp = fopen(filename.c_str(), "rb");
                //read size of file
                fseek(fp, 0, SEEK_END);
                //size_t size = ftell(fp);
                size_t size = _ftelli64(fp);

                fseek(fp, 0, SEEK_SET);
                data.resize(size / sizeof(Record));
                fread(&data[0], size, 1, fp);
                fclose(fp);
                clock_t t3 = clock();
                std::cout << filename + " 읽기시간 : " << t3 - t2 << std::endl;
            }

            size_t rows = data.size();

            int progress_step = 1000000;

            SPDLOG_INFO("Estimated rows is {}", rows);

            //내부변수인 bucket과 multiplier는 아래 단계에서 생성
            int buckets = find_prime_number(rows / LOAD_FACTOR); //
            int multiplier = network.get_node_count();
            std::shared_ptr<UBODT> table = std::make_shared<UBODT>(buckets, multiplier);
            long NUM_ROWS = 0;

            for (Record& rc : data)
            {
                NUM_ROWS++;
                Record* r = (Record*)malloc(sizeof(Record));
                r->source = rc.source;
                r->target = rc.target;
                r->first_n = rc.first_n;
                r->prev_n = rc.prev_n;
                r->next_e = rc.next_e;
                r->cost = rc.cost;
                r->next = nullptr;
                table->insert(r);
                if (NUM_ROWS % progress_step == 0) {
                    SPDLOG_INFO("Read rows {}", NUM_ROWS);
                }
            }
            double lf = NUM_ROWS / (double)buckets;
            SPDLOG_INFO("Estimated load factor #elements/#tablebuckets {}", lf);
            if (lf > 10) {
                SPDLOG_WARN("Load factor is too large.");
            }
            SPDLOG_INFO("Finish reading UBODT with rows {}", NUM_ROWS);
            return table;
        }

        /**
       * Find a large prime number according to input value
       * @param  value input value
       * @return  a large prime number
       */
        static int find_prime_number(double value) {
            std::vector<int> prime_numbers = {
                5003, 10039, 20029, 50047, 100669, 200003, 500009,
                1000039, 2000083, 5000101, 10000103, 20000033, 5000047 ,
                  100000037, 200000051, 500000069,
                  1000000021 };
            int N = prime_numbers.size();
            for (int i = 0; i < N; ++i) {
                if (value <= prime_numbers[i]) {
                    return prime_numbers[i];
                }
            }
            return prime_numbers[N - 1];
        }



        static std::string generate_ubodt_new(
            const std::string& filename, double delta,
            const Network& network_,
            int threads = 4)
        {
            using std::cout;
            using std::endl;

            //cout << "generate_ubodt_new graph.addr : " << (Network*)&network_ << endl;

            std::ostringstream oss;
            std::chrono::steady_clock::time_point begin =
                std::chrono::steady_clock::now();

            precompute_ubodt_omp(filename, delta, network_, threads);

            std::chrono::steady_clock::time_point end =
                std::chrono::steady_clock::now();
            double time_spent =
                std::chrono::duration_cast<std::chrono::milliseconds>
                (end - begin).count() / 1000.;
            oss << "Status: success\n";
            oss << "Time takes " << time_spent << " seconds\n";
            return oss.str();
        };


        constexpr static double LOAD_FACTOR = 1.0; /**< factor measuring the
                                                    average number of elements in
                                                    a bucket. */
        static const int BUFFER_LINE = 1024; /**< Number of characters to store in
                                                  a line */
    private:


        static void precompute_ubodt_omp(
            const std::string& filename, double delta,
            const Network& network_,
            int threads = 4)
        {
            int num_vertices = network_.get_num_vertices();
            int step_size = num_vertices / 10;
            if (step_size < 10) step_size = 10;

            SPDLOG_INFO("Start to generate UBODT with delta {}, number of Vertices : {}", delta, num_vertices);
            clock_t timeT = clock();

            size_t searchCount = 0;

            omp_set_num_threads(threads);
            std::vector<std::vector<Record>> source_map(threads);
            int progress = 0;
#pragma omp parallel
            {
#pragma omp for schedule(dynamic, 1) nowait
                for (int source = 0; source < num_vertices; ++source) {
                    int threadNum = omp_get_thread_num();
                    ++progress;
                    if (progress % step_size == 0) {
                        SPDLOG_INFO("Progress {} / {}", progress, num_vertices);
                    }

                    PredecessorMap predecessor;
                    DistanceMap distance;
                    std::stringstream node_output_buf;

                    network_.single_source_upperbound_dijkstra_new(source, delta, predecessor, distance);

                    for (auto& pd : predecessor)
                    { //기록된 무작위 순서로 시작

                        NodeIndex cur_node = pd.first; //현재 노드  --> 이번 기록에서 target이 된다.
                        if (cur_node != source) { //현재 노드가 source가 아니면,
                            NodeIndex prev_node = pd.second; //직전 노드
                            NodeIndex v = cur_node; // 현재 노드
                            NodeIndex u; //임시 노드
                            // When u=s, v is the next node visited
                            while ((u = predecessor[v]) != source) {
                                v = u; //현재 노드를 target으로 할 때, node visited after source를 구하는 과정. 거슬러 올라간다.
                            }
                            NodeIndex successor = v;
                            // Write the result
                            double cost = distance[successor];

                            EdgeIndex edge_index = network_.get_edge_index(source, successor);

                            source_map[threadNum].push_back(
                                { (NodeIndex)source,               //source: 이번 기록에서는 모두 같다.
                                 cur_node,        //target: for 구문으로 하나하나씩 돌아가는 중이다.
                                 successor,       //next_n: //위에서 while문으로 추적해가면서 구함
                                 prev_node,       //prev_n: //처음에 바로 구함
                                 edge_index,      //next_e:  //마지막에 source와 next_n으로 구함
                                 distance[cur_node],  //distance: //target에 대해서 기록해 놓았음
                                 nullptr });  /**< the next record stored in hashtable */
                        }
                    } //for 기록용

                }
            }
            SPDLOG_INFO("계산 완료. 벡터 통합");
            SPDLOG_INFO("UBODT 순수 계산시간 : {}ms", clock() - timeT);
            SPDLOG_INFO("UBODT 탐색 회수 : {}", searchCount);

            //새로운 공간 만들기
            size_t size = 0;
            for (auto& subvec : source_map) size += subvec.size();
            std::vector<Record> finaldata(size);

            size_t idx = 0;
            for (auto& subvec : source_map)
            {
                for (Record& r : subvec)
                {
                    finaldata[idx] = r;
                    idx++;
                }
            }
            SPDLOG_INFO("통합 완료. bin 저장");

            std::string fileNameBin = filename + ".dat";
            FILE* fp = fopen(fileNameBin.c_str(), "wb");
            fwrite(&finaldata[0], finaldata.size() * sizeof(Record), 1, fp);
            fclose(fp);

            if (true)
            {
                SPDLOG_INFO("bin 저장 완료. tsv 저장");

                std::string fileNameTsv = filename + ".tsv";
                std::ofstream myfile(fileNameTsv);
                myfile << "source\ttarget\tnext_n\tprev_n\tnext_e\tdistance\n";
                for (Record& r : finaldata)
                {
                    myfile << r.source << "\t"
                        << r.target << "\t"
                        << r.first_n << "\t"
                        << r.prev_n << "\t"
                        << r.next_e << "\t"
                        << r.cost << "\n";
                }
                myfile.close();
                SPDLOG_INFO("tsv 저장 완료");
            }
        }


        long long multiplier;   // multiplier to get a unique ID
        int buckets;   // number of buckets
        long long num_rows = 0;   // multiplier to get a unique ID
        double delta = 0.0;
        Record** hashtable;



    };

}

#endif //FMM_SRC_FMM_FFMM_UBODT_H_