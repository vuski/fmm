
#define _CRT_SECURE_NO_WARNINGS

//이 hpp 하나에 연쇄적으로 모두 물려있음
#include "fmm/fmm_algorithm.hpp"

//아래가 의존 라이브러리들
//include spdlog
//include rapidjson
//include THST
//include cxxopts


using namespace FMM;

using std::cout;
using std::endl;

size_t Network::calcCount = 0;
size_t FastMapMatch::ubodtCount = 0;
size_t Network::calcDuration = 0;
size_t FastMapMatch::calcDuration = 0;
size_t FastMapMatch::loopCount = 0; //루프 수
size_t FastMapMatch::updateTGDuration = 0; //루프 수
size_t Network::calcDuration_bg = 0; //루프 수
size_t Network::calcDuration_thst = 0; //루프 수


void main()
{



	spdlog::set_level(spdlog::level::level_enum::debug);

	spdlog::set_pattern("[%l][%s:%-3#] %v");

	//네트워크를 읽는다.

	//geojson으로 교체(2021.11.24). GDAL 라이브러리 떼어버림


	bool timeCalc = true;

	int srid = 5179;
	Network network("..\\testData\\2019_network.geojson",
		"id", "source", "target",
		(timeCalc ? "timeFT" : ""),
		srid);

	//네트워크 계산을 위한 부스트 그래프 형식
	network.initNetworkGraph(); //특별한 것은 없고, 앞에서 설정한 거리가 그대로 들어간다.


	//ubodt 파일을 만든다. 저장하고 끝난다.
	if (true)
	{		
		int threads = 30;

		std::string status = UBODT::generate_ubodt_new(
			(timeCalc ? "..\\testData\\ubodt1_time" :
				"..\\testData\\ubodt1_dist"),
			(timeCalc ? 5.0 : 1000.0),// <-- 각각의 점에서 어느 반경까지 ubodt를 저장할 것인가. 시간이 될 수도 있음
			network,
			threads);

		std::cout << status << std::endl;
	}


	if (true)
	{
		//다시 저장된 ubodt읽는다.   
		clock_t timeT = clock();
		std::shared_ptr<UBODT> ubodt = UBODT::read_ubodt_binary_new(
			(timeCalc ? "..\\testData\\ubodt1_time.dat" :
				"..\\testData\\ubodt1_dist.dat"),
			network);

		SPDLOG_INFO("time : {}ms", clock() - timeT);


		//여기까지 거리/비용은 네트워크 최단거리 계산 딱 하나임

		////모델 구성
		timeT = clock();
		//역시 앞에서 읽은 세 형식을 저장한다. 그대로 옮기기만 함
		FastMapMatch model;
		model.init(network, ubodt);



		//emission prob. 는
		//double a = dist / error;
		//exp(-0.5 * a * a); //표준편차 공식에서 지수값만. 앞의 상수는 안쓴다.
		//	  0.1   0.995012479
		//    0.2   0.980198673
		//    0.3	0.955997482
		//    0.4	0.923116346
		//    0.5	0.882496903
		//    0.6	0.835270211
		//    0.7	0.782704538
		//    0.8	0.726149037
		//    0.9	0.666976811
		//    1.0   0.60653066
		//    1.1	0.546074427
		//    1.2	0.486752256
		//    1.3	0.429557358
		//    1.4	0.375311099
		//    1.5	0.324652467
		//    1.6	0.2780373
		//    1.7	0.235746077
		//    1.8	0.197898699
		//    1.9	0.164474457
		//    2.0   0.135335283
		//    2.1	0.110250525
		//    2.2	0.088921617
		//    2.3	0.071005354
		//    2.4	0.056134763
		//    2.5	0.043936934
		//    2.6	0.034047455

		//reverse tolerance는 0~1사이가 되어야 함. fmm_algorithm 의 validate 함수에 0~1 사이임을 체크하고 있음
		FastMapMatchConfig config{
			48,      // 후보군 숫자, 궤적에 근접한 네트워크 선분들 <-- 이게 모자라면 복잡하게 꼬일 수 있다.
			1000,    // search Radius --> cost가 무엇이든, 실제 거리로 찾는다.
			2000,     //gps_error -> gps error --> emission 확률을 계산하는데 쓰이며, gps 점과 네트워크의 실제 거리로 비교하는게 맞다.
			0.0//마지막 reverse_tolerence 여기서는 생략함
		}; 
		SPDLOG_INFO("model loaded : {}ms", clock() - timeT);




		////사용자 로그 읽기
		//오리지널 데이터를 저장한다.
		std::vector<Trajectory> trajectories;
		if (true) {
			SignalReader::convertAndreadData("..\\testData\\testdata.geojson",
				trajectories, SignalReader::GEOJSON, SignalReader::YMDHMS);
		}			
		

		///사용되는 거리는 두 가지 개념이 있다.
		//하나는 rtree등에 넣거나 네트워크 선으로부터 떨어진 거리를 측정하는 진짜 물리적 거리
		//다른 하나는 비용거리. 시간으로 넣어도 된다. offset이나 최단거리 계산용 length 등에 쓰인다.


		clock_t startT = clock();
		//결과 쓰기
		size_t calcCnt = 0;
		size_t calcDur = 0;
		ResultConfig result_config;
		result_config.output_config.write_cpath = true;

		CSVMatchResultWriter writer("..\\testData\\result.tsv",	result_config.output_config);
		int numThreads = 30;
	
		int cnt = 0;
		omp_set_num_threads(numThreads);
#pragma omp parallel
		{
			int threadNum = omp_get_thread_num();

#pragma omp for schedule(dynamic, 1) nowait
			for (int i = 0; i < trajectories.size(); i++)
			{
				Trajectory& trajectory = trajectories[i];
				Network::calcCount = 0;
				FastMapMatch::ubodtCount = 0;
				Network::calcDuration = 0;
				FastMapMatch::calcDuration = 0;
				FastMapMatch::loopCount = 0;
				FastMapMatch::updateTGDuration = 0;

				////맵 매칭 계산
				auto startTime = std::chrono::high_resolution_clock::now();
				MatchResult result = model.match_traj(trajectory, config);  // <-------- 맵 매칭 계산
				auto endTime = std::chrono::high_resolution_clock::now();
	

				double durTotal = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count() / 1000.0;
				double countSP = Network::calcCount;
				double countUbodt = FastMapMatch::ubodtCount;
				double durSP = (double)Network::calcDuration / 1000.0;
				double durSearchCddt = FastMapMatch::calcDuration / 1000.0;
				double durUpdateTg = FastMapMatch::updateTGDuration / 1000.0;

#pragma omp critical
				writer.write_result(trajectory, result);
#pragma omp critical
				cnt++;				
#pragma omp critical
				SPDLOG_INFO("no: {} fmm done: {}ms   {}번째 / {}s / {}%  ubodt사용 : {}회  길찾기 : {}회  길찾기 시간 : {}ms   후보찾기 : {}ms  비터비 루프 : {}회   updateTg : {}ms  나머지 시간 {}ms",
					i,
					durTotal,
					cnt,
					(int)((clock() - startT) / 1000),
					(int)((double)i / trajectories.size() * 100000) / 1000.0,
					countUbodt,
					countSP,
					durSP,
					durSearchCddt,
					FastMapMatch::loopCount,
					durUpdateTg,
					durTotal - durSearchCddt - durUpdateTg);
				
			}
		} //#pragma omp parallel

	}

	std::cout << "완료! " << std::endl;
};


