
#define _CRT_SECURE_NO_WARNINGS

//�� hpp �ϳ��� ���������� ��� ��������
#include "fmm/fmm_algorithm.hpp"

//�Ʒ��� ���� ���̺귯����
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
size_t FastMapMatch::loopCount = 0; //���� ��
size_t FastMapMatch::updateTGDuration = 0; //���� ��
size_t Network::calcDuration_bg = 0; //���� ��
size_t Network::calcDuration_thst = 0; //���� ��


void main()
{



	spdlog::set_level(spdlog::level::level_enum::debug);

	spdlog::set_pattern("[%l][%s:%-3#] %v");

	//��Ʈ��ũ�� �д´�.

	//geojson���� ��ü(2021.11.24). GDAL ���̺귯�� �������


	bool timeCalc = true;

	int srid = 5179;
	Network network("..\\testData\\2019_network.geojson",
		"id", "source", "target",
		(timeCalc ? "timeFT" : ""),
		srid);

	//��Ʈ��ũ ����� ���� �ν�Ʈ �׷��� ����
	network.initNetworkGraph(); //Ư���� ���� ����, �տ��� ������ �Ÿ��� �״�� ����.


	//ubodt ������ �����. �����ϰ� ������.
	if (true)
	{		
		int threads = 30;

		std::string status = UBODT::generate_ubodt_new(
			(timeCalc ? "..\\testData\\ubodt1_time" :
				"..\\testData\\ubodt1_dist"),
			(timeCalc ? 5.0 : 1000.0),// <-- ������ ������ ��� �ݰ���� ubodt�� ������ ���ΰ�. �ð��� �� ���� ����
			network,
			threads);

		std::cout << status << std::endl;
	}


	if (true)
	{
		//�ٽ� ����� ubodt�д´�.   
		clock_t timeT = clock();
		std::shared_ptr<UBODT> ubodt = UBODT::read_ubodt_binary_new(
			(timeCalc ? "..\\testData\\ubodt1_time.dat" :
				"..\\testData\\ubodt1_dist.dat"),
			network);

		SPDLOG_INFO("time : {}ms", clock() - timeT);


		//������� �Ÿ�/����� ��Ʈ��ũ �ִܰŸ� ��� �� �ϳ���

		////�� ����
		timeT = clock();
		//���� �տ��� ���� �� ������ �����Ѵ�. �״�� �ű�⸸ ��
		FastMapMatch model;
		model.init(network, ubodt);



		//emission prob. ��
		//double a = dist / error;
		//exp(-0.5 * a * a); //ǥ������ ���Ŀ��� ��������. ���� ����� �Ⱦ���.
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

		//reverse tolerance�� 0~1���̰� �Ǿ�� ��. fmm_algorithm �� validate �Լ��� 0~1 �������� üũ�ϰ� ����
		FastMapMatchConfig config{
			48,      // �ĺ��� ����, ������ ������ ��Ʈ��ũ ���е� <-- �̰� ���ڶ�� �����ϰ� ���� �� �ִ�.
			1000,    // search Radius --> cost�� �����̵�, ���� �Ÿ��� ã�´�.
			2000,     //gps_error -> gps error --> emission Ȯ���� ����ϴµ� ���̸�, gps ���� ��Ʈ��ũ�� ���� �Ÿ��� ���ϴ°� �´�.
			0.0//������ reverse_tolerence ���⼭�� ������
		}; 
		SPDLOG_INFO("model loaded : {}ms", clock() - timeT);




		////����� �α� �б�
		//�������� �����͸� �����Ѵ�.
		std::vector<Trajectory> trajectories;
		if (true) {
			SignalReader::convertAndreadData("..\\testData\\testdata.geojson",
				trajectories, SignalReader::GEOJSON, SignalReader::YMDHMS);
		}			
		

		///���Ǵ� �Ÿ��� �� ���� ������ �ִ�.
		//�ϳ��� rtree� �ְų� ��Ʈ��ũ �����κ��� ������ �Ÿ��� �����ϴ� ��¥ ������ �Ÿ�
		//�ٸ� �ϳ��� ���Ÿ�. �ð����� �־ �ȴ�. offset�̳� �ִܰŸ� ���� length � ���δ�.


		clock_t startT = clock();
		//��� ����
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

				////�� ��Ī ���
				auto startTime = std::chrono::high_resolution_clock::now();
				MatchResult result = model.match_traj(trajectory, config);  // <-------- �� ��Ī ���
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
				SPDLOG_INFO("no: {} fmm done: {}ms   {}��° / {}s / {}%  ubodt��� : {}ȸ  ��ã�� : {}ȸ  ��ã�� �ð� : {}ms   �ĺ�ã�� : {}ms  ���ͺ� ���� : {}ȸ   updateTg : {}ms  ������ �ð� {}ms",
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

	std::cout << "�Ϸ�! " << std::endl;
};


