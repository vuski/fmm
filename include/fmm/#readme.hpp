/***
util.hpp				여러가지 파일 입출력 관련 함수						<-- spdlog.h

geometry.hpp			모든 파일의 기초									<-- 의존 없음

type.hpp				여러가지 기본 타입								<-- geometry.hpp

geom_algorithm.hpp		지오메트리 관련 여러가지 알고리즘					<-- type.hpp

transition_graph.hpp	HMM과 비터비 알고리즘 관련. 기본 타입과 함수		<-- type.hpp

gps_reader.hpp			궤적 파일 읽는데 필요. Trajectory 관리			<-- type.hpp

network.hpp				네트워크, 그래프 관리								#include "fmm/type.hpp"
																		#include "fmm/util.hpp"
																		#include "fmm/geom_algorithm.hpp"
																		#include "fmm/transition_graph.hpp"
																		#include "fmm/gps_reader.hpp"
																		rapidjson.hpp
																		THST/RTree.hpp

ubodt.hpp				네트워크를 ubodt 형식으로 저장					#include "fmm/util.hpp"
																		#include "fmm/type.hpp"
																		#include "fmm/network.hpp"
																		#include "fmm/transition_graph.hpp"

mm_writer.hpp			맵매칭 결과를 쓰는데 필요한 함수들					#include "fmm/util.hpp"
																		#include "fmm/type.hpp"
																		#include "fmm/network.hpp"

fmm_algorithm.hpp		주요 함수										#include "fmm/util.hpp"
																		#include "fmm/type.hpp"				<--geometry.hpp
																		#include "fmm/geom_algorithm.hpp"
																		#include "fmm/transition_graph.hpp"
																		#include "fmm/network.hpp"
																		#include "fmm/ubodt.hpp"
																		#include "fmm/mm_writer.hpp"
																		#include "fmm/gps_reader.hpp"
																		cxxopts.hpp
***/

///아래 세 부분 수정
//transition_graph.hpp
static double calc_ep(double dist, double error) {
	double a = std::max(dist - 5000.0, 0.0) / error;
	//double a = dist / error;
	return exp(-0.5 * a * a); //표준편차 공식에서 지수값만. 앞의 상수는 안쓴다
}

//fmm_flgorithm.hpp
double tp = TransitionGraph::calc_tp(sp_cost, trj_cost);
if (tp < 0.6) continue;
//double tp =


//network.hpp
double calc_heuristic_dist(
	const Point & p1, const Point & p2) const {
	double hdist = sqrt((p2.x - p1.x) * (p2.x - p1.x) +
		(p2.y - p1.y) * (p2.y - p1.y));
	//SPDLOG_INFO("calc_heuristic_dist network.costDefault :{}", costDefault);
	return 0;
	//return costDefault ? hdist : hdist / 2000;
	//pow(hdist / 10, 0.8) * 0.03; //가까우면 속도 느리고, 멀면 좀 빨라짐
//return 0;
}