/***
util.hpp				�������� ���� ����� ���� �Լ�						<-- spdlog.h

geometry.hpp			��� ������ ����									<-- ���� ����

type.hpp				�������� �⺻ Ÿ��								<-- geometry.hpp

geom_algorithm.hpp		������Ʈ�� ���� �������� �˰���					<-- type.hpp

transition_graph.hpp	HMM�� ���ͺ� �˰��� ����. �⺻ Ÿ�԰� �Լ�		<-- type.hpp

gps_reader.hpp			���� ���� �дµ� �ʿ�. Trajectory ����			<-- type.hpp

network.hpp				��Ʈ��ũ, �׷��� ����								#include "fmm/type.hpp"
																		#include "fmm/util.hpp"
																		#include "fmm/geom_algorithm.hpp"
																		#include "fmm/transition_graph.hpp"
																		#include "fmm/gps_reader.hpp"
																		rapidjson.hpp
																		THST/RTree.hpp

ubodt.hpp				��Ʈ��ũ�� ubodt �������� ����					#include "fmm/util.hpp"
																		#include "fmm/type.hpp"
																		#include "fmm/network.hpp"
																		#include "fmm/transition_graph.hpp"

mm_writer.hpp			�ʸ�Ī ����� ���µ� �ʿ��� �Լ���					#include "fmm/util.hpp"
																		#include "fmm/type.hpp"
																		#include "fmm/network.hpp"

fmm_algorithm.hpp		�ֿ� �Լ�										#include "fmm/util.hpp"
																		#include "fmm/type.hpp"				<--geometry.hpp
																		#include "fmm/geom_algorithm.hpp"
																		#include "fmm/transition_graph.hpp"
																		#include "fmm/network.hpp"
																		#include "fmm/ubodt.hpp"
																		#include "fmm/mm_writer.hpp"
																		#include "fmm/gps_reader.hpp"
																		cxxopts.hpp
***/

///�Ʒ� �� �κ� ����
//transition_graph.hpp
static double calc_ep(double dist, double error) {
	double a = std::max(dist - 5000.0, 0.0) / error;
	//double a = dist / error;
	return exp(-0.5 * a * a); //ǥ������ ���Ŀ��� ��������. ���� ����� �Ⱦ���
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
	//pow(hdist / 10, 0.8) * 0.03; //������ �ӵ� ������, �ָ� �� ������
//return 0;
}