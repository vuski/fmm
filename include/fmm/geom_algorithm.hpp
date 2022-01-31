/**
 * Fast map matching.
 * 
 * Algorithms for geoprocessing.
 *
 * Updated by Can Yang 2020-04-05
 * Updated by Diaolin 18.01.17
 * @author: Can Yang
 * @version: 2017.11.11
 * 
 * @revision : vuski
 * @version: 2021.11.29
 */

#ifndef FMM_GEOM_ALGORITHM_HPP
#define FMM_GEOM_ALGORITHM_HPP

#include <cmath>
#include <cstdlib>


#include "fmm/type.hpp"



namespace FMM {
/**
 * Algorithms for geoprocessing.
 */

/**
 * Calculate segment length of a trajectory
 * @param trajectory a trajectory as input
 * @return A vector of double values with size of N-1.
 * N is the number of points in trajectory.
 */
std::vector<double> cal_eu_dist(const FMM::LineString &trajectory) {
    int N = trajectory.get_num_points();
    std::vector<double> lengths(N - 1);
    double x0 = trajectory.get_x(0);
    double y0 = trajectory.get_y(0);
    for (int i = 1; i < N; ++i) {
        double x1 = trajectory.get_x(i);
        double y1 = trajectory.get_y(i);
        double dx = x1 - x0;
        double dy = y1 - y0;
        lengths[i - 1] = sqrt(dx * dx + dy * dy);
        x0 = x1;
        y0 = y1;
    }
    return lengths;
}


std::vector<double> cal_timestamp_cost(
    const FMM::LineString& trajectory,
    const std::vector<double>& timestamp) {

    //그냥 모두 시속 100km 로 간주.

    int N = trajectory.get_num_points();
    std::vector<double> times(N - 1);
    double x0 = trajectory.get_x(0);
    double y0 = trajectory.get_y(0);
    for (int i = 1; i < N; ++i) {
        double x1 = trajectory.get_x(i);
        double y1 = trajectory.get_y(i);
        double dx = x1 - x0;
        double dy = y1 - y0;
        times[i - 1] = sqrt(dx * dx + dy * dy) / (100.0 / 60.0 * 1000);
        x0 = x1;
        y0 = y1;
    }
    return times;


    if (false)
    {
        int N = timestamp.size();
        std::vector<double> duration(N - 1);

        for (int i = 0; i < N - 1; i++) {
            duration[i] = (timestamp[i + 1] - timestamp[i]) / 60.0; //분단위
        }

        return duration;
    }
}

/**
 * Concatenate a linestring segs to a linestring line, used in the
 * function network.complete_path_to_geometry
 *
 * @param line: linestring which will be updated
 * @param segs: segs that will be appended to line
 * @param offset: the number of points skipped in segs.
 */
void append_segs_to_line(FMM::LineString *line,
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
 * Create a new linestring as the reverse of an input linestring
 *
 * @param rhs the input linestring
 * @return A linestring containing points in the reverse order.
 */
FMM::LineString reverse_geometry(const FMM::LineString &rhs) {
    FMM::LineString line;
    int Npoints = rhs.get_num_points();
    for (int i = Npoints - 1; i >= 0; --i) {
        line.add_point(rhs.get_x(i), rhs.get_y(i));
    }
    return line;
}

/**
 * Interpolate a linestring at a fixed distance threshold and return the result
 * in as a vector of linestring
 *
 * @param line input line
 * @param delta the distance threshold to split the line
 * @return A vector of linestring
 */
std::vector<FMM::LineString> split_line(
    const FMM::LineString &line, double delta) {
    // SPDLOG_INFO("FMM::LineString {}", line.exportToWkt());
    std::vector<FMM::LineString> segments;
    if (line.get_length() <= delta) {
        segments.push_back(line);
        return segments;
    }
    int Npoints = line.get_num_points();
    if (Npoints < 2) {
        return segments;
    }
    double length_visited = 0;
    FMM::LineString interpolated_line;
    int i = 1;       // Point in line
    double cx = line.get_x(0);
    double cy = line.get_y(0);
    interpolated_line.add_point(cx, cy);
    while (i < Npoints) {
        double nx = line.get_x(i);
        double ny = line.get_y(i);
        double temp = std::sqrt((nx - cx) * (nx - cx) + (ny - cy) * (ny - cy));
        if (length_visited <= delta
            && length_visited + temp > delta) {
            // Split the edge
            double ratio = (delta - length_visited) / temp;
            // SPDLOG_INFO("Ratio is {}",ratio);
            cx = ratio * (nx - cx) + cx;
            cy = ratio * (ny - cy) + cy;
            interpolated_line.add_point(cx, cy);
            segments.push_back(interpolated_line);
            interpolated_line.clear();
            interpolated_line.add_point(cx, cy);
            length_visited = 0;
        }
        else {
            cx = nx;
            cy = ny;
            interpolated_line.add_point(cx, cy);
            length_visited += temp;
            ++i;
        }
    }
    segments.push_back(interpolated_line);
    // for (int i = 0;i<segments.size();++i){
    //   SPDLOG_INFO(" Segments {} {}", i, segments[i].exportToWkt());
    // }
    return segments;
}

/**
 * Interpolate a linestring according to a vector of distances to
 * the start point of a line
 * @param line
 * @param distances a vector of distance values to the start point
 * @return a linestring containing the interpolated points
 */
FMM::LineString interpolate_line_distances(
    const FMM::LineString &line, const std::vector<double> &distances) {
    FMM::LineString interpolated_line;
    int Npoints = line.get_num_points();
    int distance_count = distances.size();
    double length_visited = 0;
    int i = 0;       // Point in line
    int j = 0;       // Element visited in distances vector
    while (i < Npoints - 1 && j < distance_count) {
        double x1 = line.get_x(i);
        double y1 = line.get_y(i);
        double x2 = line.get_x(i + 1);
        double y2 = line.get_y(i + 1);
        double temp = std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
        while (j < distance_count && length_visited <= distances[j]
            && length_visited + temp > distances[j]) {
            double ratio = (distances[j] - length_visited) / temp;
            // SPDLOG_INFO("Ratio is {}",ratio);
            interpolated_line.add_point(ratio * (x2 - x1) + x1,
                ratio * (y2 - y1) + y1);
            ++j;
        }
        length_visited += temp;
        ++i;
    }
    if (length_visited < distances.back()) {
        interpolated_line.add_point(line.get_x(Npoints - 1),
            line.get_y(Npoints - 1));
    }
    return interpolated_line;
}

/**
 * Interpolate a linestring according to a distance step value
 * @param line the input linestring
 * @param distance the distance step value
 * @return a linestring containing the interpolated points
 */
FMM::LineString interpolate_line_distance(
    const FMM::LineString &line, double distance) {
    double length = line.get_length();
    int k = (int)ceil(length / distance);
    std::vector<double> distances;
    for (int i = 0; i < k + 1; ++i) {
        distances.push_back(distance * i);
    }
    return interpolate_line_distances(line, distances);
}

/**
 * Interpolate k points in a linestring with equal distance
 * @param line the input linestring
 * @param k interpolate k points in equal distances
 * @return a linestring containing the interpolated points
 */
FMM::LineString interpolate_line_kpoints(
    const FMM::LineString &line, int k) {
    double length = line.get_length();
    std::vector<double> distances;
    for (int i = 0; i < k; ++i) {
        distances.push_back((double)rand() / RAND_MAX * length);
    }
    std::sort(distances.begin(), distances.end());
    // SPDLOG_DEBUG("Length is {}",length);
    // SPDLOG_DEBUG("Sorted length is {}",UTIL::vec2string<double>(distances));
    return interpolate_line_distances(line, distances);
}

/**
 * Compute the boundary of an FMM::LineString and returns the result in
 * the passed x1,y1,x2,y2 variables.
 *
 * @param linestring input line
 * @param x1 the bottom left point x coordinate
 * @param y1 the bottom left point y coordinate
 * @param x2 the top right point x coordinate
 * @param y2 the top right point y coordinate
 */
void boundingbox_geometry(const FMM::LineString &linestring,
                          double &x1, double &y1, double &x2, double &y2) {
    int Npoints = linestring.get_num_points();
    x1 = DBL_MAX;
    y1 = DBL_MAX;
    x2 = DBL_MIN;
    y2 = DBL_MIN;
    double x, y;
    for (int i = 0; i < Npoints; ++i) {
        x = linestring.get_x(i);
        y = linestring.get_y(i);
        if (x < x1) x1 = x;
        if (y < y1) y1 = y;
        if (x > x2) x2 = x;
        if (y > y2) y2 = y;
    }
}

/**
 * Calculate the distance from each point in a linestring to the end point
 * of a linestring
 * @param geom input linestring
 * @return a vector of distance values
 */
std::vector<double> calc_length_to_end_vec(const FMM::LineString &geom) {
    int N = geom.get_num_points();
    if (N < 2) return std::vector<double>();
    std::vector<double> result(N - 1);
    double x0 = geom.get_x(0);
    double y0 = geom.get_y(0);
    for (int i = 1; i < N; ++i) {
        double x1 = geom.get_x(i);
        double y1 = geom.get_y(i);
        double dx = x1 - x0;
        double dy = y1 - y0;
        result[i - 1] = sqrt(dx * dx + dy * dy);
        x0 = x1;
        y0 = y1;
    }
    double temp = 0;
    for (int i = N - 2; i >= 0; --i) {
        result[i] = temp + result[i];
        temp = result[i];
    }
    return result;
}

/**
 * Calculate the closest point p' on a segment p1 (x1,y1) p2 (x2,y2) to a
 * specific point p (x,y)
 *
 * @param x
 * @param y
 * @param x1
 * @param y1
 * @param x2
 * @param y2
 * @param dist the distance from p to p'
 * @param offset the distance from p1 to p'
 */
void closest_point_on_segment(double x, double y, double x1, double y1,
                              double x2, double y2,
                              double *dist, double *offset) {
    double L2 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
    if (L2 == 0.0) {
        *dist = std::sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1));
        *offset = 0.0;
        return;
    }
    double x1_x = x - x1;
    double y1_y = y - y1;
    double x1_x2 = x2 - x1;
    double y1_y2 = y2 - y1;
    double ratio = (x1_x * x1_x2 + y1_y * y1_y2) / L2;
    ratio = (ratio > 1) ? 1 : ratio;
    ratio = (ratio < 0) ? 0 : ratio;
    double prj_x = x1 + ratio * (x1_x2);
    double prj_y = y1 + ratio * (y1_y2);
    *offset = std::sqrt((prj_x - x1) * (prj_x - x1) +
        (prj_y - y1) * (prj_y - y1));
    *dist = std::sqrt((prj_x - x) * (prj_x - x) + (prj_y - y) * (prj_y - y));
} // closest_point_on_segment

/**
 *
 * Calculate the closest point p' on a segment p1 (x1,y1) p2 (x2,y2) to a
 * specific point p (x,y)
 *
 * @param x
 * @param y
 * @param x1
 * @param y1
 * @param x2
 * @param y2
 * @param dist the distance from p to p'
 * @param offset the distance from p1 to p'
 * @param closest_x the x coordinate of p'
 * @param closest_y the y coordinate of p'
 */
void closest_point_on_segment(double x, double y, double x1, double y1,
                              double x2, double y2, double *dist,
                              double *offset, double *closest_x,
                              double *closest_y) {
    double L2 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
    if (L2 == 0.0) {
        *dist = std::sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1));
        *offset = 0.0;
        *closest_x = x1;
        *closest_y = y1;
        return;
    }
    double x1_x = x - x1;
    double y1_y = y - y1;
    double x1_x2 = x2 - x1;
    double y1_y2 = y2 - y1;
    double ratio = (x1_x * x1_x2 + y1_y * y1_y2) / L2;
    ratio = (ratio > 1) ? 1 : ratio;
    ratio = (ratio < 0) ? 0 : ratio;
    double prj_x = x1 + ratio * (x1_x2);
    double prj_y = y1 + ratio * (y1_y2);
    *offset = std::sqrt((prj_x - x1) * (prj_x - x1) +
        (prj_y - y1) * (prj_y - y1));
    *dist = std::sqrt((prj_x - x) * (prj_x - x) + (prj_y - y) * (prj_y - y));
    *closest_x = prj_x;
    *closest_y = prj_y;
} // closest_point_on_segment

/**
 * A linear referencing function
 *
 * Given a point p (px,py) and a polyline, find the closest point p' on line
 * and return the projected distance (p to p') and offset distance
 * (the distance along the polyline from its start to ')
 * @param px x coordinate of p
 * @param py y coordinate of p
 * @param linestring input line
 * @param result_dist the projected distance, the pointer will be updated
 * @param result_offset the offset distance, the pointer will be updated
 */
//나중에 지오메트리를 추출하기 위함이므로 물리적인 거리로 offset을 구한다.
void linear_referencing_physically(double px, double py,
                        const FMM::LineString &linestring,
                        double& result_dist,
                        double& result_offset) {
    int Npoints = linestring.get_num_points();
    double min_dist = DBL_MAX;
    double final_offset = DBL_MAX;
    double length_parsed = 0;
    int i = 0;
    // Iterating to check p(i) == p(i+2)
    // int seg_idx=0;
    while (i < Npoints - 1) {
        double x1 = linestring.get_x(i);
        double y1 = linestring.get_y(i);
        double x2 = linestring.get_x(i + 1);
        double y2 = linestring.get_y(i + 1);
        double temp_min_dist;
        double temp_min_offset;
        closest_point_on_segment(px, py, x1, y1, x2, y2,
            &temp_min_dist, &temp_min_offset);
        if (temp_min_dist < min_dist) {
            min_dist = temp_min_dist;
            final_offset = length_parsed + temp_min_offset;
        }
        length_parsed += std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
        ++i;
    }
    result_dist = min_dist;
    result_offset = final_offset; //라인을 따라 누적한 길이.
} // linear_referencing

/**
 * A linear referencing function
 *
 * Given a point p (px,py) and a polyline, find the closest point p' on line
 * and return the projected distance (p to p') and offset distance
 * (the distance along the polyline from its start to ')
 * @param px x coordinate of p
 * @param py y coordinate of p
 * @param linestring input line
 * @param result_dist the projected distance
 * @param result_offset the offset distance
 * @param proj_x x coordinate of p'
 * @param proj_y y coordinate of p'
 */
void linear_referencing(
    double px, double py, const FMM::LineString& linestring,
    double& result_dist, double& result_offset,
    double& proj_x, double& proj_y, double& total_length) {
    int Npoints = linestring.get_num_points();
    double min_dist = DBL_MAX;
    double temp_x = 0, temp_y = 0;
    double final_offset = DBL_MAX;
    double length_parsed = 0;
    int i = 0;
    // Iterating to check p(i) == p(i+2)
    // int seg_idx=0;
    while (i < Npoints - 1) {
        double x1 = linestring.get_x(i);
        double y1 = linestring.get_y(i);
        double x2 = linestring.get_x(i + 1);
        double y2 = linestring.get_y(i + 1);
        double temp_min_dist;
        double temp_min_offset;
        closest_point_on_segment(px, py, x1, y1, x2, y2,
            &temp_min_dist, &temp_min_offset,
            &temp_x, &temp_y);
        if (temp_min_dist < min_dist) {
            min_dist = temp_min_dist;
            final_offset = length_parsed + temp_min_offset;
            proj_x = temp_x;
            proj_y = temp_y;
        }
        length_parsed += std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
        ++i;
    }
    result_dist = min_dist;
    result_offset = final_offset;
    total_length = length_parsed;
} // linear_referencing

/**
 * Locate the point on a linestring according to the input of offset
 * The two pointer's target value will be updated.
 */
/**
 * Locate a point p on a linestring according to an offset value
 * @param linestring
 * @param offset the distance from p to start point of linestring
 * @param x the x coordinate of p
 * @param y the y coordinate of p
 */
void locate_point_by_offset(const FMM::LineString &linestring,
                            double offset, double *x, double *y) {
    int Npoints = linestring.get_num_points();
    if (offset <= 0.0) {
        *x = linestring.get_x(0);
        *y = linestring.get_y(0);
        return;
    }
    double L_processed = 0;       // length parsed
    int i = 0;
    double px = 0;
    double py = 0;
    bool found = false;
    // Find the idx of the point to be exported close to p
    while (i < Npoints - 1) {
        double x1 = linestring.get_x(i);
        double y1 = linestring.get_y(i);
        double x2 = linestring.get_x(i + 1);
        double y2 = linestring.get_y(i + 1);
        double deltaL = std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
        double ratio = (offset - L_processed) / deltaL;
        if (offset >= L_processed && offset <= L_processed + deltaL) {
            px = x1 + ratio * (x2 - x1);
            py = y1 + ratio * (y2 - y1);
            found = true;
            break;
        }
        ++i;
        L_processed += deltaL;
    }
    if (found) {
        *x = px;
        *y = py;
    }
    else {
        *x = linestring.get_x(Npoints - 1);
        *y = linestring.get_y(Npoints - 1);
    }
} // calculate_offset_point

/**
 * Cut a linestring at two offset values
 * @param linestring input line
 * @param offset1 starting offset, distance to the start point of linestring
 * @param offset2 ending offset, distance to the start point of linestring
 * @return a linestring containing only the part covering starting offset to
 * ending offset
 */
FMM::LineString cutoffseg_unique(
    const FMM::LineString &linestring, double offset1, double offset2) {
    FMM::LineString cutoffline;
    // SPDLOG_TRACE("Offset1 {} Offset2 {}", offset1, offset2);
    int Npoints = linestring.get_num_points();
    if (Npoints == 2) {
        // A single segment
        double x1 = linestring.get_x(0);
        double y1 = linestring.get_y(0);
        double x2 = linestring.get_x(1);
        double y2 = linestring.get_y(1);
        double L = std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
        double ratio1 = offset1 / L;
        double new_x1 = x1 + ratio1 * (x2 - x1);
        double new_y1 = y1 + ratio1 * (y2 - y1);
        double ratio2 = offset2 / L;
        double new_x2 = x1 + ratio2 * (x2 - x1);
        double new_y2 = y1 + ratio2 * (y2 - y1);
        cutoffline.add_point(new_x1, new_y1);
        cutoffline.add_point(new_x2, new_y2);
    }
    else {
        // Multiple segments
        double l1 = 0;
        double l2 = 0;
        int i = 0;
        while (i < Npoints - 1) {
            double x1 = linestring.get_x(i);
            double y1 = linestring.get_y(i);
            double x2 = linestring.get_x(i + 1);
            double y2 = linestring.get_y(i + 1);
            double deltaL = std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
            l2 = l1 + deltaL;
            // Insert p1
            // SPDLOG_TRACE("  L1 {} L2 {} ", l1, l2);
            if (l1 >= offset1 && l1 <= offset2) {
                cutoffline.add_point(x1, y1);
                // SPDLOG_TRACE("  add p1 {} {}", x1, y1);
            }

            // Insert p between p1 and p2
            if (offset1 > l1 && offset1 < l2) {
                double ratio1 = (offset1 - l1) / deltaL;
                double px = x1 + ratio1 * (x2 - x1);
                double py = y1 + ratio1 * (y2 - y1);
                cutoffline.add_point(px, py);
                // SPDLOG_TRACE("  add p {} {} between p1 p2", px, py);
            }

            if (offset2 > l1 && offset2 < l2) {
                double ratio2 = (offset2 - l1) / deltaL;
                double px = x1 + ratio2 * (x2 - x1);
                double py = y1 + ratio2 * (y2 - y1);
                cutoffline.add_point(px, py);
                // SPDLOG_TRACE("  add p {} {} between p1 p2", px, py);
            }

            // last point
            if (i == Npoints - 2 && offset2 >= l2) {
                cutoffline.add_point(x2, y2);
                // SPDLOG_TRACE("  add p2 {} {} for last point", x2, y2);
            }
            l1 = l2;
            ++i;
        }
    }
    return cutoffline;
} //cutoffseg_twoparameters

/**
 * Added by Diao 18.01.17
 * modified by Can 18.01.19
 * modified by Can 18.03.14
 *
 * Cut a linestring at an offset value according to a mode value.
 * If mode is 0, cutting from offset to the end of linestring,
 * otherwise cutting from the starting point to the offset value
 * @param linestring input linestring
 * @param offset offset value, distance to the start point of linestring
 * @param mode cutting mode value
 * @return a linestring that is cut from input linestring
 */
FMM::LineString cutoffseg(
    const FMM::LineString &linestring, double offset, int mode) {
    double L = linestring.get_length();
    double offset1, offset2;
    if (mode == 0) {
        offset1 = offset;
        offset2 = L;
    }
    else {
        offset1 = 0;
        offset2 = offset;
    }
    return FMM::cutoffseg_unique(linestring, offset1, offset2);
}


} // FMM
#endif /* FMM_ALGORITHM_HPP */