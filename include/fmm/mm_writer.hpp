/**
 * Fast map matching.
 *
 * Definition of MatchResultWriter Class, which contains functions for
 * writing the results.
 *
 * @author: Can Yang
 * @version: 2017.11.11
 * 
 * @revision : vuski
 * @version: 2021.11.29
 */


#ifndef FMM_MM_WRITER_HPP
#define FMM_MM_WRITER_HPP


#define _CRT_SECURE_NO_WARNINGS

#include <sstream>
#include <iostream>
#include <fstream>
#include <omp.h>

#include "fmm/util.hpp"
#include "fmm/type.hpp"
#include "fmm/network.hpp"



namespace FMM {

    /**
    *  Output configuration class defining the fields exported for map matching
    */
    struct OutputConfig {
        bool write_opath = false; /**< if true, opath (edge id matched for
                                        each point) will be exported */
        bool write_offset = false; /**< if true, offset (distance to the start point
                                        of a matched edge) will be exported */
        bool write_error = false; /**< if true, gps error (distance from GPS point to
                                        the matched point) will be exported */
        bool write_cpath = true; /**< if true, cpath (a list of edge ID representing
                                        the matched path) will be exported */
        bool write_tpath = false; /**< if true, tpath (the path traversed between
                                        each two consecutive observations)
                                        will be exported */
        bool write_mgeom = true; /**< if true, mgeom (the geometry of
                                        the matched path) will be exported */
        bool write_spdist = false; /**< if true, spdist (the distance traversed
                                        between each two consecutive points)
                                        will be exported */
        bool write_pgeom = false; /**< if true, pgeom (a linestring connecting
                                        the matched point) will be exported */
        bool write_ep = false; /**< if true, ep (emission proability of each point)
                                    will be exported */
        bool write_tp = false; /**< if true, tp (transition probability) will be
                                    exported */
        bool write_length = false; /**< if true, length (length of each matched edge)
                                        will be exported */
        bool write_duration = false; /**< if true, duration (time difference between
                                        two points) will be exported */
        bool write_speed = false; /**< if true, speed (sp_dist/duration)
                                        will be exported */
    };

    /**
        * Result Configuration class, defining output file and output fields
        */
    struct ResultConfig {
        std::string file; /**< Output file to write the result */
        OutputConfig output_config; /**< Output fields to export */
        /**
            * Check the validation of the configuration
            * @return true if valid otherwise false
            */

        bool validate() const {
            if (FMM::file_exists(file))
            {
                SPDLOG_WARN("Overwrite existing result file {}", file);
            };
            std::string output_folder = get_file_directory(file);
            if (!FMM::folder_exist(output_folder)) {
                SPDLOG_CRITICAL("Output folder {} not exists", output_folder);
                return false;
            }
            return true;
        };


        /**
            * Print the configuration information
            */
        void print() const {
            std::ostringstream ss;
            if (output_config.write_opath)
                ss << "opath ";
            if (output_config.write_pgeom)
                ss << "pgeom ";
            if (output_config.write_offset)
                ss << "offset ";
            if (output_config.write_error)
                ss << "error ";
            if (output_config.write_spdist)
                ss << "spdist ";
            if (output_config.write_cpath)
                ss << "cpath ";
            if (output_config.write_tpath)
                ss << "tpath ";
            if (output_config.write_mgeom)
                ss << "mgeom ";
            if (output_config.write_ep)
                ss << "ep ";
            if (output_config.write_tp)
                ss << "tp ";
            if (output_config.write_length)
                ss << "length ";
            if (output_config.write_duration)
                ss << "duration ";
            if (output_config.write_speed)
                ss << "speed ";
            SPDLOG_INFO("ResultConfig");
            SPDLOG_INFO("File: {}", file);
            SPDLOG_INFO("Fields: {}", ss.str());
        };

        std::string to_string() const {
            std::ostringstream oss;
            oss << "Result file : " << file << "\n";
            oss << "Output fields: ";
            if (output_config.write_opath)
                oss << "opath ";
            if (output_config.write_pgeom)
                oss << "pgeom ";
            if (output_config.write_offset)
                oss << "offset ";
            if (output_config.write_error)
                oss << "error ";
            if (output_config.write_spdist)
                oss << "spdist ";
            if (output_config.write_cpath)
                oss << "cpath ";
            if (output_config.write_tpath)
                oss << "tpath ";
            if (output_config.write_mgeom)
                oss << "mgeom ";
            if (output_config.write_ep)
                oss << "ep ";
            if (output_config.write_tp)
                oss << "tp ";
            if (output_config.write_length)
                oss << "length ";
            if (output_config.write_duration)
                oss << "duration ";
            if (output_config.write_speed)
                oss << "speed ";
            return oss.str();
        };
        /**
            * Parse a string separated by , into a set of strings.
            *
            * This is used for extracting the output fields from argument
            *
            * @param s input string
            * @return a set of strings
            */

        std::set<std::string> string2set(
            const std::string& s) {
            char delim = ',';
            std::set<std::string> result;
            std::stringstream ss(s);
            std::string intermediate;
            while (getline(ss, intermediate, delim)) {
                result.insert(intermediate);
            }
            return result;
        };

    };
    /**
     * An interface defined for writing the map match result
     */
    class MatchResultWriter {
    public:
        /**
         * Write the match result to a file
         * @param result the match result of a trajectory
         */
        virtual void write_result(
            const FMM::Trajectory& traj,
            const FMM::MatchResult& result) = 0;
    };

    /**
     * A writer class for writing matche result to a CSV file.
     */
    class CSVMatchResultWriter : public MatchResultWriter {
    public:
        /**
         * Constructor
         *
         * The output fields are defined only once and later all the match result
         * will be exported according to that configuration.
         *
         * @param result_file the filename to write result
         * @param config_arg the fields that will be exported
         *
         */
        CSVMatchResultWriter(const std::string& result_file,
            const FMM::OutputConfig& config_arg) :
            m_fstream(result_file), config_(config_arg) {
            write_header();
        }
        /**
         * Write a header line for the fields exported
         */
        void write_header() {
            std::string header = "id";
            if (config_.write_opath) header += ";opath";
            if (config_.write_error) header += ";error";
            if (config_.write_offset) header += ";offset";
            if (config_.write_spdist) header += ";spdist";
            if (config_.write_pgeom) header += ";pgeom";
            if (config_.write_cpath) header += ";cpath";
            if (config_.write_tpath) header += ";tpath";
            if (config_.write_mgeom) header += ";mgeom";
            if (config_.write_ep) header += ";ep";
            if (config_.write_tp) header += ";tp";
            if (config_.write_length) header += ";length";
            if (config_.write_duration) header += ";duration";
            if (config_.write_speed) header += ";speed";
            m_fstream << header << '\n';
        }
        /**
         * Write match result
         * @param traj Input trajectory
         * @param result Map match result
         */
        void write_result(const FMM::Trajectory& traj,
            const FMM::MatchResult& result) {
            
            if (result.cpath.size() == 0) return;

            std::stringstream buf;
            buf << result.id;
            if (config_.write_opath) {
                buf << ";" << result.opath;
            }
            if (config_.write_error) {
                buf << ";";
                if (!result.opt_candidate_path.empty()) {
                    int N = result.opt_candidate_path.size();
                    for (int i = 0; i < N - 1; ++i) {
                        buf << result.opt_candidate_path[i].c.dist << ",";
                    }
                    buf << result.opt_candidate_path[N - 1].c.dist;
                }
            }
            if (config_.write_offset) {
                buf << ";";
                if (!result.opt_candidate_path.empty()) {
                    int N = result.opt_candidate_path.size();
                    for (int i = 0; i < N - 1; ++i) {
                        buf << result.opt_candidate_path[i].c.offsetLength << ",";
                    }
                    buf << result.opt_candidate_path[N - 1].c.offsetLength;
                }
            }
            if (config_.write_spdist) {
                buf << ";";
                if (!result.opt_candidate_path.empty()) {
                    int N = result.opt_candidate_path.size();
                    for (int i = 1; i < N; ++i) {
                        buf << result.opt_candidate_path[i].sp_dist
                            << (i == N - 1 ? "" : ",");
                    }
                }
            }
            if (config_.write_pgeom) {
                buf << ";";
                if (!result.opt_candidate_path.empty()) {
                    int N = result.opt_candidate_path.size();
                    FMM::LineString pline;
                    for (int i = 0; i < N; ++i) {
                        const FMM::Point& point = result.opt_candidate_path[i].c.point;
                        pline.add_point(point);
                    }
                    buf << pline;
                }
            }
            // Write fields related with cpath
            if (config_.write_cpath) {
                buf << ";" << result.cpath;
            }
            if (config_.write_tpath) {
                buf << ";";
                if (!result.cpath.empty()) {
                    // Iterate through consecutive indexes and write the traversed path
                    int J = result.indices.size();
                    for (int j = 0; j < J - 1; ++j) {
                        int a = result.indices[j];
                        int b = result.indices[j + 1];
                        for (int i = a; i < b; ++i) {
                            buf << result.cpath[i];
                            buf << ",";
                        }
                        buf << result.cpath[b];
                        if (j < J - 2) {
                            // Last element should not have a bar
                            buf << "|";
                        }
                    }
                }
            }
            if (config_.write_mgeom) {                
                buf << ";" << result.mgeom;
            }
            if (config_.write_ep) {
                buf << ";";
                if (!result.opt_candidate_path.empty()) {
                    int N = result.opt_candidate_path.size();
                    for (int i = 0; i < N - 1; ++i) {
                        buf << result.opt_candidate_path[i].ep << ",";
                    }
                    buf << result.opt_candidate_path[N - 1].ep;
                }
            }
            if (config_.write_tp) {
                buf << ";";
                if (!result.opt_candidate_path.empty()) {
                    int N = result.opt_candidate_path.size();
                    for (int i = 0; i < N - 1; ++i) {
                        buf << result.opt_candidate_path[i].tp << ",";
                    }
                    buf << result.opt_candidate_path[N - 1].tp;
                }
            }
            if (config_.write_length) {
                buf << ";";
                if (!result.opt_candidate_path.empty()) {
                    int N = result.opt_candidate_path.size();
                    SPDLOG_TRACE("Write length for {} edges", N);
                    for (int i = 0; i < N - 1; ++i) {
                        // SPDLOG_TRACE("Write length {}",i);
                        buf << result.opt_candidate_path[i].c.edge->cost << ",";
                    }
                    // SPDLOG_TRACE("Write length {}",N-1);
                    buf << result.opt_candidate_path[N - 1].c.edge->cost;
                }
            }
            if (config_.write_duration) {
                buf << ";";
                if (!traj.timestamps.empty()) {
                    int N = traj.timestamps.size();
                    SPDLOG_TRACE("Write duration for {} points", N);
                    for (int i = 1; i < N; ++i) {
                        // SPDLOG_TRACE("Write length {}",i);
                        buf << traj.timestamps[i] - traj.timestamps[i - 1]
                            << (i == N - 1 ? "" : ",");
                    }
                }
            }
            if (config_.write_speed) {
                buf << ";";
                if (!result.opt_candidate_path.empty() && !traj.timestamps.empty()) {
                    int N = traj.timestamps.size();
                    for (int i = 1; i < N; ++i) {
                        double duration = traj.timestamps[i] - traj.timestamps[i - 1];
                        buf << (duration > 0 ? (result.opt_candidate_path[i].sp_dist / duration) : 0)
                            << (i == N - 1 ? "" : ",");
                    }
                }
            }
            buf << '\n';
            // Ensure that fstream is called corrected in OpenMP
#pragma omp critical
            m_fstream << buf.rdbuf();
        }

    private:
        std::ofstream m_fstream;
        const FMM::OutputConfig& config_;
    }; // CSVMatchResultWriter


} //FMM
#endif // FMM_MM_WRITER_HPP