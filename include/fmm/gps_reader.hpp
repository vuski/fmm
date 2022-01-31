/**
 * Fast map matching.
 *
 * Definition of GPS reader classes
 *
 * @author: Can Yang
 * @version: 2017.11.11
 * 
 *  * @revision : vuski
 * @version: 2021.11.29
 */



#ifndef FMM_GPS_READER_HPP
#define FMM_GPS_READER_HPP

#define _CRT_SECURE_NO_WARNINGS


#include <iostream>
#include <fstream>
#include <string>

#include "fmm/type.hpp"
#include "csv.h"
#include "rapidjson/document.h"

namespace FMM {
/**
 * Classes related with input and output
 */


/**
 * Trajectory Reader Interface.
 */


 /**
    * %Trajectory class
    *
    * A GPS trajectory represented with id, geometry and timestamps
    */
struct Trajectory {
    Trajectory() {};
    Trajectory(int id_arg,
        const LineString& geom_arg) :id(id_arg), geom(geom_arg) {};
    Trajectory(int id_arg,
        const LineString& geom_arg, const std::vector<double>& timestamps_arg)
        :id(id_arg), geom(geom_arg), timestamps(timestamps_arg) {};
    int id; /**< Id of the trajectory */
    LineString geom; /**< Geometry of the trajectory */
    std::vector<double> timestamps; /**< Timestamps of the trajectory */
};

struct GPS { //24
    int id;
    double x;
    double y;
    int time;
};

using std::string;
using std::cout;
using std::endl;
using std::vector;
using std::ifstream; 


class SignalReader {

public :

    

    static enum RAWTYPE { CSV, TSV, GEOJSON, GPX, DAT_SIGNAL, DAT_GPS };
    static enum TIMEFORMAT { UNIXTIME, YMDHMS };

    SignalReader()
    {

    }

    //텍스트 형식을 dat로 변환한다.
    static void convertAndreadData(std::string fileNameData,
        std::vector<Trajectory>& trajectories, int rawtype, int timeformat)
    {
        if (rawtype == CSV || rawtype == TSV)
        {
            clock_t time = clock();
            cout << fileNameData << "변환시작" << endl;

            io::CSVReader<5, io::trim_chars<' '>, io::no_quote_escape<'\t'>> table(fileNameData);
            try {
                table.read_header(io::ignore_extra_column, "serial", "id", "x", "y", "time");
            }
            catch (int e)
            {
                table.set_header("id0", "id", "x", "y", "time");
            }
            size_t count = 0;

            bool isFirst = true;
            char* timestr = (char*)malloc(sizeof(char) * 21);
            size_t idPre = -1;
            size_t id;
            double x;
            double y;
            int time0 = 0; 
			int serial;
			int idCount = 0;
            FMM::LineString geom;
            std::vector<double> timestamps;

            while (timeformat == UNIXTIME ? 
                table.read_row(serial, id, x, y, time0) :
                table.read_row(serial, id, x, y, timestr))
            {
				id += serial * 100ULL;

                if (id != idPre)
                {
                    if (isFirst)
                    {
                        isFirst = false;
                        idPre = id;
                        geom.clear();
                        timestamps.clear();
                    }
                    else {
                        trajectories.push_back({ idCount, geom, timestamps });
						idCount++;

                        idPre = id;
                        geom.clear();
                        timestamps.clear();
                    }
                }                

                geom.add_point(x, y);

                time0 = timeformat == UNIXTIME ?
                    time0 : (int)getTimeT(timestr);
                timestamps.push_back(time0); //체류시작시각까지 기록

                if (count % 1000000 == 0) cout << ".";
                count++;
            }
            trajectories.push_back({ idCount, geom, timestamps }); //마지막 추가

            cout << fileNameData << "입력 완료" << clock() - time << "ms" << endl;  

        }
        else if (rawtype == GEOJSON) 
        {
            clock_t time = clock();

            using std::cout;
            using std::endl;

            rapidjson::Document doc;

            {
                using namespace rapidjson;
                char* buffer;
                {
                    clock_t t2 = clock();
                    int count;
                    FILE* fp = NULL;
                    fopen_s(&fp, fileNameData.c_str(), "r");

                    //read size of file
                    fseek(fp, 0, SEEK_END);
                    int64_t size = ftell(fp);

                    buffer = (char*)malloc(size + 1);
                    memset(buffer, 0, size + 1);

                    fseek(fp, 0, SEEK_SET);
                    count = fread(buffer, size, 1, fp);
                    //처리	
                    fclose(fp);
                    clock_t t3 = clock();
                    std::cout << " 읽기시간 : " << t3 - t2 << std::endl;
                }

                std::cout << "문자열 길이:" << strlen(buffer) << std::endl;
                doc.Parse(buffer);

                std::cout << "parsed" << std::endl;
                free(buffer);
            }


            using namespace rapidjson;
            const Value& features = doc["features"];
            cout << "geojson 멤버 수:" << features.Size() << endl;
            unsigned int featureSize = features.Size();


            std::unordered_map<size_t, Edge> keymap;

            for (unsigned int i = 0; i < featureSize; i++)//도형 하나하나 루프
            {
                const Value& feature = features[i];
                const Value& properties = feature["properties"];
                std::string linetype = feature["geometry"]["type"].GetString();

                const Value& linestring = linetype == "MultiLineString" ?
                    feature["geometry"]["coordinates"].GetArray()[0] :
                    feature["geometry"]["coordinates"];//폴리곤들의 배열    


                //속성을 읽는다.
                int id = properties["id"].GetInt();


                //좌표값을 읽는다.
                LineString geom;
                unsigned int coordSize = linestring.Size();

                if (coordSize > 1) { //점이 하나로 이루어져 있을 때는 버린다.                
                    for (unsigned int k = 0; k < coordSize; k++) {
                        //cout << linestring[k].GetArray()[0].GetDouble() << "\t"
                        //    << linestring[k].GetArray()[1].GetDouble() << endl;
                        geom.add_point(linestring[k].GetArray()[0].GetDouble(),
                            linestring[k].GetArray()[1].GetDouble());
                    } //k                 
                }

                vector<double> timestamps;

                trajectories.push_back({ id, geom, timestamps }); //마지막 추가

            } //i 도형 하나하나

            cout << fileNameData << "입력 완료" << clock() - time << "ms" << endl;
        }
        else if (rawtype == GPX)
        {

           
        }
        

    }




    static void writeOriginalGeojson(const std::string fileName,
        std::vector<FMM::Trajectory> trjVec) {

        using std::cout;
        using std::endl;
        using std::string;

        auto startTime = std::chrono::high_resolution_clock::now();
        cout << "link geojson 파일 쓰기 시작.....";
        std::stringstream header(std::stringstream::out | std::stringstream::binary);
        header << "{\"type\": \"FeatureCollection\"," << endl;
        header << "\"crs\": { \"type\": \"name\", \"properties\": { \"name\": \"urn:ogc:def:crs:EPSG::5179\" } }," << endl;
        header << "\"features\":[" << endl;
        string headerStr = header.str();

        std::stringstream str(std::stringstream::out | std::stringstream::binary);
        str.precision(10);


        auto myfile = std::fstream(fileName, std::ios::out | std::ios::binary);
        myfile.write(header.str().c_str(), header.str().length());

        const int writeLimit = 100000;

        for (int idx = 0 ; (idx < trjVec.size() && idx <= writeLimit ); idx++)
        {
            FMM::Trajectory & trj = trjVec[idx];
            int size = trj.geom.get_num_points();
            //SPDLOG_INFO("id : {}  size: {}", trj.id, size);
            int i = 0;
            str.str("");

            for (; i < size - 1; i++)
            {
                double x0 = trj.geom.get_x(i);
                double y0 = trj.geom.get_y(i);
                double x1 = trj.geom.get_x(i + 1);
                double y1 = trj.geom.get_y(i + 1);


                str << "{ \"type\": \"Feature\", \"properties\": {";
                str << "\"trid\": " << trj.id << ", \"time\":";
                if (trj.timestamps.size() > i + 1) str << (trj.timestamps[i + 1] - trj.timestamps[i]) / 60;
                else str << 0;
                str << "}, \"geometry\": { \"type\": \"LineString\", \"coordinates\": [ ";
                str << "[ " << x0 << ", " << y0 << " ], ";
                str << "[ " << x1 << ", " << y1 << " ]]}}";// , " <<endl;

                if ( (idx == writeLimit ||idx == trjVec.size() -1) && i == size - 2) //마지막줄에는
                {
                    str << endl;
                }
                else {
                    str << "," << endl;
                }

            }

            myfile.write(str.str().c_str(), str.str().length());


        }

        str.str("]}");
        myfile.write(str.str().c_str(), str.str().length());
        myfile.close();


        auto endTime = std::chrono::high_resolution_clock::now();
        cout << "파일 쓰기 소요시간 :" << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << endl;

    }



private:

    inline static time_t getTimeT(const char* ymdhms)
    {
        //sscanf_s(datetime, "%4d-%2d-%2dT%2d:%2d:%2d", &yyyy, &MM, &dd, &hh, &mm, &ss);
        //sscanf_s 는 아래처럼 직접 파싱하는 것에 비해 6배 이상 시간이 걸린다.

        struct tm* t = (struct tm*)malloc(sizeof(struct tm));
        //t->tm_year = 2018 - 1900;
        t->tm_year = (ymdhms[0] * 1000 + ymdhms[1] * 100 + ymdhms[2] * 10 + ymdhms[3] - '0' * 1111) - 1900;
        t->tm_mon = ymdhms[5] * 10 + ymdhms[6] - '0' * 11 - 1;
        t->tm_mday = ymdhms[8] * 10 + ymdhms[9] - '0' * 11;
        t->tm_hour = ymdhms[11] * 10 + ymdhms[12] - '0' * 11;
        t->tm_min = ymdhms[14] * 10 + ymdhms[15] - '0' * 11;
        t->tm_sec = ymdhms[17] * 10 + ymdhms[18] - '0' * 11;
        //t->tm_sec = datetime[17] * 10 + datetime[18] - '0' * 11;
        time_t result = mktime(t);
        free(t);
        return result;
    }

    template <typename DATATYPE>
    static size_t readVectorBuffer(std::string fileName, std::vector<DATATYPE>& data) {

        clock_t t2 = clock();
        FILE* fp = NULL;
        fopen_s(&fp, fileName.c_str(), "rb");
        if (fp == NULL) {
            fclose(fp);

            std::cout << "오픈에러!" << std::endl;
            return 0;
        }

        //read size of file
        fseek(fp, 0, SEEK_END);
        size_t size = _ftelli64(fp);
        fseek(fp, 0, SEEK_SET);
        data.resize(size / sizeof(DATATYPE));
        fread(&data[0], size, 1, fp);
        fclose(fp);
        clock_t t3 = clock();
        std::cout << fileName + " 읽기시간 : " << t3 - t2 << std::endl;

        return size;
    }


    static void getFileNameList(string folder, vector<string>& fileNameList, string extension)
    {

        char search_path[400];
        sprintf(search_path, "%s*.%s", folder.c_str(), extension.c_str());
        WIN32_FIND_DATA fd;
        HANDLE hFind = ::FindFirstFile(search_path, &fd);
        if (hFind != INVALID_HANDLE_VALUE) {
            do {
                // read all (real) files in current folder
                // , delete '!' read other 2 default folder . and ..
                if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
                    fileNameList.push_back(fd.cFileName);
                    //cout << "파일목록:" << fd.cFileName << endl;
                }
            } while (::FindNextFile(hFind, &fd));
            ::FindClose(hFind);
        }
    }

};


} // FMM

#endif // FMM_READER_HPP