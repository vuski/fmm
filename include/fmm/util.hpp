/**
 * Fast map matching.
 *
 * Utility functions
 *
 * @author: Can Yang
 * @version: 2017.11.11
 
  * @revision : vuski
 * @version: 2021.11.29
 */

#define _CRT_SECURE_NO_WARNINGS


#ifndef FMM_UTIL_HPP
#define FMM_UTIL_HPP


#include <iomanip>
#include <cfloat>
#include <iostream>
#include <string>
#include <sstream>
#include <sys/stat.h>
#include <cmath>
#include <chrono>
#include <vector>
#include <ctime>


#include "spdlog/spdlog.h"

namespace FMM {

    /**
* Log level strings for printing the log level information
*/
    static const std::vector<std::string>
        LOG_LEVESLS{ "0-trace","1-debug","2-info",
                     "3-warn","4-err","5-critical","6-off" };
/**
 * Utility functions
 */

 /**
* Time point
*/
    typedef std::chrono::system_clock::time_point TimePoint;


/**
 * Get the current timestamp
 */
TimePoint get_current_time() {
    return std::chrono::system_clock::now();
};

/**
 * Print a timestamp
 * @param timestamp timestamp point
 */
void print_time(
    const TimePoint& timestamp) {
    std::time_t timestamp_t = std::chrono::system_clock::to_time_t(timestamp);

    struct tm* t;
    localtime_s(t, &timestamp_t);
    char s[18];
    sprintf_s(s, "%04d%02d%02d_%02d%02d%02d",
        t->tm_year + 1900, t->tm_mon + 1, t->tm_mday,
        t->tm_hour, t->tm_min, t->tm_sec);

    std::cout << "Time " << s << '\n';
};

/**
 * Get the duration between two timestamps
 * @param  the starting timestamp
 * @param  the ending timestamp
 * @return duration in seconds
 */
double get_duration(const TimePoint& t1, const TimePoint& t2) {
    return std::chrono::duration_cast<
        std::chrono::milliseconds>(t2 - t1).count() / 1000.;
};

/**
 * Check if file exist or not
 * @param  filename file name
 * @return true if file exists
 */
bool file_exists(const char* filename) {
    struct stat buf;
    if (stat(filename, &buf) != -1) {
        return true;
    }
    return false;
}
/**
 * Check if file exist or not
 * @param  filename file name
 * @return true if file exists
 */
bool file_exists(const std::string& filename) {
    return file_exists(filename.c_str());
}
/**
 * Check if folder exists or not
 * @param  folder_name folder name
 * @return true if folder exists
 */
bool folder_exist(const std::string& folder_name) {
    if (folder_name.empty()) return true;
    struct stat sb;
    if (stat(folder_name.c_str(), &sb) == 0 && (((sb.st_mode) & _S_IFDIR) == _S_IFDIR)) {
        return true;
    }
    return false;
}
/**
 * Get folder path from file path
 * @param  fn File path
 * @return Folder path to a file
 */
std::string get_file_directory(const std::string& fn) {
    std::size_t found = fn.find_last_of("/");
    if (found != std::string::npos) {
        return fn.substr(0, found);
    }
    return {};
};

/**
 * Convert string to bool
 * @param  str input string
 * @return true if str equals "true", "t" or "1";
 */
bool string2bool(const std::string& str) {
    return str == "true" || str == "t" || str == "1";
}

/**
 * Convert bool to string
 * @param  value bool value
 * @return  "true" if value is true
 */
inline std::string bool2string(bool value) {
  return (value ? "true" : "false");
}


/**
 * Split a string containing string separated by , into a vector of string
 * @param str input string
 * @return a vector of strings
 */
std::vector<std::string> split_string(const std::string& str) {
    char delim = ',';
    std::vector<std::string> result;
    std::stringstream ss(str);
    std::string intermediate;
    while (getline(ss, intermediate, delim)) {
        result.push_back(intermediate);
    }
    return result;
}


/**
 * Check if the filename has an extension in the list
 * @param filename
 * @param extension_list_str a list of string separated by ,
 * @return true if file extension is in the list
 */
bool check_file_extension(const std::string& filename,
    const std::string& extension_list_str) {
    bool result = false;
    std::stringstream ss;
    std::string fn_extension = filename.substr(
        filename.find_last_of('.') + 1);
    std::vector<std::string> extensions =
        split_string(extension_list_str);
    for (const auto& extension : extensions) {
        if (fn_extension == extension)
            result = true;
    }
    return result;
}

/**
 * Convert a vector of type into a string with delimiter of ,
 * @param  vec input vector
 * @return  a string of values in the vector delimited by ,
 */
template<typename T>
std::string vec2string(
  const std::vector<T> &vec) {
  std::ostringstream vts;
  if (!vec.empty()) {
    std::copy(vec.begin(), vec.end() - 1,
              std::ostream_iterator<T>(vts, ","));
    vts << vec.back();
  }
  return vts.str();
}

/**
 * Convert string to vector
 * @param  str a string containing a list of values separated by ,
 * @return  a vector of values
 */
template<typename T>
std::vector<T> string2vec(
  const std::string &str) {
  std::vector<T> vec;
  std::stringstream ss(str);
  T i;
  while (ss >> i) {
    vec.push_back(i);
    if (ss.peek() == ',')
      ss.ignore();
  }
  return vec;
}

/**
 * Get line correctly by handling `\r\n` line endings
 */
std::istream& safe_get_line(std::istream& is, std::string& t, char delim)
{
    t.clear();
    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();
    for (;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if (sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case std::streambuf::traits_type::eof():
            // Also handle the case when the last line has no line ending
            if (t.empty())
                is.setstate(std::ios::eofbit);
            return is;
        default:
            if (c == delim) {
                return is;
            }
            t += (char)c;
        }
    }
}

} // FMM
#endif /* FMM_UTIL_HPP */