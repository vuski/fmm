/**
 * Fast map matching.
 *
 * Definition of geometry
 * 2020-11-25 Remove linestring2ogr
 *
 * @author: Can Yang
 * @version: 2017.11.11
 * 
 * @revision : vuski
 * @version: 2021.11.29
 * 
 * //가장 우선순위 파일. 모든 파일에 등장하는 타입의 기초
 */

#ifndef FMM_GEOMTYPES_HPP
#define FMM_GEOMTYPES_HPP

#include <string>
#include <sstream>
#include <vector>
#include <iomanip> 

namespace FMM {
/**
 * Core data types
 */


/**
 *  Point class
 */
//typedef boost::geometry::model::point<double, 2,
//    boost::geometry::cs::cartesian> Point; // Point for rtree box

struct Point {

    double x;
    double y;

    Point() {
        x = 0.0;
        y = 0.0;
    }

    Point(double x_, double y_) {
        x = x_;
        y = y_;
    }
};

                                           
                                           
                                           /**
 *  Linestring geometry class
 *
 *  This class wraps a boost linestring geometry.
 */
class LineString {
public:
  /**
   * This is the boost geometry linestring class, stored inside the
   * LineString class.
   */
  typedef std::vector<Point> linestring_t;
  /**
   * Get the x coordinate of i-th point in the line
   * @param i point index
   * @return x coordinate
   */
  inline double get_x(int i) const{
    return line[i].x;
  };
  /**
   * Get the y coordinate of i-th point in the line
   * @param i point index starting from 0 to N-1, where N is the
   * number of points in the line
   * @return y coordinate
   */
  inline double get_y(int i) const{
    return line[i].y;
  };
  /**
   * Set x coordinate of i-th point in the line
   * @param i point index
   * @param v the value to update the old coordinate
   */
  inline void set_x(int i, double v){
    line[i].x = v;
  };
  /**
   * Set y coordinate of i-th point in the line
   * @param i point index
   * @param v the value to update the old coordinate
   */
  inline void set_y(int i, double v){
     line[i].y = v;
  };
  /**
   * Add a point to the end of the current line
   * @param x x coordinate of the point to add
   * @param y y coordinate of the point to add
   */
  inline void add_point(double x,double y){
    line.push_back(Point(x,y));
  };
  /**
   * Add a point to the end of the current line
   * @param point the point to be added
   */
  inline void add_point(const Point& point){
      line.push_back(point);
  };
  /**
   * Get the i-th point in the line
   * @param  i point index starting from 0 to N-1
   * @return The i-th point of the line.
   *
   * Note that the point is a copy of the original point.
   * Manipulating the returned point will not change the
   * original line.
   */
  inline Point get_point(int i) const{
    return line[i];
  };
  /**
   * Get a constance reference of the i-th point in the line
   * @param  i point index
   * @return  A constant reference to the ith point of line, which
   * avoids create a new point.
   */
  inline const Point &at(int i) const{
    return line[i];
  }
  /**
   * Get the number of points in a line
   * @return the point number
   */
  inline int get_num_points() const{
    return (int)line.size();
  };
  /**
   * Check if the line is empty or not
   * @return true if the line is empty, otherwise false
   */
  inline bool is_empty(){
    return line.empty();
  };
  /**
   * Remove all points in the current line.
   */
  inline void clear(){
    line.clear();
  };
  /**
   * Get the length of the line
   * @return the length value
   */
  inline double get_length() const {

      double length = 0.0;
      for (int i = 0; i < line.size() - 1; i++)
      {
          length += std::sqrt((line[i].x - line[i+1].x) * (line[i].x - line[i + 1].x)
              + (line[i].y - line[i + 1].y) * (line[i].y - line[i + 1].y));          
      }
    return length;
  };

  /**
   * Export a string containing GeoJSON representation of the line.
   * @return The GeoJSON of the line
   */
  inline std::string export_json() const{
    std::ostringstream ss;
    int N = get_num_points();
    if (N>0){
      ss << "{\"type\":\"LineString\",\"coordinates\": [";
      for (int i=0;i<N;++i){
        ss << "[" << get_x(i) << "," << get_y(i) <<"]"
           << (i==N-1 ? "": ",");
      }
      ss << "]}";
    }
    return ss.str();
  };
  /**
   * Get a const reference to the inner boost geometry linestring
   * @return const reference to the inner boost geometry linestring
   */
  inline const linestring_t &get_geometry_const() const{
    return line;
  };
  /**
   * Get a reference to the inner boost geometry linestring
   * @return a reference to the inner boost geometry linestring
   */
  linestring_t &get_geometry(){
    return line;
  };

  /**
   * Compare if two linestring are the same.
   *
   * It the two lines overlap with each other within a threshold of 1e-6,
   * they are considered equal. This is used in the test class.
   *
   * @param rhs the linestring
   * @return true if the two lines are equal.
   */
  inline bool operator==(const LineString& rhs) const {
    int N = get_num_points();
    if (rhs.get_num_points()!=N)
      return false;
    bool result = true;
    for (int i=0;i<N;++i){

        Point a = get_point(i);
        Point b = rhs.get_point(i);
        double dist = sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) + (a.y - b.y));
      if (dist >1e-6) result = false;
    }
    return result;
  };
  /**
   * Overwrite the operator of << of linestring
   * @param os an input stream
   * @param rhs a linestring object
   * @return the wkt representation of the line will be written to the stream.
   */
  friend std::ostream& operator<<(std::ostream& os, const LineString& rhs)
  {     
      os << std::setprecision(12); //boost::geometry::wkt(rhs.line);
      os << "LINESTRING(";
      int N = rhs.get_num_points();
      for (int i = 0; i < N - 1; i++)
      {
          os << rhs.get_x(i) << " " << rhs.get_y(i) << ",";
      }
      os << rhs.get_x(N - 1) << " " << rhs.get_y(N - 1) << ")";

      return os;
};

private:
  linestring_t line;  
}; // LineString

  std::ostream& operator<<(std::ostream& os, const FMM::LineString& rhs);


}; // FMM


#endif // FMM_GEOMTYPES_HPP