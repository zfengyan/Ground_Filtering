/*
  GEO1015.2021
  hw03 
  --
  [YOUR NAME] 
  [YOUR STUDENT NUMBER] 
  [YOUR NAME] 
  [YOUR STUDENT NUMBER] 
*/

#include <string>
#include <vector>

// -- json parsing
#include <json.hpp>

// -- CGAL kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

using json = nlohmann::json;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
//typedef Kernel::Point_3 Point;


struct Point {
	double x;
	double y;
	double z;

	Point():x(0),y(0),z(0){}
	Point(const double& p_x, const double& p_y, const double& p_z):
		x(p_x),y(p_y),z(p_z){}

};

std::vector<Point> read_lasfile(const json& jparams);
void write_lasfile(const std::string filename, const std::vector<Point>& pointcloud, const std::vector<int>& class_labels);

void groundfilter_tin(const std::vector<Point>& pointcloud, const json& jparams);
void groundfilter_csf(const std::vector<Point>& pointcloud, const json& jparams);