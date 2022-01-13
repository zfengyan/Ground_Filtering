/********************************************************************
* GEO1015.2021
* hw03
* Yitong Xia
* 5445825
* Fengyan Zhang
* 5462150
*
* The idea of CSF is from this paper:
* http://www.mdpi.com/2072-4292/8/6/501/htm
* The basic architectures and ideas are inspired by this article (its related code is open source).
*
*********************************************************************/


#ifndef _GROUND_FILTER_H_
#define _GROUND_FILTER_H_

#include <string>
#include <vector>

// -- json parsing
#include <json.hpp>

// -- CGAL kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

using json = nlohmann::json;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;

std::vector<Point> read_lasfile(const json& jparams);
void write_lasfile(const std::string filename, const std::vector<Point>& pointcloud, const std::vector<int>& class_labels);

void groundfilter_tin(const std::vector<Point>& pointcloud, const json& jparams);
void groundfilter_csf(const std::vector<Point>& pointcloud, const json& jparams);


#endif