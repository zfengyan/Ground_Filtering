/*
  GEO1015.2021
  hw03
  --
  Ravi Peters
  r.y.peters@tudelft.nl
*/

/*
!!! DO NOT MODIFY this file !!!
*/

#include <iostream>
#include <fstream>
#include <chrono>
#include <filesystem>
namespace fs = std::filesystem;

//-- our GroundFilter functions
#include "GroundFilter.h"

int main(int argc, char** argv)
{
  //-- setup path to params.json file
  std::string json_path = JSON_PARAMS_PATH;
  //-- take the json path from command line if supplied so
  if (argc == 2) {
    json_path = argv[1];
  }
  
  //-- set current working directory to the parent directory of json_path. As a result the
  //-- filepaths in the json file will be read relative to the location of the json file
  fs::path working_directory = fs::path(json_path).parent_path();
  fs::current_path(working_directory);
  std::cout << "Active working directory: " << working_directory << std::endl;

  //-- read the params.json file
  std::ifstream json_file(json_path);
  if (!json_file) {
    std::cerr << "JSON file " << json_path << " not found.\n";
    return 1;
  }
  json j; json_file >> j;

  //-- read pointcloud from input file
  std::vector<Point> pointcloud = read_lasfile(j["input"]);

  //-- groundfilter if in the params  
  std::cout << std::fixed << std::setprecision(3);

  if (j.contains("groundfilter_tin")) {
    std::cout << "=== TIN refinement groundfilter ===\n";

    auto start = std::chrono::high_resolution_clock::now();
    groundfilter_tin(pointcloud, j["groundfilter_tin"]);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "-->" << elapsed.count() << "s" << std::endl;
  }
  
  if (j.contains("groundfilter_csf")) {
    std::cout << "=== CSF groundfilter ===\n";

    auto start = std::chrono::high_resolution_clock::now();
    groundfilter_csf(pointcloud, j["groundfilter_csf"]);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "-->" << elapsed.count() << "s" << std::endl;
  }

  //-- we're done, return 0 to say all went fine
  return 0;
}
