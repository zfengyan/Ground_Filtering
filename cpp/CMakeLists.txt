cmake_minimum_required(VERSION 3.17)
project(hw03)

# get point cloud data
# if(NOT EXISTS ${PROJECT_SOURCE_DIR}/data/ahn3.ply)
#   message(STATUS "Downloading ahn3.ply file...")
#   file(DOWNLOAD http://godzilla.bk.tudelft.nl/tmp/geo1015_2020_hw03/ahn3.ply ${PROJECT_SOURCE_DIR}/data/ahn3.ply
#     TIMEOUT 60)
# endif()

# -- Dependencies
# nlohmann-json
include_directories(${PROJECT_SOURCE_DIR}/external/nlohmann-json)
# cgal
find_package(CGAL QUIET COMPONENTS Core REQUIRED)
# laslib
add_subdirectory(external/LAStools)

add_definitions(
  -DJSON_PARAMS_PATH=\"${PROJECT_SOURCE_DIR}/../data/params.json\"
)

add_executable(hw03 
  src/main.cpp
  src/GroundFilter.cpp
 "src/Vector3d.h" "src/Vector3d.cpp" "src/Particle.h"  "src/Particle.cpp" "src/Cloth.h")
target_link_libraries(hw03 PRIVATE
  CGAL::CGAL CGAL::CGAL_Core
  LASlib
)
set_target_properties(hw03 PROPERTIES 
  CXX_STANDARD 17
)

# Make the executable a default target to build & run in Visual Studio
set_property(DIRECTORY ${PROJECT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT hw03)