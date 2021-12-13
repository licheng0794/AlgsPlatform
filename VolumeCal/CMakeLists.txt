cmake_minimum_required(VERSION 3.8.2)
project(VolumeCal)
find_package(PDAL 1.0.0 REQUIRED CONFIG)
find_package(Open3D REQUIRED CONFIG)
add_compile_definitions(NOMINMAX)
include_directories(${PDAL_INCLUDE_DIRS} ${Open3D_INCLUDE_DIRS})
link_directories(${PDAL_LIBRARY_DIRS} ${Open3D_LIBRARY_DIRS})

add_executable(VolumeCal Main.cpp VolumeCal.cpp VolumeCal.hpp ../shared/FileLAS.cpp ../shared/FileLAS.hpp ../shared/PointCloudAI.cpp ../shared/PointCloudAI.hpp)
target_link_libraries(VolumeCal ${PDAL_LIBRARIES} ${Open3D_LIBRARIES}) 



