cmake_minimum_required(VERSION 3.8.2)
project(M3C2EP)
find_package(PDAL 1.0.0 REQUIRED CONFIG)
find_package(Open3D REQUIRED CONFIG)
add_compile_definitions(NOMINMAX)
include_directories(${PDAL_INCLUDE_DIRS} ${Open3D_INCLUDE_DIRS} ${Armadillo_INCLUDE_DIRS})
link_directories(${PDAL_LIBRARY_DIRS} ${Open3D_LIBRARY_DIRS})

add_executable(M3C2EP Main.cpp M3C2EP.cpp M3C2EP.hpp M3C2.cpp M3C2.hpp ../shared/FileLAS.cpp ../shared/FileLAS.hpp ../shared/PointCloudAI.cpp ../shared/PointCloudAI.hpp)
target_link_libraries(M3C2EP ${PDAL_LIBRARIES} ${Open3D_LIBRARIES} ${Armadillo_LIBRARY_DIRS}/libopenblas.lib) 

