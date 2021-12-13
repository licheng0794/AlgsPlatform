cmake_minimum_required(VERSION 3.8.2)
project(CSF)
find_package(PDAL 1.0.0 REQUIRED CONFIG)
find_package(Open3D REQUIRED CONFIG)
add_compile_definitions(NOMINMAX)
include_directories(${PDAL_INCLUDE_DIRS} ${Open3D_INCLUDE_DIRS})
link_directories(${PDAL_LIBRARY_DIRS} ${Open3D_LIBRARY_DIRS})

set(CSF_SOURCES
    c2cdist.cpp
    Cloth.cpp
    Constraint.cpp
    CSF.cpp
    Particle.cpp
    point_cloud.cpp
    Rasterization.cpp
    XYZReader.cpp
    CSFMain.cpp
)

set(CSF_HEADERS
    c2cdist.h
    Cloth.h
    Constraint.h
    CSF.h
    Particle.h
    point_cloud.h
    Rasterization.h
    Vec3.h
    XYZReader.h
)

add_executable(CSF ${CSF_SOURCES} ${CSF_HEADERS} ../shared/FileLAS.cpp ../shared/FileLAS.hpp ../shared/PointCloudAI.cpp ../shared/PointCloudAI.hpp)
target_link_libraries(CSF ${PDAL_LIBRARIES} ${Open3D_LIBRARIES}) 



