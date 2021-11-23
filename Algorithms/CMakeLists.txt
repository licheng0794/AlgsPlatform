cmake_minimum_required(VERSION 3.8.2)
project(M3C2)
find_package(PDAL 1.0.0 REQUIRED CONFIG)
find_package(Open3D REQUIRED CONFIG)
add_compile_definitions(NOMINMAX)
include_directories(${PDAL_INCLUDE_DIRS} ${Open3D_INCLUDE_DIRS} ${Armadillo_INCLUDE_DIRS})
link_directories(${PDAL_LIBRARY_DIRS} ${Open3D_LIBRARY_DIRS})

add_executable(M3C2 Main.cpp M3C2.cpp M3C2.hpp FileLAS.cpp FileLAS.hpp PointCloudAI.cpp PointCloudAI.hpp)
target_link_libraries(M3C2 ${PDAL_LIBRARIES} ${Open3D_LIBRARIES} ${Armadillo_LIBRARY_DIRS}/libopenblas.lib) 


if (MSVC)
  file(GLOB Open3D_DLLS "../open3d_dlls/*.dll")
  add_custom_command(TARGET M3C2
                     POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E copy_if_different
                     ${Open3D_DLLS}
                     $<TARGET_FILE_DIR:M3C2>)

  file(GLOB PDAL_DLLS "../pdal_dlls/*.dll")
  add_custom_command(TARGET M3C2
                     POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E copy_if_different
                     ${PDAL_DLLS}
                     $<TARGET_FILE_DIR:M3C2>)
  file(GLOB Armadillo_DLLS "../armadillo_dlls/*.dll")
  add_custom_command(TARGET M3C2
                     POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E copy_if_different
                     ${Armadillo_DLLS}
                     $<TARGET_FILE_DIR:M3C2>)
endif (MSVC)
