cmake_minimum_required(VERSION 3.8.2) # 3.8.0 is minimum version that allow Linux remote building and debugging
project(AIalgorithms)

set(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
add_subdirectory(M3C2)
add_subdirectory(M3C2EP)
add_subdirectory(CSF)
add_subdirectory(VolumeCal)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

install(TARGETS M3C2 M3C2EP CSF VolumeCal RUNTIME DESTINATION bin)

file(GLOB open3d "open3d_dlls/*.dll")
file(COPY ${open3d} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

file(GLOB pdal "PDAL_DLLS/*.dll")
file(COPY ${pdal} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

file(GLOB Armadillo "armadillo_dlls/*.dll")
file(COPY ${Armadillo} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

file(GLOB winfiles "windows_dlls/*.dll")
file(COPY ${winfiles} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

file(COPY "Readme.txt" DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)