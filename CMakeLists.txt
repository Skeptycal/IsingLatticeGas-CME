cmake_minimum_required (VERSION 2.8.8)

project (CME_Model)

set (CME_MDL_SRC
     ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
	 ${CMAKE_CURRENT_SOURCE_DIR}/src/Lattice.cpp)

set (CME_MDL_HDR
     ${CMAKE_CURRENT_SOURCE_DIR}/include/Lattice.h)
	 
add_executable(CmeModel ${CME_MDL_SRC} ${CME_MDL_HDR})
