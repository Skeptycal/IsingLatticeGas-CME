cmake_minimum_required (VERSION 2.8.8)

project (CME_Model)

set(BOOST_DIR
	${CMAKE_CURRENT_SOURCE_DIR}/boost-1_64_0)

include_directories(${BOOST_DIR}/include)

link_directories(${BOOST_DIR}/lib)

file(COPY ${BOOST_DIR}/lib/boost_filesystem-vc120-mt-gd-1_64.dll 
	 DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/build/Debug/bin)
	 
file(COPY ${BOOST_DIR}/lib/boost_filesystem-vc120-mt-1_64.dll 
	 DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/build/RelWithDebInfo/bin)
	
file(COPY ${BOOST_DIR}/lib/boost_filesystem-vc120-mt-1_64.dll 
	 DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/build/Release/bin)
	 
file(COPY ${BOOST_DIR}/lib/boost_thread-vc120-mt-gd-1_64.dll 
	 DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/build/Debug/bin)
	 
file(COPY ${BOOST_DIR}/lib/boost_thread-vc120-mt-1_64.dll 
	 DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/build/RelWithDebInfo/bin)
	
file(COPY ${BOOST_DIR}/lib/boost_thread-vc120-mt-1_64.dll 
	 DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/build/Release/bin)

set (CME_MDL_SRC
     ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
	 ${CMAKE_CURRENT_SOURCE_DIR}/src/Lattice.cpp
	 ${CMAKE_CURRENT_SOURCE_DIR}/src/Protein.cpp
	 ${CMAKE_CURRENT_SOURCE_DIR}/src/Clump.cpp
	 ${CMAKE_CURRENT_SOURCE_DIR}/src/Export.cpp
	 ${CMAKE_CURRENT_SOURCE_DIR}/src/SizingProcess.cpp)

set (CME_MDL_HDR
     ${CMAKE_CURRENT_SOURCE_DIR}/include/Lattice.h
	 ${CMAKE_CURRENT_SOURCE_DIR}/include/Protein.h
	 ${CMAKE_CURRENT_SOURCE_DIR}/include/Clump.h
	 ${CMAKE_CURRENT_SOURCE_DIR}/include/Export.h
	 ${CMAKE_CURRENT_SOURCE_DIR}/include/SizingProcess.h)


add_executable(CmeModel ${CME_MDL_SRC} ${CME_MDL_HDR})

target_link_libraries(CmeModel boost_filesystem*)
