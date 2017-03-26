# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(kTracker CXX C Fortran)

# Put executables into a bin
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_MACOSX_RPATH 1)  # only needed to suppress cmake warning

if(APPLE)
  set(CMAKE_Fortran_COMPILER gfortran)
endif()

if(UNIX AND NOT APPLE)
  set(CMAKE_Fortran_COMPILER g++)
endif()

# Command line options
# Additional compile targets
set(ADDTIONAL "" CACHE STRING "Additional compilation targets to build")
set(ADDTGTLIST "")
if(NOT("${ADDTIONAL}" STREQUAL ""))
  string(REPLACE "," ";" ADDTGTLIST ${ADDTIONAL})
endif()

# Add version number from git if exists
if(EXISTS "${PROJECT_SOURCE_DIR}/.git")
  execute_process(COMMAND git rev-parse --abbrev-ref HEAD OUTPUT_VARIABLE GIT_BRANCH OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND git rev-parse --short HEAD      OUTPUT_VARIABLE GIT_HASH   OUTPUT_STRIP_TRAILING_WHITESPACE)
  add_definitions("-DGIT_VERSION=\"${GIT_BRANCH}-${GIT_HASH}\"")
else()
  add_definitions("-DGIT_VERSION=\"1.7.0\"")
endif()

# ROOT dict generation
add_custom_command (
  OUTPUT ${PROJECT_SOURCE_DIR}/src/FastTrackletDict.cxx ${PROJECT_SOURCE_DIR}/src/FastTrackletDict.h
  COMMAND rootcint
  ARGS -f ${PROJECT_SOURCE_DIR}/src/FastTrackletDict.cxx -c ${PROJECT_SOURCE_DIR}/inc/FastTracklet.h ${PROJECT_SOURCE_DIR}/inc/FastTrackletLinkDef.h
  )

add_custom_command (
  OUTPUT ${PROJECT_SOURCE_DIR}/src/SRawEventDict.cxx ${PROJECT_SOURCE_DIR}/src/SRawEventDict.h
  COMMAND rootcint
  ARGS -f ${PROJECT_SOURCE_DIR}/src/SRawEventDict.cxx -c ${PROJECT_SOURCE_DIR}/inc/SRawEvent.h ${PROJECT_SOURCE_DIR}/inc/SRawEventLinkDef.h
  )

add_custom_command (
  OUTPUT ${PROJECT_SOURCE_DIR}/src/SRecEventDict.cxx ${PROJECT_SOURCE_DIR}/src/SRecEventDict.h
  COMMAND rootcint
  ARGS -f ${PROJECT_SOURCE_DIR}/src/SRecEventDict.cxx -c ${PROJECT_SOURCE_DIR}/inc/SRecEvent.h ${PROJECT_SOURCE_DIR}/inc/SRecEventLinkDef.h
  )

add_custom_command (
  OUTPUT ${PROJECT_SOURCE_DIR}/src/TriggerRoadDict.cxx ${PROJECT_SOURCE_DIR}/src/TriggerRoadDict.h
  COMMAND rootcint
  ARGS -f ${PROJECT_SOURCE_DIR}/src/TriggerRoadDict.cxx -c ${PROJECT_SOURCE_DIR}/inc/TriggerRoad.h ${PROJECT_SOURCE_DIR}/inc/TriggerRoadLinkDef.h
  )

add_custom_command (
  OUTPUT ${PROJECT_SOURCE_DIR}/src/SMillepedeUtilDict.cxx ${PROJECT_SOURCE_DIR}/src/SMillepedeUtilDict.h
  COMMAND rootcint
  ARGS -f ${PROJECT_SOURCE_DIR}/src/SMillepedeUtilDict.cxx -c ${PROJECT_SOURCE_DIR}/inc/SMillepedeUtil.h ${PROJECT_SOURCE_DIR}/inc/SMillepedeUtilLinkDef.h
  )

# DPSim source code
include_directories(${PROJECT_SOURCE_DIR}/inc)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cxx   ${PROJECT_SOURCE_DIR}/src/*.cc     ${PROJECT_SOURCE_DIR}/src/*.f)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h ${PROJECT_SOURCE_DIR}/include/*.hh)

# ROOT
find_program(ROOTCONF "root-config")
if(ROOTCONF)
  message("-- Detecting ROOT:    found at ${ROOTCONF}")
else()
  message(FATAL_ERROR "-- Detecting ROOT:    not found")
endif()
execute_process(COMMAND root-config --prefix OUTPUT_VARIABLE ROOT_PREFIX  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CFLAGS  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND root-config --libs   OUTPUT_VARIABLE ROOT_LINK    OUTPUT_STRIP_TRAILING_WHITESPACE)

# the rest packages are not needed for libonly mode
# Geant4
find_program(G4CONF "geant4-config")
if(G4CONF)
  message("-- Detecting Geant4:  found at ${G4CONF}")
else()
  message(FATAL_ERROR "-- Detecting Geant4:  not found")
endif()
find_package(Geant4 REQUIRED gdml vis_all)
execute_process(COMMAND geant4-config --cflags OUTPUT_VARIABLE GEANT4_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND geant4-config --libs   OUTPUT_VARIABLE GEANT4_LINK   OUTPUT_STRIP_TRAILING_WHITESPACE)

# MySQL
find_program(MYSQLCONF "mysql_config")
if(MYSQLCONF)
  message("-- Detecting MySQL:   found at ${MYSQLCONF}")
else()
  message(FATAL_ERROR "-- Detecting MySQL:   not found")
endif()
execute_process(COMMAND mysql_config --cflags OUTPUT_VARIABLE MYSQL_CFLAGS  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND mysql_config --libs   OUTPUT_VARIABLE MYSQL_LINK    OUTPUT_STRIP_TRAILING_WHITESPACE)

# Warning suppression
string(REGEX REPLACE "-W[a-z-]*" "" GEANT4_CFLAGS ${GEANT4_CFLAGS})
string(REPLACE "-pedantic" "" GEANT4_CFLAGS ${GEANT4_CFLAGS})

#find BOOST
find_package(Boost)
if(Boost_FOUND)
  message("-- Detecting BOOST:   found at ${Boost_INCLUDE_DIRS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${Boost_INCLUDE_DIRS}")
elseif(EXISTS $ENV{BOOST_ROOT})
  message("-- Detecting BOOST:   using Boost installation at $ENV{BOOST_ROOT}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I$ENV{BOOST_ROOT}/include")
else()
  message(FATAL_ERROR "-- Detecting BOOST:   not found")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GEANT4_CFLAGS} ${ROOT_CFLAGS} ${MYSQL_CFLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GEANT4_LINK} ${ROOT_LINK} ${MYSQL_LINK}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${GEANT4_LINK} ${ROOT_LINK} ${MYSQL_LINK}" )

set(TARGETLIST "kFastTracking;kVertex;milleAlign;sqlDataReader;sqlMCReader;sqlResWriter")
foreach(TGT ${TARGETLIST})
  message(STATUS "Adding targets to the compile list: ${TGT}")

  # add executables
  add_executable(${TGT} ${TGT}.cxx ${sources} ${headers} ${PROJECT_SOURCE_DIR}/src/FastTrackletDict.cxx ${PROJECT_SOURCE_DIR}/src/SRawEventDict.cxx ${PROJECT_SOURCE_DIR}/src/SRecEventDict.cxx ${PROJECT_SOURCE_DIR}/src/SMillepedeUtilDict.cxx ${PROJECT_SOURCE_DIR}/src/TriggerRoadDict.cxx)
  #target_link_libraries(${TGTNAME} ${ROOT_LINK} ${Geant4_LIBRARIES})
endforeach(TGT)
add_library(kTracker SHARED ${sources} ${headers} ${PROJECT_SOURCE_DIR}/src/FastTrackletDict.cxx ${PROJECT_SOURCE_DIR}/src/SRawEventDict.cxx ${PROJECT_SOURCE_DIR}/src/SRecEventDict.cxx ${PROJECT_SOURCE_DIR}/src/SMillepedeUtilDict.cxx ${PROJECT_SOURCE_DIR}/src/TriggerRoadDict.cxx)

# add additional targets
foreach(TGT ${ADDTGTLIST})
  message(STATUS "Adding additional targets to the compile list: ${TGT}")

  # add executables
  add_executable(${TGT} ${TGT}.cxx ${sources} ${headers} ${PROJECT_SOURCE_DIR}/src/FastTrackletDict.cxx ${PROJECT_SOURCE_DIR}/src/SRawEventDict.cxx ${PROJECT_SOURCE_DIR}/src/SRecEventDict.cxx ${PROJECT_SOURCE_DIR}/src/SMillepedeUtilDict.cxx ${PROJECT_SOURCE_DIR}/src/TriggerRoadDict.cxx)
endforeach(TGT)
