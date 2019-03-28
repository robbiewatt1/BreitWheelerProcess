# $Id: CMakeLists.txt 86065 2014-11-07 08:51:15Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(BreitWheeler)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

FIND_PACKAGE(Eigen3 3.0.1 REQUIRED)

# ----- Add include directories -----

include(${Geant4_USE_FILE})
include_directories( ${EIGEN3_INCLUDE_DIR} /home/robbie/CppLibs/libgp/include)



file(GLOB sources ${PROJECT_SOURCE_DIR}/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/*.hh)

#
add_executable(BWT BreitWheelerTest.cpp ${sources} ${headers})
target_link_libraries(BWT ${Geant4_LIBRARIES} /home/robbie/CppLibs/libgp/libgp.a)
