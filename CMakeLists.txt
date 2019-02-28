#########################################################
#                  Garfield Simulation                  #
#               Author : LAGARDE Francois               #
#########################################################

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(GarfieldSimulation VERSION 0.1 DESCRIPTION "Garfield Simulation" HOMEPAGE_URL "https://github.com/RPClab/Garfield" LANGUAGES CXX C)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard used for compiling" FORCE)
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "C++ standard required" FORCE)
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Allow C++ extensions" FORCE)

# Set the C standard
set(CMAKE_C_STANDARD 11 CACHE STRING "C standard used for compiling" FORCE)
set(CMAKE_C_STANDARD_REQUIRED ON CACHE BOOL "C standard required" FORCE)
set(CMAKE_C_EXTENSIONS OFF CACHE BOOL "Allow C extensions" FORCE)

option( BUILD_GARFIELDPP "Set to OFF to disable CLHEP build" ON )
set(GARFIELDPP_version    "master")
set(GARFIELDPP_repository "https://github.com/RPClab/Garfieldpp.git")

# Set default install prefix to project root directory if is not changed by the user
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/bin" CACHE PATH "Default install path for ${PROJECT_NAME}" FORCE)
endif()

# Set default cmake build type to RelWithDebInfo
# Possible options are: None Debug Release RelWithDebInfo MinSizeRel
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE )
endif()

# Default path for cmake modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" CACHE INTERNAL "Default path for cmake modules" FORCE)

# Append link pathes to rpath list
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE CACHE INTERNAL "Append directories in the linker search path and outside the project to the INSTALL_RPATH" FORCE)
# Add library install path to the rpath list
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" CACHE INTERNAL "Add library install path to the rpath list" FORCE)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS Geom Gdml)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include(ExternalProject)

# ----- GARFIELDPP package -----
ExternalProject_Add( 
                        GARFIELDPP
                        GIT_REPOSITORY ${GARFIELDPP_repository}
                        GIT_TAG ${GARFIELDPP_version}
                        GIT_PROGRESS TRUE
                        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} -DROOT_CMAKE_DIR=ROOT_USE_FILE -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
                        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
                        LIST_SEPARATOR %
                        UPDATE_COMMAND ""
)

find_package(Garfield)
message("LLLLL ${GARFIELD_INCLUDE_DIRS}")
include_directories(${ROOT_INCLUDE_DIRS} ${GARFIELD_INCLUDE_DIRS})
add_subdirectory(source)
