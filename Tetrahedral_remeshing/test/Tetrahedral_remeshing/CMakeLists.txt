# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.14)

project(Tetrahedral_remeshing_Tests)

# CGAL and its components
find_package(CGAL REQUIRED COMPONENTS ImageIO)

# Boost and its components
find_package(Boost REQUIRED)
if(NOT Boost_FOUND)
  message(
    STATUS "This project requires the Boost library, and will not be compiled.")
  return()
endif()

# Creating entries for all C++ files with "main" routine
# ##########################################################
create_single_source_cgal_program("test_tetrahedral_remeshing.cpp")
create_single_source_cgal_program(
  "test_tetrahedral_remeshing_with_features.cpp")
create_single_source_cgal_program(
  "test_tetrahedral_remeshing_of_one_subdomain.cpp")
create_single_source_cgal_program("test_tetrahedral_remeshing_io.cpp")
create_single_source_cgal_program(
  "test_tetrahedral_remeshing_from_mesh_file.cpp")

# Tests using Mesh_3 require Eigen
find_package(Eigen3 3.1.0 REQUIRED) #(3.1.0 or greater)
include( ${EIGEN3_USE_FILE} )

create_single_source_cgal_program(
  "test_mesh_and_remesh_polyhedral_domain_with_features.cpp")

if(CGAL_ImageIO_USE_ZLIB)
  create_single_source_cgal_program("test_mesh_and_remesh_image.cpp")
endif()
