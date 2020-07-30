# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Code_ )

cmake_minimum_required(VERSION 2.8.10)

set(CMAKE_CXX_STANDARD 14)

set(CGAL_DIR "~/CGAL-4.14")

set(CGAL_Boost_USE_STATIC_LIBS ON)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_BUILD_TYPE Release)

find_package(CGAL QUIET COMPONENTS Core )


if ( CGAL_FOUND )

  add_executable(points_inside_sphere "points_inside_sphere.cpp")
  target_link_libraries(points_inside_sphere ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

  add_executable(points_on_sphere "points_on_sphere.cpp")
  target_link_libraries(points_on_sphere ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

