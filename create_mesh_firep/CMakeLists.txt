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

#include_directories(AFTER /home/mkerber/gudhi-2.1.0/include)

#find_package(Boost REQUIRED Timer)

#find_package(GUDHI)



if ( CGAL_FOUND )

  add_executable(main "main.cpp")
  target_link_libraries(main ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

