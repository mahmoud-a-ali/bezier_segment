cmake_minimum_required(VERSION 3.0.2)
project(bezier_segment)

## Compile as C++11, supported in ROS Kinetic and newer
# add_compile_options(-std=c++11)


find_package(PythonLibs 2.7 REQUIRED)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  trajectory_msgs
  joint_trajectory_controller
)
catkin_package( )


## Your package locations should be listed before other locations
include_directories(
# include
  ${catkin_INCLUDE_DIRS}
)

include_directories(include ${catkin_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS}) 

add_library(bezier_segment
    include/bezier_segment/bezier_quintic_segment.h
    include/bezier_segment/bezier_quintic_segment.cpp
)

add_executable(segment_example              src/segment_example.cpp)

target_link_libraries(segment_example                ${catkin_LIBRARIES} bezier_segment ) 

