cmake_minimum_required(VERSION 3.15)
project(opencvTest)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib/)

set(CMAKE_CXX_STANDARD 17)

# if necessary, set OpenCV_DIR variable, like below
# set(OpenCV_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../build)

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
add_executable(main
  ${PROJECT_SOURCE_DIR}/src/main.cpp
)
target_link_libraries(main ${OpenCV_LIBS})
