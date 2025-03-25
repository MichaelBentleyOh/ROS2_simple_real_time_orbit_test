cmake_minimum_required(VERSION 3.5)
project(orbit_test)
# done for now
# Default to C++14
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14)
endif()

# =====================Find dependencies=============================
find_package(ament_cmake REQUIRED)
# find_package(ament_python REQUIRED)

# ROS2 package
find_package(rclcpp REQUIRED)
find_package(std_msgs REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
# Eigen
find_package(eigen3_cmake_module REQUIRED)
find_package(Eigen3 REQUIRED)

set(dependencies "rclcpp" "std_msgs")

# launch
# find_package(launch REQUIRED)
# find_package(launch_ros REQUIRED)
# ===================================================================

# ====================Add source libraries===========================
add_library(orbit_test src/orbit_test.cpp)

# ===================Add execution symbols===========================
add_executable(execution src/orbit_test.cpp src/execution.cpp)

# Link your executable to the necessary dependencies
ament_target_dependencies(
  orbit_test
  rclcpp
  std_msgs
  Eigen3
)

ament_target_dependencies(
  execution
  rclcpp
  std_msgs
  Eigen3
)

# =====================Include headers===============================
include_directories(include)
include_directories(include/${PROJECT_NAME})

# =====================Include headers===============================
# Install the executable
install(TARGETS execution DESTINATION lib/${PROJECT_NAME})

# Install launch files
install(DIRECTORY launch
  DESTINATION share/${PROJECT_NAME}
)

# Install parameters file
install(DIRECTORY config
  DESTINATION share/${PROJECT_NAME}
)

# Install header and source files
install(DIRECTORY include
  DESTINATION include/${PROJECT_NAME}
)
install(DIRECTORY src
  DESTINATION src/${PROJECT_NAME}
)

ament_package()
