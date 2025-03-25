#ifndef ORBIT_TEST_HPP_
#define ORBIT_TEST_HPP_

#include <chrono>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <rclcpp/rclcpp.hpp>
#include <std_msgs/msg/float32_multi_array.hpp>
#include <std_msgs/msg/float32.hpp>

#include "orbit_test/Kinematics.hpp"
#include "orbit_test/OrbitMechanics.hpp"

using namespace std;
using namespace std::chrono_literals;

class OrbitTestNode : public rclcpp::Node
{
public:
  OrbitTestNode();
  void orbitPropCallback();
private:
  // Classes for orbit calculation
  OrbitMechanics orbit_;

  // Time step for integration
  double dt_;
  std::vector<float> state_;
  // ROS2 part

  rclcpp::TimerBase::SharedPtr timer_;
  rclcpp::Publisher<std_msgs::msg::Float32MultiArray>::SharedPtr publisher_;
}; //??;??

#endif // ORBIT_TEST_HPP_
// TBD
