#include "orbit_test/orbit_test.hpp"

OrbitTestNode::OrbitTestNode()
    : Node("orbit_test_node"),
      orbit_(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0){

  // Declare parameters first
  this->declare_parameter<double>("semi_major", 0.0);
  this->declare_parameter<double>("eccentricity", 0.0);
  this->declare_parameter<double>("inclination", 0.0);
  this->declare_parameter<double>("ascending", 0.0);
  this->declare_parameter<double>("perigee", 0.0);
  this->declare_parameter<double>("mean_anomaly", 0.0);

  // Retrieve parameters
  double semi_major = this->get_parameter("semi_major").as_double();
  double eccentricity = this->get_parameter("eccentricity").as_double();
  double inclination = this->get_parameter("inclination").as_double();
  double ascending = this->get_parameter("ascending").as_double();
  double perigee = this->get_parameter("perigee").as_double();
  double mean_anomaly = this->get_parameter("mean_anomaly").as_double();

  // Initialize OrbitMechanics AFTER parameter declaration
  this->orbit_ = OrbitMechanics(
      semi_major, eccentricity, inclination,
      ascending, perigee, mean_anomaly,
      0, 0, 0, 0, 0, 0  // Provide all 12 expected arguments
  );

  this->publisher_ = this->create_publisher<std_msgs::msg::Float32MultiArray>("/orbit_data", 10);
  this->dt_ = 0.1; // you have to match dt to timer_
  this->timer_ = this->create_wall_timer(
      std::chrono::milliseconds(100), // 10 Hz
      std::bind(&OrbitTestNode::orbitPropCallback, this));

  RCLCPP_INFO(this->get_logger(), "Orbit Test Node initialized.");
}

void OrbitTestNode::orbitPropCallback() {
  auto msg = std_msgs::msg::Float32MultiArray();
  this->state_.clear();
  Vector3d zero_accel = Vector3d::Zero();
  Vector3d zero_force = Vector3d::Zero();
  Vector6d cartesian_state = this->orbit_.getCartesian();
  auto state = this->orbit_.integrateStateRK4(cartesian_state,
                                              this->dt_,
                                              100.0,
                                              zero_accel,
                                              zero_force);
  this->state_.push_back(state(0));
  this->state_.push_back(state(1));
  this->state_.push_back(state(2));
  this->state_.push_back(state(3));
  this->state_.push_back(state(4));
  this->state_.push_back(state(5));

  msg.layout.dim.resize(1);
  msg.layout.dim[0].label = "orbit propagation test";
  msg.layout.dim[0].size = 6; // x y z vx vy vz
  msg.layout.dim[0].stride = 1;
  msg.layout.data_offset = 0;
  msg.data = this->state_;

  this->publisher_->publish(msg);
}
