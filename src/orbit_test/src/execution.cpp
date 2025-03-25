#include "orbit_test/orbit_test.hpp"

int main(int argc, char **argv)
{
    rclcpp::init(argc, argv);
    rclcpp::spin(std::make_shared<OrbitTestNode>());
    rclcpp::shutdown();
    return 0;
}
// done for now
