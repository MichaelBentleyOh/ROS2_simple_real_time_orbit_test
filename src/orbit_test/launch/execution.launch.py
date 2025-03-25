import os
from launch import LaunchDescription
from launch_ros.actions import Node
from ament_index_python.packages import get_package_share_directory

def generate_launch_description():

    ld = LaunchDescription()

    config = os.path.join(
            get_package_share_directory('orbit_test'),
            'config',
            'initial_value.yaml'
    )

    orbit_node = Node(
                package="orbit_test",
                executable="execution",
                name="orbit_test_node",
                parameters=[config]
    )

    plot_node = Node(
                package="orbit_viewer",
                executable="orbit_viewer",
                name="orbit_plot_node",
    )

    ld.add_action(orbit_node)
    ld.add_action(plot_node)

    return ld
# done for now
