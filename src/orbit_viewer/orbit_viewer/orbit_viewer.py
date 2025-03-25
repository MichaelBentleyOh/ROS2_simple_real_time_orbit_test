import rclpy
from rclpy.node import Node
from std_msgs.msg import Float32MultiArray
import multiprocessing
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import csv
from mpl_toolkits.mplot3d import Axes3D

def plot_process(data_queue):
    """Separate process for real-time 3D plotting."""
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.set_xlim([-10000000, 10000000])
    ax.set_ylim([-10000000, 10000000])
    ax.set_zlim([-10000000, 10000000])

    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('z [m]')

    x_data, y_data, z_data = [], [], []

    def update_plot(frame):
        """Update function for animation."""
        while not data_queue.empty():
            try:
                new_data = data_queue.get_nowait()
                x_data.append(new_data[0])
                y_data.append(new_data[1])
                z_data.append(new_data[2])

                # Limit stored data to prevent memory overflow
                if len(x_data) > 100:
                    x_data.pop(0)
                    y_data.pop(0)
                    z_data.pop(0)
            except Exception:  # Avoid queue.Empty as multiprocessing.Queue does not use it
                pass

        ax.set_xlim([-10000000, 10000000])
        ax.set_ylim([-10000000, 10000000])
        ax.set_zlim([-10000000, 10000000])

        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_zlabel('z [m]')

        if x_data:
            ax.plot(x_data, y_data, z_data, marker='o', linestyle='-', color='b')

    ani = FuncAnimation(fig, update_plot, interval=100)
    plt.show()

class OrbitPlotNode(Node):
    def __init__(self, data_queue):
        super().__init__("orbit_plot_node")
        self.subscription = self.create_subscription(
            Float32MultiArray,
            'orbit_data',
            self.listener_callback,
            10)
        self.data_queue = data_queue

        # Open CSV file once and append data
        self.csv_file = open('data.csv', 'a', newline='')
        self.csv_writer = csv.writer(self.csv_file)

        self.pos_x = 0
        self.pos_y = 0
        self.pos_z = 0

    def listener_callback(self, msg):
        """Callback function for ROS 2 subscriber."""
        self.get_logger().info(f"Data received: {msg.data}")
        if len(msg.data) >= 3:
            self.data_queue.put([msg.data[0], msg.data[1], msg.data[2]])

        self.pos_x = msg.data[0]
        self.pos_y = msg.data[1]
        self.pos_z = msg.data[2]

        # Write to CSV
        self.csv_writer.writerow([self.pos_x, self.pos_y, self.pos_z])

    def destroy_node(self):
        """Ensure the CSV file is properly closed on shutdown."""
        self.csv_file.close()
        super().destroy_node()

def main():
    rclpy.init()

    data_queue = multiprocessing.Queue()

    # Start plotting process
    plotter_process = multiprocessing.Process(target=plot_process, args=(data_queue,))
    plotter_process.start()

    # Start ROS 2 node
    orbit_node = OrbitPlotNode(data_queue)

    try:
        rclpy.spin(orbit_node)
    except KeyboardInterrupt:
        pass
    finally:
        orbit_node.destroy_node()
        rclpy.shutdown()
        plotter_process.terminate()
        plotter_process.join()  # Ensure the process terminates cleanly

if __name__ == '__main__':
    main()

# # Reload the CSV with proper column names
# import pandas as pd
# df = pd.read_csv(file_path, header=None, names=["x", "y", "z"])

# # Display the first few rows again to confirm structure
# df.head()
