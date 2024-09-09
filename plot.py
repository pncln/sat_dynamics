"""
This module provides a GUI application for visualizing satellite data. It reads a CSV file containing the satellite data, and creates a tabbed interface with a plot for each parameter in the data.

The `create_plot` function takes a parameter name as input and returns a matplotlib figure object that plots the data for that parameter over time.

The `MainWindow` class is the main entry point for the application. It creates the tab widget and adds a tab for each parameter, with the corresponding plot displayed in each tab.
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from PyQt6.QtWidgets import QApplication, QMainWindow, QTabWidget, QVBoxLayout, QWidget
import sys

# Read the CSV file
data = pd.read_csv("satellite_data.csv")

# Get the parameter names from the column names
parameter_names = data.columns

# Create a function to create a plot for a given parameter
def create_plot(parameter):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(data['Time'], data[parameter])
    ax.set_xlabel('Time (seconds)')
    ax.set_ylabel(parameter)
    ax.set_title(f'{parameter} over Time')
    ax.grid(True)
    return fig

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Satellite Data Plots")
        self.setGeometry(100, 100, 800, 600)

        # Create a tab widget
        self.tab_widget = QTabWidget()
        self.setCentralWidget(self.tab_widget)

        # Create a tab for each parameter (excluding the Time column)
        for parameter in data.columns[1:]:
            tab = QWidget()
            layout = QVBoxLayout()
            tab.setLayout(layout)

            fig = create_plot(parameter)
            canvas = FigureCanvasQTAgg(fig)
            layout.addWidget(canvas)

            self.tab_widget.addTab(tab, parameter)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
