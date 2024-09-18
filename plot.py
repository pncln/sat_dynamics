"""
This module provides a GUI application for visualizing satellite data. It reads a CSV file containing the satellite data, and creates a tabbed interface with a plot for each parameter in the data.

The `create_plot` function takes a parameter name as input and returns a matplotlib figure object that plots the data for that parameter over time.

The `MainWindow` class is the main entry point for the application. It creates the tab widget and adds a tab for each parameter, with the corresponding plot displayed in each tab.
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from PyQt6.QtWidgets import QApplication, QMainWindow, QTabWidget, QVBoxLayout, QWidget, QScrollBar
from PyQt6.QtCore import Qt
import sys
import numpy as np
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT

# Read the CSV file
data = pd.read_csv("satellite_data.csv")

# Get the parameter names from the column names
parameter_names = data.columns

# Create a function to create a plot for a given parameter
import numpy as np

from matplotlib.widgets import RangeSlider



def create_plot(parameter):
    fig, ax = plt.subplots(figsize=(10, 6))
    
    valid_mask = ~(np.isnan(data['Time']) | np.isinf(data['Time']) | np.isnan(data[parameter]) | np.isinf(data[parameter]))
    valid_data = data[valid_mask]
    
    if valid_data.empty:
        ax.text(0.5, 0.5, f"No valid data for {parameter}", ha='center', va='center')
        return fig
    
    ax.plot(valid_data['Time'], valid_data[parameter])
    ax.set_xlabel('Time (seconds)')
    ax.set_ylabel(parameter)
    ax.set_title(f'{parameter} over Time')
    ax.grid(True)

    x_min, x_max = valid_data['Time'].min(), valid_data['Time'].max()
    y_min, y_max = valid_data[parameter].min(), valid_data[parameter].max()
    
    if x_min == x_max:
        x_min, x_max = x_min - 1, x_max + 1
    if y_min == y_max:
        y_min, y_max = y_min - 1, y_max + 1
    
    x_range = x_max - x_min
    ax.set_xlim(x_min, x_min + x_range * 0.1)
    ax.set_ylim(y_min, y_max)

    return fig

class ScrollableCanvas(QWidget):
    def __init__(self, fig):
        super().__init__()
        self.canvas = FigureCanvasQTAgg(fig)
        self.scroll_bar = QScrollBar(Qt.Orientation.Horizontal)
        
        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.addWidget(self.scroll_bar)
        self.setLayout(layout)
        
        self.ax = self.canvas.figure.axes[0]
        self.full_xlim = self.ax.get_xlim()
        self.visible_range = (self.full_xlim[1] - self.full_xlim[0]) * 0.1
        
        self.scroll_bar.setRange(0, 1000)
        self.scroll_bar.valueChanged.connect(self.update_plot)
        
        self.canvas.mpl_connect('draw_event', self.on_draw)
        self.canvas.mpl_connect('motion_notify_event', self.on_pan)
        
        self.is_panning = False
        
    def update_plot(self, value):
        if not self.is_panning:
            total_range = self.full_xlim[1] - self.full_xlim[0]
            new_min = self.full_xlim[0] + (total_range - self.visible_range) * (value / 10000)
            if np.isfinite(new_min):
                self.ax.set_xlim(new_min, new_min + self.visible_range)
                self.canvas.draw()
        
    def on_draw(self, event):
        visible_range = self.ax.get_xlim()[1] - self.ax.get_xlim()[0]
        if visible_range != self.visible_range:
            self.visible_range = visible_range
        self.update_scroll_bar()
        
    def on_pan(self, event):
        if event.button == 1:  # Left mouse button
            self.is_panning = True
            self.update_scroll_bar()
            self.is_panning = False
    
    def update_scroll_bar(self):
        total_range = self.full_xlim[1] - self.full_xlim[0]
        visible_start = self.ax.get_xlim()[0] - self.full_xlim[0]
        position = (visible_start / (total_range - self.visible_range)) * 10000
        if np.isfinite(position):
            self.scroll_bar.setValue(max(min(int(position), 10000), 0))

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Satellite Data Plots")
        self.setGeometry(100, 100, 800, 600)

        self.tab_widget = QTabWidget()
        self.setCentralWidget(self.tab_widget)

        for parameter in data.columns[1:]:
            tab = QWidget()
            layout = QVBoxLayout()
            tab.setLayout(layout)

            fig = create_plot(parameter)
            scrollable_canvas = ScrollableCanvas(fig)
            layout.addWidget(scrollable_canvas)

            toolbar = NavigationToolbar2QT(scrollable_canvas.canvas, self)
            layout.addWidget(toolbar)

            self.tab_widget.addTab(tab, parameter)



if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
