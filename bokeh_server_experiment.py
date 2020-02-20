# myapp.py

import numpy as np

from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Slider, TextInput
from bokeh.plotting import figure

# create a plot and style its properties
plot=figure(plot_height=400, plot_width=400, title="my sine wave",
              tools="crosshair,pan,reset,save,wheel_zoom",
              x_range=[0, 5], y_range=[0, 5])
time = np.array([1,2,3])
intensity = np.array([1,2,3])
source_dict = {'x':time, 'y':intensity}
source = ColumnDataSource(data=source_dict)
plot.line('x','y',source=source)







# put the button and plot in a layout and add to the document
curdoc().add_root(row(plot, width = 800))
