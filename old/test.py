
import numpy as np
from pyteomics import mzml, auxiliary
import pdb
from copy import copy

from bokeh.io import curdoc
from bokeh.layouts import column, row, layout
from bokeh.models import ColumnDataSource, Slider, TextInput
from bokeh.plotting import figure
time = np.array([1,2,3,0,0])
intensity = np.array([1,2,3,0,0])
source_dict = {'x':time, 'y':intensity}
source = ColumnDataSource(data=source_dict)

plot=figure(title='total ion count vs. retention time',plot_height=300, plot_width=1000)
plot.line('x','y',source=source)


l = layout([
  [plot],
], sizing_mode='fixed')

curdoc().add_root(l)
