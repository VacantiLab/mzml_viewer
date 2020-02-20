from pyteomics import mzml, auxiliary
import numpy as np

from bokeh.layouts import column, row, widgetbox, layout
from bokeh.models import CustomJS, ColumnDataSource, Slider, TextInput, Select
from bokeh.plotting import Figure, output_file, show, reset_output
from bokeh.io import curdoc


MZML = mzml.read('1051x_peptides_50nmol.mzML',dtype=dict)

intensity = np.zeros(12950)
time = np.zeros(12950)

MyDict = {}
i = 0
for key in MZML:
#    MyDict[i] = key
    intensity[i] = sum(key['intensity array'])
    time[i] = float(key['scanList']['scan'][0]['scan start time'])
    i = i+1
    print(i)

n_keys = i

source_dict = {'x':time, 'y':intensity}
source = ColumnDataSource(data=source_dict)

output_file('mzml_plot.html', title='peptides', mode='inline')

plot=Figure(title='ion counts vs. time', x_axis_label='retention time',y_axis_label='ion counts',plot_width=950,plot_height=300)

plot.line('x','y',source=source)
l = layout([
    [plot],
],sizing_mode='fixed')

show(l)
