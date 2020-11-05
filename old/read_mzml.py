from pyteomics import mzml, auxiliary
import numpy as np
from pdb import set_trace

from bokeh.layouts import column, row, widgetbox, layout
from bokeh.models import CustomJS, ColumnDataSource, Slider, TextInput, Select
from bokeh.plotting import Figure, output_file, show, reset_output
from bokeh.io import curdoc

name = '20200226_1522_HeLa_NoFilters'
directory = '/Users/nate/Dropbox/Research/Vacanti_Laboratory/mzml_files/'

MZML = mzml.read(directory + name + '.mzML',dtype=dict)
i = 0
for key in MZML:
    i = i+1
    print(i)
n_keys = i

intensity = np.zeros(n_keys)
time = np.zeros(n_keys)

MZML = mzml.read(directory + name + '.mzML',dtype=dict)
MyDict = {}
i = 0
for key in MZML:
#    MyDict[i] = key
    intensity[i] = sum(key['intensity array'])
    time[i] = float(key['scanList']['scan'][0]['scan start time'])
    i = i+1
    print(i)

source_dict = {'x':time, 'y':intensity}
source = ColumnDataSource(data=source_dict)

output_file(name+'.html', title=name, mode='inline')

plot=Figure(title='ion counts vs. time', x_axis_label='retention time',y_axis_label='ion counts',plot_width=950,plot_height=300)

plot.line('x','y',source=source)
l = layout([
    [plot],
],sizing_mode='fixed')

show(l)
