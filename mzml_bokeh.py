import numpy as np
from pyteomics import mzml, auxiliary
import pdb
from copy import copy

from bokeh.io import curdoc
from bokeh.layouts import column, row, layout
from bokeh.models import ColumnDataSource, Slider, TextInput
from bokeh.plotting import figure

MZML = mzml.read('1051x_peptides_50nmol.mzML',dtype=dict)

intensity = np.zeros(12950)
time = np.zeros(12950)

MyDict = {}
xvalues = np.arange(100,3000,0.0001)
n_mzs = len(xvalues)
MyDict['x'] = xvalues
i = 0
for key in MZML:
    time[i] = float(key['scanList']['scan'][0]['scan start time'])
    mzs = np.round(key['m/z array'],decimals=4)
    mz_indices = np.isin(xvalues,mzs)
    intensities = np.zeros(n_mzs)
    intensities[mz_indices] = key['intensity array']
    if i==0:
        MyDict['y'] = copy(intensities[mz_indices])
    MyDict1['x'+'%.3f'%(time[i])] = np.round(key['m/z array'],decimals=4)
    MyDict1['y'+'%.3f'%(time[i])] = key['intensity array']
    intensity[i] = sum(key['intensity array'])
    i = i+1
    print(i)

n_keys = i

source_dict = {'x':time, 'y':intensity}
source = ColumnDataSource(data=source_dict)

source2 = ColumnDataSource(data=MyDict)
plot2 = figure(title='intensity vs. mz', x_axis_label='m/z',y_axis_label='ion counts',plot_width=1000,plot_height=300)
plot2.line('x','y',source=source2,color='firebrick')
#plot2.vbar(x='x', bottom=0, width=0.01, top='y',color='firebrick',source=source2)


# time = np.array([1,2,3])
# intensity = np.array([1,2,3])
# source_dict = {'x':time, 'y':intensity}
# source = ColumnDataSource(data=source_dict)

plot=figure(title='total ion count vs. retention time',plot_height=300, plot_width=1000)
plot.line('x','y',source=source)

l = layout([
  [plot],
  [plot2],
], sizing_mode='fixed')

curdoc().add_root(l)
