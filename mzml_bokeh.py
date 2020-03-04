import numpy as np
from pyteomics import mzml, auxiliary
import pdb
from copy import copy

from bokeh.io import curdoc
from bokeh.layouts import column, row, layout
from bokeh.models import ColumnDataSource, Slider, TextInput
from bokeh.plotting import figure

MZML = mzml.read('/Users/nate/Dropbox/Research/Vacanti_Laboratory/mzml_files/20200226_1522_HeLa.mzML',dtype=dict)

intensity = np.zeros(18440)
time = np.zeros(18440)
n_intensities = np.zeros(18440)
most_intensities = 28463

MyDict = {}
#xvalues = np.arange(100,3000,0.0001)
#n_mzs = len(xvalues)
#MyDict['x'] = xvalues
i = 0
for key in MZML:
    time[i] = float(key['scanList']['scan'][0]['scan start time'])
    #mzs = np.round(key['m/z array'],decimals=4)
    #mz_indices = np.isin(xvalues,mzs)
    #intensities = np.zeros(n_mzs)
    #intensities[mz_indices] = key['intensity array']
    #if i==0:
    #    MyDict['y'] = copy(intensities[mz_indices])
    xarray = np.array(key['m/z array'],dtype='float') #these must be redefined as type 'float' from 'object' so that the 'maximum' attribute works in np.pad
    yarray = np.array(key['intensity array'],dtype='float') #these must be redefined as type 'float' from 'object' so that the 'maximum' attribute works in np.pad

    n_zeros = most_intensities - len(xarray)

    xarray = np.pad(xarray,[0,n_zeros],'maximum')
    yarray = np.pad(yarray,[0,n_zeros],'constant')


    if i==0:
        MyDict['x'] = xarray
        MyDict['y'] = yarray
        source_dict2 = {}
        source_dict2['x'] = xarray
        source_dict2['y'] = yarray

    #MyDict seems to get too big and the function ColumnDataSource seems to time out
    #    the nunber of scans seems to be the problem, but the number of mz values may also be problematic
    if i < 4:
        MyDict['x'+'%.3f'%(time[i])] = xarray
        MyDict['y'+'%.3f'%(time[i])] = yarray
        n_intensities[i] = len(MyDict['y'+'%.3f'%(time[i])])

    intensity[i] = sum(key['intensity array'])
    i = i+1
    print(i)

n_keys = i
#most_intensities = max(n_intensities)

source_dict = {'x':time, 'y':intensity}
source = ColumnDataSource(data=source_dict)

plot=figure(title='total ion count vs. retention time',plot_height=300, plot_width=1000)
plot.line('x','y',source=source)


print('sourcing 2')
source2 = ColumnDataSource(data=MyDict)
print('done sourcing 2')

plot2 = figure(title='intensity vs. mz', x_axis_label='m/z',y_axis_label='ion counts',plot_width=1000,plot_height=300)
plot2.line('x','y',source=source2,color='firebrick')
#plot2.vbar(x='x', bottom=0, width=0.01, top='y',color='firebrick',source=source2)


# time = np.array([1,2,3])
# intensity = np.array([1,2,3])
# source_dict = {'x':time, 'y':intensity}
# source = ColumnDataSource(data=source_dict)


l = layout([
  [plot2],
], sizing_mode='fixed')

curdoc().add_root(l)
