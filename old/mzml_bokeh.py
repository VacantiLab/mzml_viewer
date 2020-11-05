#This script plots the mzml files from proteomics experiments

import numpy as np
from pyteomics import mzml, auxiliary
from pdb import set_trace
from copy import copy

from bokeh.io import curdoc
from bokeh.layouts import column, row, layout
from bokeh.models import ColumnDataSource, Slider, TextInput
from bokeh.plotting import figure
from bokeh.events import DoubleTap

# bokeh serve --show mzml_bokeh.py
#mzml_file_directory = '/Users/nate/Dropbox/Research/Vacanti_Laboratory/mzml_files/QE1_QC_HeLa_20200225_r2.mzML'
mzml_file_directory = '/Users/nate/Dropbox/Research/Vacanti_Laboratory/projects/PolyMID/correction_program/references/mzml_files/2018_1016_02.mzML'

MZML = mzml.read(mzml_file_directory,dtype=dict)
n_intensities_array = np.array([])
i=0
for key in MZML:
    n_intensities_key = len(np.array(key['m/z array'],dtype='float'))
    n_intensities_array = np.append(n_intensities_array,n_intensities_key)
    i = i+1
    print(i)
most_intensities = int(max(n_intensities_array))

MZML = mzml.read(mzml_file_directory,dtype=dict)
time = np.array([])
tic = np.array([])
ic_mz_dict = {}
i = 0
for key in MZML:
    time_point = float(key['scanList']['scan'][0]['scan start time'])
    time = np.append(time,time_point)
    xarray = np.array(key['m/z array'],dtype='float') #these must be redefined as type 'float' from 'object' so that the 'maximum' attribute works in np.pad
    yarray = np.array(key['intensity array'],dtype='float') #these must be redefined as type 'float' from 'object' so that the 'maximum' attribute works in np.pad
    n_zeros = most_intensities - len(xarray)
    xarray = np.pad(xarray,[0,n_zeros],'maximum')
    yarray = np.pad(yarray,[0,n_zeros],'constant')

    if i==0:
        ic_mz_plot_dict = {}
        ic_mz_plot_dict['x'] = xarray
        ic_mz_plot_dict['y'] = yarray

    ic_mz_dict['x'+'%.3f'%(time[i])] = xarray
    ic_mz_dict['y'+'%.3f'%(time[i])] = yarray

    tic = np.append(tic,sum(key['intensity array']))

    i = i+1
    print(i)

# tic = np.zeros(n_keys)
# time = np.zeros(n_keys)

# MZML = mzml.read('/Users/nate/Dropbox/Research/Vacanti_Laboratory/mzml_files/20200226_1522_HeLa_NoFilters.mzML',dtype=dict)
# ic_mz_dict = {}
# i = 0
# for key in MZML:
#     time[i] = float(key['scanList']['scan'][0]['scan start time'])
#     #mzs = np.round(key['m/z array'],decimals=4)
#     #mz_indices = np.isin(xvalues,mzs)
#     #intensities = np.zeros(n_mzs)
#     #intensities[mz_indices] = key['intensity array']
#     #if i==0:
#     #    MyDict['y'] = copy(intensities[mz_indices])
#     xarray = np.array(key['m/z array'],dtype='float') #these must be redefined as type 'float' from 'object' so that the 'maximum' attribute works in np.pad
#     yarray = np.array(key['intensity array'],dtype='float') #these must be redefined as type 'float' from 'object' so that the 'maximum' attribute works in np.pad
#
#     n_zeros = most_intensities - len(xarray)
#
#     xarray = np.pad(xarray,[0,n_zeros],'maximum')
#     yarray = np.pad(yarray,[0,n_zeros],'constant')
#
#
#     if i==0:
#         ic_mz_plot_dict = {}
#         ic_mz_plot_dict['x'] = xarray
#         ic_mz_plot_dict['y'] = yarray
#
#     #MyDict seems to get too big and the function ColumnDataSource seems to time out
#     #    the nunber of scans seems to be the problem, but the number of mz values may also be problematic
#
#     ic_mz_dict['x'+'%.3f'%(time[i])] = xarray
#     ic_mz_dict['y'+'%.3f'%(time[i])] = yarray
#
#     tic[i] = sum(key['intensity array'])
#     i = i+1
#     print(i)
#
# n_keys = i
# #most_intensities = max(n_intensities)

tic_dict = {'x':time, 'y':tic}
tic_source = ColumnDataSource(data=tic_dict)

tic_plot=figure(title='total ion count vs. retention time',plot_height=300, plot_width=1000)
tic_plot.line('x','y',source=tic_source)


print('sourcing 2')
ic_mz_plot_source = ColumnDataSource(data=ic_mz_plot_dict)
print('done sourcing 2')

ic_mz_plot = figure(title='intensity vs. mz', x_axis_label='m/z',y_axis_label='ion counts',plot_width=1000,plot_height=300)
ic_mz_plot.line('x','y',source=ic_mz_plot_source,color='firebrick')
#ic_mz_plot.vbar(x='x', bottom=0, width=0.1, top='y',color='firebrick',source=ic_mz_plot_source)

# #add a dot where the click happened
# coordList=[]
# def callback(event):
#     Coords=(event.x,event.y)
#     coordList.append(Coords)
#     tic_source.data = dict(x=[i[0] for i in coordList], y=[i[1] for i in coordList])
# tic_plot.on_event(DoubleTap, callback)

#add a dot where the click happened
def callback(event):
    rt_click = event.x
    rt_index = np.argmin(abs(time-rt_click))
    rt = time[rt_index]
    x_index = 'x'+'%.3f'%(rt)
    y_index = 'y'+'%.3f'%(rt)
    ic_mz_plot_source.data = dict(x=ic_mz_dict[x_index], y=ic_mz_dict[y_index])
tic_plot.on_event(DoubleTap, callback)


# time = np.array([1,2,3])
# intensity = np.array([1,2,3])
# source_dict = {'x':time, 'y':intensity}
# source = ColumnDataSource(data=source_dict)


l = layout([
  [tic_plot],
  [ic_mz_plot],
], sizing_mode='fixed')

curdoc().add_root(l)
