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
mzml_file_directory = '/Users/nate/Dropbox/Research/Vacanti_Laboratory/projects/PolyMID/correction_program/references/mzml_files/2018_1016_02_filtered.mzML'

# Set the lower and upper bounds on the mzs for which a time trace will be available
#     This is necessary because the number of mzs in high resolution data is over 1 million and the script will take too long to create the necessary dictionary
mzLower = 105
mzUpper = 110

# Read the mzML file as an iterable object
MZML = mzml.read(mzml_file_directory,dtype=dict)

# Initialize array of all time points
time_array = np.array([])

# Initialize array of unique mz values measured
unique_mzs = np.array([])

# Initialize dictionaries with time points as the keys
time_mz_dict = {}
time_intensity_dict = {}

# Initialize dictionaries with time points as keys
#     these only include the mz values specified by mzLower and mzUpper
time_mz_dict_limited = {}
time_intensity_dict_limited = {}

# Iterate through the scnas of the iterable object of the mzML file
i=0
for key in MZML:
    # get the current time associated with the current scan
    time_point = float(key['scanList']['scan'][0]['scan start time'])
    time_array = np.append(time_array,time_point)

    # get all of the mz values scanned in the current scan
    mzs = np.array(key['m/z array'],dtype='float')

    # get the intensities associated with those mz values
    intensities = np.array(key['intensity array'],dtype='float')

    # store the mz values and intensity values in their corresponding dictionaries
    time_mz_dict[time_point] = mzs
    time_intensity_dict[time_point] = intensities

    # parse out the values assiciated with mz values that are between the lower and upper bounds
    indices = (mzs > mzLower) & (mzs < mzUpper)
    mzsToKeep = mzs[indices]
    IntesitiesToKeep = intensities[indices]

    # store those values in separate "limited" dictionaries
    time_mz_dict_limited[time_point] = mzsToKeep
    time_intensity_dict_limited[time_point] = IntesitiesToKeep

    # track the unique mz values across all scans
    unique_mzs = np.append(unique_mzs,mzsToKeep)
    unique_mzs = np.unique(unique_mzs)

    i = i+1
    print(i)

# sort the array of unique mz values from smallest to largest
unique_mzs = np.sort(unique_mzs)

# Initialize a dictionary where the keys are the unique mz values
#     Each entry will be an array of intensities
#         Each position of the array is associated with the time point at the same position in the time_array
n_scans = i
mz_time_dict = {key:np.zeros(n_scans) for key in unique_mzs}

# Fill the dictionary where unique mz values are the keys
# Iterate over the time points
time_point_index = 0
for time_point in time_array:
    # Iterate over the mz values scanned at the current time point
    #     This is performed in the limited time_mz_dict because it was created to only hold the mz values of interest
    #       All mz values cannot be considered because there are too many
    mz_index = 0
    for mz in time_mz_dict_limited[time_point]:
        # Fill the dictionary with unique mz values as the keys
        #     Each unique mz value key is associated with an array with intensity entries for each time point
        #
        CurrentIntensity = time_intensity_dict_limited[time_point][mz_index]
        mz_time_dict[mz][time_point_index] = CurrentIntensity
        mz_index = mz_index + 1

    # Iterate to the next time point
    time_point_index = time_point_index + 1
    print('time point index = '+ str(time_point_index))

set_trace()



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
    #xarray = np.pad(xarray,[0,n_zeros],'maximum')
    #yarray = np.pad(yarray,[0,n_zeros],'constant')

    if i==0:
        ic_mz_plot_dict = {}
        ic_mz_plot_dict['x'] = xarray
        ic_mz_plot_dict['y'] = yarray

    ic_mz_dict['x'+'%.3f'%(time[i])] = xarray
    ic_mz_dict['y'+'%.3f'%(time[i])] = yarray

    tic = np.append(tic,sum(key['intensity array']))

    i = i+1
    print(i)

tic_dict = {'x':time, 'y':tic}
tic_source = ColumnDataSource(data=tic_dict)

tic_plot=figure(title='total ion count vs. retention time',plot_height=300, plot_width=1000)
tic_plot.line('x','y',source=tic_source)


print('sourcing 2')
ic_mz_plot_source = ColumnDataSource(data=ic_mz_plot_dict)
print('done sourcing 2')

ic_mz_plot = figure(title='intensity vs. mz', x_axis_label='m/z',y_axis_label='ion counts',plot_width=1000,plot_height=300)
ic_mz_plot.line('x','y',source=ic_mz_plot_source,color='firebrick')

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
