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

# Initialize the total ion chromatograph array
tic_array = np.array([])

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

    # get the value of the sum of all ions at the current time point
    tic = sum(intensities)
    tic_array = np.append(tic_array,tic)

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
mz_intensity_dict = {key:np.zeros(n_scans) for key in unique_mzs}

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
        mz_intensity_dict[mz][time_point_index] = CurrentIntensity
        mz_index = mz_index + 1

    # Iterate to the next time point
    time_point_index = time_point_index + 1
    print('time point index = '+ str(time_point_index))

# place the total ion chromatograph in the dictionary with mz values as keys
mz_intensity_dict['tic'] = tic_array

# Initialize the dictionary containing the data for the intensity vs. time plot
intensity_vs_time_plot_dict = {'x':time_array,'y':mz_intensity_dict['tic']}
intensity_vs_time_plot_source = ColumnDataSource(data=intensity_vs_time_plot_dict)
# Create the intensity vs. time plot
intensity_vs_time_plot=figure(title='total ion count vs. retention time',plot_height=300, plot_width=1000)
intensity_vs_time_plot.line('x','y',source=intensity_vs_time_plot_source)

# Initialize the dictionary containing data for the intensity vs. mz plot
InitialTimePoint = time_array[0]
intesity_vs_mz_plot_dict = {}
intesity_vs_mz_plot_dict['x'] = time_mz_dict[InitialTimePoint]
intesity_vs_mz_plot_dict['y'] = time_intensity_dict[InitialTimePoint]
intesity_vs_mz_plot_source = ColumnDataSource(data=intesity_vs_mz_plot_dict)
# Create the intensity vs. mz plot
intesity_vs_mz_plot = figure(title='intensity vs. mz', x_axis_label='m/z',y_axis_label='ion counts',plot_width=1000,plot_height=300)
intesity_vs_mz_plot.line('x','y',source=intesity_vs_mz_plot_source,color='firebrick')

# Create the callback to change the time point from which the data for the intensity vs. mz plot is taken
#     Double clicking on the intensity vs. time plot selects this new time point
def callback(event):
    rt_click = event.x
    rt_index = np.argmin(abs(time_array-rt_click))
    rt = time_array[rt_index]
    x_index = rt
    y_index = rt
    intesity_vs_mz_plot_source.data = dict(x=time_mz_dict[x_index], y=time_intensity_dict[y_index])
intensity_vs_time_plot.on_event(DoubleTap, callback)






l = layout([
  [intensity_vs_time_plot],
  [intesity_vs_mz_plot],
], sizing_mode='fixed')

curdoc().add_root(l)
