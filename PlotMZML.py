# This script plots the mzml files
# It supplies boxes to plot mz-time traces within a certain window
#     The window is hard-coded for now
# This script requires the files created with Convert()
# Call with: # bokeh serve --show PlotMZML.py --port 2001

import numpy as np
from pdb import set_trace
from copy import copy
import pickle
import hdfdict
import h5py

from bokeh.io import curdoc
from bokeh.layouts import column, row, layout, widgetbox
from bokeh.models import ColumnDataSource, Slider, TextInput, Label
from bokeh.plotting import figure
from bokeh.events import DoubleTap

print('Loading Dictionaries')
TimeMZFile = '/Users/nate/Desktop/temporary/time_mz.hdf5'
TimeIntensityFile = '/Users/nate/Desktop/temporary/time_intensity.hdf5'
TimeMZDict = hdfdict.load(TimeMZFile)
TimeIntensityDict = hdfdict.load(TimeIntensityFile)

print('Extracting Data')
nScans = len(TimeMZDict.keys())
TimeArray = np.array([])
TimeArray = [np.append(TimeArray,float(key)) for key in TimeMZDict.keys()]

for i in np.arange(0,len(TimeArray)):
    TimeArray[i] = TimeArray[i][0]
TimeArray = np.sort(TimeArray)


TicArray = np.array([])
TicArray = [np.append(TicArray,sum(TimeIntensityDict[str(key)])) for key in TimeArray]
for i in np.arange(0,len(TicArray)):
    TicArray[i] = TicArray[i][0]


# create a blank array in the dictionary with mz values as keys
BlankData = np.zeros(nScans)
for i in range(0,len(BlankData)):
    BlankData[i] = np.nan


# Initialize the dictionary containing the data for the intensity vs. time plot
#     Do so for both lines plotted
intensity_vs_time_plot_dict = {'x':TimeArray,'y':TicArray}
intensity_vs_time_plot_source = ColumnDataSource(data=intensity_vs_time_plot_dict)

intensity_vs_time_plot_dict2 = {'x':TimeArray,'y':BlankData}
intensity_vs_time_plot_source2 = ColumnDataSource(data=intensity_vs_time_plot_dict2)

# Create the intensity vs. time plot
#     There are two line objects because there are two lines plotted
intensity_vs_time_plot=figure(title='total ion count vs. retention time',plot_height=300, plot_width=1000)
intensity_vs_time_plot.line('x','y',source=intensity_vs_time_plot_source,color='blue')
intensity_vs_time_plot.line('x','y',source=intensity_vs_time_plot_source2,color='red')

# Initialize the dictionary containing data for the intensity vs. mz plot
InitialTimePoint = TimeArray[0]
intesity_vs_mz_plot_dict = {}
intesity_vs_mz_plot_dict['x'] = TimeMZDict[str(InitialTimePoint)]
intesity_vs_mz_plot_dict['y'] = TimeIntensityDict[str(InitialTimePoint)]
intesity_vs_mz_plot_source = ColumnDataSource(data=intesity_vs_mz_plot_dict)
# Create the intensity vs. mz plot
intesity_vs_mz_plot = figure(title='intensity vs. mz', x_axis_label='m/z',y_axis_label='ion counts',plot_width=1000,plot_height=300)
intesity_vs_mz_plot.line('x','y',source=intesity_vs_mz_plot_source,color='firebrick')

print('Defining callbacks')

# Create the callback to change the time point from which the data for the intensity vs. mz plot is taken
#     Double clicking on the intensity vs. time plot selects this new time point
def callback(event):
    rt_click = event.x
    rt_index = np.argmin(abs(TimeArray-rt_click))
    rt = TimeArray[rt_index]
    x_index = rt
    y_index = rt
    intesity_vs_mz_plot_source.data = dict(x=TimeMZDict[str(x_index)], y=TimeIntensityDict[str(y_index)])
intensity_vs_time_plot.on_event(DoubleTap, callback)

#####Need to Fix UpdateMZ callbacks!!


# Create the callback to change the mz plotted for the 1st line
def UpdateMZ(attrname, old, new):
    BoxValue = MZ1.value
    # Obtain the mz range that is plotted
    #     The range is centered on BoxValue
    BinValue = float(MZ1Bin.value)
    if (BoxValue != 'tic') & (BoxValue != 'blank'):
        mz = float(BoxValue)
        mzRange = np.array([mz-BinValue,mz+BinValue])
        MZsInRangeIndices = (unique_mzs >= mzRange[0]) & (unique_mzs <= mzRange[1])
        MZsInRange = unique_mzs[MZsInRangeIndices]
        NewY = np.zeros(n_scans)
        # Sum the intensities for all mz values within the range for all time points
        for mz in MZsInRange:
            # Initialize the counter for the position of the array corresponding to the mz key of the dictionary containing intensities
            TimeIndexIterator = 0
            # Iterate through the indices of time_array where an intensity for the current mz is recoreded
            for TimeIndex in mz_TimePointIndex_dict[mz]:
                TimeIndex = TimeIndex.astype('int')
                NewY[TimeIndex] = NewY[TimeIndex] + mz_intensity_dict[mz][TimeIndexIterator]
                TimeIndexIterator = TimeIndexIterator + 1
        y = NewY
    if BoxValue == 'tic':
        y = tic_array
    if BoxValue == 'blank':
        y = blank_data
    intensity_vs_time_plot_source.data = dict(x=time_array,y=y)

# Create the callback to change the mz plotted for the 2nd line
def UpdateMZ2(attrname, old, new):
    BoxValue = MZ2.value
    BinValue = float(MZ2Bin.value)
    if (BoxValue != 'tic') & (BoxValue != 'blank'):
        mz = float(BoxValue)
        mzRange = np.array([mz-BinValue,mz+BinValue])
        MZsInRangeIndices = (unique_mzs >= mzRange[0]) & (unique_mzs <= mzRange[1])
        MZsInRange = unique_mzs[MZsInRangeIndices]
        NewY = np.zeros(n_scans)
        for mz in MZsInRange:
            TimeIndexIterator = 0
            for TimeIndex in mz_TimePointIndex_dict[mz]:
                TimeIndex = TimeIndex.astype('int')
                NewY[TimeIndex] = NewY[TimeIndex] + mz_intensity_dict[mz][TimeIndexIterator]
                TimeIndexIterator = TimeIndexIterator + 1
        y = NewY
    if BoxValue == 'tic':
        y = tic_array
    if BoxValue == 'blank':
        y = blank_data
    intensity_vs_time_plot_source2.data = dict(x=time_array,y=y)

# Create the first text box set
#     The Bin is the box that sets the range of the mz values included
MZ1 = TextInput(title = 'mz range', value=str('tic')) #the textbox widget, the value must be a string
MZ1.on_change('value',UpdateMZ)
MZ1Bin = TextInput(title='MZ1 Bin', value=str(0))
MZ1Bin.on_change('value',UpdateMZ)

# Create the second text box set
#     The Bin is the box that sets the range of the mz values included
MZ2 = TextInput(title = 'mz value range', value=str('tic')) #the textbox widget, the value must be a string
MZ2.on_change('value',UpdateMZ2)
MZ2Bin = TextInput(title='MZ2 Bin', value=str(0))
MZ2Bin.on_change('value',UpdateMZ2)


l = layout([
  [MZ1,MZ1Bin],
  [MZ2,MZ2Bin],
  [intensity_vs_time_plot],
  [intesity_vs_mz_plot],
], sizing_mode='fixed')

curdoc().add_root(l)
