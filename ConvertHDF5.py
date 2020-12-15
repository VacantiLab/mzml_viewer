def Convert():
# This function converts an mzML file to a HDF5 file in a form readable by PlotMZML.py

    from pyteomics import mzml, auxiliary
    import numpy as np
    from pdb import set_trace
    import hdfdict

    time_mz_HDF5_file_directory = '/Users/nate/Desktop/temporary/time_mz.hdf5'
    time_intensity_HDF5_file_directory = '/Users/nate/Desktop/temporary/time_intensity.hdf5'
    mz_intensity_HDF5_file_directory = '/Users/nate/Desktop/temporary/mz_intensity.hdf5'
    mz_TimePointIndex_dict_HDF5_file_directory = '/Users/nate/Desktop/temporary/mz_TimePointIndex.hdf5'

    # Specify the mzml file
    mzml_file_directory = '/Users/nate/Dropbox/Research/Vacanti_Laboratory/projects/PolyMID/correction_program/references/mzml_files/Tryptophan/YY2019060827.mzML'

    # Read the mzML file as an iterable object
    MZML = mzml.read(mzml_file_directory,dtype=dict)


    # Initialize dictionaries with time points as the keys
    time_mz_dict = {}
    time_intensity_dict = {}
    time_mz_dict['tic'] = np.array([])
    time_mz_dict['time'] = np.array([])

    # Initialize dictionaries with time points as keys
    #     these only include the mz values specified by mzLower and mzUpper
    time_mz_dict_limited = {}
    time_intensity_dict_limited = {}


    # Initialize the total ion chromatograph array
    tic_array = np.array([])

    n_total_mzs = 0

    # Iterate through the scnas of the iterable object of the mzML file
    i=0
    for key in MZML:
        # get the current time associated with the current scan
        time_point = float(key['scanList']['scan'][0]['scan start time'])

        # get all of the mz values scanned in the current scan
        mzs = np.array(key['m/z array'],dtype='float')

        n_total_mzs = n_total_mzs + len(mzs)

        # get the intensities associated with those mz values
        intensities = np.array(key['intensity array'],dtype='float')

        # get the value of the sum of all ions at the current time point
        tic = sum(intensities)
        tic_array = np.append(tic_array,tic)

        # store the mz values and intensity values in their corresponding dictionaries
        time_mz_dict[str(time_point)] = mzs
        time_mz_dict['tic'] = np.append(time_mz_dict['tic'],tic)
        time_mz_dict['time'] = np.append(time_mz_dict['time'],time_point)
        time_intensity_dict[str(time_point)] = intensities

        # Find another way to get unique MZs - from dictionary already made! This takes too much time
        # track the unique mz values across all scans
        #unique_mzs = np.append(unique_mzs,mzs)
        #unique_mzs = np.unique(unique_mzs)

        i = i+1
        if i%200 == 0:
            print('Reading mzML scan number: ' + str(i))

        #if i >= 200:
        #    break

    n_scans=i


    print('writing files')
    hdfdict.dump(time_mz_dict,time_mz_HDF5_file_directory)
    hdfdict.dump(time_intensity_dict,time_intensity_HDF5_file_directory)
