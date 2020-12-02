def Convert():
# Inputs:
#   mzLower: This is the lower bound of the mz range that can be plotted vs. time
#   mzUpper: This is the upper bound of the mz range that can be plotted vs. time
# This function opens a mzML file that is hard-coded for now

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

    n_total_mzs = 0

    # Iterate through the scnas of the iterable object of the mzML file
    i=0
    for key in MZML:
        # get the current time associated with the current scan
        time_point = float(key['scanList']['scan'][0]['scan start time'])
        time_array = np.append(time_array,time_point)

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
        time_intensity_dict[str(time_point)] = intensities

        # Find another way to get unique MZs - from dictionary already made! This takes too much time
        # track the unique mz values across all scans
        #unique_mzs = np.append(unique_mzs,mzs)
        #unique_mzs = np.unique(unique_mzs)

        i = i+1
        if i%100 == 0:
            print(i)

        #if i >= 200:
        #    break

    hdfdict.dump(time_mz_dict,time_mz_HDF5_file_directory)
    hdfdict.dump(time_intensity_dict,time_intensity_HDF5_file_directory)


    AllMZs = np.zeros(n_total_mzs)
    n_MZsTallied = 0
    i = 0

    # Read the mzML file as an iterable object
    MZML = mzml.read(mzml_file_directory,dtype=dict)
    for key in MZML:
        # get all of the mz values scanned in the current scan
        mzs = np.array(key['m/z array'],dtype='float')
        n_mzs_CurrentScan = len(mzs)
        AllMZs[n_MZsTallied:n_MZsTallied+n_mzs_CurrentScan] = mzs
        n_MZsTallied = n_MZsTallied+n_mzs_CurrentScan

        i = i+1
        if i%100 == 0:
            print('TotalMZ Scan: ' + str(i) + ', Total MZs Tallied: ' + str(n_MZsTallied))

    unique_mzs = np.unique(AllMZs)
    unique_mzs = np.sort(unique_mzs)


    n_scans=i

    # variables to save: time_mz_dict, time_intensity_dict, time_array, unique_mzs, tic_array, n_scans

    # Initialize a dictionary where the keys are the unique mz values
    #     Each entry will be an array of intensities
    #         Each position of the array is associated with the time point at the same position in the time_array
    print('initializing mz_intensity_dict')
    mz_intensity_dict = {str(key):np.array([]) for key in unique_mzs}
    print('initializing mz_TimePointIndex_dict')
    mz_TimePointIndex_dict = {str(key):np.array([]) for key in unique_mzs}

    # Fill the dictionary where unique mz values are the keys
    # Iterate over the time points
    time_point_index = 0
    for time_point in time_array:
        # Iterate over the mz values scanned at the current time point
        mz_index = 0
        for mz in time_mz_dict[str(time_point)]:
            # Fill the dictionaries with unique mz values as the keys
            #     One dictionary will hold intensities for every mz value
            #     The other dictionary will hold the corresponding indices of the time_array for those intensities
            CurrentIntensity = time_intensity_dict[str(time_point)][mz_index]
            # Make an array holding intensities corresponding to the mz value that specifies the key to this dictionary
            mz_intensity_dict[str(mz)] = np.append(mz_intensity_dict[str(mz)],CurrentIntensity)
            # Make an array holding time point indices corresponding to the above intensities
            mz_TimePointIndex_dict[str(mz)] = np.append(mz_TimePointIndex_dict[str(mz)],time_point_index)
            mz_index = mz_index + 1

        # Iterate to the next time point
        time_point_index = time_point_index + 1
        if time_point_index%100 == 0:
            print('time point index = '+ str(time_point_index))

    hdfdict.dump(mz_intensity_dict,mz_intensity_HDF5_file_directory)
    hdfdict.dump(mz_TimePointIndex_dict,mz_TimePointIndex_dict_HDF5_file_directory)
