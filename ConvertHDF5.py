def Convert(mzLower,mzUpper):
# Inputs:
#   mzLower: This is the lower bound of the mz range that can be plotted vs. time
#   mzUpper: This is the upper bound of the mz range that can be plotted vs. time
# This function opens a mzML file that is hard-coded for now

    from pyteomics import mzml, auxiliary
    import numpy as np
    from pdb import set_trace
    import hdfdict

    HDF5_file_directory = '/Users/nate/Desktop/temporary/test.hdf5'

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
        time_mz_dict[str(time_point)] = mzs
        time_intensity_dict[str(time_point)] = intensities

        # track the unique mz values across all scans
        unique_mzs = np.append(unique_mzs,mzs)
        unique_mzs = np.unique(unique_mzs)

        i = i+1
        if i%100 == 0:
            print(i)

        #if i >= 200:
        #    break

    hdfdict.dump(time_mz_dict,HDF5_file_directory)
    set_trace()

    n_scans=i

    # save the desired outputs as a pickle file
    output = [time_array,unique_mzs,time_mz_dict,time_intensity_dict,time_mz_dict_limited,time_intensity_dict_limited,tic_array,n_scans,mzLower,mzUpper]

    # sort the array of unique mz values from smallest to largest
    unique_mzs = np.sort(unique_mzs)

    # Initialize a dictionary where the keys are the unique mz values
    #     Each entry will be an array of intensities
    #         Each position of the array is associated with the time point at the same position in the time_array
    print('initializing mz_intensity_dict')
    mz_intensity_dict = {key:np.array([]) for key in unique_mzs}
    print('initializing mz_TimePointIndex_dict')
    mz_TimePointIndex_dict = {key:np.array([]) for key in unique_mzs}

    # Fill the dictionary where unique mz values are the keys
    # Iterate over the time points
    time_point_index = 0
    for time_point in time_array:
        # Iterate over the mz values scanned at the current time point
        #     This is performed in the limited time_mz_dict because it was created to only hold the mz values of interest
        #       All mz values cannot be considered because there are too many
        mz_index = 0
        for mz in time_mz_dict[time_point]:
            # Fill the dictionaries with unique mz values as the keys
            #     One dictionary will hold intensities for every mz value
            #     The other dictionary will hold the corresponding indices of the time_array for those intensities
            CurrentIntensity = time_intensity_dict[time_point][mz_index]
            # Make an array holding intensities corresponding to the mz value that specifies the key to this dictionary
            mz_intensity_dict[mz] = np.append(mz_intensity_dict[mz],CurrentIntensity)
            # Make an array holding time point indices corresponding to the above intensities
            mz_TimePointIndex_dict[mz] = np.append(mz_TimePointIndex_dict[mz],time_point_index)
            mz_index = mz_index + 1

        # Iterate to the next time point
        time_point_index = time_point_index + 1
        if time_point_index%100 == 0:
            print('time point index = '+ str(time_point_index))


    #return(time_array,unique_mzs,time_mz_dict,time_intensity_dict,time_mz_dict_limited,time_intensity_dict_limited,tic_array,n_scans)
