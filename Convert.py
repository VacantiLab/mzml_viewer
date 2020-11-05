def Convert(mzLower,mzUpper):
# Inputs:
    # Set the lower and upper bounds on the mzs for which a time trace will be available
    #     This is necessary because the number of mzs in high resolution data is over 1 million and the script will take too long to create the necessary dictionary

    from pyteomics import mzml, auxiliary
    import numpy as np
    import pickle
    from pdb import set_trace

    #mzml_file_directory = '/Users/nate/Dropbox/Research/Vacanti_Laboratory/mzml_files/QE1_QC_HeLa_20200225_r2.mzML'
    mzml_file_directory = '/Users/nate/Dropbox/Research/Vacanti_Laboratory/projects/PolyMID/correction_program/references/mzml_files/2018_1016_10.mzML'

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
        if i%100 == 0:
            print(i)

    n_scans=i

    output = [time_array,unique_mzs,time_mz_dict,time_intensity_dict,time_mz_dict_limited,time_intensity_dict_limited,tic_array,n_scans]
    StorageFile = '/Users/nate/Desktop/temporary/2018_1016_10_file1.p'
    with open(StorageFile, 'wb') as PickleFile:
        pickle.dump(output,PickleFile)


    #return(time_array,unique_mzs,time_mz_dict,time_intensity_dict,time_mz_dict_limited,time_intensity_dict_limited,tic_array,n_scans)
