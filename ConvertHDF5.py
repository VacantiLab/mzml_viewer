def Convert():
# This function converts an mzML file to a HDF5 file in a form readable by PlotMZML.py


    from pyteomics import mzml, auxiliary
    import numpy as np
    from pdb import set_trace
    import hdfdict
    import h5py

    time_mz_HDF5_file_directory = '/Users/nate/Desktop/temporary/time_mz.hdf5'
    time_intensity_HDF5_file_directory = '/Users/nate/Desktop/temporary/time_intensity.hdf5'
    mz_intensity_HDF5_file_directory = '/Users/nate/Desktop/temporary/mz_intensity.hdf5'
    mz_TimePointIndex_dict_HDF5_file_directory = '/Users/nate/Desktop/temporary/mz_TimePointIndex.hdf5'

    # Specify the mzml file
    #mzml_file_directory = '/Users/nate/Dropbox/Research/Vacanti_Laboratory/projects/PolyMID/correction_program/references/mzml_files/Tryptophan/YY2019060827.mzML'
    mzml_file_directory = '/Users/nate/Desktop/temporary/YY2019060827.mzML'


    # Read the mzML file as an iterable object
    MZML = mzml.read(mzml_file_directory,dtype=dict)


    # Initialize dictionaries with time points as the keys
    tics = np.array([])
    times = np.array([])


    # Initialize the total ion chromatograph array
    tic_array = np.array([])

    n_total_mzs = 0

    # Iterate through the scnas of the iterable object of the mzML file
    with h5py.File(time_mz_HDF5_file_directory, "w") as hf1, h5py.File(time_intensity_HDF5_file_directory, "w") as hf2:
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
            hf1.create_dataset(str(time_point), data=mzs)
            hf2.create_dataset(str(time_point), data=intensities)
            tics = np.append(tics,tic)
            times = np.append(times,time_point)

            # Find another way to get unique MZs - from dictionary already made! This takes too much time
            # track the unique mz values across all scans
            #unique_mzs = np.append(unique_mzs,mzs)
            #unique_mzs = np.unique(unique_mzs)

            i = i+1
            if i%200 == 0:
                print('Reading mzML scan number: ' + str(i))

            if i >= 200:
                break

        hf1.create_dataset('tic', data=tics)
        hf1.create_dataset('time', data=times)

    n_scans=i



    print('writing files')

    # create the files
    #   "a" specifies read write permissions and to create the file if it does not exist


    #hdfdict.dump(time_mz_dict,time_mz_HDF5_file_directory)
    #hdfdict.dump(time_intensity_dict,time_intensity_HDF5_file_directory)
