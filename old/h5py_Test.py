import hdfdict
import h5py
import numpy as np

filename = '/Users/nate/Desktop/temporary/test.hdf5'
d = {
    'a': np.array([1,2,3]),
    'b': np.array([2,3,4]),
}

hdfdict.dump(d,filename)
res = hdfdict.load(filename)

# with h5py.File(filename, "w") as f:
#     dset = f.create_dataset("mydataset", data)
#
# f2 = h5py.File(filename, 'r')
