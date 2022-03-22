import h5py
import numpy as np
import pandas as pd

HDF5_File = '/Users/nate/Desktop/temporary/JohanssonProteome.hdf5'

quantities = np.array([[1,2,3],
                [4,5,6],
                [7,8,9]])
tumors = ['tumor1','tumor2','tumor3']
genes = ['gene1','gene2','gene3']

TumorProteomeDF = pd.DataFrame(data=quantities,index=genes,columns=tumors)

with h5py.File(HDF5_File, "w") as hdf5:
    for gene in genes:
        data_to_write = np.array(TumorProteomeDF.loc[gene,:])
        hdf5.create_dataset(str(gene), data=data_to_write)

    hdf5.create_dataset('tumors', data=tumors)


with h5py.File(HDF5_File, "r") as hdf5:
    tumor_list = np.array(hdf5.get('tumors'))
    gene1_quants = np.array(hdf5.get('gene1'))
