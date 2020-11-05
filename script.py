import pickle
StorageFile = '/Users/nate/Desktop/temporary/2018_1016_10.p'
with open(StorageFile, 'rb') as PickleFile:
    input = pickle.load(PickleFile)
