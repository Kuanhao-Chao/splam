import h5py
filename = "datafile_test.h5"

with h5py.File(filename, "r") as f:
    # List all groups
    print(("Keys: %s" % list(f.keys())))
    a_group_key = list(f.keys())[3]
    print(("a_group_key: ", a_group_key))

    # Get the data
    data = list(f[a_group_key])
    print(("data: ", data))
    print(("data: ", len(data)))
