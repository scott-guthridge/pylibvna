#!/usr/bin/python3
import argparse
from libvna.data import NPData, PType
import numpy as np
import sys


try:
    argp = argparse.ArgumentParser()
    argp.add_argument('file1')
    argp.add_argument('file2')
    args = argp.parse_args()
except Exception as e:
    print(e)
    exit(1)

npd1 = NPData(filename=args.file1, ptype=PType.S)
npd2 = NPData(filename=args.file2, ptype=PType.S)
assert npd1.rows == npd2.rows
assert npd1.columns == npd2.columns
assert npd1.frequencies == npd2.frequencies
assert np.allclose(npd1.frequency_vector[...], npd2.frequency_vector[...],
                   rtol=1e-4, atol=1e-4)
assert np.allclose(npd1.data_array[...], npd2.data_array[...],
                   rtol=1e-4, atol=1e-4)
