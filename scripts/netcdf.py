'''
module pixel_cloud.py
module author Capgemini
Copyright (c) 2018 CNES.All rights reserved
'''
from __future__ import print_function
import os
import netCDF4 as nc
import numpy as np

FILL_VALUE = -9.99e9

def set_variable(dataset, key, array, dimensions):
    '''Set the NetCDF variable, dealing with complex numbers.

    If array is complex, it is stored in the dataset with a third dimension,
    'depth', so that:
        variable.shape == (lines, pixels, 2)
        variable[:, :, 0] == array.real
        variable[:, :, 1] == array.imag
    '''
    # np.ma.mask_array has fill_value attr, else use FILL_VALUE
    fill_value = getattr(array, 'fill_value', FILL_VALUE)

    def _make_variable(key, data, dimensions):
        dataset.createVariable(key, data.dtype, dimensions, fill_value=fill_value)
        dataset[key][:] = data

    if 'complex' in array.dtype.name:
        # Add the depth dimension
        if 'depth' not in dataset.dimensions:
            dataset.createDimension('depth', 2)
        shape = array.shape
        n_bytes = int(array.dtype.itemsize/2)
        float_type = np.dtype('f'+str(n_bytes))
        if isinstance(array, np.ma.core.MaskedArray):
            # Somehow MaskedArray.view() doesn't work, so convert to a normal
            # array.
            array = array.filled()
        tmp = array.view(dtype=float_type).reshape(shape+(2,))
        _make_variable(key, tmp, dimensions+['depth'])
    else:
        _make_variable(key, array, dimensions)
