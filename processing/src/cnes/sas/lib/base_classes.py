from __future__ import print_function
from __future__ import absolute_import
import collections
import os
import warnings
#import shapefile
import numpy as np
import netCDF4 as nc
import itertools

import lib.netcdf
import copy

FIELD_WARNING = '"{}" is not a valid field for "{}", ignoring'

# This probably should be handled somewhere else...
def get_filepath(filename, out_dir=os.path.curdir, suffix_replace='', suffix=''):

    out_dir = [os.path.abspath(os.path.dirname(filename)), out_dir][out_dir is not None]
    
    # Smash up the filename
    name, extension = os.path.splitext(filename.replace(suffix_replace, ''))
    out_name = os.path.join(out_dir, os.path.basename(name + suffix + extension))
    
    return out_name


class Variable(object):
    '''Wrapper for arrays with defined dimension names.'''
    def __init__(self, data, dimensions):
        self.data = data
        self.dimensions = dimensions


class PixelCloud(object):
    def __init__(self, dimensions=['record']):
        attributes = [
            'wavelength',
            'near_range',
            'range_spacing',
            'azimuth_spacing',
            'noise_power_left',
            'noise_power_right',
            # For PODAAC ingestion
            'start_time',
            'stop_time',
            'pass_number',
            'cycle_number',
            'tile_ref',
            # Swath bounding box
            'inner_first_lat',
            'inner_first_lon',
            'inner_last_lat',
            'inner_last_lon',
            'outer_first_lat',
            'outer_first_lon',
            'outer_last_lat',
            'outer_last_lon',
            'description',
            'nr_pixels',
            'nr_lines',
            'looks_to_efflooks', # ?
        ]
        # TODO: variables to handle bigger sensor file than pixc azimuth extent
        variables = [
            'azimuth_index',
            'range_index',
            'pixel_area',
        # rare
            'classification',
            'continuous_classification',
            'num_rare_looks',
        # medium
            'latitude',
            'longitude',
            'height',
            'cross_track',# can get look vector from lat/lon/height and illumination time 'look_unit_x','look_unit_y','look_unit_z',
            'illumination_time', # was sensor_s, time coresponding to s_profile, proxy for doppler
        ]
        self._attributes = collections.OrderedDict()
        self._dimensions = collections.OrderedDict()
        self._variables = collections.OrderedDict()
        self._groups = collections.OrderedDict()
        # Fill attributes, dimensions and variables atribute with None values
        self._attributes = dict(map(lambda attr : (attr, None), attributes))
        self._dimensions = dict(map(lambda dim : (dim, None), dimensions))
        self._variables = dict(map(lambda var : (var, Variable(None, dimensions)), variables))

    @property
    def attributes(self):
        dictionary = self._attributes.copy()
        dictionary.update(self._dimensions)
        return dictionary
    
    def copy_attributes(self, other_product):
        '''Copy the attributes of given product into this one'''
        for key, value in other_product.attributes.items():
            self[key] = value

    def _copy(self, new_product, with_variables=True):
        '''Copy all of self into new_product'''
        new_product._attributes = dict(map(lambda attr_items : (attr_items[0], attr_items[1]), self._attributes.items()))
        new_product._dimensions = dict(map(lambda dim_items : (dim_items[0], dim_items[1]), self._dimensions.items()))
        if with_variables:
            new_product._variables = dict(map(lambda var_items : (var_items[0], var_items[1]), self._variables.items()))

        return new_product
    
    def copy(self, **kwargs):
        '''Return a deep copy of self'''
        # Make a new blank object of this type
        new_product = type(self)(dimensions=self._dimensions)
        return self._copy(new_product, **kwargs)
    
    def __setattr__(self, key, item):
        if key in ['_attributes', '_dimensions', '_variables', '_groups']:
            # key is present in NetCDF
            super(PixelCloud, self).__setattr__(key, item)
        elif isinstance(item, PixelCloud):
            # Try to set a new group
            if key in self._groups:
                self._groups[key] = item
            else:
                # Group is not defined in the product
                warnings.warn(FIELD_WARNING.format(key, type(self).__name__))
        elif isinstance(item, np.ndarray):
            # Try to set a new variable
            if key in self._variables:
                variable = self._variables[key]
                assert len(item.shape) == len(variable.dimensions)
                for i, dimension in enumerate(variable.dimensions):
                    if self._dimensions[dimension] is None:
                        self._dimensions[dimension] = item.shape[i]
                    assert self._dimensions[dimension] == item.shape[i]
                variable.data = item
            else:
                warnings.warn(FIELD_WARNING.format(key, type(self).__name__))
        else:
        # Try to set a new attribute
            if key in self._attributes:
                self._attributes[key] = item
            else:
                # Attributes is not defined
                warnings.warn(FIELD_WARNING.format(key, type(self).__name__))

    def __getattr__(self, key):
        if key in self._attributes:
            return self._attributes[key]
        if key in self._variables:
            return self._variables[key].data
        if key in self._groups:
            return self._groups[key]
        raise AttributeError('{} not in product'.format(key))
    
    def __setitem__(self, key, item):
        return setattr(self, key, item)

    def __getitem__(self, key):
        return getattr(self, key)

    def to_ncfile(self, *args, **kwargs):
        '''Write self to a netCDF file.'''
        outfile = get_filepath(*args, **kwargs)
        dataset = nc.Dataset(outfile, 'w')
        self.to_dataset(dataset)
        dataset.sync()
        dataset.close()
        return outfile

    def to_dataset(self, dataset):
        '''Convert self data in NetCDF dataset/group.'''
        # Loop over groups
        for key, value in self._groups.items():
            netcdf_group = dataset.createGroup(key)
            # Recursively add the group members
            value.to_dataset(netcdf_group)
        # Loop over attributes, dimensions and variables
        for key_attr, value_attr, key_dim, value_dim, key_var, value_var in itertools.zip_longest(self._attributes.keys(), self._attributes.values(),\
                                                                        self._dimensions.keys(), self._dimensions.values(),\
                                                                        self._variables.keys(), self._variables.values(), fillvalue=None):
            if key_attr is not None:
                dataset.setncattr(key_attr, [value_attr, lib.netcdf.FILL_VALUE][value_attr is None])
            if key_dim is not None:
                dataset.createDimension(key_dim, [value_dim, lib.netcdf.FILL_VALUE][value_dim is None])
            if value_var is not None:
                lib.netcdf.set_variable(dataset, key_var, value_var.data, list(value_var.dimensions))

    @classmethod
    def from_ncfile(cls, filename, force=False, variables=None):
        '''Generate directly an output File from
        an input NetCDF file.                '''
        product = cls()
        dataset = nc.Dataset(filename, 'r')
        product.from_dataset(dataset, force, variables)
        return product

    def from_dataset(self, dataset, force=False, variables=None):
        '''Import NetCDF dataset/group into self'''
        # Loop over attributes
        for key in dataset.ncattrs():
            if force:
                # Force attributes definition
                self._attributes[key] = None
            # Try to find attribute in input file
            try:
                getattr(self, key)
            except AttributeError:
                pass
            else:
                setattr(self, key, dataset.getncattr(key))
        # Loop over keys from groups
        for key in dataset.groups.keys():
            if force:
                # Force groups key definition
                self._groups[key] = PixelCloud()
            # Try to find key in input file
            try:
                getattr(self, key)
            except AttributeError:
                pass
            # Recursively load NetCDF groups into self groups
            self[key].from_dataset(dataset[key], force, variables)
