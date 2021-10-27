# -*- coding: utf8 -*-
'''
Ply writer

Copyright (c) 2018 CNES. All rights reserved.
'''

from plyfile import PlyData, PlyElement, _data_type_reverse, _lookup_type, PlyProperty
import numpy as np
import utm
import geopandas as gpd
from shapely.geometry import Polygon
from typing import List

class MyPlyElement(PlyElement):
    '''
    Derived from PlyElement
    '''
    def __init__(self, name, properties, count, comments=[]):
        super().__init__(name, properties, count, comments)
        
    @staticmethod
    def describe(data, name, len_types={}, val_types={},
                 comments=[]):
        '''
        Construct a PlyElement from an array's metadata.

        len_types and val_types can be given as mappings from list
        property names to type strings (like 'u1', 'f4', etc., or
        'int8', 'float32', etc.). These can be used to define the length
        and value types of list properties.  List property lengths
        always default to type 'u1' (8-bit unsigned integer), and value
        types default to 'i4' (32-bit integer).

        '''
        if not isinstance(data, np.ndarray):
            raise TypeError("only numpy arrays are supported")

        if len(data.shape) != 1:
            raise ValueError("only one-dimensional arrays are "
                             "supported")

        count = len(data)

        properties = []
        descr = data.dtype.descr

        for t in descr:
            if not isinstance(t[1], str):
                raise ValueError("nested records not supported")

            if not t[0]:
                raise ValueError("field with empty name")

            if len(t) != 2 or t[1][1] == 'O':
                # non-scalar field, which corresponds to a list
                # property in PLY.

                if t[1][1] == 'O':
                    if len(t) != 2:
                        raise ValueError("non-scalar object fields not "
                                         "supported")

                len_str = _data_type_reverse[len_types.get(t[0], 'u1')]
                if t[1][1] == 'O':
                    val_type = val_types.get(t[0], 'i4')
                    val_str = _lookup_type(val_type)
                else:
                    val_str = _lookup_type(t[1][1:])

                prop = PlyListProperty(t[0], len_str, val_str)
            else:
                val_str = _lookup_type(t[1][1:])
                prop = PlyProperty(t[0], val_str)

            properties.append(prop)

        elt = MyPlyElement(name, properties, count, comments)
        elt.data = data

        return elt


    def _write_txt(self, stream):
        '''
        Save a PLY element to an ASCII-format PLY file.  The element may
        contain list properties.

        '''
        for rec in self.data:
            fields = []
            for prop in self.properties:
                fields.extend(prop._to_fields(rec[prop.name]))

            np.savetxt(stream, [fields], '%.6f', newline='\n')
            
    @property
    def header(self):
        '''
        Format this element's metadata as it would appear in a PLY
        header.
        '''
        # Some information is lost here, since all comments are placed
        # between the 'element' line and the first property definition.
        lines = []
        
        for c in self.comments:
            lines.append('comment ' + c)
        lines.append('element %s %d' % (self.name, self.count))

        lines.extend(list(map(str, self.properties)))

        return '\n'.join(lines)


            
def gdf_to_file(filename: str, points: gpd.GeoDataFrame, mode: str ="text"):
    '''
    Write to Ply file (UTM coordinates)

    :param filename: File path
    :param points: GeoDataFrame ("latitude","longitude","height","x","y","z")
    :type IN_inputvecfiles: GeoPandas Dataframe
    ''' 
    # Extract informations for header
    nb = points.shape[0]
    zone_number = utm.latlon_to_zone_number(points.iloc[0]['latitude'],
                                                              points.iloc[0]['longitude'])
    pos = 'N' if points.iloc[0]['latitude'] > 0 else 'S'
    # Utm coordinates to numpy array
    vertex = np.array([tuple(value) for value in points[['x','y','z', 'elevation']].values],
                 dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4'), ('elevation', 'f4')])
    el = MyPlyElement.describe(vertex, 'vertex',
                    comments=['projection: UTM {}{}'.format(zone_number, pos)])
    if mode == "text":
        PlyData([el], text=True).write(filename)
    elif mode == "binary":
        PlyData([el]).write(filename)
    else:
        raise Exception("Mode unknown")


def polygons_to_file(filename: str, polygons: List[Polygon], mode: str ="text"):
    '''
    Write to Ply file (UTM coordinates)

    :param filename: File path
    :param points: GeoDataFrame ("latitude","longitude","height","x","y","z")
    :type IN_inputvecfiles: GeoPandas Dataframe
    ''' 
    raise Exception("Not implemented")
