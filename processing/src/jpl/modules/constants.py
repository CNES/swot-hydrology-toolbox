"""
Module for swot constants
"""

PIXC_CLASSES = {
    'invalid': -1,
    'land': 1,
    'land_near_water': 2,
    'water_near_land': 3,
    'open_water': 4,
    'dark_water': 5,
    'low_coh_water_near_land': 6,
    'low_coh_water': 7,
    'land_near_dark_water': 22,  # legacy/depreciated
    'dark_water_near_land': 23,  # legacy/depreciated
    'dark_water_legacy': 24      # legacy/depreciated
    }

# TODO: figure out how we want to aggregate low-coherence water near shore
#       When coherence is low, but it is bright it is likely due to layover
#       in which case the water fraction may be bad, but so is the total area
#       since there is water from two areas on the ground...
#       for now just aggregate the pixel area so it will only map to the one
#       ground-plane area that it is assigned to, while the other laid-over
#       water will just be missing 
AGG_CLASSES = {
    'interior_water_klasses':[
        PIXC_CLASSES['open_water'],
        PIXC_CLASSES['low_coh_water'],
        PIXC_CLASSES['low_coh_water_near_land'], # TODO: do we want this?
        ],
    'water_edge_klasses':[
        PIXC_CLASSES['water_near_land'],
        #PIXC_CLASSES['low_coh_water_near_land'], # TODO: do we want this?
        ],
    'land_edge_klasses':[
        PIXC_CLASSES['land_near_water'],
        ],
    'dark_water_klasses':[
        PIXC_CLASSES['dark_water'],
        PIXC_CLASSES['land_near_dark_water'],
        PIXC_CLASSES['dark_water_near_land'],
        PIXC_CLASSES['dark_water_legacy'],
        ]
    }

GDEM_PIXC_CLASSES = {
    'open_water': 4, 'dark_water': 5, #
    'open_water_lake': 34, 'dark_water_lake': 5} # 54
