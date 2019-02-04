# LOCNES - Lake Observation Cover aNd Extent from Swot
This project is dedicated to the production of lake shapefile from pixel cloud (L2_HR_PIXC products -> L2_HR_LakeSP products). It is composed by two main processes:
* LakeTile processing: use a tile of PixC and produces L2_HR_LakeTile products. Only lakes entirely within a tile are processed. Pixels of lakes at the edge of two tiles are set apart in dedicated files (L2_HR_LakeTile_edge).
* LakeSP processing: use outputs of LakeTile processing (L2_HR_LakeTile.shp and L2_HR_LakeTile_edge files) and gather all lakes within a single-pass in a L2_HR_LakeSP product.

Specific documentation is available in PGE/lake_tile and PGE/lake_sp packages.

# Workflow

![Alt text](20180312_AlgosLakes_v0_12.png?raw=true "Workflow diagram")