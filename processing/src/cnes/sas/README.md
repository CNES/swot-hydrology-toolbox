# LOCNES - Lake Observation Cover aNd Extent from Swot
This project is dedicated to the production of lake shapefile from pixel cloud (L2_HR_PIXC products -> L2_HR_LakeSP products). It is composed by two main processes:
* LakeTile processing: uses a tile of PixC and produces L2_HR_LakeTile products. Only lakes entirely within a tile are processed. Pixels of lakes at the edge of two tiles are set apart in dedicated files (L2_HR_LakeTile_Edge).
* LakeSP processing: uses outputs of LakeTile processing (L2_HR_LakeTile.shp and L2_HR_LakeTile_Edge files) and gather all lakes within a single continent-pass in a L2_HR_LakeSP product.
* LakeAvg processing: uses outputs of LakeSP processing (L2_HR_LakeSP_Prior files) over a CBB basin and a cycle, and produces LakeAvg products.

Specific documentation is available in PGE/lake_tile, PGE/lake_sp and PGE/lake_avg packages.

# Workflow

![Alt text](20180312_AlgosLakes_v0_12.png?raw=true "Workflow diagram")