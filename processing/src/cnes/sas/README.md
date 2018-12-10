# locnes
This project is dedicated to the production of lake shapefile from pixel clouds by (L2_HR_PIXC products -> L2_HR_LAKE_SP products). It is composed by two main processings :
* Lake_tile processing : use a tile of PixC and produces L2_HR_LAKE_TILE products. Only lakes entierly within a tile are processed. Lakes at the edge of two tiles are written to specific files L2_HR_LAKE_TILE_edgefile.
* Lake_sp processing : use outputs of lake_tile processing (L2_HR_LAKE_TILE and L2_HR_LAKE_TILE_edgefile) and gather all lakes within a L2_HR_LAKE_SP product (Single Pass).

Specific documentation is available in lake_tile and lake_sp projects.

# Workflow

![Alt text](20180312_AlgosLakes_v0_12.png?raw=true "Workflow diagram")