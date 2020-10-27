# Aquifer_Risk_Map
The data processing steps used to create the aquifer risk map. 
Map available at: https://gispublic.waterboards.ca.gov/portal/apps/webappviewer/index.html?id=b5488c0911314ef6a4796245dd584d0b

The data processing methodology is split into two R files. The first one ("waterquality_datadownload")
handles the data collection and data standardization. The second one ("aquiferriskmap_method") handles
the de-clustering and aggregation of the water quality data into census block group "risk" scores.

If you would like to re-run this data processing, or are interested in any of the dependant reference files
(shapefile boundaries, chemical comparison concentration tables, etc.), please contact Emily Houlihan 
(Emily.Houlihan@Waterboards.ca.gov). These reference files are too large to host on GitHub.
The R files are provided here to document the methodology steps.
