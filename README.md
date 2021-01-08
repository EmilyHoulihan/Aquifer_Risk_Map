# Aquifer_Risk_Map
The data processing steps used to create the aquifer risk map. 
Map available at: https://gispublic.waterboards.ca.gov/portal/apps/webappviewer/index.html?id=17825b2b791d4004b547d316af7ac5cb

The data processing documentation is available in the "Methodology" folder. "waterquality_datadownload.R"
handles the data collection and data standardization. "aquiferriskmap_method.R" handles
the de-clustering and aggregation of the water quality data into census block group "risk" scores.

The final tables containing the risk scores in .csv format are available in the "Tables" folder. Data is available for census block groups, square mile sections (PLSS) and well points.

If you would like to re-run this data processing, or are interested in any of the dependant reference files
(shapefile boundaries, chemical comparison concentration tables, etc.), please contact Emily Houlihan 
(Emily.Houlihan@Waterboards.ca.gov). These reference files are too large to host on GitHub.
The R files are provided here to document the methodology steps.
