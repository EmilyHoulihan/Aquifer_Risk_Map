# Aquifer_Risk_Map
The data processing steps used to create the aquifer risk map. 
Map available at: https://gispublic.waterboards.ca.gov/portal/apps/webappviewer/index.html?id=17825b2b791d4004b547d316af7ac5cb

Methodology for the "2021 Aquifer Risk Map" is available in the "2021" folder, and includes:
1. "Aquifer Risk Map Methodology.docx" (methodology write-up)
2. "waterquality_datadownload.R" (water quality data download steps)
3. "waterquality_datadownload2.R" (water quality data download for tables downloaded after 3/1/2021 - table structure/headings have changed)
4. "aquiferriskmap_method.R" (de-clustering and aggregation of water quality data to section and block group scores)
5. "NA_allchems.xlsx" (list of chemicals and comparison concentrations used)

Methodology for the "2022 Aquifer Risk Map" is available in the "2022" folder, and includes:
1. "ARM_2022.R" (aquifer risk map data download and data processing steps)


If you would like to re-run this data processing, or are interested in any of the dependant reference files
(shapefile boundaries, etc.), please contact Emily Houlihan (Emily.Houlihan@Waterboards.ca.gov).
