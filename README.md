# Aquifer_Risk_Map
The data processing steps used to create the aquifer risk map.

2021 map: https://gispublic.waterboards.ca.gov/portal/apps/webappviewer/index.html?id=d11cd558dd4945729ae4f222034bd9c9

2022 map: [archived link coming soon]

2023 map: https://gispublic.waterboards.ca.gov/portal/apps/webappviewer/index.html?id=17825b2b791d4004b547d316af7ac5cb [updating soon]

Methodology for the "2023 Aquifer Risk Map" is available in the "2023" folder, and includes:
1. "ARM_2023.R" (aquifer risk map data download and data processing steps)
2. "NA_allchems_updated.xlsx" (list of chemicals and comparison concentrations)

Methodology for the "2022 Aquifer Risk Map" is available in the "2022" folder, and includes:
1. "ARM_2022.R" (aquifer risk map data download and data processing steps)
3. "oswcr_filtering.R" (steps for identifying domestic well records)

Methodology for the "2021 Aquifer Risk Map" is available in the "2021" folder, and includes:
1. "Aquifer Risk Map Methodology.docx" (methodology write-up)
2. "waterquality_datadownload.R" (water quality data download steps)
3. "waterquality_datadownload2.R" (water quality data download for tables downloaded after 3/1/2021 - table structure/headings have changed)
4. "aquiferriskmap_method.R" (de-clustering and aggregation of water quality data to section and block group scores)
5. "NA_allchems.xlsx" (list of chemicals and comparison concentrations used)

If you would like to re-run this data processing, or are interested in any of the dependant reference files
(shapefile boundaries, etc.), please contact Emily Houlihan (Emily.Houlihan@Waterboards.ca.gov).
