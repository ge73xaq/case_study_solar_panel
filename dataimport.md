### Data Import

## Conversion from .oms to .geojson
For the conversion you have to install the gdal package which causes some trouble. Click on the link to get a
description how to import it. https://stackoverflow.com/questions/32066828/install-gdal-in-virtualenvwrapper-environment

In this link you can find a detailed description how you can convert it. https://www.compose.com/articles/how-to-transform-and-use-openstreetmap-data-into-geojson-using-gdal/

Run this in the command line:
`ogr2ogr -f GEOjson test.geojson test.osm multipolygons`

For download of DEM information see the following link:
http://geologyandpython.com/dem-processing.html
Attention: For elevation you have to install all required packages (especially curl ist important)

For elevation stuff:
https://stackoverflow.com/questions/24956653/read-elevation-using-gdal-python-from-geotiff

You can run the dataflow.py file on the command line with: `python3 dataflow.py -lat 10.9776 11.0550 10.9776 11.0550 -long 47.4078 47.4078 47.450514 47.450514`

Good documentation to read something on handling projections of the data
http://geopandas.org/projections.html

Information on splitting polygon:
https://gis.stackexchange.com/questions/232771/splitting-polygon-by-linestring-in-geodjango

Information on plotting
https://sgillies.net/2010/04/07/descartes.html
https://deparkes.co.uk/2015/03/11/how-to-plot-polygons-in-python/

## Projection stuff
We always used 3879 which is the projection for Finland.
https://automating-gis-processes.github.io/2016/Lesson3-projections.html
https://de.wikipedia.org/wiki/Gau%C3%9F-Kr%C3%BCger-Koordinatensystem
https://spatialreference.org/ref/epsg/dhdn-gauss-kruger-zone-3/
