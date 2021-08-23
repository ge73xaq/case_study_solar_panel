"""
Authors: Daniela Schlager (daniela.schlager@online.de), Yichen Lou (yichen.lw99@gmail.com), Jessica Reichelt
(j.s.reichelt@gmail.com), Lukas Dreier (lukas.dreier@gmx.de)

Contains the preprocessing class
"""
from shapely.geometry import Polygon
import geopandas as gpd
import elevation
import os
from osgeo import gdal
import matplotlib.pyplot as plt


class PreProcessor:
    """
    This class preprocesses the data. One part is responsible for transforming the given vertices into a plain polgon.
    The other part is responsible for the generation of a tiff file to get topographic data and the plotting of the
    topographic data.
    """
    def __init__(self, coordinates, path_for_tiff, obstacles):
        """
        :param coordinates: String copied from the API with information on the corners of our region
        :param path_for_tiff: Path where tiff file should be saved
        :param obstacles: List of strings copied from the API with information on the corners of different obstacles
        """
        self.coordinates = coordinates
        self.path_for_tiff = path_for_tiff
        try:
            if len(obstacles) == 0:
                self.obstacles = []
            else:
                self.obstacles = obstacles
        except:
            self.obstacles = []

    def manipulate_string(self, coordinates):
        """
        This function manipulates the string we get from the internet and transform it to a list of coordinates which
        can be used to transform it to a proper polygon.
        :param coordinates: String copied from the API with information on the corners of a polygon
        :return: List of coordinate tuples
        """
        coord_split = coordinates.split(',')
        list_of_coordinates = []
        for i in coord_split:
            corner = i.split()
            list_of_coordinates.append((float(corner[0]), float(corner[1])))

        return list_of_coordinates

    def create_polygon_4326(self, list_of_coordinates):
        """
        This function created a polygon in ESPG4326 representation. This representation is used by Open Street Map to
        represent their coordinates.
        :param list_of_coordinates: List of coordinates in tuples with (lat, long)
        :return: Geoseries which contains the polygon in ESPG4326 representation
        """
        polygon = Polygon(list_of_coordinates)
        geoseries = gpd.GeoSeries(polygon)
        geoseries.crs = {'init': 'epsg:4326'}

        return geoseries

    def transform_polygon_to_31467(self, geoseries):
        """
        This function transforms the polygon from EPSG4326 to ESPG31467. This is the so called Gauß-Krüger transformation
        of longitude and latitude data to get distances. It is used in parts of Germany. Since the differences are not
        very big we used it as default configuration. In the following article you can adjust the projection based on
        the region: https://spatialreference.org/ref/epsg/dhdn-gauss-kruger-zone-3/
        :param geoseries: Geoseries which stores the polygon
        :return: Shapely polygon
        """
        geoseries_proj = geoseries.to_crs(epsg=31467)
        polygon = geoseries_proj.iloc[0]

        return polygon

    def generate_tiff_file(self):
        """
        This function generates a tiff file and saves it to predetermined path.
        :return: Saved tiff file in predetermined path
        """
        output = os.getcwd() + self.path_for_tiff
        elevation.clip(bounds=self.geoseries.bounds.iloc[0], output=output, product='SRTM1')

    def add_topographic_information(self):
        """
        This function extracts the topographic information from the tiff file. The topographic information is r
        epresented through the matrix elevation. Rows and columns discretize the region and the respective value is the
        altitude. The resolution is about 20 metres.
        :return: Extraction of topographic information
        """
        gdal.UseExceptions()
        output = os.getcwd() + self.path_for_tiff
        ds = gdal.Open(output)
        band = ds.GetRasterBand(1)
        elevation = band.ReadAsArray()

        self.ds = ds
        self.band = band
        self.elevation = elevation

    def plot_topographic_information(self):
        """
        This function plots the topographic information and gives a feeling how the topology of the region looks like.
        :return: Plot of topographic information
        """
        nrows, ncols = self.elevation.shape
        x0, dx, dxdy, y0, dydx, dy = self.ds.GetGeoTransform()

        x1 = x0 + dx * ncols
        y1 = y0 + dy * nrows

        plt.imshow(self.elevation, cmap='gist_earth', extent=[x0, x1, y1, y0])
        plt.show()
        plt.close()

    def get_polygon(self):
        """
        This function generates needed plain polygon.
        :return: Polygon in 3879 representation
        """
        list_of_coordinates = self.manipulate_string(self.coordinates)
        geoseries = self.create_polygon_4326(list_of_coordinates)
        polygon = self.transform_polygon_to_31467(geoseries)

        self.geoseries = geoseries
        self.polygon = polygon
        self.list_of_coordinates = list_of_coordinates

    def get_obstacles(self):
        """
        This function takes the list of obstacles and creates a shapely polygon for every obstacle such that it can be
        used in further calculations.
        :return: List of shapely polygons representing obstacles in the region
        """
        list_of_obstacles = []
        try:
            for i in self.obstacles.keys():
                obstacle_coords = self.manipulate_string(self.obstacles[i])
                geoseries = self.create_polygon_4326(obstacle_coords)
                polygon = self.transform_polygon_to_3879(geoseries)
                list_of_obstacles.append(polygon)

            return list_of_obstacles
        except:
            return []
