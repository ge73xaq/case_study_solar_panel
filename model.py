"""
Authors: Daniela Schlager (daniela.schlager@online.de), Yichen Lou (yichen.lw99@gmail.com), Jessica Reichelt
(j.s.reichelt@gmail.com), Lukas Dreier (lukas.dreier@gmx.de)

Contains the model classes; one for the shifting and another one for the clustering and cabeling
"""
from math import ceil
from shapely.geometry import LineString, box, Point
import shapely
import numpy as np
import matplotlib.pyplot as plt
import os
import imageio
import random
import math
from dataclasses import dataclass
from typing import List
from datetime import timedelta
from matplotlib import cm
import sys


@dataclass
class Result:
    description: str
    path_to_cfg_file: str
    path_to_tiff: str
    coordinates: List[float]
    number_of_placed_tables: int
    number_of_placed_half_tables: int
    number_of_placed_tables_in_total: int
    number_of_inverters: int
    length_of_cables: int
    coverage: float
    coverage_of_tables: float
    coverage_of_paths: float
    coverage_of_obstacles: float
    ratio_of_tables: float
    ratio_of_cable_length: float
    preprocessing_time: timedelta
    shifting_time: timedelta
    clustering_time: timedelta
    overall_time: timedelta

    def write_result(self):
        header = f'description: {self.description} \n' \
            f'path_to_cfg_file: {self.path_to_cfg_file} \n' \
            f'path_to_tiff: {self.path_to_tiff} \n' \
            f'coordinates: {self.coordinates} \n' \
            f'number_of_placed_tables: {self.number_of_placed_tables} \n' \
            f'number_of_placed_half_tables: {self.number_of_placed_half_tables} \n' \
            f'number_of_placed_tables_in_total: {self.number_of_placed_tables_in_total} \n' \
            f'number_of_inverters: {self.number_of_inverters} \n' \
            f'length_of_cables: {self.length_of_cables} \n' \
            f'coverage: {self.coverage} \n' \
            f'coverage_of_tables: {self.coverage_of_tables} \n' \
            f'coverage_of_paths: {self.coverage_of_paths} \n' \
            f'coverage_of_obstacles: {self.coverage_of_obstacles} \n' \
            f'ratio_of_tables: {self.ratio_of_tables} \n' \
            f'ratio_of_cable_length: {self.ratio_of_cable_length} \n' \
            f'preprocessing_time: {self.preprocessing_time} \n' \
            f'shifting_time: {self.shifting_time} \n' \
            f'clustering_time: {self.clustering_time} \n' \
            f'overall_time: {self.overall_time} \n\n\n'

        return header

class Shifter:
    """
    This class contains the whole shifting logic. The idea is to fix one maintenance path and place the tables and other
    maintenance paths according to the first maintenance path. Afterwards we shift the maintenance path and the
    smallest fitting rectangles in discretized shiftings. The shifting checks the possibilities how to place the
    maintenance paths and tables and saves the best one. Furthermore, the topographic information is taken into account
    when we calculate the distance between two rows. In regions where obstacles are no tables are placed.
    """
    def __init__(self, polygon, elevation, width_between_maint_paths, width_of_maint_path, width_table, length_table,
                 alpha, gamma, list_of_obstacles, shifting_stepsx, shifting_stepsy, path_for_plots, path_for_gif,
                 save=False, plot=False):
        """
        :param polygon: Shapely polygon of the region
        :param elevation: Matrix which contains the topographic information
        :param width_between_maint_paths: Width between two maintenance paths (PV tables can be placed here)
        :param width_of_maint_path: Width of one maintenance path (no PV tables can be placed here)
        :param width_table: Width of PV module table
        :param length_table: Length of PV module table
        :param alpha: Angle between PV module and zero line in DEG
        :param gamma: Angle between sun and zero line in DEG
        :param list_of_obstacles: List of obstacles represented as Shapely polygons
        :param shifting_stepsx: Number of shifting steps in x-direction (East-West)
        :param shifting_stepsy: Number of shifting steps in y-direction (North-South)
        :param path_for_plots: Relative path where to save plots
        :param path_for_gif: Relative path where to save gif
        :param save: True/False whether plots should be saved
        :param plot: True/False whether shifting steps in x direction should be plotted
        """
        self.polygon = polygon
        self.elevation = elevation
        if width_between_maint_paths % width_table == 0:
            self.width_between_maint_paths = width_between_maint_paths
        else:
            self.width_between_maint_paths = math.floor(width_between_maint_paths/width_table)*width_table
        self.width_of_maint_path = width_of_maint_path
        self.width_table = width_table
        self.length_table_unprojected = length_table
        self.length_table = length_table * math.cos(math.pi*alpha/180)
        self.list_of_obstacles = list_of_obstacles
        self.shifting_stepsx = shifting_stepsx
        self.shifting_stepsy = shifting_stepsy
        self.path_for_plots = path_for_plots
        self.path_for_gif = path_for_gif
        self.save = save
        self.plot = plot
        self.alpha = math.pi*alpha/180
        self.gamma = math.pi*gamma/180
        self.height = self.length_table_unprojected * math.sin(self.alpha)

    def get_lines_for_paths(self, shift_path=0):
        """
        This function places the maintenance paths based on the shifting input.
        :param shift_path: Distance of first maintenance path to left border
        :return: List of lines representing the maintenance paths
        """
        minx, miny, maxx, maxy = self.polygon.bounds
        dist = minx + shift_path
        counter = 0
        lines = {}
        while dist < maxx:
            counter += 1
            if counter % 2 == 1:
                lines[f'left_{ceil(counter / 2)}'] = LineString([(dist, maxy), (dist, miny)])
                dist += self.width_of_maint_path
            else:
                lines[f'right_{ceil(counter / 2)}'] = LineString([(dist, maxy), (dist, miny)])
                dist += self.width_between_maint_paths

        lines = list(lines.values())

        return lines

    def create_smallest_fitting_rectangles(self, list_lines, paths=False):
        """
        This function fits smallest fitting rectangle to region between two maintenance paths
        :param list_lines: List of lines representing the maintenance paths
        :param paths: True/False whether boxes representing maintenance paths should be plotted
        :return: List of smallest fitting rectangles as Shapely box
        """
        merge_lines = [self.polygon.boundary]
        for i in list_lines:
            merge_lines.append(i)
        merged_lines = shapely.ops.linemerge(merge_lines)
        border_lines = shapely.ops.unary_union(merged_lines)
        decomposition = shapely.ops.polygonize(border_lines)
        decomposition_list = []
        checker_list = []

        if not paths:
            for i in decomposition:
                len_list_lines = len(list_lines)
                maint_path = False
                for j in range(math.ceil(len_list_lines/2)):
                    maint_path = False
                    if (((i.bounds[0] == list_lines[2*j].bounds[0]) &
                         (i.bounds[2] == list_lines[min(2*j+1, len_list_lines-1)].bounds[0])) or
                            (i.bounds[0] == list_lines[len_list_lines-1].bounds[0]) & (len_list_lines%2 != 0)):
                        maint_path = True
                        break
                if not maint_path:
                    if (i.bounds[2] - i.bounds[0] == self.width_between_maint_paths) & (self.polygon.contains(i)):
                        decomposition_list.append(i)
                    else:
                        checker_list.append(i)

            for i in checker_list:
                check = True
                for j in decomposition_list:
                    if not j.disjoint(i):
                        check = False
                        break
                if check:
                    decomposition_list.append(i)

            return decomposition_list
        else:
            for i in decomposition:
                if i.bounds[2] - i.bounds[0] == self.width_of_maint_path:
                    decomposition_list.append(i)

            return decomposition_list

    def fit_tables_with_shifting(self, decomposition_list):
        """
        This function fits the PV tables to every smallest fitting rectangle. Furthermore, it shifts the smallest
        fitting rectangle in y-direction to check which position is best for a specific position in x-direction in terms
        of our discretization. The optimal distance regarding the shadowing of the PV tables is calculated to ensure
        that the output is optimal in terms of shadowing.
        :param decomposition_list: List of smallest fitting rectangles as Shapely boxes
        :return: List of placed PV tables, list of placed half PV tables, corresponding list of smallest fitting
        rectangles
        """
        list_of_tables = []
        list_of_half_tables = []
        for i in decomposition_list:
            bminx, bminy, bmaxx, bmaxy = i.bounds
            max_tabs = 0
            list_of_tabs = []
            list_of_half_tabs = []
            opt_shift = 0
            height_table = math.sin(self.alpha) * self.length_table_unprojected
            init_rowspace = height_table / math.tan(self.gamma)
            for j in np.arange(-self.length_table / 2 - init_rowspace, self.length_table / 2,
                               (self.length_table + init_rowspace) / self.shifting_stepsy):
                ls = []
                ls_half = []
                b = box(bminx, bminy + j, bmaxx, bmaxy + j)
                if bminx == self.polygon.bounds[0]:
                    counterx = bmaxx
                    countery = bmaxy + j
                    while (countery - self.length_table) >= bminy:
                        rowspace = max(self.calculate_rowspace(bminx, bmaxx, countery - self.length_table/2), 0)
                        if rowspace > 7*self.length_table:
                            print('Area may be to steep for placing PV tables')
                        while (counterx - self.width_table / 2) >= bminx:
                            table = box(counterx - self.width_table, countery - self.length_table, counterx, countery)
                            if self.polygon.contains(table) & i.contains(table):
                                in_obst = False
                                for obst in self.list_of_obstacles:
                                    if obst.intersects(table):
                                        in_obst = True
                                        break
                                if not in_obst:
                                    ls.append(table)
                            else:
                                table = box(counterx - self.width_table/2, countery - self.length_table, counterx,
                                            countery)
                                if self.polygon.contains(table) & i.contains(table):
                                    in_obst = False
                                    for obst in self.list_of_obstacles:
                                        if obst.intersects(table):
                                            in_obst = True
                                            break
                                    if not in_obst:
                                        ls_half.append(table)
                                else:
                                    table = box(counterx - self.width_table, countery - self.length_table,
                                                counterx - self.width_table/2, countery)
                                    if self.polygon.contains(table) & i.contains(table):
                                        in_obst = False
                                        for obst in self.list_of_obstacles:
                                            if obst.intersects(table):
                                                in_obst = True
                                                break
                                        if not in_obst:
                                            ls_half.append(table)
                            counterx -= self.width_table
                        countery = countery - self.length_table - rowspace
                        counterx = bmaxx
                else:
                    counterx = bminx
                    countery = bmaxy + j
                    while (countery - self.length_table) >= bminy:
                        rowspace = max(self.calculate_rowspace(bminx, bmaxx, countery - self.length_table/2), 0)
                        if rowspace > 7*self.length_table:
                            print('Area may be to steep for placing PV tables')
                        while (counterx + self.width_table / 2) <= bmaxx:
                            table = box(counterx, countery - self.length_table, counterx + self.width_table, countery)
                            if self.polygon.contains(table) & i.contains(table):
                                in_obst = False
                                for obst in self.list_of_obstacles:
                                    if obst.intersects(table):
                                        in_obst = True
                                        break
                                if not in_obst:
                                    ls.append(table)
                            else:
                                table = box(counterx, countery - self.length_table, counterx + self.width_table/2,
                                            countery)
                                if self.polygon.contains(table) & i.contains(table):
                                    in_obst = False
                                    for obst in self.list_of_obstacles:
                                        if obst.intersects(table):
                                            in_obst = True
                                            break
                                    if not in_obst:
                                        ls_half.append(table)
                                else:
                                    table = box(counterx + self.width_table/2, countery - self.length_table,
                                                counterx + self.width_table, countery)
                                    if self.polygon.contains(table) & i.contains(table):
                                        in_obst = False
                                        for obst in self.list_of_obstacles:
                                            if obst.intersects(table):
                                                in_obst = True
                                                break
                                        if not in_obst:
                                            ls_half.append(table)
                            counterx += self.width_table
                        countery = countery - self.length_table - rowspace
                        counterx = bminx

                if (len(ls) + len(ls_half)/2) > max_tabs:
                    list_of_tabs = ls
                    list_of_half_tabs = ls_half
                    max_tabs = len(ls) + len(ls_half)/2
                    opt_shift = b
            if max_tabs != 0:
                for k in list_of_tabs:
                    list_of_tables.append(k)
                for k in list_of_half_tabs:
                    list_of_half_tables.append(k)
                index = decomposition_list.index(i)
                decomposition_list[index] = opt_shift

        if len(list_of_half_tables) % 2 == 1:
            del list_of_half_tables[0]

        return list_of_tables, list_of_half_tables, decomposition_list

    def shifting_x_direction(self):
        """
        This function shifts the maintenance paths in x direction. It starts with the maintenance path aligned to the
        left border and shifts all maintenance paths successively shifting_stepsx times through the width between two
        maintenance paths.
        :return: List of maximal number of photovoltaic tables,  list of maximal number of half photovoltaic tables,
        corresponding list of smallest fitting rectangles, corresponding list of lines representing maintenance paths
        """
        max_number_of_tables = 0
        max_list_of_tables = []
        max_list_of_half_tables = []
        max_decomposition_list = []
        max_list_lines = []
        counter = 1
        min_coverage = 1
        max_coverage = 0
        for i in np.arange(0, self.width_between_maint_paths, self.width_between_maint_paths / self.shifting_stepsx):
            list_lines = self.get_lines_for_paths(i)
            decomposition_list = self.create_smallest_fitting_rectangles(list_lines)
            list_of_tables, list_of_half_tables, decomposition_list = self.fit_tables_with_shifting(decomposition_list)

            if (len(list_of_tables) + len(list_of_half_tables)/2) > max_number_of_tables:
                max_number_of_tables = len(list_of_tables) + len(list_of_half_tables)/2
                max_list_of_tables = list_of_tables
                max_list_of_half_tables = list_of_half_tables
                max_decomposition_list = decomposition_list
                max_list_lines = list_lines

            if self.plot & self.save:
                self.plot_region(decomposition_list, list_of_tables, list_of_half_tables, list_lines, counter)
                counter += 1
            elif self.plot:
                self.plot_region(decomposition_list, list_of_tables, list_of_half_tables, list_lines)

            coverage, coverage_tables, coverage_paths, coverage_obstacles = self.get_used_area(list_of_tables,
                                                                                               list_of_half_tables,
                                                                                               list_lines)
            if coverage > max_coverage:
                max_coverage = coverage
                max_coverage_tables = coverage_tables
                max_coverage_paths = coverage_paths
            if coverage < min_coverage:
                min_coverage = coverage

        #print(f'Maximum coverage for x-shifting: {round(max_coverage*100, 2)}')
        #print(f'Minimum coverage for x-shifting: {round(min_coverage*100, 2)}')

        #print(f'Coverage of tables: {round(max_coverage_tables*100, 2)}')
        #print(f'Coverage of paths: {round(max_coverage_paths*100, 2)}')

        return max_list_of_tables, max_list_of_half_tables, max_decomposition_list, max_list_lines

    def calculate_rowspace_alternative(self, minx, maxx, y):
        """
        This function calculates the space between to rows of PV tables for a given y point and a range of x points
        sucht that the shadowing is minimal. First the respective index in the elevation file is calculated. Afterwards
        the optimal space is calculated via some properties of sin, cos and tan.
        :param minx: Minimal x point
        :param maxx: Maximal x point
        :param y: Y point
        :return:
        """
        pminx, pminy, pmaxx, pmaxy = self.polygon.bounds
        pdiffx = pmaxx - pminx
        pdiffy = pmaxy - pminy
        diffx = maxx - minx
        nrows, ncols = self.elevation.shape
        index_y = math.floor((pmaxy - y) / (pdiffy / nrows))
        max_index_x = math.ceil((maxx - pminx) / (pdiffx / ncols))
        min_index_x = math.floor((minx - pminx) / (pdiffx / ncols))
        max_dist = 0
        for i in range(min_index_x, max_index_x):
            slope = (self.elevation[index_y, i] - self.elevation[min(index_y+1, nrows-1), i])/(pdiffy / nrows)
            beta = math.atan(slope)
            if abs(beta) > abs(self.alpha):
                sys.exit('Choose other alpha in cfg file')
            height_region = self.length_table * math.tan(beta)
            height_adjusted = self.height - height_region
            question_mark = height_adjusted * (math.sin(math.pi/2 - self.gamma) /
                                               math.sin(math.pi-(math.pi/2-self.gamma) - (math.pi/2-beta)))
            distance = math.sin(math.pi/2 - beta) * question_mark
            if distance >= max_dist:
                max_dist = distance

        return max_dist


    def calculate_rowspace(self, minx, maxx, y):
        """
        This function calculates the space between to rows of PV tables for a given y point and a range of x points
        sucht that the shadowing is minimal. First the respective index in the elevation file is calculated. Afterwards
        the optimal space is calculated via some properties of sin, cos and tan.
        :param minx: Minimal x point
        :param maxx: Maximal x point
        :param y: Y point
        :return:
        """
        pminx, pminy, pmaxx, pmaxy = self.polygon.bounds
        pdiffx = pmaxx - pminx
        pdiffy = pmaxy - pminy
        diffx = maxx - minx
        nrows, ncols = self.elevation.shape
        index_y = math.floor((pmaxy - y) / (pdiffy / nrows))
        max_index_x = math.ceil((maxx - pminx) / (pdiffx / ncols))
        min_index_x = math.floor((minx - pminx) / (pdiffx / ncols))
        max_dist = 0
        for i in range(min_index_x, max_index_x):
            slope = (self.elevation[index_y, i] - self.elevation[min(index_y+1, nrows-1), i])/(pdiffy / nrows)
            beta = math.atan(slope)
            if beta > self.alpha:
                sys.exit('Choose other alpha in cfg file')
            height_table = math.sin(self.alpha - beta) * self.length_table_unprojected
            length_table_projected = math.cos (self.alpha - beta) * self.length_table_unprojected
            help = height_table / math.tan(beta + self.gamma)
            distance_sum = math.cos(beta) * (length_table_projected + help)
            distance1 = math.cos(self.alpha) * self.length_table_unprojected
            distance = distance_sum - distance1

            if distance >= max_dist:
                max_dist = distance

        return max_dist

    def get_used_area(self, list_of_tables, list_of_half_tables, list_lines):
        """
        This function calculates the area which is used. It should give an overview on how good the shifting works and
        information on the area which is covered.
        :param list_of_tables: List of placed PV tables
        :param list_of_half_tables: List of placed half PV tables
        :param list_lines: List of lines
        :return:
        """
        area_polygon = self.polygon.area
        area_tables = 0
        area_half_tables = 0

        if len(list_of_tables) > 0:
            area_tables = len(list_of_tables)*list_of_tables[0].area
        if len(list_of_half_tables) > 0:
            area_half_tables = len(list_of_half_tables)*list_of_half_tables[0].area
        area_tables_combined = area_tables + area_half_tables

        area_paths = 0
        decomposition_list = self.create_smallest_fitting_rectangles(list_lines, True)
        for i in decomposition_list:
            area_paths = area_paths + i.area

        area_obstacles = 0
        for i in self.list_of_obstacles:
            area_obstacles = area_obstacles + i.area

        coverage = (area_tables_combined + area_paths)/area_polygon
        coverage_tables = area_tables_combined/area_polygon
        coverage_paths = area_paths/area_polygon
        coverage_obstacles = area_obstacles/area_polygon

        return coverage, coverage_tables, coverage_paths, coverage_obstacles

    def get_ratio_of_placed_tables(self, list_of_tables, list_of_half_tables):
        """
        This function calculates the ratio of placed tables. First one calculates the number of PV tables which could
        be placed in theory. You divide the number of placed tables by this number and get a ratio for the placed
        tables.
        :param list_of_tables: List of placed tables
        :param list_of_half_tables: List of placed half tables
        :return:
        """
        area_polygon = self.polygon.area
        area_table = list_of_tables[0].area
        possible_tables = math.floor(area_polygon/area_table)
        number_of_placed_tables = len(list_of_tables) + len(list_of_half_tables)
        ratio = number_of_placed_tables/possible_tables

        return ratio

    def plot_region(self, decomposition_list, list_of_tables, list_of_half_tables, lines, save=''):
        """
        This function plots the area with smallest fitting rectangles, placed PV tables, obstacles and maintenance
        paths.
        :param decomposition_list: List of smallest fitting rectangles
        :param list_of_tables: List of placed PV tables
        :param list_of_half_tables: List of placed half PV tables
        :param lines: List of lines representing maintenance paths
        :param save: Relativ path where to save plot
        :return: Plot of area with smallest fitting rectangles, placed PV tables (full and half), obstacles and
        maintenance paths
        """
        minx, miny, maxx, maxy = self.polygon.bounds
        distx = maxx - minx
        disty = maxy - miny
        x, y = self.polygon.exterior.xy
        plt.plot(x, y)

        for i in decomposition_list:
            x, y = box(i.bounds[0], i.bounds[1], i.bounds[2], i.bounds[3]).exterior.xy
            plt.plot(x, y, color='green', lw=2)

        for i in lines:
            x1, y1 = i.xy
            plt.plot(x1, y1, color='red')

        for i in list_of_tables:
            x, y = i.exterior.xy
            plt.plot(x, y, color='black', lw=.5)

        for i in list_of_half_tables:
            x, y = i.exterior.xy
            plt.plot(x, y, color='black', lw=.5)

        for i in self.list_of_obstacles:
            x, y = i.exterior.xy
            plt.plot(x, y, color='black')

        plt.xlim([minx - 0.05 * distx, maxx + 0.05 * distx])
        plt.ylim([miny - 0.05 * disty, maxy + 0.05 * disty])
        plt.suptitle(f'Number of tables: {int(len(list_of_tables) + len(list_of_half_tables)/2)}')
        plt.axis('off')
        if save != '':
            plt.savefig(f'{os.getcwd() + self.path_for_plots}plot_{save}.png')
        plt.show()
        plt.close()

    def get_filenames_for_gif(self, include_final=False):
        """
        This function gets the filenames which are needed to generate a gif.
        :param include_final: True/False whether to include the final plot with converters and cables
        :return: List of filenames
        """
        filenames = []
        for i in range(1, self.shifting_stepsx+1):
            filenames.append(f'{os.getcwd() + self.path_for_plots}plot_{i}.png')

        if include_final:
            filenames.append(f'{os.getcwd() + self.path_for_plots}plot_final.png')

        return filenames

    def convert_images_to_gif(self, filenames):
        """
        This function converts several images to a gif.
        :param filenames: List of filenames
        :return: Gif
        """
        images = []

        for i in filenames:
            images.append(imageio.imread(i))

        imageio.mimsave(os.getcwd() + self.path_for_gif + 'shifting.gif', images)

    def create_gif(self, include_final=False):
        """
        This function creates a gif file.
        :param include_final: True/False whether to include the final plot with converters and cables
        :return: Gif
        """
        filenames = self.get_filenames_for_gif(include_final)
        self.convert_images_to_gif(filenames)


class KMeans:
    """
    This class contains the logic for clustering the inverters and assigning the cables to the corresponding inverters.
    The general idea is to use a modification of the k-means-algorithm to assign the tables to the inverters. A few
    modifications have to be done to fulfill several properties, e.g. the capacity of the inverters.
    """
    def __init__(self, polygon, list_of_obstacles, list_of_tables, list_of_half_tables, capacity_inverter,
                 number_of_inverters, width_table, length_table, decomposition_list, maintenance_paths, list_lines,
                 path_for_plots, path_for_gif, save=False, plot=False):
        """
        :param polygon: Shapely polygon of the region
        :param list_of_obstacles: List of obstacles
        :param list_of_tables: List of placed PV tables
        :param list_of_half_tables: List of placed half PV tables
        :param capacity_inverter: Maximum capacity of inverters
        :param number_of_inverters: Number of tables assigned to one inverter
        :param decomposition_list: List of smallest fitting rectangles as Shapely boxes
        :param list_lines: List of lines representing maintenance paths
        :param path_for_plots: Relative paths where plots should be saved
        """
        self.polygon = polygon
        self.list_of_obstacles = list_of_obstacles
        self.list_of_tables = list_of_tables
        self.list_of_half_tables = list_of_half_tables
        self.capacity_inverter = capacity_inverter
        self.number_of_inverters = number_of_inverters
        self.width_table = width_table
        self.length_table = length_table
        self.decomposition_list = decomposition_list
        self.maintenance_paths = maintenance_paths
        self.list_lines = list_lines
        self.path_for_plots = path_for_plots
        self.path_for_gif = path_for_gif
        self.save = save
        self.plot = plot

    def get_mid_points_of_tables(self):
        """
        This functions calculates the mid points of the placed PV tables.
        :return: List of tuples representing the mid points of the placed PV tables
        """
        mid_points_of_tables = []
        mid_points_of_half_tables = []

        for i in self.list_of_tables:
            bminx, bminy, bmaxx, bmaxy = i.bounds
            mid_points_of_tables.append(((bmaxx-bminx)*0.5+bminx, (bmaxy-bminy)*0.5+bminy))

        for i in self.list_of_half_tables:
            bminx, bminy, bmaxx, bmaxy = i.bounds
            mid_points_of_half_tables.append(((bmaxx-bminx)*0.5+bminx, (bmaxy-bminy)*0.5+bminy))

        self.list_of_mid_points_of_tables = mid_points_of_tables
        self.list_of_mid_points_of_half_tables = mid_points_of_half_tables

    def L1(self, p1, p2):
        """
        This function calculates the L1 (Manhattan) norm for two points.
        :param p1: Point 1
        :param p2: Point 2
        :return: Distance between Point 1 and Point 2
        """
        if len(p1) != len(p2):
            print('Error')
        else:
            return sum([abs(p1[i] - p2[i]) for i in range(len(p1))])

    def L1_pen_outside(self, p1, p2):
        """
        This function is a modified calculation of the L1 norm. It penalizes connections in which one part of the
        connection is outside the feasible region.
        :param p1: Point 1
        :param p2: Point 2
        :return: Distance between Point 1 and Point 2
        """
        if (not self.polygon.contains(LineString([p1, (p1[0], p2[1])]))
                or not self.polygon.contains(LineString([p2, (p1[0], p2[1])]))) \
                and (self.polygon.contains(LineString([p1, (p2[0], p1[1])]))
                     or not self.polygon.contains(LineString([p2, (p2[0], p1[1])]))):
            return self.L1(p1, p2)*1000
        else:
            return self.L1(p1, p2)

    def L1_pen_cross(self, p1, p2):
        """
        This function is a modified calculation of the L1 norm. It penalizes connections in which a service path is
        crossed.
        :param p1: Point 1
        :param p2: Point 2
        :return: Distance between Point 1 and Point 2
        """
        horizontal_line = LineString([p1, (p2[0], p1[1])])
        inter = False
        for i in self.list_lines:
            if i.intersects(horizontal_line):
                inter = True
                break
        if inter:
            return self.L1(p1, p2)*1000
        else:
            return self.L1(p1, p2)

    def L1_pen_outside_soph(self, p1, p2):
        """
        This function is a modified calculation of the L1 norm. It penalizes connections in which one part of the
        connection is outside the feasible region. The penalty depends on the length of the cable outside the area.
        :param p1: Point 1
        :param p2: Point 2
        :return: Distance between Point 1 and Point 2
        """
        way11 = LineString([p1, (p1[0], p2[1])])
        way12 = LineString([p2, (p1[0], p2[1])])
        int11 = self.polygon.intersection(way11)
        int12 = self.polygon.intersection(way12)
        obst11 = 0
        obst12 = 0
        for i in self.list_of_obstacles:
            obst11 += i.intersection(way11).length
            obst12 += i.intersection(way12).length

        dist1 = (way11.length - int11.length)*1000 + obst11*10000 + (int11.length - obst11) + \
                (way12.length - int12.length)*1000 + obst12*10000 + (int12.length - obst12)

        way21 = LineString([p1, (p2[0], p1[1])])
        way22 = LineString([p2, (p2[0], p1[1])])
        int21 = self.polygon.intersection(way21)
        int22 = self.polygon.intersection(way22)
        obst21 = 0
        obst22 = 0
        for i in self.list_of_obstacles:
            obst21 += i.intersection(way21).length
            obst22 += i.intersection(way22).length

        dist2 = (way21.length - int21.length)*1000 + obst21*10000 + (int21.length - obst21) + \
                (way22.length - int22.length)*1000 + obst22*10000 + (int22.length - obst22)

        return min(dist1, dist2)

    def L1_penalized(self, p1, p2):
        """
        This function is a modified calculation of the L1 norm. It penalizes connections in which one part of the
        connection is outside the feasible region. The penalty depends on the length of the cable outside the area and
        the length of cables crossing service paths.
        :param p1: Point 1
        :param p2: Point 2
        :return: Distance between Point 1 and Point 2
        """
        penalty_outside = 10000
        penalty_crossing = 1000
        way11 = LineString([p1, (p1[0], p2[1])])
        way12 = LineString([p2, (p1[0], p2[1])])
        int11 = self.polygon.intersection(way11)
        int12 = self.polygon.intersection(way12)
        obst11 = 0
        obst12 = 0
        for i in self.list_of_obstacles:
            obst11 += i.intersection(way11).length
            obst12 += i.intersection(way12).length
        cross11 = 0
        cross12 = 0
        for i in self.maintenance_paths:
            cross11 += i.intersection(way11).length
            cross12 += i.intersection(way12).length

        dist1 = (way11.length - int11.length)*penalty_outside + obst11*penalty_outside + \
                (int11.length - obst11 - cross11) + cross11*penalty_crossing + \
                (way12.length - int12.length)*penalty_outside + obst12*penalty_outside + \
                (int12.length - obst12 - cross12) + cross12*penalty_crossing

        way21 = LineString([p1, (p2[0], p1[1])])
        way22 = LineString([p2, (p2[0], p1[1])])
        int21 = self.polygon.intersection(way21)
        int22 = self.polygon.intersection(way22)
        obst21 = 0
        obst22 = 0
        for i in self.list_of_obstacles:
            obst21 += i.intersection(way21).length
            obst22 += i.intersection(way22).length
        cross21 = 0
        cross22 = 0
        for i in self.maintenance_paths:
            cross21 += i.intersection(way21).length
            cross22 += i.intersection(way22).length

        dist2 = (way21.length - int21.length)*penalty_outside + obst21*penalty_outside + \
                (int21.length - obst21 - cross21) + cross21*penalty_crossing + \
                (way22.length - int22.length)*penalty_outside + obst22*penalty_outside + \
                (int22.length - obst22 - cross22) + cross22*penalty_crossing

        return min(dist1, dist2)

    def L1_penalized_combined_simple(self, p1, p2):
        """
        This function is a modified calculation of the L1 norm. It penalizes connections in which a service path is
        crossed or one part of the connection is outside the feasible region.
        :param p1: Point 1
        :param p2: Point 2
        :return: Distance between Point 1 and Point 2
        """
        horizontal_line = LineString([p1, (p2[0], p1[1])])
        inter = False
        outside_area = ((not self.polygon.contains(LineString([p1, (p1[0], p2[1])]))
                         or not self.polygon.contains(LineString([p2, (p1[0], p2[1])])))
                        and (self.polygon.contains(LineString([p1, (p2[0], p1[1])]))
                             or not self.polygon.contains(LineString([p2, (p2[0], p1[1])]))))
        for i in self.list_lines:
            if i.intersects(horizontal_line):
                inter = True
                break
        if outside_area and inter:
            return self.L1(p1, p2)*10000*1000
        elif outside_area:
            return self.L1(p1, p2)*10000
        elif inter:
            return self.L1(p1, p2)*1000
        else:
            return self.L1(p1, p2)

    def L1_penalized_combined_soph(self, p1, p2):
        """
        This function is a modified calculation of the L1 norm. It penalizes connections in which a service path is
        crossed or one part of the connection is outside the feasible region. The penalty depends on the length of the
        cable outside the area.
        :param p1: Point 1
        :param p2: Point 2
        :return: Distance between Point 1 and Point 2
        """
        horizontal_line = LineString([p1, (p2[0], p1[1])])
        inter = False
        for i in self.list_lines:
            if i.intersects(horizontal_line):
                inter = True
                break
        if inter:
            return self.L1_penalized3(p1, p2) + 20 * 1000
        else:
            return self.L1_penalized3(p1, p2)

    def dist_change_clusters(self, clusters, oldclusters):
        """
        This function calculates the distange which changes in one iteration of the k-means-algorithm.
        :param clusters: Dictionary contaning the new cluster
        :param oldclusters: Dictionary containing the old cluster
        :return:
        """
        dist = 0
        for i in range(len(clusters)):
            dist = dist + self.L1(clusters[i], oldclusters[i])

        return dist

    def get_k(self, number_of_tables):
        """
        This function calculates the optimal k, number of inverters, for the k-means-algorithm.
        :param number_of_tables: Integer number of tables
        :return: Optimal k
        """
        k = math.ceil(number_of_tables/self.capacity_inverter)
        counter = 0
        while (k-counter)*self.capacity_inverter >= number_of_tables - (k-counter):
            counter += 1

        return (k-counter+1)

    def get_clusters(self, number_of_iterations, tol, conv_criterion=True, med=False):
        """
        This function assigns every PV module to a inverter. The basic idea is the k-means-algorithm. This algorithm had
        to be adjusted to our specific problem since we are limited by capacity. Furthermore, the inverter has to be
        in the feasible region and the cables should also be in the feasible region.
        :param number_of_iterations: Number of iterations in k-means-algorithm
        :param tol: Tolerance for the convergence criterion
        :param conv_criterion: True/False which convergence criterion should be used (True = change in cable length,
        False = change in inverter location)
        :param med: True/False whether you want to place the inverter in a new iteration at the median or mean position
        (median ensures that it is in the feasible region but in general the results are worse)
        :return: List of position of the inverters, dictionary of assigned PV tables to inverters
        """
        minx, miny, maxx, maxy = self.polygon.bounds
        tables = self.list_of_tables + self.list_of_half_tables
        mid_points_tables = self.list_of_mid_points_of_tables + self.list_of_mid_points_of_half_tables
        number_of_tables = len(self.list_of_tables) + len(self.list_of_half_tables)/2
        k = self.get_k(number_of_tables)
        if type(self.number_of_inverters) != int:
            self.number_of_inverters = k
        elif type(self.number_of_inverters) == int \
                and self.number_of_inverters*self.capacity_inverter < k:
            sys.exit('Not enough inverters')

        clusters = [random.choice(mid_points_tables) for i in range(self.number_of_inverters)]
        best_clusters = clusters.copy()

        lastmatches = None

        counter = 0
        cable_length = np.inf
        cycle_check = []
        cluster_list = []
        while counter < number_of_iterations or conv > tol:
            dict_capacity = {i: 0 for i in range(self.number_of_inverters)}
            bestmatches = {i: {} for i in range(self.number_of_inverters)}

            mid_points_tables_copy = mid_points_tables.copy()

            while len(mid_points_tables_copy) != 0:
                point = mid_points_tables_copy.pop()
                bestmatch = 0
                best_d = np.inf
                if point in self.list_of_mid_points_of_half_tables:
                    capa = 0.5
                else:
                    capa = 1
                for k in range(self.number_of_inverters):
                    d = self.L1_penalized(point, clusters[k]) #, counter)
                    if (d < best_d) & (dict_capacity[k] + capa <= self.capacity_inverter + 1):
                        bestmatch = k
                        best_d = d
                        delete = None
                    elif (d < best_d) & (dict_capacity[k] + capa > self.capacity_inverter + 1):
                        ind = None
                        key = list(bestmatches[k].keys())
                        val = list(bestmatches[k].values())
                        val_check = val.copy()
                        check = True
                        while check:
                            val_max = max(val)
                            if key[val_check.index(val_max)] in self.list_of_mid_points_of_half_tables:
                                capa_swap = 0.5
                            else:
                                capa_swap = 1

                            if val_max > d and capa == capa_swap:
                                ind = key[val_check.index(val_max)]
                                break
                            elif val_max > d and capa != capa_swap:
                                val.remove(val_max)
                            else:
                                break

                        if ind is not None:
                            bestmatch = k
                            best_d = d
                            delete = ind

                bestmatches[bestmatch][point] = best_d
                if delete is None:
                    dict_capacity[bestmatch] += capa
                else:
                    try:
                        del bestmatches[bestmatch][delete]
                        mid_points_tables_copy.append(delete)
                        if (point in self.list_of_mid_points_of_half_tables) \
                                and (delete in self.list_of_mid_points_of_tables):
                            dict_capacity[bestmatch] -= 0.5
                        elif (point in self.list_of_mid_points_of_tables) \
                                and (delete in self.list_of_mid_points_of_half_tables):
                            dict_capacity[bestmatch] += 0.5
                    except:
                        dict_capacity[bestmatch] += capa

            if bestmatches == lastmatches:
                break

            if med:
                for i in range(self.number_of_inverters):
                    if len(bestmatches[i]) != 0:
                        key = list(bestmatches[i].keys())
                        val = list(bestmatches[i].values())
                        median = np.median(val)
                        arg = np.argmin(abs(np.array(val)-median))
                    clusters[i] = key[arg]
            else:
                for i in range(self.number_of_inverters):
                    average = [0, 0]
                    if len(bestmatches[i]) != 0:
                        key = list(bestmatches[i].keys())
                        for j in range(len(bestmatches[i])):
                            average[0] = average[0] + key[j][0]
                            average[1] = average[1] + key[j][1]
                        average[0] = average[0] / len(bestmatches[i])
                        average[1] = average[1] / len(bestmatches[i])
                    clusters[i] = average

            if self.plot & self.save:
                self.cluster = bestmatches
                horizontal_lines, vertical_lines = self.get_cabling(clusters)
                self.plot_region_with_inverters(clusters, horizontal_lines, vertical_lines, option='color+inverter',
                                                save=counter)
            elif self.plot:
                self.cluster = bestmatches
                horizontal_lines, vertical_lines = self.get_cabling(clusters)
                self.plot_region_with_inverters(clusters, horizontal_lines, vertical_lines, option='color+inverter')

            if conv_criterion:
                new_cable_length = self.get_cable_length_kmean(bestmatches)
                conv = 1-new_cable_length/cable_length
                if abs(conv) < tol and tol < 0:
                    break
                elif abs(conv) < tol and tol >= 0:
                    lastmatches = bestmatches
                    best_clusters = clusters
                    break
                else:
                    lastmatches = bestmatches
                    best_clusters = clusters
                    cable_length = new_cable_length
                    if lastmatches in cycle_check:
                        index = cycle_check.index(lastmatches)
                        best_val = np.inf
                        best_ind = None
                        for checker in range(index, len(cycle_check)):
                            if self.get_cable_length_kmean(cycle_check[checker]) < best_val:
                                best_val = self.get_cable_length_kmean(cycle_check[checker])
                                best_ind = checker
                            lastmatches = cycle_check[best_ind]
                            best_clusters = cluster_list[checker]
                        break
                    else:
                        cycle_check.append(lastmatches)
                        cluster_list.append(best_clusters)
                    counter += 1
            else:
                new_dist_change = self.dist_change_clusters(clusters, best_clusters)
                if counter == 0:
                    dist_change = new_dist_change
                conv = new_dist_change/dist_change
                lastmatches = bestmatches
                best_clusters = clusters.copy()
                if conv < tol:
                    break
                else:
                    dist_change = new_dist_change
                    if lastmatches in cycle_check:
                        index = cycle_check.index(lastmatches)
                        best_val = np.inf
                        best_ind = None
                        for checker in range(index, len(cycle_check)):
                            if self.get_cable_length_kmean(cycle_check[checker]) < best_val:
                                best_val = self.get_cable_length_kmean(cycle_check[checker])
                                best_ind = checker
                            lastmatches = cycle_check[best_ind]
                            best_clusters = cluster_list[checker]
                        break
                    else:
                        cycle_check.append(lastmatches)
                        cluster_list.append(best_clusters)
                    counter += 1

        position_inverters = best_clusters
        self.cluster = lastmatches

        unconnected_tables = []
        for i in dict_capacity.keys():
            if dict_capacity[i] % 1 != 0:
                rem_table = None
                max_dist = 0
                for j in self.cluster[i].keys():
                    if j in self.list_of_mid_points_of_half_tables and self.cluster[i][j] > max_dist:
                        rem_table = j
                        max_dist = self.cluster[i][j]
                del self.cluster[i][rem_table]
                dict_capacity[i] = dict_capacity[i] - 0.5
                self.remove_table(rem_table)
                table = self.get_table(rem_table, half=True)
                unconnected_tables.append(table)

        self.unconnected_tables = unconnected_tables

        return position_inverters

    def transform_cluster_to_point(self, position_inverters):
        """
        This function transforms inverters to points.
        :param position_inverters: List with position of the inverters
        :return: List with position of inverters as points
        """
        points = []
        for i in position_inverters:
            points.append(Point(i[0], i[1]))

        return position_inverters

    def remove_table(self, point):
        """
        This function removes a table for a given midpoint.
        :param point: Midpoint of a tables as tuple
        :return: Updated list of (half) PV tables
        """
        p = Point(point[0], point[1])
        for i in self.list_of_tables + self.list_of_half_tables:
            if i.contains(p):
                if (p.x, p.y) in self.list_of_mid_points_of_tables:
                    self.list_of_tables.remove(i)
                    break
                else:
                    self.list_of_half_tables.remove(i)
                    break

    def get_table(self, midpoint, half):
        """
        This functions transforms a midpoint to a table.
        :param midpoint: Midpoint of a table
        :param half: True/False whether it is a half table or not
        :return: Shapely box representing the table
        """
        if half:
            return box(midpoint[0] - self.width_table / 4, midpoint[1] - self.length_table / 2,
                       midpoint[0] + self.width_table / 4, midpoint[1] + self.length_table / 2)
        else:
            return box(midpoint[0] - self.width_table / 2, midpoint[1] - self.length_table / 2,
                       midpoint[0] + self.width_table / 2, midpoint[1] + self.length_table / 2)

    def get_cabling(self, position_inverters):
        """
        This function returns the cabling of the PV tables to the inverters.
        :param position_inverters:
        :return: Two lists for cables (one for horizontal, one for vertical)
        """
        horizontal_lines = []
        vertical_lines = []
        for i in range(self.number_of_inverters):
            for j in self.cluster[i].keys():
                if position_inverters[i] is not None:
                    way11 = LineString([position_inverters[i], (position_inverters[i][0], j[1])])
                    way12 = LineString([j, (position_inverters[i][0], j[1])])

                    int11 = self.polygon.intersection(way11)
                    int12 = self.polygon.intersection(way12)

                    obst11 = 0
                    obst12 = 0
                    for k in self.list_of_obstacles:
                        obst11 += k.intersection(way11).length
                        obst12 += k.intersection(way12).length

                    dist1 = (way11.length - int11.length)*1000 + obst11*10000 + (int11.length - obst11) + \
                            (way12.length - int12.length)*1000 + obst12*10000 + (int12.length - obst12)

                    way21 = LineString([j, (j[0], position_inverters[i][1])])
                    way22 = LineString([position_inverters[i], (j[0], position_inverters[i][1])])

                    int21 = self.polygon.intersection(way21)
                    int22 = self.polygon.intersection(way22)

                    obst21 = 0
                    obst22 = 0
                    for k in self.list_of_obstacles:
                        obst21 += k.intersection(way21).length
                        obst22 += k.intersection(way22).length

                    dist2 = (way21.length - int21.length)*1000 + obst21*10000 + (int21.length - obst21) + \
                            (way22.length - int22.length)*1000 + obst22*10000 + (int22.length - obst22)

                    if dist1 < dist2:
                        horizontal_lines.append(way11)
                        vertical_lines.append(way12)
                    else:
                        horizontal_lines.append(way21)
                        vertical_lines.append(way22)

        return horizontal_lines, vertical_lines

    def get_cable_length(self, horizontal_lines, vertical_lines):
        """
        This function gets the length of the placed cables. This is the real length of the cables.
        :param horizontal_lines: List of horizontal cables
        :param vertical_lines: List of vertical cables
        :return: Length of cables
        """
        counter = 0
        for i in (horizontal_lines + vertical_lines):
            counter += i.length

        return counter

    def get_cable_length_kmean(self, inverter_table_assignment):
        """
        This function gets the length of the placed cables. This is the length based on the L1 metric.
        :param inverter_table_assignment: Dictionary for the inverters with the information on which tables are assigned
        and the corresponding L1 norm.
        :return: Length of cables
        """
        counter = 0
        for i in range(len(inverter_table_assignment)):
            for j in inverter_table_assignment[i].values():
                counter += j

        return counter

    def get_nearest_table(self, position_inverters):
        """
        This functions looks for the nearest table for a certain point. It is used to get the inverter positions. A list
        with the current positions is inputed and will be transformed to a list with midpoints of closest tables.
        :param position_inverters: List of inverter positions
        :return: List of new inverter positions representing midpoints of tables
        """
        new_position_inverters = []
        for i in range(len(position_inverters)):
            min = np.inf
            new_position = (0, 0)
            for j in self.cluster[i]:
                if self.L1(position_inverters[i], j) < min and j in self.list_of_mid_points_of_tables:
                    min = self.L1(position_inverters[i], j)
                    new_position = j

            if new_position in self.list_of_mid_points_of_tables:
                new_position_inverters.append(new_position)
            else:
                new_position_inverters.append(None)

        return new_position_inverters

    def get_inverters_instead_of_tables(self, position_inverters):
        """
        This function removes the tables and places inverters instead. It returns a list with the new positions of the
        inverters. In addition to the called function it removes the tables.
        :param position_inverters: List of inverter positions
        :return: List of new inverter positions representing midpoints of tables
        """
        new_position_inverters = self.get_nearest_table(position_inverters)
        for i in new_position_inverters:
            if i is not None:
                self.remove_table(i)
                del self.list_of_mid_points_of_tables[self.list_of_mid_points_of_tables.index(i)]

        return new_position_inverters

    def get_lower_bound(self, alpha, gamma):
        """
        This function should give a lower bound for the length of the cables. The idea is to have one inverter fixed in
        the middle. Around the inverter tables are placed and this function finds the closest tables to inverter until
        capacity is full. Then this lower bound for one inverter gets multiplied by number of tables divided by capacity
        of inverters to get an overall lower bound.
        :param alpha: Angle for PV table
        :param gamma: Angle for sun
        :return: Overall lower bound for cabling
        """
        alpha = alpha/180 * math.pi
        gamma = gamma/180 * math.pi
        height_table = math.sin(alpha) * self.length_table/math.cos(alpha)
        rowspace = height_table / math.tan(gamma)
        up = self.length_table + rowspace
        side = self.width_table
        counter = 1
        while counter*counter-1 < self.capacity_inverter:
            counter += 2
        mat = [[0 for i in range(counter)] for j in range(counter)]
        for i in range(counter):
            for j in range(counter):
                mat[i][j] = ((i-(counter-1)/2)*side, (j-(counter-1)/2)*up)

        dict = {}
        for i in range(counter):
            for j in range(counter):
                if (i != 0 or j != 0) and len(dict) == self.capacity_inverter:
                    if self.L1(mat[i][j], mat[int((counter-1)/2)][int((counter-1)/2)]) < max(list(dict.values())):
                        del dict[list(dict.keys())[list(dict.values()).index(max(list(dict.values())))]]
                        dict[mat[i][j]] = self.L1(mat[i][j], mat[int((counter-1)/2)][int((counter-1)/2)])
                elif (i != 0 or j != 0) and len(dict) < self.capacity_inverter:
                    dict[mat[i][j]] = self.L1(mat[i][j], mat[int((counter-1)/2)][int((counter-1)/2)])

        length_one_cluster = sum(list(dict.values()))
        lower_bound = (len(self.list_of_tables) + len(self.list_of_half_tables)/2)/self.capacity_inverter * \
                      length_one_cluster

        return lower_bound

    def plot_region_with_inverters(self, position_inverters, horizontal_lines, vertical_lines, option='color', save=''):
        """
        This function plots the feasible regions with the placed inverters and the cabling.
        :param position_inverters: List of the position of inverters
        :param horizontal_lines: List of horizontal cables
        :param vertical_lines: List of vertical cables
        :param option: Option for plotting; either color or color+inverter or nothing
        :param save: True/False whether to save the plot or not
        :return: Plot of feasible regions with placed tables, service paths and inverters
        """
        minx, miny, maxx, maxy = self.polygon.bounds
        distx = maxx - minx
        disty = maxy - miny
        x, y = self.polygon.exterior.xy
        plt.plot(x, y)

        for i in self.decomposition_list:
            x, y = box(i.bounds[0], i.bounds[1], i.bounds[2], i.bounds[3]).exterior.xy
            plt.plot(x, y, color='green', lw=2)

        #for i in self.list_lines:
        #    x1, y1 = i.xy
        #    plt.plot(x1, y1, color='red')
        for i in self.maintenance_paths:
            x, y = i.exterior.xy
            plt.fill(x, y, color='grey')

        for i in self.list_of_tables:
            x, y = i.exterior.xy
            plt.plot(x, y, color='black', lw=.5)

        for i in self.list_of_half_tables:
            x, y = i.exterior.xy
            plt.plot(x, y, color='black', lw=.5)

        try:
            for i in self.unconnected_tables:
                x, y = i.exterior.xy
                plt.plot(x, y, color='red', lw=.5)
        except:
            pass

        for i in self.list_of_obstacles:
            x, y = i.exterior.xy
            plt.plot(x, y, color='black')

        if option == 'color':
            viridis = cm.get_cmap('twilight')
            for i in range(self.number_of_inverters):
                key = list(self.cluster[i].keys())
                for j in key:
                    if j in self.list_of_mid_points_of_half_tables:
                        table = self.get_table(j, half=True)
                        x, y = table.exterior.xy
                        plt.fill(x, y, c=viridis((i + 1) / self.number_of_inverters), zorder=2)
                    elif j in self.list_of_mid_points_of_tables:
                        table = self.get_table(j, half=False)
                        x, y = table.exterior.xy
                        plt.fill(x, y, c=viridis((i + 1) / self.number_of_inverters), zorder=2)
        elif option == 'color+inverter':
            viridis = cm.get_cmap('RdYlBu')
            for i in range(self.number_of_inverters):
                key = list(self.cluster[i].keys())
                for j in key:
                    if j in self.list_of_mid_points_of_half_tables:
                        table = self.get_table(j, half=True)
                        x, y = table.exterior.xy
                        plt.fill(x, y, c=viridis((i + 1) / self.number_of_inverters), zorder=1)
                    elif j in self.list_of_mid_points_of_tables:
                        table = self.get_table(j, half=False)
                        x, y = table.exterior.xy
                        plt.fill(x, y, c=viridis((i + 1) / self.number_of_inverters), zorder=1)

            for i in horizontal_lines:
                x1, y1 = i.xy
                plt.plot(x1, y1, color='olive', lw=.4, zorder=5)

            for i in vertical_lines:
                x1, y1 = i.xy
                plt.plot(x1, y1, color='olive', lw=.4, zorder=5)

            for i in position_inverters:
                if i is not None:
                    plt.scatter(i[0], i[1], color='black', lw=2, zorder=5)
        else:
            for i in horizontal_lines:
                x1, y1 = i.xy
                plt.plot(x1, y1, color='olive', lw=.4, zorder=5)

            for i in vertical_lines:
                x1, y1 = i.xy
                plt.plot(x1, y1, color='olive', lw=.4, zorder=5)

            for i in position_inverters:
                plt.scatter(i[0], i[1], color='orange', lw=2)

        plt.xlim([minx - 0.05 * distx, maxx + 0.05 * distx])
        plt.ylim([miny - 0.05 * disty, maxy + 0.05 * disty])
        plt.suptitle(f'Number of tables: {int(len(self.list_of_tables) + len(self.list_of_half_tables) / 2)}')
        plt.axis('off')
        if save != '':
            plt.savefig(f'{os.getcwd() + self.path_for_plots}plot_inverter_{save}.png')
        plt.show()
        plt.close()

    def get_filenames_for_gif(self, number_of_iterations, include_final=False):
        """
        This function gets the filenames which are needed to generate a gif.
        :param number_of_iterations: Number of iterations in clustering step
        :param include_final: True/False whether to include the final plot with converters and cables
        :return: List of filenames
        """
        filenames = []
        for i in range(number_of_iterations+1):
            filenames.append(f'{os.getcwd() + self.path_for_plots}plot_inverter_{i}.png')

        if include_final:
            filenames.append(f'{os.getcwd() + self.path_for_plots}plot_inverter_final.png')

        return filenames

    def convert_images_to_gif(self, filenames):
        """
        This function converts several images to a gif.
        :param filenames: List of filenames
        :return: Gif
        """
        images = []

        for i in filenames:
            try:
                images.append(imageio.imread(i))
            except:
                break

        kargs = {'duration': 2}
        imageio.mimsave(os.getcwd() + self.path_for_gif + 'clustering.gif', images, **kargs)

    def create_gif(self, number_of_iterations, include_final=False):
        """
        This function creates a gif file.
        :param number_of_iterations: Number of iterations in clustering step
        :param include_final: True/False whether to include the final plot with converters and cables
        :return: Gif
        """
        filenames = self.get_filenames_for_gif(number_of_iterations, include_final)
        self.convert_images_to_gif(filenames)
