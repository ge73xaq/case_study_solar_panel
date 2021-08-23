"""
Authors: Daniela Schlager (daniela.schlager@online.de), Yichen Lou (yichen.lw99@gmail.com), Jessica Reichelt
(j.s.reichelt@gmail.com), Lukas Dreier (lukas.dreier@gmx.de)

Runs the main logic
"""
import yaml
import argparse
import os
import preprocess
import model
from model import Result
from time import time
from datetime import timedelta as td
import logging


# create logger with 'spam_application'
logger = logging.getLogger('discrete-optimization-siemens.main')
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
fh = logging.FileHandler('siemens.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)

topo = False
save = False
gif = False
plot = False
savekmean = False
plotkmean = False
gifkmean = False


def parse_cfg(path):
    if not os.path.isfile(path):
        exit(1)
    with open(path, 'r') as stream:
        try:
            return yaml.load(stream, Loader=yaml.SafeLoader)
        except yaml.YAMLError as exc:
            exit(1)


_tick = time()
parser = argparse.ArgumentParser('CaseStudiesDiscreteOptimization')
parser.add_argument('-path', '--path', action='store', dest='path', required=True)
parser.add_argument('-output', '--output', action='store', dest='output', required=True)
args = parser.parse_args()

cfg = parse_cfg(args.path)


coordinates = cfg['region'][1]
desc = cfg['region'][0]
obstacles = cfg['obstacles']
width_between_maint_paths = cfg['parameters']['width_between_maint_paths']
width_of_maint_path = cfg['parameters']['width_of_maint_path']
width_table = cfg['parameters']['width_table']
length_table = cfg['parameters']['length_table']
alpha = cfg['parameters']['alpha']
gamma = cfg['parameters']['gamma']
capacity_inverter = cfg['parameters']['capacity_inverter']
number_of_inverters = cfg['parameters']['number_of_inverters']
shifting_stepsx = cfg['hyperparameter']['x']
shifting_stepsy = cfg['hyperparameter']['y']
number_of_iterations = cfg['hyperparameter']['number_of_iterations']
tol = cfg['hyperparameter']['tol']
path_for_tiff = cfg['input']['path_for_tiff']
path_for_plots = cfg['input']['path_for_plots']
path_for_gif = cfg['input']['path_for_gif']

logger.info(f'Start running logic for {desc}')
_tick1 = time()
preprocessor = preprocess.PreProcessor(coordinates, path_for_tiff, obstacles)
preprocessor.get_polygon()
if topo:
    preprocessor.generate_tiff_file()
preprocessor.add_topographic_information()
preprocessor.plot_topographic_information()

polygon = preprocessor.polygon
elevation = preprocessor.elevation
list_of_obstacles = preprocessor.get_obstacles()
_tock1 = time()
preprocessing_time = td(seconds=_tock1-_tick1)
logger.info(f'Preprocessing Time: {preprocessing_time}')

_tick1 = time()
shifter = model.Shifter(polygon, elevation, width_between_maint_paths, width_of_maint_path, width_table, length_table,
                        alpha, gamma, list_of_obstacles, shifting_stepsx, shifting_stepsy, path_for_plots, path_for_gif,
                        save, plot)

list_of_tables, list_of_half_tables, decomposition_list, list_lines = shifter.shifting_x_direction()
maintenance_paths = shifter.create_smallest_fitting_rectangles(list_lines, paths=True)

if len(list_of_tables) > 0 and (number_of_inverters == 'Calc' or number_of_inverters <= len(list_of_tables)):
    ratio_of_tables = round(shifter.get_ratio_of_placed_tables(list_of_tables, list_of_half_tables), 2)
    number_of_placed_tables = int(len(list_of_tables))
    number_of_placed_half_tables = int(len(list_of_half_tables))
    number_of_placed_tables_in_total = int(len(list_of_tables) + len(list_of_half_tables)/2)
    coverage, coverage_of_tables, coverage_of_paths, coverage_of_obstacles = shifter.get_used_area(list_of_tables,
                                                                                                   list_of_half_tables,
                                                                                                   list_lines)

    shifter.plot_region(decomposition_list, list_of_tables, list_of_half_tables, list_lines, save='final')

    if gif:
        shifter.create_gif(include_final=True)
    _tock1 = time()
    shifting_time = td(seconds=_tock1-_tick1)
    logger.info(f'Shifting Time: {shifting_time}')

    _tick1 = time()
    clustering = model.KMeans(polygon, list_of_obstacles, list_of_tables, list_of_half_tables, capacity_inverter,
                              number_of_inverters, width_table, shifter.length_table, decomposition_list,
                              maintenance_paths, list_lines, path_for_plots, path_for_gif, savekmean, plotkmean)

    clustering.get_mid_points_of_tables()
    position_inverters = clustering.get_clusters(number_of_iterations, tol, conv_criterion=True, med=False)

    position_inverters = clustering.get_inverters_instead_of_tables(position_inverters)
    horizontal_lines, vertical_lines = clustering.get_cabling(position_inverters)
    clustering.plot_region_with_inverters(position_inverters, horizontal_lines, vertical_lines, option='color+inverter',
                                          save='final_color_inverter')
    clustering.plot_region_with_inverters(position_inverters, horizontal_lines, vertical_lines, option='test',
                                          save='final')
    if gifkmean:
        clustering.create_gif(number_of_iterations)

    cable_length = clustering.get_cable_length(horizontal_lines, vertical_lines)
    length_of_cables = round(cable_length/1000, 2)

    lower_bound = clustering.get_lower_bound(alpha, gamma)
    ratio_of_cable_length = lower_bound/cable_length

    _tock1 = time()
    clustering_time = td(seconds=_tock1-_tick1)
    logger.info(f'Clustering Time: {td(seconds=_tock1-_tick1)}')

    _tock = time()
    overall_time = td(seconds=_tock-_tick)
    logger.info(f'Overall Time: {overall_time}')
    logger.info('End of run')

    r = Result(description=desc,
               path_to_cfg_file=args.output,
               path_to_tiff=f'{os.getcwd()}{path_for_tiff}',
               coordinates=preprocessor.list_of_coordinates,
               number_of_placed_tables=len(clustering.list_of_tables),
               number_of_placed_half_tables=len(clustering.list_of_half_tables),
               number_of_placed_tables_in_total=int(len(clustering.list_of_tables)
                                                    + len(clustering.list_of_half_tables)/2),
               number_of_inverters=clustering.number_of_inverters,
               length_of_cables=length_of_cables,
               coverage=round(coverage*100, 2),
               coverage_of_tables=round(coverage_of_tables*100, 2),
               coverage_of_paths=round(coverage_of_paths*100, 2),
               coverage_of_obstacles=round(coverage_of_obstacles*100, 2),
               ratio_of_tables=ratio_of_tables,
               ratio_of_cable_length=round(ratio_of_cable_length, 2),
               preprocessing_time=preprocessing_time,
               shifting_time=shifting_time,
               clustering_time=clustering_time,
               overall_time=overall_time)

    output = r.write_result()

    with open(args.output, 'a') as f:
        f.write(output)
elif number_of_inverters > len(list_of_tables):
    logger.info('Too many inverters for modules')
    print('Too many inverters for PV modules')
else:
    logger.info('No module placed')
    print('No module placed')
