# Discrete Optimization SIEMENS
This project was proposed by Siemens Corporate Technology. The idea is to define a feasible region with obstacles and
get information on how many tables could be placed in that region and how many cables and inverters do you need under
several constraints. The constraints are the topography of the region, the capacity of inverters, the shape of the 
feasible region or the distance between maintenance paths.

## Implementation
This project was built in Python 3.6. The implementation contains the `main.py`, `preprocess.py` and `model.py` file. 
Every function is described with input and output and the general functionality.
- `main.py` file contains the main logic of the project
- `preprocess.py` file contains the data input and the preprocessing of the data in general 
- `model.py` file the mathematical implementation of the shifting and the clustering is done 

## Required packages
All packages used in the three scripts (main.py, preprocess.py and model.py) are contained in the requirements.txt file.
Most of the packages are straight forward to pip install. One causing troubles is the osgeo/gdal package. If you use a 
virtual environment the description in this article 
(https://stackoverflow.com/questions/51934231/gdal-virtualenv-python-3-6-installation) is pretty good. 
If someone needs help then please let us know. In the best case installing all packages works with 
```pip install requirements.txt```.

## Configuration file
The configuration file specifies the feasible regions. For this you can copy and paste the output from the Copernicus 
Open Access Hub (https://scihub.copernicus.eu/dhus/#/home). There you can select a polygon and paste the information in 
the following form to the config file: __'13.188942243765236 48.76164346213886,13.192288620520241 
48.76600913233696,13.203164344974017 48.76270044868767,13.198005347476714 48.755530886306786,13.189778837953986 
48.75874812420179,13.188942243765236 48.76164346213886,13.188942243765236 48.76164346213886'__. Inputing the information 
on obstacles works analogous to the polygon. Furthermore, it contains some parameters which need to be set in order to 
run the main logic. For saving plots correctly the right path has to be specified. The relative path from your current 
working directory is meant. Also the data for the PV modules can be specified in the configuration file.

## Run main logic
As soon as the config file is specified correctly the main logic can be ran. This is rather simple since you only have 
to specify the input config file and the output file. An example for this is:
```
python3 main.py -path /home/lukas/Dokumente/CaseStudies/discrete-optimization-siemens/cfg.yml -output result/test.result
```

## Generation of TIFF file
Obviously the generation of the TIFF file does not work for Windows user. In our default implementation the generation
of the TIFF file is switched off. For the delivered example shapes it is no problem since we generated these tiff files.
They are saved in the data folder. If you want to test the tool on new data you have to set the variable __topo__ to 
__True__ in the `main.py` file. 

## Output
The output file contains some statistics on the model. Besides the placed number of PV modules the ratio of possible 
modules is calculated as well as the duration of several parts of the model is contained.

## Contact
TUM advisors: Katie Fitch, Ulf Friedrich <br>
Project team: Daniela Schlager (daniela.schlager@online.de), Yichen Lou (yichen.lw99@gmail.com), Jessica Reichelt 
(j.s.reichelt@gmail.com), Lukas Dreier (lukas.dreier@gmx.de)
