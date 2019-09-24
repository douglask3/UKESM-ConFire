# fire_limitation
Bayesian parameterisation of the Kelley et al. global fire limitation model.

## How to optimize
Optimization is performed in two jupter notebooks in "notebooks" dir:
* [notebooks/prepare_data.ipynb](https://github.com/rhyswhitley/fire_limitation/blob/master/notebooks/prepare_data.ipynb) prepares data ready for optimization.
* [notebooks/bayesian_inference.ipynb](https://github.com/rhyswhitley/fire_limitation/blob/master/notebooks/bayesian_inference.ipynb) performs the optimzation, and can be used without prepare_data if using csv input data.

## What to install
To run the optimation, you first need jupyter installed ([see here](https://jupyter.org/install.html)), along withthe following libraries:
* numpy as
* pandas
* pymc3 
* scipy
* theano

netCDF4 will also be required to run prepare_data, and matplotlib to perform plotting.

## Reproducing "shifts in global fire regimes due to emerging trends in bioclimatic and human drive"

This is run using the [paper1_bayes](https://github.com/rhyswhitley/fire_limitation/tree/paper1_bayes/notebooks) branch. 
To run with the original data, please contact [douglas.i.kelley@gmail.com](mailto:douglas.i.kelley@gmail.com).


## How to extract data
Added [24/09/19] by Megan Brown

_Note: iris is needed in order to run these scripts_

To extract data from UK Earth System Model (UKESM), use the [multi-year_extract_inputs](https://github.com/douglask3/UKESM-ConFire/tree/mb_Tree/notebooks/multi-year_extract_inputs.ipynb) and set the years of interest in the file (notebooks folder).
Both the .ipynb and .py are identical wrt code - .py is designed for running in JASMIN, while .ipynb gives a description of what the code is doing.
Also see: [retrieve_stash.ipynb](https://github.com/douglask3/UKESM-ConFire/tree/mb_Tree/notebooks/retrieve_stash.ipynb) which gives a more detailed breakdown of the code.

The outputs of the extracted inputs from UKESM can be plotted up using the [plot_stash_outputs](https://github.com/douglask3/UKESM-ConFire/tree/mb_Tree/notebooks/plot_stash_outputs.ipynb) to check it has been retrieved properly.

## Running the model

Use [make_model_output](https://github.com/douglask3/UKESM-ConFire/tree/mb_Tree/notebooks/make_model_output2.ipynb) to run the model.
Both .py and .ipynb can be used (.py designed for JASMIN use).
Here are some guidelines and useful tips for make_model_output script:


* Depending on which inputs you want to use for the model, the dictionary at the top of the file must be changed.
* Comment/comment out the files you want to use (e.g. all UKESM inputs with the exception of observed vegcover).
* All input files must be in the same location.
* Change the ``title_output`` to whatever you'd like to call the run. This will be the title of all figures as well as the .nc output
* You must mannually create a folder in the figures directory with the same title as ``title_output``.
* Toggle figures on/off by setting True/False at the beginning.

## Regridding and formatting inputs

There are a few scripts in the notebooks folder which are helpful for regridding certain data and filtering out unwanted years, as well as reformatting cubes.

* [Regridding JULES](https://github.com/douglask3/UKESM-ConFire/tree/mb_Tree/notebooks/Regridding_JULES_data.ipynb): JULES data is in a funny format and needs condensing into one .nc file.
* [Regridding to n96-e](https://github.com/douglask3/UKESM-ConFire/tree/mb_Tree/notebooks/regrid_nc.ipynb): If you need to regrid the observations, this is the script for you.
* [Regridding and condensing population density](https://github.com/douglask3/UKESM-ConFire/tree/mb_Tree/notebooks/pop_dens_regrid.ipynb):
This script is also useful for changing coordinates on cubes and various other bits and pieces.
* [UKESM cube format](https://github.com/douglask3/UKESM-ConFire/tree/mb_Tree/notebooks/formatting_files.ipynb): All file inputs must have the same mapping coordinates and time frames. This script will take the data from the file of interest and put it into a UKESM-framed cube.
* [Filtering months](https://github.com/douglask3/UKESM-ConFire/tree/mb_Tree/notebooks/filter_months.ipynb): The observation inputs run from July 2000 to June 2014, so this script will adapt the UKESM inputs from Jan 2000-Dec 2014 to match this.
The script isn't written that well and is quite clunky - ye hath been warn&eacute;d.

## Plot model outputs

[Plot burnt area](https://github.com/douglask3/UKESM-ConFire/tree/mb_Tree/notebooks/plot_burnt_area.ipynb) takes the outputs from the model run and will plot them - both spatially and temporally if one so desires.
I've also written some code in there which will plot the standard controls and can take the difference between two runs. Again, this script is poorly written and cringeworthy. 

If there are any other strange scripts I've forgotten to mention, just let me (Megan) know.
