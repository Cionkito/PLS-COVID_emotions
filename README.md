# PLS COVID-19 emotions

Code to compute partial least square (PLS) correlation results and figure for the publication : Thomasson, M., et al. (2023) _Markers of limbic system damage following SARS-CoV-2 infection_


The code includes two parts: the first part runs in Matlab (folder MATLAB) and is taken from the _myPLS_ toolbox (available here: https://github.com/MIPLabCH/myPLS); the second part (generation of figures) is implemented in Python.

## System requirement & installation

The following softwares need to be installed:

- For the matlab code you (obviously) need MATLAB (originally created on version R2017a).

- For python code, you will need python
  - Recommended to follow the steps from [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation) to install python.
  - You may install the environment from the provided file as described [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file).

## How to run

### Part 1 (MATLAB):

Instructions: You may specifiy your parameters in the `myPLS_inputs_CSTM.m` file, then run the main script by typing `myPLS_main` in the Matlab command window or launching the `myPLS_main.m` file.

Expected output: all PLS computations described in the PLS part of the paper are performed and stored in the `outputs/PLS_results/Emo-Rerun-wcov3-wmemo1_GRP` folder.

Expected run time: around 2 minutes.


### Part 2 (Python):

Instructions: open and run the script `PLSC_figures.ipynb` (in the `python` folder) in the Jupyter Notebook. 

Expected output: The notebook loads the PLS results and produces the A and B subplots of figure 4. 

Expected run time: around 1 minutes.
