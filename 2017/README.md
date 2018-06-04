aciccomp2017
============

Data and simulations from the 2017 Atlantic Causal Inference competition.

## Installation

Steps to install from source:

    install.packages("devtools")
    devtools::install_github("vdorie/aciccomp/2017")

## Usage

The raw data given to contestants in the 2017 competition is exported by the package in the object `input_2017`. The main function used to generate data is `dgp_2017`. The 32 simulation settings and their 250 replications can be accessed by running:

    dgp_2017(parameterNum, simulationNum)

where `parameterNum` is in `1:32` and `simulationNum` is in `1:250`. The mapping of parameter numbers to settings is provided by the data set `parameters_2017`.
