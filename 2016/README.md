aciccomp2016
============

Data and simulations from the 2016 Atlantic Causal Inference competition.

## Installation

Steps to install from source:

    install.packages("devtools")
    devtools::install_github("vdorie/aciccomp/2016")

## Usage

The raw data given to contestants in the 2016 competition is exported by the package in the object `input_2016`. The main function used to generate data is `dgp_2016`. The 77 simulation settings and their 100 replications can be accessed by running:

    dgp_2016(input_2016, parameterNum, simulationNum)

where `s parameterNum ` is in `1:77` and `simulationNum` is in `1:100`.
