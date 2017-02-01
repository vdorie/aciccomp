# Overview

This repo contains the code to run causal inference methods against a set of data files and track the results.

# Directory Structure

Because the contest contents are a secret, most of the important files go in sub-directories that are currently empty. The required subdirectories are

  * `data`
  * `log`
  * `methods`
  * `results`
  * `job`
  * `src`

They can be located anywhere, but in general the code expects them to defined before running. The easiest way to do this to create a file `site_setup.R` at the top-level that will be `source`d and contains something like:

    dirs <- list(data    = "data",
                 log     = "log",
                 methods = "methods",
                 results = "results",
                 src     = "src",
                 job     = "job")

## `data`

At the top-level of the data folder needs to be a file `x.csv`, that is used as the covariates for all simulations.

Below that, files can be organized however desired within the constraint that at leaf folders, individual simulations are named `iteration_number.csv`. For example:

    data/dgp_1/sigma_1/1.csv
    data/dgp_1/sigma_1/2.csv
    ...
    data/dgp_1/sigma_1/100.csv
    data/dgp_1/sigma_2/1.csv

and so forth.

## `methods`

Each file in the `methods` folder should be an executable (or more-likely script) that accepts 2 to 3 arguments. The first argument should be the input csv with columns `y`, `z`, and `x_1` through `x_d`. The second argument is the name of the output file, where the method should write a csv with columns `est`, `ci_lower`, and `ci_upper`. The third file, if present, is for individual effects and has the same format as the second but has as many rows as the data file.

The entries in the methods should correspond to a file `methods.csv` that is detailed below. If methods require any additional files, they should be placed in `methods/src`. Help can be provided for commands that cannot find their own executable path.

## `log` and `results`

These are created automatically.

# General Procedure

Since the code needs to track submissions that come in over time and the results of those submissions that come in and may or may not fail, the procedure involves a number of intermediate steps. In general:

  1. Drop datasets of interest in the `data` folder
  2. Create a `runCases` csv which scans the contents of the data folder and enumerates all the iterations (`bin/create_runcases`). These should probably be set before the competition is run.
  3. Synchronize the executables/scripts in the `methods` folder with the contents of `methods.csv`.
  3. Create a `runStatus` Rdata object (`bin/update_runstatus`), which holds the status of every method run on every data set (success, failure, hang, etc.).
  4. Create a `results` Rdata object (`bin/create_results bias coverage ...`), which holds values derived from the input files and estimates computed by each method.

The above steps need to be run prior to spawning off any processes. Once they are in place, methods can be fit by running:

  5. Queue jobs (`bin/queue_jobs`) or run them locally (`bin/run_locally`). This proceeds through each method, simulation setting/DGP, and iteration, creating a temp file that contains the input and storing the output in the `results` folder.
  6. Update the `runStatus` (`bin/update_runstatus`) - this scans the results and log folders to see what has completed.
  7. Update the `results` (`bin/update_results evaluation_functions.R`) - this uses the runStatus and current value of results to see what needs to be computed to bring the derived results in line with raw computations.

# `methods.csv`

This needs to be a csv with columns:

    name,display,abbreviation,individual_effects,language

  * `name` is the same as that as the executable
  * `display` long name to be used in plotting
  * `abbreviation` used in plotting
  * `individual_effects` `0/1` if the method is capable of computing individual effects
  * `language` `R`/`matlab`/`stata`/`python`/etc; Each needs to be supported
