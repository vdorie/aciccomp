# Overview

This repo contains the code to run causal inference methods against a set of data files and track the results.

# Directory Structure

Because the contest contents are a secret, most of the important files go in sub-directories that are currently empty. Some of the empty subdirectories should be created automatically, but the totality of subdirectories includes:

  * `data`
  * `log`
  * `methods`
  * `results`
  * `job`
  * `src`

They can be located anywhere, but in general the code expects them to defined before running. The easiest way to do this to create a file `site_setup.R` at the top-level that will be `source`d and containing something like:

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

Each file in the `methods` folder should be an executable that accepts 2 to 3 arguments. This can be a script so long as it can be run from the command line. The first argument should be the input csv with columns `y`, `z`, and `x_1` through `x_d`. The second argument is the name of the output file, where the method should write a csv with columns `est`, `ci_lower`, and `ci_upper`. The third file, if present, is for individual effects and has the same format as the second but has as many rows as the data file.

The entries in the methods should correspond to a file `methods.csv` that is detailed below. If methods require any additional files, they should be placed in `methods/src`. Help can be provided for commands to find heir own executable path.

## `log`, `results`, `jobs`

These are created automatically.

# General Procedure

Since the code needs to track submissions that come in over time and the results of those submissions that come in and may or may not fail, the procedure involves a number of intermediate steps. In general:

  1. Drop datasets of interest in the `data` folder
  2. Create a `runCases` csv which scans the contents of the data folder and enumerates all the iterations. These should probably be set before the competition is run. `bin/create_runcases` performs this action.
  3. Synchronize the executables/scripts in the `methods` folder with the contents of `methods.csv`.
  3. Create a `runStatus` Rdata object (`bin/update_runstatus`), which holds the status of every method run on every data set (success, failure, hang, etc.).
  4. Create a `results` Rdata object (`bin/create_results bias coverage ...`), which holds values derived from the input files and estimates computed by each method.

The above steps need to be run prior to spawning off any processes. Once they are in place, methods can be fit by:

  5. Queueing jobs (`bin/queue_jobs`) or runing them locally (`bin/run_locally`). This proceeds through each method, simulation setting/DGP, and iteration, creating a temp file that contains the input and storing the output in the `results` folder.
  6. Updating the `runStatus` (`bin/update_runstatus`) - this scans the results and log folders to see what has completed.
  7. Updating the `results` (`bin/update_results evaluation_functions.R`) - this uses the runStatus and current value of results to see what needs to be computed to bring the derived results in line with raw computations.

# `methods.csv`

This needs to be a csv with at least the columns:

    name,individual_effects,language

  * `name` is the same as that as the executable
  * `individual_effects` `0/1` if the method is capable of computing individual effects
  * `language` `R`/`matlab`/`stata`/`python`/etc; Each needs to be supported

Additional columns can be added to help with plotting/display.

# Executeables

As indicated above, a number of helper commands have been placed in the `bin` folder.

  * `sync_methods` - examines `methods.csv` for any added/dropped rows and modifies the `runStatus` and `results` files to match
  * `create_runcases` - scans the `data` directory and tracks the number of unique simulation settings and iterations under each
  * `update_runstatus [method ...]` - creates or updates the `runStatus` file by examining the `results` and `log` directories; optionally updates only for methods listed
  * `create_results result [result ...]` - creates a `results` object with named elements corresponding to the results listed on the command line
  * `update_results evaluation_functions.R` - updates the `results` object by using the functions defined in its first argument to collapse the data and estimate down to a single-value summary; `evaluation_functions.R` must define an object `evaluation_functions` that is a named list of functions with the prototype `function(data, results, results.ind)`, where good `results.ind` is supplied only if the method provides individual treatment effect estimates
  * `run_locally [method ...]` - runs methods against datasets on the local machine; optionally does so only for the methods listed
  * `queue_jobs [method ...]` - copies and `template.pbs` to the `job` folder, edits it, and spawns off a process for each method and each simulation setting; optionally does so only for the methods listed

# Examples

## Getting Started

    mkdir ~/acic
    cd ~/acic
    
    git clone https://github.com/vdorie/aciccomp.git
    
    ## transfer data, including x.csv and simulation subdirectories
    scp https://somewhere.edu/data.tar.gz .
    tar xzf data.tar.gz
    rm data.tar.gz
    
    ## transfer methods, including methods.csv and methods folder
    scp https://somewhere.edu/methods.tar.gz .
    tar xzf methods.csv
    rm methods.tar.gz
    
    ## configure directories
    echo 'dirs <- list(
      data    = "data",
      log     = "log",
      methods = "methods",
      results = "results",
      src     = "src",
      job     = "job")' > site_setup.R
    
    ## create resource files
    bin/create_runcases
    bin/update_runstatus
    bin/create_results bias coverage
    
    ## define evaluation procedure
    echo 'evaluation_functions <- list(
      bias = function(data, results, results.ind) {
        mean(data$mu.1[data$z == 1] - data$mu.0[data$z == 1]) - results$est
      }
      coverage = function(data, results, results.ind) {
        te <- mean(data$mu.1[data$z == 1] - data$mu.0[data$z == 1])
        results$ci_lower <= te && te <= resutls$ci_upper
      })' > evaluation_functions.R
    
    ## run methods
    bin/queue_jobs
    
    ## wait 2 days to harvest the results
    bin/update_runstatus
    bin/update_results
    
    ## transfer run status and results off of the cluster

## Changing the Data

    ## transfer in new data
    rm -rf data
    scp https://somewhere.edu/data.tar.gz .
    tar xzf data.tar.gz
    rm data.tar.gz
    
    bin/create_runcases
    rm runStatus.Rdata
    bin/update_runstatus
    bin/create_results bias coverage
    
    ## run as before

## Adding Methods
   
    ## add a method file to the methods fold, line to methods.csv
    bin/sync_methods
    
    ## run as before

## Other Issues

It can be helpful to setup R on a new, shared environment to use a local library folder:

    echo 'R_LIBS_USER="~/.R/%p-library/%v"' > ~/.Renviron
    R_PLATFORM=$(R --no-save --slave -e 'cat(R.version$platform)')
    R_VERSION=$(R --no-save --slave -e 'cat(R.version$major, ".", sub("([0-9]+)\\.[0-9]+", "\\1", R.version$minor), sep = "")')
    mkdir -p ~/.R/${R_PLATFORM}-library/${R_VERSION}
    unset R_PLATFORM
    unset R_VERSION
