"
Author/maintainer: Can Zhou [eidotog@gmail.com]
Date: Oct 4 2021
Version: 0.2
"

if (!require("snow")) {
    install.packages("snow")
}
if (!require("snowfall")) {
    install.packages("snowfall")
}

## Data
rm(list = ls())
MULTI <- FALSE
source("prep.r")

## Specs
PAR1 <- 0 # 0: von Bertalanffy, 1: G. West, 2: Gompertz, 3: Logistic, 4: Generalized logistic, 5+: General von Bertalanffy
pickone <- 1

## Warm-up run
TMB::compile("v5_3.cpp")
DEBUG <- FALSE
source("wu.r")

# Verification: here, it should spit out 14255.44, which is the AIC for the baseline VB model.

## Jobs
source("job.r")
nrep <- 20

## Snow
# 5 jobs will be spawned here.
# At least 32 GB of RAM and 6 cpu cores are recommended.
# The run time for the following code ranges from 1h to 12hs depending on the complexity of the model on a multicore desktop computer. This code is both CPU and RAM intensive. If the following code soaks up all the CPU and RAM resources on your machine, the run time will be substantialy longer, and a bigger machine is recommended.
# The following code is configured to run on your local machine.
# To offload the job to a cluster, you need to change the cluster type below.
# Please refer to the documentationn of the snowfall package and the configuration of your cluster.

sfInit(parallel = T, cpus = 5, type = "SOCK")
sfExportAll()

tri <- sfLapply(1:5, fun = job)

sfStop()

## CV score
x <- unlist(tri)
cat(mean(x), "\n")