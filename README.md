![cover image](https://hvoltbb.github.io/pics/cover_pic3.jpg)
# Tracking the effects of climate change on marine life - case study 2


![version: v0.5.3](https://img.shields.io/badge/version-v0.5.3-blue.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Platform: win linux ios](https://img.shields.io/badge/platform-win%20%7C%20linux%20%7C%20ios-lightgrey)

[4/23/21 update: article submitted. Text for each section will appear after being finalized]

[6/15/21 update: major revision]

[9/23/21 update: accepted]

## Quick View
**Background**: 
**Problem statement**: 
**Contributions**: 

## Full Story

### Issues
1.
2.
### Data
[Tagging database](https://iccat.int/Data/Tag/_tagBFT.7z)

[NAO index](https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/nao.shtml)

[AO index](https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/ao.shtml)

[PNA index](https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/pna.shtml)

[AMO index](https://psl.noaa.gov/data/timeseries/AMO/)

Disclaimer: None of the datasets are owned/maintained by me. Please go to the original source to obtain your own copy of the file. Questions solely related to the dataset should be directed to the data maintainers found in the links above. However, if the links provided above are broken, let me know. The format of the current version of the data files might be different from my copy, and you may need to change some lines in the [preprocessing code](src/prep.r) to take that into account. 

My own copy of the data files will be conditionally provided through PM (latency: hours) or through an automatic email server by clicking [here](mailto:eidotog@gmail.com?subject=XxCLIMATE02xX&body=Do%20not%20modify%20the%20subject%20line.%20Not%20monitored.) (latency: secs. Currently offline. Outstanding requests will be fullfilled once online). You may want to check your spam folder for the reply because replying an email in milliseconds isn't humanly possible and it will be flagged as spam the majority of the times. Note that I don't monitor these data requests. Once a request is fullfilled, the message will be permanantly removed from the server. No personal information will be collected by me.

### Source code
The source files contain four R scripts [the main program](src/CV.r), [data preprocessor](src/prep.r), [single threaded code to be run on your local machine](src/wu.r), [jobs that you may run on your local machine but preferably offloaded to a cluster](src/job.r), two header files [1](src/growth.h) and [2](src/growth_imp.h), and a [TMB script](src/v5_3.cpp).

This code has been tested on both Windows and Linux systems with 32 GB of RAM. 

Warning: it is recommended to run the full program on a system with at least 32 GB of RAM. Some of the large models require a substantial amount of memory to execute. Your system may freeze if you are low on available memory. Your R may crash if the TMB program is given a wrong set of inputs. The way RStudio handles R crashes is not elegant. You will lose all unsaved changes. Preferably, put a light wrapper around base R such that whenever base R crashes it doesn't bring your IDE with it. You can thank me later.

It is unfortunate that the name of this repo reads "climate-change-02". As one of the reviewers pointed out, correctly, that this is not "climate change" but "climate variability". So just treat this name as a string of characters, just as there's no dolphin in a MySQL server, nor is there a penguin sitting on your linux desktop.

Bug reports, feature requests, and colabs are welcome. 
