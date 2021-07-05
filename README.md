![cover image](https://hvoltbb.github.io/pics/cover_pic3.jpg)
# Tracking the effects of climate change on marine life - Case study 2
[4/23/21 update: article submitted. Text for each section will appear after being finalized]
[6/15/21 update: major revision]

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

Disclaimer: Neither of the datasets are owned/maintained by me. Please go to the original source to obtain your own copy of the file. Questions solely related to the dataset should be directed to the data maintainers found in the links above. However, if the links provided above are broken, let me know. The format of the current version of the data files might be different from my copy, and you may need to change some lines in the [preprocessing code](src/prep.R) to take that into account. 

My own copy of the data files will be conditionally provided through PM (latency: hours) or through an automatic email server by clicking [here](mailto:eidotog@gmail.com?subject=XxCLIMATE02xX&body=Do%20not%20modify%20the%20subject%20line.%20Not%20monitored.) (latency: secs. Currently offline. Outstanding requests will be fullfilled once online). You may want to check your spam folder for the reply because replying an email in milliseconds isn't humanly possible and it will be flagged as spam the majority of the times. Note that I don't monitor these data requests. Once a request is fullfilled, the message will be permanantly removed from the server. No personal information will be collected by me.

### Source code
The source files contain R scripts [1](src/v3.R) and [2](src/prep.R), two header files [1](src/growth.h) and [2](src/growth_imp.h), and a [TMB script](src/v5.cpp).

The R code has been annotated and sectioned according to RStudio style. It is recommended to use RStudio to view the code.

This code has been tested on both Windows and Linux systems with 8+ GB of RAM. 

Warning: it is recommended to run the full program on a system with 8+ GB of RAM. Some of the large models require a substantial amount of memory to execute. Your system may freeze if you are low on available memory. Your R may crash if the TMB program is given a wrong set of inputs. You have been warned. Do save your work first.

Bug reports, feature requests, and colabs are welcome. 

