![cover image](https://hvoltbb.github.io/pics/cover_pic3.jpg)

# Somatic growth of Atlantic bluefin tuna _Thunnus thynnus_ under global climate variability

![version: v0.6](https://img.shields.io/badge/version-v0.6-blue.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Platform: win linux ios](https://img.shields.io/badge/platform-win%20%7C%20linux%20%7C%20ios-lightgrey)

> Ultimately, we cannot interpret the growth of an organism except with reference to its particular environment. - R. Beverton and S. Holt

## Quick View

**Background**: Somatic growth is integral to fishery stock productivity. Under climate variability, omitting growth variability renders fishery management strategies non-optimal.

**Problem statement**: Traditionally, somatic growth is not directly measured but assessed through bio-chronologies based on proxies. This dependency presents obstacles for its general application to detect climate signals. In addition, the coarse temporal resolution of those studies are ill-adapted for studying fast processes operating on a much finer scale, such as somatic growth, and may fail to identify important climate signals due to excessive averaging.

**Contributions**: Based on a multidecadal tag-recapture database, a case study is presented to investigate the potential growth response of the Atlantic bluefin tuna (BFT) to three regionally relevant large-scale climate patterns, i.e., the North Atlantic Oscillation, Arctic Oscillation, and Pacific North America pattern. An additional simulation study is conducted to explore the effect of the overall scale and the distribution of measurement error on the detection probability of extrinsic effects and the estimation of growth parameters.

**Key takeaways**:

1. Climate variability substantially affects the growth of the Atlantic BFT
2. Complex growth responses to large-scale teleconnection patterns
3. Substantial bias in growth parameters when observation error at tag release is high[^1]
4. Integrating different growth datasets in a naive fashion is not recommended

### Data

[Tagging database](https://iccat.int/Data/Tag/_tagBFT.7z)

[NAO index](https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/nao.shtml)

[AO index](https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/ao.shtml)

[PNA index](https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/pna.shtml)

[AMO index](https://psl.noaa.gov/data/timeseries/AMO/)

Disclaimer: None of the datasets are maintained by me. Please go to the original source to obtain your own copy of the file. Questions solely related to the dataset should be directed to the data maintainers found in the links above. However, if the links provided above are broken, let me know. The format of the current version of the data files might be different from my copy, and you may need to change some lines in the [preprocessing code](src/prep.r) to take that into account.

### Source code

The source files contain four R scripts [the main program](src/CV.r), [data preprocessor](src/prep.r), [single threaded code to be run on your local machine](src/wu.r), [jobs that you may run on your local machine but preferably offloaded to a cluster](src/job.r), two header files [1](src/growth.h) and [2](src/growth_imp.h), and a [TMB script](src/v5_3.cpp)[^2].

This code has been tested on both Windows[^3] and Linux systems with 32 GB of RAM.

Warning: it is recommended to run the full program on a system with at least 32 GB of RAM. Some of the large models require a substantial amount of memory to execute. Your system may freeze if you are low on available memory. Your R may crash if the TMB program is given a wrong set of inputs. The way RStudio handles R crashes is not elegant. You will lose all unsaved changes. Preferably, put a light wrapper around base R such that whenever base R crashes it doesn't bring your IDE down with it. You can thank me later.

It is unfortunate that the name of this repo reads "climate-change-02". As one of the reviewers pointed out, correctly, that this is not "climate change" but "climate variability". So just treat this name as a string of characters, just as there's no dolphin in a MySQL server, nor is there a penguin sitting on your linux desktop.

Bug reports, feature requests, and colabs are welcome.

[^1]: Note that biases are properties of estimators. Points 3 and 4 are specific about just the estimation method used in this study. Building a robust climate-ready growth increment model is an area of current research. I am hopeful that these restrictions can be lifted shortly (probably next year depending on progress).
[^2]: The state of this repo is currently frozen such that it provides a snapshot view of the model development stage as presented in the publication. The code, however, has been continuously developed with new features added and tested. The next version of the code will appear in a different repo. Future updates of this repo are limited to minor changes, such as typo or bug fixes.
[^3]: Although the code provided here has been tested on `R 4.1.2` with `Rtools40`, I have encountered some incompatibilities of the `Eigen` package (the matrix backend of the `TMB` package) with the build tools provided by `Rtools`, which is essential for code compilation on Windows. It is not recommended to further develop the code on Windows. Meanwhile, the Windows Linux Subsystem (wsl) works fine for me and may be used instead if you do not have a standalone Linux machine.
