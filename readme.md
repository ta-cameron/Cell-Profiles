### Introduction

This R package will take a series of fluorescent profiles acquired in ImageJ/Fiji (or MicrobeTracker), arrange them from shortest to longest, then use Hadley Wickam’s `ggplot()` to plot the profiles out as a vertical stack. Each horizontal line represents the fluorescence intensity profile of a single cell. By ordering cells according to length, general localization patterns over the cell cycle can be easily visualized.

### Installation

Install directly from Github using the devtools package:
```R
install.packages(c("devtools","reshape2","scales"))
library(devtools)

devtools::install_github("ta-cameron/Cell-Profiles")
library(cellProfiles)
```

**Vignette update**

Apparently the install_github recently stopped compiling vignettes, for reasons unknown to me. To compile the vignette, you can clone the Git repository and build the package yourself. Alternatively, the vignette can also be viewed here in this project's [wiki page](https://github.com/ta-cameron/Cell-Profiles/wiki/Vignette).  

**Old vignette instructions**

Note, if you want to install with the vignette (**recommended**), you'll need to load several additional packages in order to compile the vignette file:
```R
packages <- c("ggplot2","RColorBrewer","devtools", "scales", "knitr")
install.packages(packages)

packages <- c("ggplot2","RColorBrewer","devtools","scales", "grid")
sapply(packages, require, character.only=TRUE)

devtools::install_github("ta-cameron/Cell-Profiles", build_vignettes = TRUE)
library(cellProfiles)
```


### Usage

This package loads two functions, `cellProfiles()` and `cellProfilesTruncate()`. `cellProfiles()` is the primary function of interest, and it will prepare the raw data for graphing with `ggplot()`. It accepts ImageJ ROI “multiplot” output data (a wide table with paired length and intensity data). `cellProfilesTruncate()` is an optional helper function to trim off data points from each end of all cells. 

See vignette / [wiki](https://github.com/ta-cameron/Cell-Profiles/wiki/Vignette) for detailed usage examples: `vignette("cellProfiles")`

### Changes

The are some incompatible differences between the output of the older function-styled versions and the newer package version. The upside is that the new data output is more intuitively named and less complicated to graph. Additionally, the method of reusing prior profile orientations has changed. Please refer to the help file & vignette / [wiki](https://github.com/ta-cameron/Cell-Profiles/wiki/Vignette) for current usage instructions.
