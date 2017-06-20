### Introduction

This R package will take a series of fluorescent profiles acquired in ImageJ/Fiji (or MicrobeTracker), arrange them from shortest to longest, then use Hadley Wickam’s `ggplot()` to plot the profiles out as a vertical stack. Each horizontal line represents the fluorescence intensity profile of a single cell. By ordering cells according to length, general localization patterns over the cell cycle can be easily visualized.

### Installation

Install directly from Github using the devtools package:
```R
install.packages("devtools")
library(devtools)
devtools::install_github("ta-cameron/Cell-Profiles")
```

### Usage

This package loads two functions, `cellProfiles()` and `cellProfilesTruncate()`. `cellProfiles()` is the primary function of interest, and it will prepare the raw data for graphing with `ggplot()`. It accepts ImageJ ROI “multiplot” output data (a wide table with paired length and intensity data). `cellProfilesTruncate()` is an optional helper function to trim off data points from each end of all cells. 

See vignette for detailed usage examples: `vignette("cellProfiles")`

### Changes

The are some incompatible differences between the output of the older function-styled versions and the newer package version. The upside is that the new data output is more intuitively named and less complicated to graph. Additionally, the method of reusing prior profile orientations has changed. Please refer to the help file & vignette for current usage instructions.