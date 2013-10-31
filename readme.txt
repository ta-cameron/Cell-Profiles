Introduction

This R script will take a series of fluorescent profiles acquired in ImageJ/Fiji (or MicrobeTracker), arrange them from shortest to longest, then use Hadley Wickam’s ggplot2 to plot the profiles out as a vertical stack. Each horizontal line represents the fluorescence intensity profile of a single cell. By ordering cells according to length, general localization patterns over the cell cycle can be easily visualized.

Usage

This script runs as a function loaded manually into your current R session. It accepts ImageJ ROI “multiplot” output data (a wide table with paired length and intensity data); see “profiles.csv” for example data. The ordered output data can be plotted with ggplot2, as detailed in “example plots.R”.