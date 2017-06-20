#==================run this each new R session==================
packages <- c("ggplot2", "reshape2", "scales", "RColorBrewer", "grid", "devtools", "Cairo")

#force install new packages
#install.packages(packages,dep=TRUE)

#ipak from stevenworthington / GitHub
# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.

ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}
ipak(packages)

#display.brewer.all()

devtools::install_github("ta-cameron/Cell-Profiles")
library(cellProfiles)
library(ggplot2)

packages <- c("ggplot2","RColorBrewer","devtools","grid", "scales", "Cairo")
sapply(packages, require, character.only=TRUE)

#================================================================
#INTENSITY HEAT MAPS OVER CELL LENGTH
#================================================================

#commands for changing to a specific directory & load data tables
#assumes there is no header... use header=TRUE if there is one
#setwd("")

dtable <- read.table(system.file("extdata", "ftsZ_profiles.csv", package = "cellProfiles"),header=FALSE,sep=",")

#optional:
#slightly truncate the poles to remove dark polar bands (artifact from width of profile line)
#dtable <- cellProfileTruncate(data=raw_table,1)

#================================================================
#the following two graphs are general templates suitable for most data & purposes (in centered, and non-centered versions)
#================================================================

#CENTERED -norm
#Note: uses interpolated data set: profileResults$interpolated
# For non-interpolated data, use: profileResults$native

profileResults <- cellProfiles(data=dtable, contrast="norm", position="center", align="orient")

g <- ggplot(data=profileResults$interpolated,
            mapping=aes(x=x, y=y, fill=intensity, color=intensity) ) +
  geom_tile() +
  #geom_raster() +
	scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.25), trans="reverse") +
	scale_x_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(0,1)) +
  theme_classic() +
  theme(axis.line.x = element_line(color="grey40", size = 0.7),
        axis.line.y = element_line(color="grey40", size = 0.7), 
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradientn(colours=rev(brewer.pal(9, "YlGnBu")), guide="colorbar", na.value=NA) +
  scale_color_gradientn(colours=rev(brewer.pal(9,"YlGnBu")), guide="colorbar", na.value=NA) +
	labs(x="distance from midcell (μm)", y="fraction of cell cycle", title="1- position=\"center\"")

dev.new(width=4.9,height=3.4, noRStudioGD=TRUE)
g


#NOT CENTERED -norm
#Note: uses actual resolution of data set: profileResults$native
# For interpolated data, use: profileResults$interpolated

profileResults<-cellProfiles(data=dtable,position="left",contrast="norm")

g2 <- ggplot(data=profileResults$native,
             mapping=aes(x=x, y=y, fill=intensity, color=intensity ) ) +
  geom_tile() +
  scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.25), trans="reverse") +
  scale_x_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(0,1)) +
  theme_classic() +
  theme(axis.line.x = element_line(color="grey40", size = 0.7),
        axis.line.y = element_line(color="grey40", size = 0.7), plot.title = element_text(hjust = 0.5) ) +
  scale_fill_gradientn(colours=rev(brewer.pal(9, "YlGnBu")), guide="colorbar", na.value=NA) +
  scale_color_gradientn(colours=rev(brewer.pal(9,"YlGnBu")), guide="colorbar", na.value=NA) +
  labs(x="length of cell (μm)", y="fraction of cell cycle", title="2- position=\"left\"")

dev.new(width=4.9,height=3.4, noRStudioGD=TRUE)
g2

#Additional useful notes:
#
#These graphs can be saved directly as png, tiff, or pdf files:
# (this will place the last drawn ggplot graph into your current working directory)
#
# ggsave("cell_profiles.png", width=4.86,height=3.4, units="in", dpi=300, type="cairo-png")
# ggsave("cell_profiles.tiff", width=4.86,height=3.4, units="in", dpi=300, type="cairo")
# ggsave("cell_profiles.pdf", width=4.86,height=3.4, units="in", dpi=300, device=cairo_pdf)
#
#
# ONLY when using interpolated data, geom_raster() can be used instead of geom_tile().
#  This will reduce file size & drawing time when dealing with many cells.
#  However, some PDF-viewing programs (such as Apple's Preview app) do not display or handle these files correctly, so I only recommend this for internal use.


#================================================================
#EXAMPLES OF GRAPHING OPTIONS
#================================================================
#The following graphs demonstrate the different graphing options available
#Select from here to the bottom and execute to view them all
raw_table <- read.table(system.file("extdata", "ftsZ_profiles.csv", package = "cellProfiles"),header=FALSE,sep=",")
dtable <- cellProfileTruncate(data=raw_table,0)
vplayout <- function(x,y) viewport(layout.pos.row=x,layout.pos.col=y)

dev.new(width=8,height=3.25, noRStudioGD=TRUE)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,2, heights = unit(c(0.25,3,3),"null"))))
grid.text("Position Options", vp=viewport(layout.pos.row=1,layout.pos.col=1:2),gp=gpar(fontsize=16))
print(g,vp=vplayout(2,1))
print(g2,vp=vplayout(2,2))

#================================================================
#3- native
#================================================================
profileResults <- cellProfiles(data=dtable)

g <- g + theme_classic(base_size = 10) + theme(axis.line.x = element_line(color="grey40", size = 0.7), axis.line.y = element_line(color="grey40", size = 0.7), plot.title = element_text(hjust = 0.5))

#the only thing changing for most of these plots is the underlying data... so I am just recycling the settings of the first plot here
pt3 <- g %+% profileResults$interpolated + labs(title="3- contrast=\"native\"\n")

#================================================================
#4- max range=c(0,1)
#================================================================
profileResults<-cellProfiles(data=dtable,contrast="max",range=c(0,1))

pt4 <- g %+% profileResults$interpolated + labs(title="4- contrast=\"max\",\n range=c(0,1)")

#================================================================
#5- norm range=c(0,1)
#================================================================
profileResults<-cellProfiles(data=dtable,contrast="norm",range=c(0,1))

pt5 <- g %+% profileResults$interpolated + labs(title="5- contrast=\"norm\",\n range=c(0,1)")

#================================================================
#6- norm range=c(0.005,0.995)
#================================================================
profileResults<-cellProfiles(data=dtable,contrast="norm",range=c(0.005,0.995))

pt6 <- g %+% profileResults$interpolated + labs(title="6- contrast=\"norm\",\n range=c(0.005,0.995)")

#================================================================
#7- norm range=c(0.02,0.98)
#================================================================
profileResults<-cellProfiles(data=dtable,contrast="norm",range=c(0.02,0.98))

pt7 <- g %+% profileResults$interpolated + labs(title="7- contrast=\"norm\",\n range=c(0.02,0.98)")

#================================================================
#8- norm range=c(0.05,0.95)
#================================================================
profileResults<-cellProfiles(data=dtable,contrast="norm",range=c(0.05,0.95))

pt8 <- g %+% profileResults$interpolated + labs(title="8- contrast=\"norm\",\n range=c(0.05,0.95)")

dev.new(width=12,height=6.4, noRStudioGD=TRUE)
grid.newpage()
pushViewport(viewport(layout=grid.layout(3,3, heights = unit(c(0.25,3,3),"null"))))
grid.text("Contrast Options", vp=viewport(layout.pos.row=1,layout.pos.col=1:3),gp=gpar(fontsize=16))
print(pt3,vp=vplayout(2,1))
print(pt5,vp=vplayout(2,2))
print(pt6,vp=vplayout(2,3))
print(pt4,vp=vplayout(3,1))
print(pt7,vp=vplayout(3,2))
print(pt8,vp=vplayout(3,3))



#================================================================
#9- native
#================================================================
profileResults <- cellProfiles(data=dtable,contrast="max")

pt9 <-g %+% profileResults$interpolated + labs(title="9- align=\"native\",\n contrast=\"max\"")

#================================================================
#10- orient
#================================================================
profileResults <- cellProfiles(data=dtable,align="orient",contrast="max")

pt10 <-g %+% profileResults$interpolated + labs(title="10- align=\"orient\",\n contrast=\"max\"")

#================================================================
#11- random
#================================================================
profileResults <- cellProfiles(data=dtable,align="random",contrast="max")

pt11 <- g %+% profileResults$interpolated + labs(title="11- align=\"random\",\n contrast=\"max\"")

#================================================================
#12- orient / reverse
#================================================================
profileResults <- cellProfiles(data=dtable,align="orient",contrast="max",reverse=TRUE)

pt12 <- g %+% profileResults$interpolated + labs(title="12- align=\"orient\", reverse=TRUE, \ncontrast=\"max\"")

#================================================================
#13- reuse / reverse #1
#================================================================
profileResults <- cellProfiles(data=dtable,align=profileResults,contrast="max",reverse=TRUE)

pt13 <- g %+% profileResults$interpolated + labs(title="13- align:(reused), reverse=TRUE, \ncontrast=\"max\"")

#================================================================
#14- reuse / reverse #2
#================================================================
profileResults <- cellProfiles(data=dtable,align=profileResults,contrast="max",reverse=TRUE)

pt14 <-g %+% profileResults$interpolated + labs(title="14- align:(reused), reverse=TRUE, \ncontrast=\"max\"")

dev.new(width=12,height=6.4, noRStudioGD=TRUE)
grid.newpage()
pushViewport(viewport(layout=grid.layout(3,3, heights = unit(c(0.25,3,3),"null"))))
grid.text("Alignment Options", vp=viewport(layout.pos.row=1,layout.pos.col=1:3),gp=gpar(fontsize=16))
print(pt9,vp=vplayout(2,1))
print(pt10,vp=vplayout(2,2))
print(pt11,vp=vplayout(2,3))
print(pt12,vp=vplayout(3,1))
print(pt13,vp=vplayout(3,2))
print(pt14,vp=vplayout(3,3))

#================================================================
#15- violin plots of length quantiles \n
# wherein I badly abuse geom_violin for something it was never intended
#================================================================

#contrast="max" yields the greatest variability, but native or norm may be more representative of the average
profileResults <- cellProfiles(data=dtable, contrast="max")

#first determine cell lengths --> cell
profileResults$interpolated$cell <- NA
profileResults$proportional$cell <- NA
for (i in unique(profileResults$interpolated$y)){
  cell <- max(profileResults$native[profileResults$native$y == i & !is.na(profileResults$native$intensity), "x"], na.rm=TRUE) -
    min(profileResults$native[profileResults$native$y == i & !is.na(profileResults$native$intensity), "x"], na.rm=TRUE)
  profileResults$interpolated[profileResults$interpolated$y == i, "cell"] <- cell
  profileResults$proportional[profileResults$proportional$y == i, "cell"] <- cell
}
rm(cell)

#next categorize into five discrete quantiles --> cell_cat
cutoffs <- quantile(profileResults$interpolated$cell, probs=c(0,0.2,0.4,0.6,0.8), na.rm=TRUE )
profileResults$interpolated$cell_cat <- NA
profileResults$proportional$cell_cat <- NA
for(i in unique(profileResults$interpolated$cell)){
  cell_cat <- cutoffs[[ which( cutoffs==max(cutoffs[(cutoffs <= i )]) ) ]]
  profileResults$interpolated[profileResults$interpolated$cell == i, "cell_cat"] <- cell_cat
  profileResults$proportional[profileResults$proportional$cell == i, "cell_cat"] <- cell_cat
}
rm(cell_cat)

pt15 <- ggplot(data=profileResults$proportional, mapping=aes(x=factor(cell_cat), y=x*cell_cat, weight=(intensity)^2 ) ) +
    #also look into using scale="count" or "area"
		geom_violin(adjust=0.5, bw="sj", scale="width") +
  	theme_classic(base_size=10) +
    theme(axis.line.x = element_line(color="grey40", size = 0.7), axis.line.y = element_line(color="grey40", size = 0.7), plot.title = element_text(hjust = 0.5)) +
		labs(x="cell length quantile\n(widths indicate relative min/max intensity)",y="cell length (um)",title="15- cell length quantile profiles")


#================================================================
#16- line plots with rank-colored lines
#================================================================
profileResults <- cellProfiles(data=dtable,contrast="norm",align="orient")

pt16 <- ggplot(data=profileResults$proportional, mapping=aes(x=x, y=intensity, group=factor(y), color=y) ) +
  geom_line(stat="identity", position="identity", alpha=0.3) +
  coord_cartesian(xlim=c(0,1)) +
  theme_classic(base_size=10) +
  theme(axis.line.x = element_line(color="grey40", size = 0.7), axis.line.y = element_line(color="grey40", size = 0.7), plot.title = element_text(hjust = 0.5), legend.position="none") +
  scale_color_gradient2(low="blue", mid="purple", high="red", midpoint=0.5) +
  labs(x="fraction of cell length", y="normalized intensity", title="16- line plots colored by\nordered length of cells")

#================================================================
#17- line plots split between short\nmedium & long cells (simple mean & confidence)
#================================================================

#the next lines add a third colum that defines how to split the cells into categories
profileResults$proportional$split <- round(profileResults$proportional$y*2)

pt17 <- ggplot(data=profileResults$proportional, mapping=aes(x=x, y=intensity, group=factor(split), color=factor(split), fill=factor(split))) +
  geom_line(mapping=aes(group=factor(y)), stat="identity", position="identity", alpha=0.8) +
  stat_summary(fun.data="mean_cl_boot", geom="smooth", alpha=0.5, color=NA, fill="black") +
  stat_summary(fun.y="mean", geom="line", size=1.25) +
  coord_cartesian(xlim=c(0,1)) +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  theme_classic(base_size=10) +
  theme(axis.line.x = element_line(color="grey40", size = 0.7), axis.line.y = element_line(color="grey40", size = 0.7), plot.title = element_text(hjust = 0.5), legend.position="none") +
  labs(x="fraction of cell length", y="normalized intensity", title="17- line plots split between short, medium,\n& long cells (simple mean & confidence)")

#note: mean_cl_boot (based on smean.cl.boot from the Hmisc package) is "a very fast implementation of the basic nonparametric bootstrap for obtaining confidence limits for the population mean without assuming normality"

#================================================================
#18- gaussian kernel density \nplot of cell lengths
#================================================================

profileResults <- cellProfiles(data=dtable, position = "left")
lengths <- rep(NA, length(profileResults$fliplist))
n <- 1
for (i in unique(profileResults$native$y)){
  lengths[n] <- max(profileResults$native[profileResults$native$y == i & !is.na(profileResults$native$intensity), "x"], na.rm=TRUE)
  n <- n+1
}
lengths <- as.data.frame(lengths)

pt18 <- ggplot(data=lengths, mapping=aes(x=lengths)) +
		layer(geom="density", stat="density", position="identity") +
    theme_classic(base_size=10) +
    theme(axis.line.x = element_line(color="grey40", size = 0.7), axis.line.y = element_line(color="grey40", size = 0.7), plot.title = element_text(hjust = 0.5), legend.position="none") +
    scale_x_continuous(expand=c(0,0)) +
		scale_y_continuous(expand=c(0,0), limits=c(0,1)) +
		labs(x="cell length (um)", y="probability density", title="18- gaussian kernel density \nplot of cell lengths")

dev.new(width=8,height=6, noRStudioGD=TRUE)
grid.newpage()
pushViewport(viewport(layout=grid.layout(3,2, heights = unit(c(0.25,3,3),"null"))))
grid.text("Additional Graph Types and Examples", vp=viewport(layout.pos.row=1,layout.pos.col=1:2),gp=gpar(fontsize=16))
print(pt15,vp=vplayout(2,1))
cat(" ---> Yes, those errors are to be expected... this isn't what violin plots were intended for. <---")
print(pt16,vp=vplayout(2,2))
cat(" ---> This one takes a moment to draw... be patient. <---")
print(pt17,vp=vplayout(3,1))
print(pt18,vp=vplayout(3,2))

#================================================================
#19 & 20 - two-color colocalization
#================================================================

# first load two data sets [NOT SHOWN]
# next run the cell profiles function for each of them
# here I am simply reusing the same data set, but oriented to the left or right
# obviously you will not want to *exactly* copy this to compare two *different* data sets

# the range / contrast of each channel will need to be carefully adjusted
profileResults1 <- cellProfiles(data=dtable,contrast="norm", align="orient", range=c(0.06,0.98))
profileResults2 <- cellProfiles(data=dtable,contrast="norm",align="orient", reverse=TRUE, range=c(0.06,0.98))

#rescale intensities to 0-1
profileResults1$interpolated$r <- rescale(profileResults1$interpolated$intensity)
profileResults1$interpolated$g <- rescale(profileResults2$interpolated$intensity)

#troubleshooting
# if your two data sets do not have exactly the same number of points, you will get issues and errors!
# MicrobeTracker occassionally produces such errors, so double check your data carefully if you use it.
# it may help to look at the final data sets, which can be exported in wide-format with the code below
# write.csv(as.matrix(profileResults1$lim_or_dtable$r), file="profile table 1 (red).csv")
# write.csv(as.matrix(profileResults2$lim_or_dtable$g), file="profile table 2 (green).csv")

#straight red/green probably works best, as straight blue too closely matches the black background

#convert values to rgb
profileResults1$interpolated$rgb[!is.na(profileResults1$interpolated$r)] <- rgb(
  na.omit(profileResults1$interpolated$r),
  na.omit(profileResults1$interpolated$g),
  0 )

pt19 <- g %+% profileResults1$interpolated + aes(fill=rgb, color=rgb) +
  scale_fill_identity() + scale_color_identity() +
  labs(x="distance from midcell (μm)", y="fraction of cell cycle", title="19- red/green colocalization")

#alternate light blue / orange colors for a more color-blind friendly palette
# in case you don't want to be evil

#convert values to rgb
profileResults1$interpolated$rgb[!is.na(profileResults1$interpolated$r)] <- rgb(
  na.omit(profileResults1$interpolated$r),
  0.5*na.omit(profileResults1$interpolated$g)+0.5*na.omit(profileResults1$interpolated$r),
  na.omit(profileResults1$interpolated$g) )

pt20 <- g %+% profileResults1$interpolated + aes(fill=rgb, color=rgb) +
  scale_fill_identity() + scale_color_identity() +
  labs(x="distance from midcell (μm)", y="fraction of cell cycle", title="20- blue/orange colocalization")

dev.new(width=7,height=3.25, noRStudioGD=TRUE)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,2, heights = unit(c(0.25,3),"null"))))
grid.text("Two-color colocalization graphs", vp=viewport(layout.pos.row=1,layout.pos.col=1:2),gp=gpar(fontsize=16))
print(pt19,vp=vplayout(2,1))
print(pt20,vp=vplayout(2,2))

#================================================================
# the following graphs are for documentation purposes only
#================================================================
# g <- g + theme(base_size=10)
# g2 <- g2 + theme(base_size=10)

# dev.new(width=8,height=3)
# grid.newpage()
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(g,vp=vplayout(1,1))
# print(g2,vp=vplayout(1,2))

# dev.new(width=8,height=3)
# grid.newpage()
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(pt3,vp=vplayout(1,1))
# print(pt4,vp=vplayout(1,2))

# dev.new(width=8,height=3)
# grid.newpage()
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(pt5,vp=vplayout(1,1))
# print(pt6,vp=vplayout(1,2))

# dev.new(width=8,height=3)
# grid.newpage()
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(pt7,vp=vplayout(1,1))
# print(pt8,vp=vplayout(1,2))

# dev.new(width=8,height=3)
# grid.newpage()
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(pt9,vp=vplayout(1,1))
# print(pt11,vp=vplayout(1,2))

# dev.new(width=8,height=3)
# grid.newpage()
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(pt10,vp=vplayout(1,1))
# print(pt12,vp=vplayout(1,2))

# dev.new(width=8,height=3)
# grid.newpage()
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(pt13,vp=vplayout(1,1))
# print(pt14,vp=vplayout(1,2))

# dev.new(width=8,height=3)
# grid.newpage()
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(pt18,vp=vplayout(1,1))
# print(pt15,vp=vplayout(1,2))

# dev.new(width=8,height=3)
# grid.newpage()
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(pt16,vp=vplayout(1,1))
# print(pt17,vp=vplayout(1,2))


# #demonstration of dark vs light background
# profileResults<-cellProfiles(data=dtable,contrast="norm")

# pt19<- ggplot() + layer(data=profileResults$lim_or_dtable, mapping=aes(x=x,y=y/profileResults$ncol, fill=intensity, color=intensity), geom="tile", width=1.01*(profileResults$data[2,1]-profileResults$data[1,1]), size=0.25) + scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.25), trans="reverse") + scale_fill_gradientn(colours=rev(brewer.pal(9, "YlGnBu"))) + scale_color_gradientn(colours=rev(brewer.pal(9,"YlGnBu"))) + labs(x="distance from midcell (um)", y="fraction of cell cycle",title="dark background") + theme(panel.background=element_rect(fill="#061542"), panel.grid.minor=element_line(color="#225EA8"), panel.grid.major=element_line(color="#225ea8"))

# pt20<- ggplot() + layer(data=profileResults$lim_or_dtable, mapping=aes(x=x,y=y/profileResults$ncol, fill=intensity, color=intensity), geom="tile", width=1.01*(profileResults$data[2,1]-profileResults$data[1,1]), size=0.25) + scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.25), trans="reverse") + theme_bw() + scale_fill_gradientn(colours=rev(brewer.pal(9, "YlGnBu"))) + scale_color_gradientn(colours=rev(brewer.pal(9,"YlGnBu"))) + labs(x="distance from midcell (um)", y="fraction of cell cycle",title="light background")

# dev.new(width=7,height=3)
# grid.newpage()
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(pt19,vp=vplayout(1,1))
# print(pt20,vp=vplayout(1,2))
