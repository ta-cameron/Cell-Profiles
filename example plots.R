#==================run this each new R session==================
packages <- c("ggplot2", "reshape2", "hexbin", "quantreg", "scales", "RColorBrewer", "grid")

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

#also load the cellProfiles function by executing the entire 'cell profiles function' file

#================================================================
#INTENSITY HEAT MAPS OVER CELL LENGTH
#be sure to load the cellProfiles function first!
#================================================================

#commands for changing to a specific directory & load data tables
#assumes there is no header... use header=TRUE if there is one
#setwd("")
raw_table<-read.table("profiles.csv",header=FALSE,sep=",")

#slightly truncate the poles to remove dark polar bands (artifact from width of profile line)
dtable<-cellProfileTruncate(data=raw_table,0)

#================================================================
#the following two graphs are general templates suitable for most data & purposes (in centered, and non-centered versions)
#================================================================

#CENTERED -norm
profileResults<-cellProfiles(data=dtable,contrast="norm")

g<-ggplot() + layer(data=profileResults$lim_or_dtable, mapping=aes(x=x,y=y/profileResults$ncol, fill=intensity, color=intensity), geom="tile", width=1.01*(profileResults$data[2,1]-profileResults$data[1,1]), size=0.25) + scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.25), trans="reverse") + theme_bw() + scale_fill_gradientn(colours=rev(brewer.pal(9, "YlGnBu"))) + scale_color_gradientn(colours=rev(brewer.pal(9,"YlGnBu"))) + labs(x="distance from midcell (um)", y="fraction of cell cycle",title="1- position=\"center\"")

dev.new(width=4.86,height=3.4)
g


#NOT CENTERED -norm
profileResults<-cellProfiles(data=dtable,position="left",contrast="norm")

g2<-ggplot() + layer(data=profileResults$lim_or_dtable, mapping=aes(x=x,y=y/profileResults$ncol, fill=intensity, color=intensity), geom="tile", width=1.01*(profileResults$data[2,1]-profileResults$data[1,1]), size=0.25) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.25), trans="reverse") + theme_bw() + scale_fill_gradientn(colours=rev(brewer.pal(9, "YlGnBu"))) + scale_color_gradientn(colours=rev(brewer.pal(9,"YlGnBu"))) + labs(x="cell length (um)",y="fraction of cell cycle",title="2- position=\"left\"")

dev.new(width=4.86,height=3.4)
g2




#================================================================
#EXAMPLES OF GRAPHING OPTIONS 
#================================================================
#The following graphs demonstrate the different graphing options available
#Select from here to the bottom and execute to view them all
raw_table<-read.table("profiles.csv",header=FALSE,sep=",")
dtable<-cellProfileTruncate(data=raw_table,0)
vplayout<-function(x,y) viewport(layout.pos.row=x,layout.pos.col=y)

standard_theme	<-theme_bw(base_size=10)
standard_labs	<-labs(x="distance from midcell (um)",y="fraction of cell cycle")
standard_y		<-scale_y_continuous(expand=c(0,0),breaks=seq(0,1,0.25),trans="reverse")
standard_fill	<-scale_fill_gradientn(colours=rev(brewer.pal(9,"YlGnBu")))
standard_color	<-scale_color_gradientn(colours=rev(brewer.pal(9,"YlGnBu")))
standard_layer	<-function(x1,x2) layer(geom="tile",width=1.01*(x1-x2),size=0.25)

dev.new(width=8,height=3.25)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,2, heights = unit(c(0.25,3,3),"null"))))
grid.text("Position Options", vp=viewport(layout.pos.row=1,layout.pos.col=1:2),gp=gpar(fontsize=16))
print(g,vp=vplayout(2,1))
print(g2,vp=vplayout(2,2))

#================================================================
#3- native
#================================================================
profileResults<-cellProfiles(data=dtable)

pt3<-ggplot(profileResults$lim_or_dtable, aes(x=x, y=y/profileResults$ncol, fill=intensity, color=intensity)) + standard_layer(profileResults$data[2,1], profileResults$data[1,1]) + standard_theme + standard_y + standard_fill + standard_color + standard_labs + labs(title="3- contrast=\"native\"\n")

#================================================================
#4- max range=c(0,1)
#================================================================
profileResults<-cellProfiles(data=dtable,contrast="max",range=c(0,1))

pt4<-ggplot(profileResults$lim_or_dtable, aes(x=x, y=y/profileResults$ncol, fill=intensity, color=intensity)) + standard_layer(profileResults$data[2,1], profileResults$data[1,1]) + standard_theme + standard_y + standard_fill + standard_color + standard_labs + standard_labs+labs(title="4- contrast=\"max\",\n range=c(0,1)")

#================================================================
#5- norm range=c(0,1)
#================================================================
profileResults<-cellProfiles(data=dtable,contrast="norm",range=c(0,1))

pt5<-ggplot(profileResults$lim_or_dtable, aes(x=x, y=y/profileResults$ncol, fill=intensity, color=intensity)) + standard_layer(profileResults$data[2,1], profileResults$data[1,1]) + standard_theme + standard_y + standard_fill + standard_color + standard_labs + labs(title="5- contrast=\"norm\",\n range=c(0,1)")

#================================================================
#6- norm range=c(0.005,0.995)
#================================================================
profileResults<-cellProfiles(data=dtable,contrast="norm",range=c(0.005,0.995))

pt6<-ggplot(profileResults$lim_or_dtable, aes(x=x, y=y/profileResults$ncol, fill=intensity, color=intensity)) + standard_layer(profileResults$data[2,1], profileResults$data[1,1]) + standard_theme + standard_y + standard_fill + standard_color + standard_labs + labs(title="6- contrast=\"norm\",\n range=c(0.005,0.995)")

#================================================================
#7- norm range=c(0.02,0.98)
#================================================================
profileResults<-cellProfiles(data=dtable,contrast="norm",range=c(0.02,0.98))

pt7<-ggplot(profileResults$lim_or_dtable, aes(x=x, y=y/profileResults$ncol, fill=intensity, color=intensity)) + standard_layer(profileResults$data[2,1], profileResults$data[1,1]) + standard_theme + standard_y + standard_fill + standard_color + standard_labs + labs(title="7- contrast=\"norm\",\n range=c(0.02,0.98)")

#================================================================
#8- norm range=c(0.05,0.95)
#================================================================
profileResults<-cellProfiles(data=dtable,contrast="norm",range=c(0.05,0.95))

pt8<-ggplot(profileResults$lim_or_dtable, aes(x=x, y=y/profileResults$ncol, fill=intensity, color=intensity)) + standard_layer(profileResults$data[2,1], profileResults$data[1,1]) + standard_theme + standard_y + standard_fill + standard_color + standard_labs + labs(title="8- contrast=\"norm\",\n range=c(0.05,0.95)")

dev.new(width=12,height=6.4)
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
profileResults<-cellProfiles(data=dtable,contrast="max")

pt9<-ggplot(profileResults$lim_or_dtable, aes(x=x, y=y/profileResults$ncol, fill=intensity, color=intensity)) + standard_layer(profileResults$data[2,1], profileResults$data[1,1]) + standard_theme + standard_y + standard_fill + standard_color + standard_labs + labs(title="9- align=\"native\",\n contrast=\"max\"")

#================================================================
#10- orient
#================================================================
profileResults<-cellProfiles(data=dtable,align="orient",contrast="max")

pt10<-ggplot(profileResults$lim_or_dtable, aes(x=x, y=y/profileResults$ncol, fill=intensity, color=intensity)) + standard_layer(profileResults$data[2,1], profileResults$data[1,1]) + standard_theme + standard_y + standard_fill + standard_color + standard_labs + labs(title="10- align=\"orient\",\n contrast=\"max\"")

#================================================================
#11- random
#================================================================
profileResults<-cellProfiles(data=dtable,align="random",contrast="max")

pt11<-ggplot(profileResults$lim_or_dtable, aes(x=x, y=y/profileResults$ncol, fill=intensity, color=intensity)) + standard_layer(profileResults$data[2,1], profileResults$data[1,1]) + standard_theme + standard_y + standard_fill + standard_color + standard_labs+labs(title="11- align=\"random\",\n contrast=\"max\"")

#================================================================
#12- orient / reverse
#================================================================
profileResults<-cellProfiles(data=dtable,align="orient",contrast="max",reverse=TRUE)

pt12<-ggplot(profileResults$lim_or_dtable, aes(x=x, y=y/profileResults$ncol, fill=intensity, color=intensity)) + standard_layer(profileResults$data[2,1], profileResults$data[1,1]) + standard_theme + standard_y + standard_fill + standard_color + standard_labs + labs(title="12- align=\"orient\", reverse=TRUE, \ncontrast=\"max\"")

#================================================================
#13- reuse / reverse #1
#================================================================
profileResults<-cellProfiles(data=dtable,align="reuse",contrast="max",reverse=TRUE)

pt13<-ggplot(profileResults$lim_or_dtable, aes(x=x, y=y/profileResults$ncol, fill=intensity, color=intensity)) + standard_layer(profileResults$data[2,1], profileResults$data[1,1]) + standard_theme + standard_y + standard_fill + standard_color + standard_labs + labs(title="13- align=\"reuse\", reverse=TRUE, \ncontrast=\"max\"")

#================================================================
#14- reuse / reverse #2
#================================================================
profileResults<-cellProfiles(data=dtable,align="reuse",contrast="max",reverse=TRUE)

pt14<-ggplot(profileResults$lim_or_dtable, aes(x=x, y=y/profileResults$ncol, fill=intensity, color=intensity)) + standard_layer(profileResults$data[2,1], profileResults$data[1,1]) + standard_theme + standard_y + standard_fill + standard_color + standard_labs + labs(title="14- align=\"reuse\", reverse=TRUE, \ncontrast=\"max\"")

dev.new(width=12,height=6.4)
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
#15- line plots of individual profiles \n
#================================================================
profileResults<-cellProfiles(data=dtable,contrast="norm",align="orient")

pt15<-ggplot()+layer(data=profileResults$prop_dtable,mapping=aes(x=x,y=intensity,group=cell),geom="line",alpha=1/sqrt(profileResults$ncol))+theme_bw(base_size=10)+coord_cartesian(xlim=c(0,1))+labs(x="fraction of cell length",y="normalized intensity",title="15- line plots of individual profiles \n")

#================================================================
#16- line plots split between \nshort & long cells
#================================================================
profileResults<-cellProfiles(data=dtable,contrast="norm",align="orient")
profileResults$prop_dtable[5]<-profileResults$prop_dtable["y"]/profileResults$ncol>0.5
names(profileResults$prop_dtable)<- c("cell","x","y","intensity","half")

min<-colQuantiles(profileResults$or_dtable["intensity"],0.05,na.rm=TRUE)
max<-profileResults$max

pt16<-ggplot() + layer(data=profileResults$prop_dtable, mapping=aes(x=x, y=intensity, group=cell, color=factor(half)), geom="line", alpha=2/sqrt(profileResults$ncol)) + theme_bw(base_size=10) + coord_cartesian(xlim=c(0,1)) + scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") + theme(legend.position="none") + labs(x="fraction of cell length", y="normalized intensity", title="16- line plots split between short & long cells \n")

#================================================================
#17- line plots split between short\nmedium & long cells (simple mean & confidence)
#================================================================
profileResults<-cellProfiles(data=dtable,contrast="norm",align="orient")

#interpolate profile results to align the data points
nrow<-profileResults$nrow
ncol<-profileResults$ncol
interp_dtable<-matrix(data=NA,nrow=0,ncol=3)
for (i in 1:ncol){
	slice<-profileResults$prop_dtable[1:nrow+nrow*(i-1),]
	len<-length(profileResults$form_dtable[!is.na(profileResults$form_dtable[i]),i])
	interp<-approx(slice[c(2,4)],n=100)
	interp_dtable<-rbind(interp_dtable,cbind(interp$x,interp$y,rep(i,100)))
}
interp_dtable<-as.data.frame(interp_dtable)

#the next lines add a third colum that defines how to split the cells into categories
#to swap to the method used for graph 16, comment out the next line, and uncomment the one after it
interp_dtable[4]<-round(interp_dtable[3]/ncol*2)
#interp_dtable[4]<-interp_dtable[3]/ncol>0.5
names(interp_dtable)<-c("x","intensity","cell","half")

min<-profileResults$min
max<-profileResults$max

pt17<-ggplot() + layer(data=interp_dtable, mapping=aes(x=x, y=intensity, group=cell, color=factor(half)), geom="line", alpha=3/sqrt(profileResults$ncol)) + stat_summary(data=interp_dtable, mapping=aes(x=x, y=intensity, group=factor(half)), fun.data="mean_cl_boot", geom="smooth", alpha=0.75, fill="black") + stat_summary(data=interp_dtable, mapping=aes(x=x, y=intensity, color=factor(half)), fun.y="mean", geom="line") + coord_cartesian(xlim=c(0,1)) + scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") + theme_bw(base_size=10) + theme(legend.position="none") + labs(x="fraction of cell length", y="normalized intensity", title="17- line plots split between short\nmedium & long cells (simple mean & confidence)")

#note: mean_cl_boot (based on smean.cl.boot from the Hmisc package) is "a very fast implementation of the basic nonparametric bootstrap for obtaining confidence limits for the population mean without assuming normality"

#================================================================
#18- gaussian kernel density \nplot of cell lengths
#================================================================
lengths<-{}
lengths_temp<-raw_table[seq(1,ncol(raw_table),by=2)]
for(i in 1:ncol(lengths_temp)){lengths[i]<-max(lengths_temp[i],na.rm=TRUE)}
lengths<-as.data.frame(lengths)

pt18<-ggplot() + layer(data=lengths, mapping=aes(x=lengths), geom="density") + theme_bw(base_size=10) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + labs(x="cell length (um)", y="probability density", title="18- gaussian kernel density \nplot of cell lengths")

dev.new(width=8,height=6)
grid.newpage()
pushViewport(viewport(layout=grid.layout(3,2, heights = unit(c(0.25,3,3),"null"))))
grid.text("Additional Graph Types and Examples", vp=viewport(layout.pos.row=1,layout.pos.col=1:2),gp=gpar(fontsize=16))
print(pt15,vp=vplayout(2,1))
print(pt16,vp=vplayout(2,2))
cat("This one takes a moment to draw... be patient.")
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
profileResults1 <- cellProfiles(data=raw_table,contrast="norm", align="orient", range=c(0.06,0.98))
profileResults2 <- cellProfiles(data=raw_table,contrast="norm",align="orient", reverse=TRUE, range=c(0.06,0.98))

#rescale intensities to 0-1
profileResults1$lim_or_dtable$r <- rescale(profileResults1$lim_or_dtable$intensity)
profileResults1$lim_or_dtable$g <- rescale(profileResults2$lim_or_dtable$intensity)

#troubleshooting
# if your two data sets do not have exactly the same number of points, you will get issues and errors!
# MicrobeTracker occassionally produces such errors, so double check your data carefully if you use it.
# it may help to look at the final data sets, which can be exported in wide-format with the code below
# write.csv(as.matrix(profileResults1$lim_or_dtable$r), file="profile table 1 (red).csv")
# write.csv(as.matrix(profileResults2$lim_or_dtable$g), file="profile table 2 (green).csv")

#straight red/green probably works best, as straight blue too closely matches the black background
pt19 <- ggplot() + layer(data=profileResults1$lim_or_dtable, mapping=aes(x=x,y=y/profileResults1$ncol, fill=rgb(r,g,0), color=rgb(r,g,0)), geom="tile", width=1.01*(profileResults1$data[2,1]-profileResults1$data[1,1]), size=0.25) + scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.25), trans="reverse") + theme_bw() + scale_fill_identity() + scale_color_identity() + labs(x="distance from midcell (μm)", y="fraction of cell cycle", title="19- red/green colocalization")

#alternate light blue / orange colors for a more color-blind friendly palette
# in case you don't want to be evil
pt20 <- ggplot() + layer(data=profileResults1$lim_or_dtable, mapping=aes(x=x,y=y/profileResults1$ncol, fill=rgb(r,(0.5*g+0.5*r),g), color=rgb(r,(0.5*g+0.5*r),g)), geom="tile", width=1.01*(profileResults1$data[2,1]-profileResults1$data[1,1]), size=0.25) + scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.25), trans="reverse") + theme_bw() + scale_fill_identity() + scale_color_identity() + labs(x="distance from midcell (μm)", y="fraction of cell cycle", title="20- blue/orange colocalization")

dev.new(width=7,height=3.25)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,2, heights = unit(c(0.25,3),"null"))))
grid.text("Two-color colocalization graphs", vp=viewport(layout.pos.row=1,layout.pos.col=1:2),gp=gpar(fontsize=16))
print(pt19,vp=vplayout(2,1))
print(pt20,vp=vplayout(2,2))

#================================================================
# the following graphs are for documentation purposes only
#================================================================
# g<-g + theme_bw(base_size=10)
# g2<-g2 + theme_bw(base_size=10)

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

# pt19<-ggplot() + layer(data=profileResults$lim_or_dtable, mapping=aes(x=x,y=y/profileResults$ncol, fill=intensity, color=intensity), geom="tile", width=1.01*(profileResults$data[2,1]-profileResults$data[1,1]), size=0.25) + scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.25), trans="reverse") + scale_fill_gradientn(colours=rev(brewer.pal(9, "YlGnBu"))) + scale_color_gradientn(colours=rev(brewer.pal(9,"YlGnBu"))) + labs(x="distance from midcell (um)", y="fraction of cell cycle",title="dark background") + theme(panel.background=element_rect(fill="#061542"), panel.grid.minor=element_line(color="#225EA8"), panel.grid.major=element_line(color="#225ea8"))

# pt20<-ggplot() + layer(data=profileResults$lim_or_dtable, mapping=aes(x=x,y=y/profileResults$ncol, fill=intensity, color=intensity), geom="tile", width=1.01*(profileResults$data[2,1]-profileResults$data[1,1]), size=0.25) + scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.25), trans="reverse") + theme_bw() + scale_fill_gradientn(colours=rev(brewer.pal(9, "YlGnBu"))) + scale_color_gradientn(colours=rev(brewer.pal(9,"YlGnBu"))) + labs(x="distance from midcell (um)", y="fraction of cell cycle",title="light background")

# dev.new(width=7,height=3)
# grid.newpage()
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(pt19,vp=vplayout(1,1))
# print(pt20,vp=vplayout(1,2))