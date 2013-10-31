cellProfiles <- function(data=NULL,position="center",align="native",reverse=FALSE,contrast="native",range=c(0.02,0.98)){
#============================
#VERSION 2.4
#USER GUIDE
#
#Command line variables
# Change the behavior of this function by defining these variables when calling it
#
# 		eg: profileResults <- cellProfiles(data=dtable,align="orient",reverse=TRUE,contrast="normx")
#		The above command would use the data "dtable" to produce a graph with:
#			1) cells automatically aligned with the brightest sides [orient]
#			2) on the *left* [reverse]
#			3) while applying the normalize function (maximum ratio enforced via truncation) [normx]
#
# data-handling options
#  position = "center" or "left"
#		[DEFAULT]->	center:	 plot with midcells aligned at point zero
#					left:	plot with left side of cells aligned at zero
#
#  align 	= "native", "orient", "random", "reuse"
#		[DEFAULT]->	native:	profiles maintain original alignment
#					orient:	automatically aligns the cell profiles so that the brightest points are on the same side
#					random:	cells are randomly flipped
#					reuse:	apply the last graphed orientations to the current graph
#
#  reverse	= FALSE or TRUE
#		[DEFAULT]->	FALSE:	no change in orientation
#					TRUE:	plot profiles in reverse orientation (can be combined with, eg, align="reuse")
#					
# contrast adjustments for graphs	
#  contrast = "native", "norm", "max"
#		[DEFAULT]->	native:	profiles maintain original fluorescence intensities
#					norm:  	divides each profile by its mean value
#					max:	maximizes contrast 
#
#  range	= c(min,max)	This vector defines the lower and upper percentile boundaries of the color scale
#		[DEFAULT]->	c(0.02,0.98)		treat these limits like changing the min/max displayed constrast in ImageJ etc
#
#Exported variables
# these variables are accessible from the script results. 
#
# 		eg, if you run: profileResults <- cellProfiles(data=dtable)
# 		then the formatted final table is accessible at profileResults$form_dtable. 
# 		write the table to a new csv file using: 
#				write.csv(as.matrix(profileResults$form_dtable),file="cellprofiles wide format.csv")
#
# or_dtable		long data table for ggplot
# lim_or_dtable	long data table for ggplot, with values truncated as defined by the range setting
# data			input data
# form_dtable	wide-formatted table for human eyes
# prop_dtable	long data table, where cell length is converted to fraction of cell length
# ncol			number of cells measured
# nrow			length of original data table
# max			max fluorescence 
# min			min fluorescence
# med			median fluorescence
# mid			midpoint between min & max fluorescence
# maxlength		longest cell length
#
#END USER GUIDE
#============================
	t0<-proc.time()

	#set up a few variables
	if (local(exists("cellProfilesFunctionEnvironment.env"),envir=.GlobalEnv) == FALSE){local(cellProfilesFunctionEnvironment.env<-new.env(),envir=.GlobalEnv)}
	if (!(is.vector(range)) || (!(length(range)==2)) || !(range[1]>=0 && range[2]<=1) || !(range[1]<range[2])){
		cat("ERROR: Range must be between 0 and 1, in the format: c(min,max).");return(NULL)
	}
	if (position=="left"){adj_tick<-5} else {adj_tick<-6}
	cmin<-range[1]
	cmax<-range[2]
	ncol<-ncol(data)/2
	nrow<-nrow(data)
	dlength<-data[seq(1,ncol*2,by=2)]
	profile<-data[seq(2,ncol*2,by=2)]
	maxlength<-ceiling(max(dlength,na.rm=TRUE))
	wtmeanleft<-log((50:5)/5)/log(10)
	wtmeanright<-log((5:50)/5)/log(10)
	
	#sets up profile fliplist, or preserves prior orientations & disables inappropriate flags
	if (align == "reuse"){
		if ((local(exists("fliplist"), envir=cellProfilesFunctionEnvironment.env) != TRUE) || !(ncol == length(local(fliplist, envir=cellProfilesFunctionEnvironment.env)))){
			cat("ERROR: Use of the align=\"reuse\" option requires prior evaluation of a dataset with the same number of cells... otherwise there is nothing to \"reuse\".")
			return(NULL)
		}
		if (reverse == TRUE){
			local(fliplist<-(fliplist-1)/-1, envir=cellProfilesFunctionEnvironment.env)
		}
	} else {
		cellProfilesFunctionEnvironment.env$fliplist<-rep(1,ncol)
		if (align == "random"){
			cellProfilesFunctionEnvironment.env$fliplist[sample(ncol,ncol/2)]<-0
		}
	}
	
	
	#perform contrast adjustments and aligns or randomizes data as set in the options above
	cat("Calculating: \n ..")
	tick<-txtProgressBar(min=1,max=ncol*adj_tick,style=3)
	tmp_profile<-profile
	for(ii in 1:ncol){
		#for orient flag, find out if the top half of the column is brighter than the bottom half
		#if not, then set fliplist[ii] to 1 so the column will be reversed
		real_rows<-colCounts(as.matrix(profile[ii]+1),na.rm=TRUE)
		half_rows<-round(real_rows/2)
		if ((align == "orient") && mean(profile[1:half_rows,ii]*approx(wtmeanleft,n=half_rows)$y) < mean(profile[(real_rows-half_rows+1):real_rows,ii]*approx(wtmeanright,n=half_rows)$y)){
			cellProfilesFunctionEnvironment.env$fliplist[ii]<-0
		}
		if ((reverse == TRUE) && (align != "reuse")){
			cellProfilesFunctionEnvironment.env$fliplist[ii]<-(cellProfilesFunctionEnvironment.env$fliplist[ii]-1)/-1
		}
	
		#normalize profile data, as appropriate
		if (contrast == "norm"){
			col_mean<-mean(as.matrix(profile[ii]),na.rm=TRUE)
			profile[ii]<-(profile[ii]/col_mean)
		}
		
		#col_min & max are used to adjust contrast when max contrast is called		
		if (contrast != "max"){
			col_min<-0
			col_max<-1
		}	else {
			col_min<-min(profile[ii],na.rm=TRUE)	
			col_max<-max(profile[ii],na.rm=TRUE)-col_min
		}
		
		#save flipped/reversed/contrast adjusted values into tmp_profile, as appropriate 
		if (cellProfilesFunctionEnvironment.env$fliplist[ii]==1) {
			tmp_profile[c(real_rows:1),ii]<-(profile[c(1:real_rows),ii]-col_min)/(col_max)
		} else {
			tmp_profile[c(1:real_rows),ii]<-(profile[c(1:real_rows),ii]-col_min)/(col_max)
		}
		setTxtProgressBar(tick,ii)
	}
	profile<-tmp_profile	
	
	

    #find maximum length of each column, and the number of real values in each column
    #then order cells by measured cell length values
	cellength<-{}
	collength<-{}
	for(i in 1:ncol){
		cellength[i]<-max(dlength[i],na.rm=TRUE)
	}
	for(i in 1:ncol){
		collength[i]<-colCounts(as.matrix(profile[i]+1),na.rm=TRUE)
	}
	setTxtProgressBar(tick,2*ncol)
	collength<-order(cellength)
	
	#setup the y values to stack each profile
	plotheight<-{}
	for(i in 1:ncol){
		temp<-rep(i,nrow)
		plotheight<-append(plotheight,temp)
	}
	setTxtProgressBar(tick,3*ncol)
	
	
	
	#shifts cells to center (if center = TRUE)
	if (position == "center"){
		for(ii in 1:ncol){
			col_adj<-(max(dlength[ii],na.rm=TRUE)+min(dlength[ii],na.rm=TRUE))/2
			dlength[ii]<-dlength[ii]-col_adj
		}
	}
	setTxtProgressBar(tick,(adj_tick-2)*ncol)
	
	
	
	#apply order info from above
	or_dlength<-melt(dlength[collength],id=NULL)
	or_profile<-melt(profile[collength],id=NULL)
	or_dtable<-cbind(or_dlength,plotheight,or_profile[2])
	names(or_dtable)<- c("cell","x","y","intensity")
	
	
	
	max<-max(or_dtable["intensity"],na.rm=TRUE)
	min<-min(or_dtable["intensity"],na.rm=TRUE)
	med<-median(as.matrix(or_dtable["intensity"]),na.rm=TRUE)
	mid<-(max-min)/2+min
	
	
	
	#simple code to breakdown the long data list into something human-readable
	fodframe<-as.data.frame(or_profile[1:nrow,2])
	colnames(fodframe)<-or_profile[1,1]
	for(c in 2:ncol){
		odframe<-as.data.frame(or_profile[(1:nrow)+(c-1)*nrow,2])
		colnames(odframe)<-or_profile[(nrow*c),1]
		fodframe<-cbind(fodframe,odframe)
	}		
	setTxtProgressBar(tick,(adj_tick-1)*ncol)
	fodframe <- round(fodframe*1000)/1000
	
	
	
	#prepared a list suitable for making x-y scatter plots (x-axis converted to proportion of cell length)
	or_dlength<-dlength[collength]
	for(i in 1:ncol){
		or_dlength[i]<-rescale(or_dlength[i])
	}
	setTxtProgressBar(tick,(adj_tick)*ncol)
	close(tick)
	or_dlength<-melt(or_dlength,id=NULL)
	por_dtable<-cbind(or_dlength,plotheight,or_profile[2])
	names(por_dtable)<- c("cell","x","y","intensity")
	
	
	
	#eliminate empty rows
	or_dtable<-or_dtable[!is.na(or_dtable[4])==TRUE,c(1:4)]
	
	
	#create contrast-truncated version
	lim_or_dtable<-or_dtable	
	lim_min<-colQuantiles(or_dtable["intensity"],cmin,na.rm=TRUE)
	lim_max<-colQuantiles(or_dtable["intensity"],cmax,na.rm=TRUE)
	cmin_count<-length(lim_or_dtable$intensity[lim_or_dtable$intensity<lim_min])
	cmax_count<-length(lim_or_dtable$intensity[lim_or_dtable$intensity>lim_max])
	lim_or_dtable[lim_or_dtable["intensity"]<lim_min,"intensity"]<-lim_min
	lim_or_dtable[lim_or_dtable["intensity"]>lim_max,"intensity"]<-lim_max
	len_count<-length(lim_or_dtable$intensity)
	cat(len_count-cmin_count-cmax_count," profile points (",round(1000*(len_count-cmin_count-cmax_count)/len_count)/10,"%) fall in the range between the lower (",round(100*lim_min[[1]])/100,") and upper (",round(100*lim_max[[1]])/100,") contrast limits.\n",sep="")
	
	
	
	profileResults <- list("or_dtable"=or_dtable, "lim_or_dtable"=lim_or_dtable, "data"=data, "form_dtable"=fodframe, "prop_dtable"=por_dtable, "ncol"=ncol, "nrow"=nrow, "max"=max, "min"=min, "med"=med, "mid"=mid, "maxlength"=maxlength)
	
	tf<-proc.time()-t0
	cat("Calculations on",ncol,"cells complete in",tf[["elapsed"]],"seconds.")
	
	return(profileResults)
}

cellProfileTruncate <- function(data=NULL,adjust){
#============================
#USER GUIDE
#
# This function will remove [adjust] rows from the beginning and end of each cell
# Use this function to correct when, for eg, wide selection lines collect too much dark background to accurately represent cell poles
# Another solution to this problem would be to simply plot the graph background to a dark color
# eg, add to the ggplot2 code: +theme(panel.background = element_rect(fill="#061542"))
#
# adjust	= a positive integer. Will normally be 1-3 rows or so, depending on how you measured the data originally
#
#END USER GUIDE
#============================
	#various sanity checks before performing the truncation
	if (!(is.numeric(adjust) && floor(adjust)==adjust) || (adjust<0)){
		cat("ERROR: adjustment value must be a positive integer.");return(NULL)
	}
	if (adjust == 0){
		return(data)
	} else if (adjust >= nrow(na.omit(data))/2){
		cat("ERROR: cannot truncate more rows than exist. Must be",floor(nrow(na.omit(data))/2),"or fewer rows.");return(NULL)
	} else if (2*adjust/nrow(na.omit(data)) > 0.25){
		cat("WARNING: truncating ",round(100*2*adjust/nrow(na.omit(data))),"% of the shortest cell. This seems unwise.",sep="")
	}
	
	tmp_dtable<-matrix(NA,nrow(data),ncol(data))
	
	#it is pretty easy to lop off the top of the table
	data<-data[1+adjust:nrow(data),]
	
	#however, taking off the bottom rows from each column (of different lengths) is rather more complicated
	#below is an optimized form of a for loop that would copy, 2 columns at a time, all but the bottom [adjust] row(s)
	for(i in seq(1,ncol(data),by=2)){
		tmp_dtable[1:(nrow(as.matrix(data[!is.na(data[,i]),c(i)]))-adjust),c(i,i+1)]<-as.matrix(data[1:(nrow(as.matrix(data[!is.na(data[,i]),c(i)]))-adjust),c(i,i+1)])
	}
	return(as.data.frame(tmp_dtable))
}

library(ggplot2)  	
library(reshape2)
library(matrixStats)
library(scales)
library(grid)
library(RColorBrewer)