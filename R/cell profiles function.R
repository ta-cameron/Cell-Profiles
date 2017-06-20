#' Cell Profiles
#'
#' \code{cellProfiles} processes fluorescence intensity profiles from, eg, ImageJ, and prepares them for graphing as length-ordered demographs by ggplot
#' @param data wide-formated input data consisting of paired columns depicting (1) fluorescent intensity and (2) cell length. Accepts output directly from, eg, ImageJ
#' @param position character string (default is "left"), matching either "left" or "center". Determines alignment of stacked profiles
#' @param align character string (default is "native"), matching either "native", "orient", or "random", OR to reuse previous alignments, the results of a prior `cellProfiles()` run. Determines orientations of individual profiles.
#' @param reverse logical (default is FALSE). Globally reverses alignment of profiles.
#' @param contrast character string (default is "native"), matching either "native", "norm", or "max". Controls contrast of individual profiles; "norm" and "max" assist in visualizing profiles with heterogeneous intensity levels
#' @param range vector of length 2 (default is c(0.02,0.98) ), defining the lower and upper percentile boundaries of the color scale. Equivalent to adjusting the min/max displayed contrast in ImageJ
#'
#' @return a list of the following values:
#' \item{interpolated}{processed interpolated data using supplied range values. Columns consist of:
#' \itemize{x: fluorescence profile x-coordinates}
#' \itemize{intensity: fluorescence profile value}
#' \itemize{y: cell number / total cells. Ordered by height for demograph plots}}
#
#' \item{native}{as above, but using native resolution}
#' \item{proportional}{as above, with absolute cell length converted to fraction of cell length}
#' \item{wide}{wide-formatted table for human eyes; NOT for ggplot2}
#' \item{fliplist}{a simple vector that records orientations of each cell profile}
#' @seealso The cellProfiles vignette for detailed usage instructions:  \code{vignette("cellProfiles")}
#' @seealso \code{\link{cellProfileTruncate}} for trimming input
#' @seealso \code{\link[ggplot2]{ggplot}} for plotting output
#' @export
#' @import reshape2
#' @import scales
#' @examples
#' #calculate cell profiles
#' data_path <- system.file("extdata", "ftsZ_profiles.csv", package = "cellProfiles")
#' dtable <- read.table(data_path, header=FALSE, sep=",")
#' profileResults <- cellProfiles(data=dtable)
#'
#' #calculate cell profiles, orienting brightest poles on the right & normalizing cell-cell contrast
#' profileResults <- cellProfiles(data=dtable, align="orient", reverse=TRUE, contrast="norm")
#'
#' #plot with ggplot
#' ggplot2::ggplot(data=profileResults$interpolated, mapping=aes(x=x, y=y, fill=intensity, color=intensity) ) + geom_tile() + theme_classic() + scale_y_continuous(expand=c(0.0,0.0), breaks=seq(0,1,0.25), trans="reverse") + scale_fill_gradient(na.value = NA) + scale_color_gradient(na.value=NA) + labs(x="distance from midcell (um)", y="fraction of cell cycle")

cellProfiles <- function(data=NULL, position="center", align="native", reverse=FALSE, contrast="native", range=c(0.02,0.98)){
#============================
#DEBUG tools
#
# data <- dtable
# position <- "center"
# position <- "left"
# align <- "native"
# reverse <- FALSE
# contrast <- "norm"
# contrast <- "max"
# range <- c(0.02,0.98)
#============================

  # 0: initialization ----------------------------------------------------------

  ##function argument validity check
  #data frame check
  cellProfileDataFail(data=data)

  #check that position is "left" or "center"
  tryCatch(position <- match.arg(position, c("left", "center")),
           error = function(...)
             stop("position may only be one of: \"left\", \"center\".", call. = FALSE))

  #check that align is a either "native", "orient", "random" or the output of a previous cellProfiles run
  tryCatch({
    if (typeof(align) == "list") {
      fliplist <- get("fliplist", align)
      align <- "reuse"
    } else {
      align <- match.arg(align, c("native", "orient", "random"))
    }
  }, error = function(...) stop( "align must be \"native\", \"orient\", \"random\" or the output of a previous cellProfiles run.", call. = FALSE)
  )

  #check that reverse is a single logical value
  if( !(is.logical(reverse) && (length(reverse) == 1) ) ){
    stop("reverse may only be set to TRUE (T) or FALSE (F).", call. = FALSE)
  }

  #check that contrast is "native", "norm", or "max"
  tryCatch(contrast <- match.arg(contrast, c("native", "norm", "max")),
           error = function(...)
             stop("contrast may only be one of: \"native\", \"norm\", \"max\".", call. = FALSE))

  #check that range is a vector of length 2, between 0 & 1
	if (!(is.vector(range)) || (!(length(range)==2)) || !(range[1]>=0 && range[2]<=1) || !(range[1]<range[2])){
	  stop("range must be between 0 and 1, in the format: c(min,max).", call. = FALSE)
	}

  ##remove any straggling rows of all NA values from the data
  data <- data[rowSums(is.na(data)) != ncol(data),]

  ##set initial values
  t0 <- proc.time()
  adj_tick <- ifelse(position=="left", 7, 8)
	cmin <- range[1]
	cmax <- range[2]
	ncol <- ncol(data)/2
	nrow <- nrow(data) #note: this is overwritten with a new value in step 4
	dlength <- data[seq(1,ncol*2,by=2)]
	profile <- data[seq(2,ncol*2,by=2)]
	maxlength <- ceiling(max(dlength,na.rm=TRUE))
	wtmeanleft <- log((50:5)/5)/log(10)
	wtmeanright <- log((5:50)/5)/log(10)
	interpol_unit <- 0.02 #this sets the resolution of the interpolation.
	                      #smaller numbers drastically reduce rendering performance & may break things

	#sets up profile fliplist, or preserves prior orientations & disables inappropriate flags
	if (align == "reuse"){
		if ( !(ncol == length(fliplist)) ){
		  stop("Datasets must have matching numbers of cells in order to reuse alignment settings. (new=",ncol,"; old=", length(fliplist),")", call. = FALSE)
		}
		if (reverse == TRUE){
		  #alternates between 0 & 1 to reverse the orientation
			fliplist <- (fliplist-1)/-1
		}
	} else {
	  #default to 0 for no modification (native)
		fliplist <- rep(0,ncol)

		if (align == "random"){
		  #randomly choose half of cells to flip the orientation
			fliplist[sample(ncol,ncol/2)] <- 1
		}
	}

	#output interpretation of options back to user
	cat("Running cellProfiles with these options:
	    position=",position,", align=",align,", reverse=",reverse,", contrast=",contrast,", range=c(",range[1],", ",range[2],")\n", sep="")

	#initialize the progress bar. everyone loves progress bars, right?
	cat("Calculating:\n")
	tick <- txtProgressBar(min=1,max=ncol*adj_tick,style=3)
	setTxtProgressBar(tick,ncol/ncol)

	# 1: orientation ----------------------------------------------------------
	#orients or randomizes data
	for(ii in 1:ncol){
	  #DEBUG: ii <- 1
	  real_rows_ii <- nrow(na.omit(profile[ii]))
	  half_rows_ii <- round(real_rows_ii/2)

		#for orient flag, find out if the first half of the column (cell) is brighter than the second half
		#if so, then set fliplist[ii] to 1 so the column will be reversed
		if ((align == "orient") &&
		    mean(profile[1:half_rows_ii,ii] * approx(wtmeanleft, n=half_rows_ii)$y) >
		    mean(profile[(real_rows_ii - half_rows_ii+1):real_rows_ii,ii] * approx(wtmeanright, n=half_rows_ii)$y)){
			fliplist[ii] <- 1
		}
		if ((reverse == TRUE) && (align != "reuse")){
		  #again, to swap between 0 & 1 for each cell
			fliplist[ii] <- (fliplist[ii]-1)/-1
		}
	}
	setTxtProgressBar(tick,1*ncol)

	# 2: contrast ----------------------------------------------------------
	#perform contrast adjustments
	tmp_profile <- profile
  for(ii in 1:ncol){
    #DEBUG: ii <- 1
    real_rows_ii <- nrow(na.omit(profile[ii]))
    half_rows_ii <- round(real_rows_ii/2)

	  #normalize profile data, as appropriate
		if (contrast == "norm"){
			col_mean <- mean(as.matrix(profile[ii]),na.rm=TRUE)
			profile[ii] <- (profile[ii]/col_mean)
		}

		#col_min & max are used to adjust contrast when max contrast is called
		if (contrast != "max"){
			col_min <- 0
			col_max <- 1
		}	else {
			col_min <- min(profile[ii],na.rm=TRUE)
			col_max <- max(profile[ii],na.rm=TRUE)-col_min
		}

		#save flipped/reversed/contrast adjusted values into tmp_profile, as appropriate
		if (fliplist[ii]==1) {
			tmp_profile[c(real_rows_ii:1),ii] <- (profile[c(1:real_rows_ii),ii]-col_min)/(col_max)
		} else {
			tmp_profile[c(1:real_rows_ii),ii] <- (profile[c(1:real_rows_ii),ii]-col_min)/(col_max)
		}
	}
	profile <- tmp_profile
	setTxtProgressBar(tick,2*ncol)

	# 3: ordering ----------------------------------------------------------
	#find maximum length of each column, and the number of real values in each column
	#then order cells by measured cell length values
	cellength <- {}
	collength <- {}
	for(i in 1:ncol){
		cellength[i] <- max(dlength[i],na.rm=TRUE)
	}
	for(i in 1:ncol){
		collength[i] <- nrow(na.omit(profile[i]))
	}
	setTxtProgressBar(tick,3*ncol)
	collength <- order(cellength)

	# 4: centering & padding ----------------------------------------------------------
	#shifts cells to center (if center = TRUE)
	# & pads ~ 5% of max length on either side of profiles (if centered)
	# or just on the right if not

	#takes the lengths from the (first) longest cell
	x_grid_o <- dlength[,!is.na(dlength[nrow,])][,1]
	x_grid <- seq(min(x_grid_o), max(x_grid_o), length.out=length(x_grid_o))
	#must round after diff due to floating point error
	x_diff <- signif(diff(x_grid), digits=10)

	#finds most common difference b/t these x-units
	x_unit <- unique(x_diff)[ which.max(tabulate(match(x_diff, unique(x_diff) ) )) ]

  x_pad <- (max(x_grid)*0.05) %/% x_unit
  x_grid <- c(x_grid_o, seq(from= max(x_grid_o)+x_unit, by= x_unit, length.out= x_pad) )

  if (position == "center"){
    x_grid <- c( rev(seq(from=min(x_grid)-x_unit, by=-x_unit, length.out=x_pad)), x_grid)
    x_grid <- x_grid - 0.5*max(x_grid_o)
  }

  #lengths
  table_pad <- dlength[1:(length(x_grid)-nrow(profile)),]
  table_pad[1:nrow(table_pad),] <- NA
  dlength <- rbind(dlength, table_pad)

  #profile
  table_pad <- profile[1:(length(x_grid)-nrow(profile)),]
  table_pad[1:nrow(table_pad),] <- NA
  profile <- rbind(profile, table_pad )

  nrow <- nrow(dlength)

  if (position == "center"){
    for(ii in 1:ncol){
      #ii <- 24
      #ii <- 26
       if(nrow(na.omit(dlength[ii])) %% 2 == 1){
         offset <- 0
       } else {
         offset <- 1
       }

      nrow(na.omit(dlength[ii]))
      before <- ceiling( (length(x_grid)-length(dlength[!is.na(dlength[ii]), ii]) ) / 2 + offset)
      after <- length(x_grid) - ceiling( (length(x_grid)-length(dlength[!is.na(dlength[ii]), ii]) ) / 2)

      #basically shifts the profile y-values around so they align with the x_grid.
      #cells that are opposite the odd/even-ness of the longest cell will be shifted left by 1/2 a unit.
      profile[ii] <- c(rep(NA,1,before-1), na.omit(profile[,ii]), rep(NA, after, nrow(profile)-after))
    }
  }
  #copies in the uniform x_grid for all cells
  dlength[1:nrow(dlength),] <- x_grid

  setTxtProgressBar(tick,(adj_tick-4)*ncol)

  # 5: stacking ----------------------------------------------------------
  #setup the y values to stack each profile
  plotheight <- {}
  for(i in 1:ncol){
    plotheight <- append(plotheight, rep(i,nrow))
  }
  setTxtProgressBar(tick,(adj_tick-3)*ncol)

	# 6: reshaping data ----------------------------------------------------------
	#apply order info from above
	or_dlength <- reshape2::melt(dlength[collength],id=NULL)
	or_profile <- reshape2::melt(profile[collength],id=NULL)
	or_dtable <- cbind(or_dlength,plotheight/ncol,or_profile[2])
	names(or_dtable) <- c("cell","x","y","intensity")
	#dropping cell identifier
	or_dtable <- or_dtable[,c("x","y","intensity")]

	#eliminate empty rows
	#or_dtable <- or_dtable[!is.na(or_dtable$intensity)==TRUE,]

	max <- max(or_dtable["intensity"],na.rm=TRUE)
	min <- min(or_dtable["intensity"],na.rm=TRUE)
	med <- median(as.matrix(or_dtable["intensity"]),na.rm=TRUE)
	mid <- (max-min)/2+min

	#simple code to breakdown the long data list into something human-readable
	fodframe <- as.data.frame(or_profile[1:nrow,2])
	colnames(fodframe) <- or_profile[1,1]
	for(c in 2:ncol){
		odframe <- as.data.frame(or_profile[(1:nrow)+(c-1)*nrow,2])
		colnames(odframe) <- or_profile[(nrow*c),1]
		fodframe <- cbind(fodframe,odframe)
	}

	#this little bit takes the original length data, indexed by cell lenth, and grabs the last (longest) column to use as the row names
	rownames(fodframe) <- round(round(x_grid, digits=3), digits=3)
	fodframe <- round(fodframe, digits=3)
	setTxtProgressBar(tick,(adj_tick-2)*ncol)

	# 7: interpolating ----------------------------------------------------------
	#prepare a data set interpolated to the same coordinates (for geom_raster)
	interp_dtable <- data.frame("x"=NA, "y"=NA, "cell"=rep(1:ncol, each=length(seq(min(x_grid), max(x_grid)+interpol_unit, interpol_unit))) )
	n <- 1
	for(i in collength){
	  interp_dtable[interp_dtable$cell==n,] <-	data.frame( approx(x = dlength[[i]], y = profile[[i]], xout=seq(min(x_grid), max(x_grid)+interpol_unit, interpol_unit), rule = 1) , "cell"=n)
	  n <- n+1
	}
	interp_dtable <- interp_dtable[,c("x","cell","y")]
	names(interp_dtable) <- c("x","y","intensity")
	setTxtProgressBar(tick,(adj_tick-1)*ncol)

	# 8: rescaling (for scatter plots) ----------------------------------------------------------
	#prepare a list suitable for making x-y scatter plots (x - axis converted to proportion of cell length)
	prop_dtable <-
	  data.frame("x" = NA,
	             "y" = NA,
	             "cell" = rep(1:ncol,each=51))
	for (i in unique(prop_dtable$cell)) {
	  prop_dtable[prop_dtable$cell == i, ] <-
	    data.frame(approx(
	      x = scales::rescale(interp_dtable[interp_dtable$y == i &
	                                          !is.na(interp_dtable$intensity), "x"],
	                          to = c(0, 1)),
	      y = interp_dtable[interp_dtable$y == i &
	                          !is.na(interp_dtable$intensity), "intensity"],
	      xout = (0:50) / 50,
	      rule = 1
	    ),
	    "cell" = i)
	}
	prop_dtable <- prop_dtable[, c("x","cell","y")]
	names(prop_dtable) <- c("x","y","intensity")

	#finally convert the cell ID to y-axis height
	prop_dtable$y <- prop_dtable$y/ncol
	interp_dtable$y <- interp_dtable$y/ncol

	#close out the progress bar
	setTxtProgressBar(tick,(adj_tick)*ncol)
	close(tick)


	#create contrast-truncated version
	lim_or_dtable <- or_dtable
	lim_min <- quantile(or_dtable["intensity"], cmin, na.rm=TRUE)[[1]]
	lim_max <- quantile(or_dtable["intensity"], cmax, na.rm=TRUE)[[1]]
	cmin_count <- length(lim_or_dtable$intensity[(lim_or_dtable$intensity < lim_min) %in% TRUE])
	cmax_count <- length(lim_or_dtable$intensity[(lim_or_dtable$intensity > lim_max) %in% TRUE])
	lim_or_dtable[(lim_or_dtable$intensity < lim_min) %in% TRUE,"intensity"] <- lim_min
	lim_or_dtable[(lim_or_dtable$intensity > lim_max) %in% TRUE,"intensity"] <- lim_max
	len_count <- length(lim_or_dtable$intensity)
	cat(len_count-cmin_count-cmax_count," profile points (",round(100*(len_count-cmin_count-cmax_count)/len_count),"%) fall in the range between the lower (",round(lim_min[[1]], digits=2),") and upper (",round(lim_max[[1]], digits=2),") contrast limits.\n",sep="")


	#create contrast-truncated version for interpolated data
	lim_interp_dtable <- interp_dtable
	lim_interp_min <- quantile(interp_dtable["intensity"], cmin, na.rm=TRUE)[[1]]
	lim_interp_max <- quantile(interp_dtable["intensity"], cmax, na.rm=TRUE)[[1]]
	cmin_interp_count <- length(lim_interp_dtable$intensity[(lim_interp_dtable$intensity < lim_interp_min)  %in% TRUE])
	cmax_interp_count <- length(lim_interp_dtable$intensity[(lim_interp_dtable$intensity > lim_interp_max)  %in% TRUE])
	lim_interp_dtable[(lim_interp_dtable$intensity < lim_interp_min) %in% TRUE,"intensity"] <- lim_interp_min
	lim_interp_dtable[(lim_interp_dtable$intensity > lim_interp_max) %in% TRUE,"intensity"] <- lim_interp_max
	len_interp_count <- length(na.omit(lim_interp_dtable$intensity))
	cat(len_interp_count-cmin_interp_count-cmax_interp_count," profile points (",round(100*(len_interp_count-cmin_interp_count-cmax_interp_count)/len_interp_count),"%) fall in the range between the lower (",round(lim_interp_min[[1]], digits=2),") and upper (",round(lim_interp_max[[1]], digits=2),") contrast limits for interpolated data.\n",sep="")

	#final output data list
	profileResults <- list("native"=lim_or_dtable, "wide"=fodframe, "proportional"=prop_dtable, "interpolated"=lim_interp_dtable, "fliplist"=fliplist)
	
  #Depreciated output
  #	profileResults <- list("native"=lim_or_dtable, "wide"=fodframe, "proportional"=prop_dtable, "interpolated"=lim_interp_dtable, "ncells"=ncol, "nrow"=nrow, "max"=max, "min"=min, "med"=med, "mid"=mid, "interpol_unit"=interpol_unit, "fliplist"=fliplist)

	tf <- proc.time()-t0
	cat("Calculations on",ncol,"cells complete in",tf[["elapsed"]],"seconds.\n")
	#cat("DEBUG version\n")

	return(profileResults)
}


#' Cell Profile Truncate
#'
#' \code{cellProfileTruncate} takes raw fluorescence intensity profiles from, eg, ImageJ, and truncates the specified number of rows from the beginning and end of each cell
#' @param data data frame, consisting of paired columns depicting (1) fluorescent intensity and (2) cell length. Accepts output directly from, eg, ImageJ
#' @param adjust integer representing the number of data points to trim from each end of the profiles
#' @return a data frame
#' @seealso \code{\link{cellProfiles}}
#' @export
#' @examples
#' data_path <- system.file("extdata", "ftsZ_profiles.csv", package = "cellProfiles")
#' dtable <- read.table(data_path, header=FALSE, sep=",")
#' dtable_trimmed <- cellProfileTruncate(data=dtable, 1)
cellProfileTruncate <- function(data=NULL,adjust=0){
	##various sanity checks before performing the truncation

  #data frame check
  cellProfileDataFail(data=data)

  #"adjust" sanity checks
  if (!(is.numeric(adjust) && floor(adjust)==adjust) || (adjust<0)){
		stop("adjustment value must be a positive integer.", call. = FALSE)
  }

  if (adjust == 0){
    warning("adjust=0 : returning unmodified input data frame.", call. = FALSE)
    return(data)
	} else if (adjust >= nrow(na.omit(data))/2){
		stop("cannot truncate more data points than exist. Must be",floor(nrow(na.omit(data))/2),"or fewer data points", call. = FALSE)
	} else if (2*adjust/nrow(na.omit(data)) > 0.2){
		warning("truncating ", round(100*2*adjust/nrow(na.omit(data))), "% of the shortest cell. This seems unwise.", call. = FALSE)
	}

	tmp_dtable <- matrix(NA,nrow(data),ncol(data))

	#it is pretty easy to lop off the top of the table
	data <- data[1+adjust:nrow(data),]

	#however, taking off the bottom rows from each column (of different lengths) is rather more complicated
	#below is an optimized form of a for loop that would copy, 2 columns at a time, all but the bottom [adjust] row(s)
	for(i in seq(1,ncol(data),by=2)){
		tmp_dtable[1:(nrow(as.matrix(data[!is.na(data[,i]),c(i)]))-adjust),c(i,i+1)] <-
		  as.matrix(data[1:(nrow(as.matrix(data[!is.na(data[,i]),c(i)]))-adjust),c(i,i+1)])
	}
	return(as.data.frame(tmp_dtable))
}

cellProfileDataFail <- function(data = NULL) {
  if (!(is.data.frame(data))) {
    stop("please supply input data as a data frame.", call. = FALSE)
  }

  #paired data check (basic check for even column number)
  x <- ncol(data) / 2
  if (abs(x - round(x)) > .Machine$double.eps ^ 0.5) {
    stop("input data does not appear to be properly paired... it has an odd number (",ncol(data),") of columns.", call. = FALSE
    )
  }
}

