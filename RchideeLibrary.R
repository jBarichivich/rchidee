# RCHIDEE function library is a collection of R functions to easily extract, plot and postprocess 
#   variables from stomate and sechiba history files of ORCHIDEE-CAN. Currently, it works only 
#   for a SINGLE PIXEL and using a nonleap (365-day) calendar for daily output.
#
# Version 1.3 -- Added the following functions:
#   * RCHIDEEsto4dim.catSiteMonthlyVars        # tree-ring width 
#   * RCHIDEEsto4dim.catSiteAnnualNCIRVars     # tree-ring width; replaced deprecated function RCHIDEEsto.catSiteAnnualNCIRVars 
#
# Version 1.2 -- Added the following functions:
#   * RCHIDEE.extractERA5halfdForcing4Point
#   * RCHIDEE.getERA5halfdLandPoint
#
# Functions available:
#   GENERAL
#   * RCHIDEE.listVariables          -- Print to the console the name and dimensions of all variables in a history file.
#   * RCHIDEE.makeWeibull3parDist    -- Write to the console a three-parameter Weibull CIRC_CLASS_DIST for ORCHIDEE-CAN run.def.
#   * RCHIDEE.day2decyy
#   * RCHIDEE.month2decyy
#   * RCHIDEE.day2decyy.cal366d
#   * RCHIDEE.extractCRUNCEPhalfdForcing4Point
#   * RCHIDEE.getCRUNCEPhalfdLandPoint
#   * RCHIDEE.extractERA5halfdForcing4Point
#   * RCHIDEE.getERA5halfdLandPoint
#
#   STOMATE
#   * RCHIDEEsto.catSiteDailyCCBA    -- CONCATENATE daily basal area, basal area increment and number of trees for each circ class.
#   * RCHIDEEsto.readSiteDailyCCBA   -- READ daily basal area, basal area increment and number of trees for each circ class from an already concatenated file.
#   * RCHIDEEsto.retrieveTRW         -- Compute annual tree-ring width time series for each circ class from annual cincumference increments.
#   * RCHIDEEsto.catSiteDailyVars.cal365  -- CONCATENATE daily time series of a list of variables for a run with a fixed 365-day calendar 
#   * RCHIDEEsto.readSiteDailyVars.cal365 -- READ daily time series of a list of variables for a run with a fixed 365-day calendar 
#   * RCHIDEEsto.catSiteDailyVars.cal366  -- CONCATENATE daily time series of a list of variables for a run with a fixed 366-day calendar 
#   * RCHIDEEsto.catSiteMonthlyVars  -- CONCATENATE monthly time series of a list of variables
#   * RCHIDEEsto.catSiteAnnualVars   -- CONCATENATE annual time series of a list of variables
#   * RCHIDEEsto.catSiteAnnualNCIRVars -- CONCATENATE annual time series for variables across NCIR classes (otherwise it is tedious to do with RCHIDEEsto.catSiteAnnualVars for large NCIR)
#   * RCHIDEEsto4dim.catSiteMonthlyVars -- CONCATENATE monthly time series of a list of variables by ncric classee from the 4dim stomate history file
#   * RCHIDEEsto.catSiteAnnualVerticalCanopyVars -- CONCATENATE annual variables for each vertical canopy level from yearly stomate history files 
#
#   SECHIBA
#   * RCHIDEEsec.catSiteMonthlyVars  -- CONCATENATE monthly time series of a list of variables from a sechiba file
#
# To do:
#   * improve general input checking and error management
#   * function to concatenate any custom variable from the list
#   * function to plot a set of variables from a single experiment
#   * function to plot a comparison of several experiments for a set of variables
#
#
# Changes
#  * v1p1: RCHIDEE.listVariables: adding case when dimension of variable is 1
#  * v1p3: 8 dec 2024, added functions to read tree-ring width from new 4dim stomate output file in ORCHIDEE4
#
#
# @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
# @   February 2015, Paris



# Load packages needed by functions
packages <- c("ncdf4", "lubridate", "lattice","RColorBrewer")

# Loop through each package, check if installed, and if not, install it
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, type = "binary")
    library(pkg, character.only = TRUE)
  }
}

####################################################################################
RCHIDEE.day2decyy = function(syear,eyear){
####################################################################################

timeStr =  seq(as.Date( paste(syear, '-1-1',sep='') ),to=as.Date( paste(eyear,'-12-31', sep='') ), by='1 day') # need to start 1 day earlier to get right array length in decimal
isnot29Feb = as.logical( -grepl('-02-29', timeStr, fixed=T) )
timeDec = decimal_date( timeStr[!isnot29Feb] ) # from lubridate package
return(timeDec)
} # end function


####################################################################################
RCHIDEE.day2decyy.cal366d = function(syear,eyear){
####################################################################################

nyears=length(c(syear:eyear))
year366=1904 # any 
timeStr =  seq(as.Date( paste(year366, '-1-1',sep='') ),to=as.Date( paste(year366,'-12-31', sep='') ), by='1 day') # need to start 1 day earlier to get right array length in decimal
timeDec = rep( decimal_date( timeStr )-1904, times=nyears) + c(syear:eyear) # from lubridate package
return(timeDec)
} # end function


####################################################################################
RCHIDEE.day2decyy.calGregorian = function(syear,eyear){
####################################################################################

timeStr =  seq(as.Date( paste(syear, '-1-1',sep='') ),to=as.Date( paste(eyear,'-12-31', sep='') ), by='1 day') # need to start 1 day earlier to get right array length in decimal
timeDec = decimal_date( timeStr ) # from lubridate package
return(timeDec)
} # end function


####################################################################################
RCHIDEE.month2decyy = function(syear,eyear){
  ####################################################################################
  
  timeStr =  seq(as.Date( paste(syear, '-1-1',sep='') ),to=as.Date( paste(eyear,'-12-31', sep='') ), by='1 month') + 15 # monthly dates with mid-month point
  timeDec = decimal_date(timeStr) # from lubridate package
  return(timeDec)
} # end function


####################################################################################
RCHIDEE.listVariables = function(infile){
####################################################################################
  # RCHIDEE.listVariables. Function to print in the console the name and dimensions 
  # of all variables in a stomate or sechiba history file.
  # Call:
  #      RCHIDEE.listVariables(infile)
  #
  # Input:
  #       infile: Name of stomate or sechiba history file including path.
  #
  # Output:
  #       No output variable.
  #       Everything is printed in the screen.
  #
  # Example call:
  #       infile='/home/scratch01/jbari/IGCM_OUT/OL2/PROD/secsto/tstrec03/SBG/Output/run_19010101_19011231_1M_stomate_history.nc'
  #       RCHIDEE.listVariables(infile)
  #
  # Dependencies:
  #      require(ncdf)
  #
  # @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
  # @   14 February 2015, Paris
  #  
  
  # open file with an identifier
  ncin = nc_open(infile)
  # print list variables in a stomate or sechiba history file []
  for (vv in 1:ncin$nvars){
    if (vv==1){cat(paste('LISTING THE ', ncin$nvars ,' VARIABLES AVAILABLE IN HISTORY FILE', sep='')); cat("\n")}

    ndimensions=ncin$var[[vv]]$ndims
    # if there are only two dimensions [only few cases]
    if (ndimensions == 2){ 
      cat(paste( vv, ncin$var[[vv]]$name,  '(', ncin$var[[vv]]$varsize[1], ncin$var[[vv]]$varsize[2], '/', ncin$var[[vv]]$dim[[1]]$name, ncin$var[[vv]]$dim[[2]]$name ,')' ,'--', ncin$var[[vv]]$longname, '[',ncin$var[[vv]]$units, ']', sep=' ')); cat("\n")

    # if there are only one dimensions [only few cases] -- modified for Curie XIOS
    }else if(ndimensions == 1){
      cat(paste( vv, ncin$var[[vv]]$name,  '(', ncin$var[[vv]]$varsize[1], ncin$var[[vv]]$varsize[1], '/', ncin$var[[vv]]$dim[[1]]$name, ncin$var[[vv]]$dim[[1]]$name ,')' ,'--', ncin$var[[vv]]$longname, '[',ncin$var[[vv]]$units, ']', sep=' ')); cat("\n")
    # otherwise assume three dimensions [most cases]
    }else{
      cat(paste( vv, ncin$var[[vv]]$name,  '(', ncin$var[[vv]]$varsize[1], ncin$var[[vv]]$varsize[2], ncin$var[[vv]]$varsize[3],'/', ncin$var[[vv]]$dim[[1]]$name, ncin$var[[vv]]$dim[[2]]$name, ncin$var[[vv]]$dim[[3]]$name ,')' ,'--', ncin$var[[vv]]$longname, '[',ncin$var[[vv]]$units, ']', sep=' ')); cat("\n")
    }
  }
  nc_close(ncin) # close file
  
}  # end function




####################################################################################
RCHIDEEsto.catSiteDailyCCBA = function(options){
####################################################################################  
  # RCHIDEEsto.catSiteDailyCCBA. Concatenate daily basal area, basal area increment and number of trees 
  # for each circ class from yearly stomate or sechiba history files of SITE/PIXEL RUNS. 
  # These variables are needed for computing annual tree-ring width.
  # 
  # Call:
  #      RCHIDEEsto.catSiteDailyCCBA(options)
  #
  # Input:
  #       options: a list object with the names of the input files and run settings.
  #         options$path2file    : Path to history files
  #         options$expFilename  : Experiment name prefix for history files, e.g., 'run' (see example below)
  #         options$years        : Sequence of years for history files, e.g., c(1901:1955)
  #         options$ncirc        : Number of circunference classes used in the run (look at your run.def or ncdump your netcdf history file)
  #         options$pftId        : Index for the PFT to extract (look at your run.def used)
  #         options$idxLat       : Index for latitude of the gridbox to extract
  #         options$idxLon       : Index for longitude of the gridbox to extract
  #
  # Output: The output is a list of the three following variables.
  #         ccba      : Variable CCBA_NCIRC (ndays, ncirc)      -- "Basal area of a trees in a circ class"
  #         ccdeltaba : Variable CCDELTABA_NCIRC (ndays, ncirc) -- "Change in basal area of a trees in a circ class"
  #         ntrees    : Variable CCN_NCIRC (ndays, ncirc)       -- "Number of trees in a circ class"
  #
  # Example call:
  #      options=list() 
  #      options$path2file    = '/home/scratch01/jbari/IGCM_OUT/OL2/PROD/secsto/tstrec03/SBG/Output/'
  #      options$expFilename  = 'run';  
  #      options$years        = c(1901:1955)
  #      options$ncirc        = 20
  #      options$pftId        = 3
  #      options$idxLat       = 1
  #      options$idxLon       = 1  
  #      BA = RCHIDEEsto.catSiteDailyCCBA(options)
  #
  # Dependencies:
  #      require(ncdf)
  #
  # @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
  # @   14 February 2015, Paris
  #  
  
  # What is the unit for CCBA and CCDELTABA, m2?
  
# set options 
nyears = length(options$years)
ndays  = nyears * 365
ncirc  = options$ncirc
pftId  = options$pftId
idxLat = options$idxLat
idxLon = options$idxLon

# preallocate output arrays [ndays x 1]
ba = matrix(nrow=ndays,ncol=ncirc)
dba = matrix(nrow=ndays,ncol=ncirc)
ntrees = matrix(nrow=ndays,ncol=ncirc)

# Loop over yearly files and read daily CCBA, CCDELTABA and CCN for each NCIRC
counter_days=0
for (yy in 1:nyears){
  infile_yy = paste(options$path2file, options$expFilename, '_', options$years[yy], '0101_', options$years[yy], '1231_1M_stomate_history.nc', sep='');
  ncin      = nc_open(infile_yy)
  
  # check if file has 365 days (no leap calendar) and if not throw error
  if( length(ncvar_get(ncin,'time_counter')) != 365 ) stop( paste('ERROR: ALL YEARLY HISTORY FILES FOT THIS FUNCTION MUST HAVE 365 days (noleap years).', 'File:', infile_yy, 'has', length(ncvar_get(ncin,'time_counter')), 'days', sep=' '))
  
  # loop over circ classes
  for (ii in 1:ncirc){
    if (ii<10){
      ba[(counter_days+1):(counter_days+365),ii] = ncvar_get(ncin, paste('CCBA_00',ii, sep=''), collapse_degen=F)[idxLon,idxLat,pftId,]               # "Basal area of a trees in a circ class"
      dba[(counter_days+1):(counter_days+365),ii] = ncvar_get(ncin, paste('CCDELTABA_00',ii, sep=''),collapse_degen=F)[idxLon,idxLat,pftId,]         # "Change in basal area of a trees in a circ class"
      ntrees[(counter_days+1):(counter_days+365),ii] = ncvar_get(ncin, paste('CCN_00',ii, sep=''), collapse_degen=F)[idxLon,idxLat,pftId,]            # "Number of trees in a circ class"
    }else{
      ba[(counter_days+1):(counter_days+365),ii] = ncvar_get(ncin, paste('CCBA_0',ii, sep=''), collapse_degen=F)[idxLon,idxLat,pftId,]                # "Basal area of a trees in a circ class"
      dba[(counter_days+1):(counter_days+365),ii] = ncvar_get(ncin, paste('CCDELTABA_0',ii, sep=''), collapse_degen=F)[idxLon,idxLat,pftId,]          # "Change in basal area of a trees in a circ class"
      ntrees[(counter_days+1):(counter_days+365),ii] = ncvar_get(ncin, paste('CCN_0',ii, sep=''), collapse_degen=F)[idxLon,idxLat,pftId,]             # "Number of trees in a circ class"
    }
  }
  counter_days=  counter_days + 365 # update day counter at each yearly pass
  nc_close(ncin)
}

# set output
return(list(ccba=ba, ccdeltaba=dba, ccn=ntrees))

} # end function




####################################################################################
RCHIDEEsto.readSiteDailyCCBA = function(options){
####################################################################################  
  # RCHIDEEsto.readSiteDailyCCBA. READ daily basal area, basal area increment and number of trees 
  # for each circ class from an already CONCATENATED stomate or sechiba history file. 
  # These variables are needed for computing annual tree-ring width.
  # 
  # Call:
  #      RCHIDEEsto.readSiteDailyCCBA(options)
  #
  # Input:
  #       options: a list object with the names of the input files and run settings.
  #         options$path2file    : Path to history files
  #         options$years        : Sequence of years for history files, e.g., c(1901:1955)
  #         options$ncirc        : Number of circunference classes used in the run (look at your run.def or ncdump your netcdf history file)
  #         options$pftId        : Index for the PFT to extract (look at your run.def used)
  #         options$idxLat       : Index for latitude of the gridbox to extract
  #         options$idxLon       : Index for longitude of the gridbox to extract
  #
  # Output: The output is a list of the three following variables.
  #         ccba      : Variable CCBA_NCIRC (ndays, ncirc)      -- "Basal area of a trees in a circ class"
  #         ccdeltaba : Variable CCDELTABA_NCIRC (ndays, ncirc) -- "Change in basal area of a trees in a circ class"
  #         ntrees    : Variable CCN_NCIRC (ndays, ncirc)       -- "Number of trees in a circ class"
  #
  # Example call:
  #      options=list() 
  #      options$path2file    = '/home/users/jjeong/ferret/sigma+.thin.o.power5_all.nc'
  #      options$years        = c(1901:2000)
  #      options$ncirc        = 20
  #      options$pftId        = 14
  #      options$idxLat       = 1
  #      options$idxLon       = 1  
  #      BA = RCHIDEEsto.readSiteDailyCCBA(options)
  #
  # Dependencies:
  #      require(ncdf4)
  #
  # @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
  # @   17 February 2015, Paris
  #  
  
  # What is the unit for CCBA and CCDELTABA, m2?
  
# set options 
nyears = length(options$years)
ndays  = nyears * 365
ncirc  = options$ncirc
pftId  = options$pftId
idxLat = options$idxLat
idxLon = options$idxLon

# preallocate output arrays [ndays x 1]
ba = matrix(nrow=ndays,ncol=ncirc)
dba = matrix(nrow=ndays,ncol=ncirc)
ntrees = matrix(nrow=ndays,ncol=ncirc)

# Read daily CCBA, CCDELTABA and CCN for each NCIRC
infile = options$infile;
ncin   = nc_open(infile)
  
# check if file has 365 days (no leap calendar) and if not throw error
#if( length(ncvar_get(ncin,'time_counter')) != 365 ) stop( paste('ERROR: ALL YEARLY HISTORY FILES FOT THIS FUNCTION MUST HAVE 365 days (noleap years).', 'File:', infile_yy, 'has', length(ncvar_get(ncin,'time_counter')), 'days', sep=' '))
  
 # loop over circ classes
  for (ii in 1:ncirc){
    if (ii<10){
      ba[,ii] = ncvar_get(ncin, paste('CCBA_00',ii, sep=''), collapse_degen=F)[idxLon,idxLat,pftId,]               # "Basal area of a trees in a circ class"
      dba[,ii] = ncvar_get(ncin, paste('CCDELTABA_00',ii, sep=''),collapse_degen=F)[idxLon,idxLat,pftId,]         # "Change in basal area of a trees in a circ class"
      ntrees[,ii] = ncvar_get(ncin, paste('CCN_00',ii, sep=''), collapse_degen=F)[idxLon,idxLat,pftId,]            # "Number of trees in a circ class"
    }else{
      ba[,ii] = ncvar_get(ncin, paste('CCBA_0',ii, sep=''), collapse_degen=F)[idxLon,idxLat,pftId,]                # "Basal area of a trees in a circ class"
      dba[,ii] = ncvar_get(ncin, paste('CCDELTABA_0',ii, sep=''), collapse_degen=F)[idxLon,idxLat,pftId,]          # "Change in basal area of a trees in a circ class"
      ntrees[,ii] = ncvar_get(ncin, paste('CCN_0',ii, sep=''), collapse_degen=F)[idxLon,idxLat,pftId,]             # "Number of trees in a circ class"
    }
  }
nc_close(ncin)

# set output
return(list(ccba=ba, ccdeltaba=dba, ccn=ntrees))

} # end function




####################################################################################
RCHIDEEsto.retrieveTRW = function(ccba, ccdeltaba){
####################################################################################  
# RCHIDEEsto.retrieveTRW. Compute annual tree-ring width from annual (stem) circumference 
#    increments for each NCIRC. It uses the output of RCHIDEEsto.catSiteDailyCCBA or RCHIDEEsto.readSiteDailyCCBA functions. 
#     
# Call:
#      RCHIDEEsto.retrieveTRW(ccba, ccdeltaba)
#
# Input:
#         ccba      : Variable CCBA_NCIRC (ndays, ncirc)      -- "Basal area of a trees in a circ class"
#         ccdeltaba : Variable CCDELTABA_NCIRC (ndays, ncirc) -- "Change in basal area of a trees in a circ class"
#
# Output: 
#         trw_mm : Annual tree-ring width time series for each NCIRC in millimeters.
#
# Example call:
#      options=list() 
#      options$path2file    = '/home/scratch01/jbari/IGCM_OUT/OL2/PROD/secsto/tstrec03/SBG/Output/'
#      options$expFilename  = 'run';  
#      options$years        = c(1901:1955)
#      options$ncirc        = 20
#      options$pftId        = 3
#      options$idxLat       = 1
#      options$idxLon       = 1
#      BA=RCHIDEEsto.catSiteDailyCCBA(options)                    # First, extract daily basal areas
#      trw_ncirc = RCHIDEEsto.retrieveTRW(BA$ccba, BA$ccdeltaba)  # Second, compute tree-ring widths
#
# Dependencies:
#      none
#
# @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
# @   14 February 2015, Paris
#  
  
  # retrieve some settings
  ndays  = dim(BA$ccba)[1]
  nyears = ndays / 365
  ncirc  = dim(BA$ccba)[2]
  
  # check that daily ccba has years with exactly 365 days
  if( nyears%%1!=0 ) stop('ERROR: CCBA input must have 365-day years!!! Check that length(ccba)/365 is a whole number.')
  
  # preallocate output matrices
  trw_mm  = matrix(nrow=nyears, ncol=ncirc)     # tree-ring width in mm

  
  
  # LOOP OVER NCIRC
  for (ii in 1:ncirc){
    # Reshape arrays to [day x year]
    ccba.dailyArray.reshape = matrix( ccba[,ii], nrow=365 ) # day x year
    ccdeltaba.dailyArray.reshape = matrix( ccdeltaba[,ii], nrow=365 )
    
    # Compute annual radial increment from basal area increment
    ccba.day1      = matrix(ccba.dailyArray.reshape[1,] - ccdeltaba.dailyArray.reshape[1,])
    radious.day1   = (ccba.day1 / pi)^0.5    # area --> radious
    radious.day365 = matrix( (ccba.dailyArray.reshape[365,]/pi)^(0.5) ) # basal area of last day
    
    trw_mm[, ii] = 1000*(radious.day365 - radious.day1)   # change unit m --> mm
  } # end loop
  
    # define output
    #return(list(trw=trw_mm))
    return(trw_mm)
  
} # end function  



####################################################################################
RCHIDEE.makeWeibull3parNCIRC_DIST = function(alpha, beta, gamma, ncirc, nha, doPlot){
####################################################################################
# RCHIDEE.makeWeibull3parDist. Write to the console a three-parameter Weibull distribution 
#     CIRC_CLASS_DIST for a given NCIRC number and a custom set of parameters (location, shape and scale).
#
# Call:
#      RCHIDEE.makeWeibull3parDist.(xi, alpha, beta, gamma, doPlot)
#
# Input:
#         alpha   : Location parameter of Weibull distribution. Ser to <0.99 to avoid zero or negative values.
#         beta    : Scale parameter of Weibull distribution  
#         gamma   : Shape parameter of Weibull distribution 
#         ncirc   : Number of CIRC (circunference) size classes
#         nha     : Density of trees per hectare (just for plotting)
#         doPlot  : Binary flag to make or not a plot of the distribution
#
# Output: 
#         CIRC_CLASS_DIST : The normalised distribution values (adding up to 1) are written to the console ready to paste into the RUN.DEF file of ORCHIDEE-CAN.
#
# Example call:
#         alpha = 0.5; beta = 5.5; gamma = 8.0; ncirc=20; nha=10000; doPlot=1
#         RCHIDEE.makeWeibull3parNCIRC_DIST(alpha, beta, gamma, ncirc, nha, doPlot)
#
# Dependencies:
#      none
#
# @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
# @   14 February 2015, Paris
  

  # Define points at which to evaluate the function
  X=seq(from=1, to=ncirc, by=1)

  # Generate Weibull distribution with shape parameter gamma (c), scale parameter beta (b), and the location parameter alpha (m)
  # alpha = 0.6525598; beta = 5.4848149; gamma = 1.0;
  # alpha=1; beta=3; gamma=2
  #alpha = 0.99 # Fix alpha to 1 in order to avoid negative values!
  wb3 = gamma/beta * ( ( (X-alpha)/beta ) ^ (gamma - 1)  * exp( -(X-alpha)/beta )^gamma )
  wb3n = wb3/sum(wb3) # normalise PDF to add up to one

  # Loop over ncirc to write
  for (nc in 1:ncirc){
    if (nc < 10){ # if there are less than 10 NCIRC
      cat(paste( 'CIRC_CLASS_DIST__0000', nc, '=' , wb3n[nc], sep='')); cat("\n");
    }else if(nc < 100){        # 
      cat(paste( 'CIRC_CLASS_DIST__000', nc, '=', wb3n[nc], sep='') ); cat("\n");
  }else{        # 
    cat(paste( 'CIRC_CLASS_DIST__00', nc, '=', wb3n[nc], sep='') ); cat("\n");
    }  
  }
  
  # plot resulting size distribution for a given forest density (N/ha)
  if(doPlot==1){
    barplot(wb3n*nha, width=1, names.arg=c(1:ncirc), xlab='NCIRC', ylab='N/ha')
  }
    
} # end function


####################################################################################  
RCHIDEEsto.catSiteDailyVars.calGregorian = function(options){
####################################################################################  
  # RCHIDEEsto.catSiteDailyVars.calGregorian. Concatenate daily time series for a list of 
  #    variables from yearly stomate or sechiba history files with gregorian calendar. 
  #  
  # 
  # Call:
  #      RCHIDEEsto.catSiteDailyVars.calGregorian(options)
  #
  # Input:
  #       options: a list object with the names of the input files and run settings.
  #         options$path2file    : Path to history files
  #         options$expFilename  : Experiment name prefix for history files, e.g., 'run' (see example below)
  #         options$years        : Sequence of years for history files, e.g., c(1901:1955)
  #         options$pftId        : Index for the PFT to extract (look at your run.def used)
  #         options$idxLat       : Index for latitude of the gridbox to extract
  #         options$idxLon       : Index for longitude of the gridbox to extract
  #         options$variables    : List of variable names to process, e.g., c('NPP','GPP','TOTAL_M')
  #         options$experimentName: String or numeric ID for the current experiment
  #
  # Output: The output is a list of the variables specified in options$variables plus a column of TIME and EXPERIMENT .
  #         dimensions(nyears*365, 1) [rows, cols+2]
  #
  # Example call:
  #      options=list() 
  #      options$path2file    = '/home/scratch01/jbari/IGCM_OUT/OL2/PROD/secsto/tstrec03/SBG/Output/'
  #      options$expFilename  = 'run';  
  #      options$years        = c(1901:1955)
  #      options$pftId        = 3
  #      options$idxLat       = 1
  #      options$idxLon       = 1
  #      options$variables    = c('NPP','GPP','TOTAL_M')
  #      options$experimentName = 'EXP1'
  #      BA = RCHIDEEsto.catSiteDailyVars.calGregorian(options)
  #
  # Dependencies:
  #      require(ncdf)
  #      require(lubridate)
  #      RCHIDEE.month2decyy
  #
  # @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
  # @   14 February 2015, Paris
  #  
  

  # Make continuous time axis
  syear=min(options$years); eyear=max(options$years)
  TIME = RCHIDEE.day2decyy.calGregorian(syear, eyear)

  # set options 
  nyears = length(options$years)
  ndays  = length(TIME)
  pftId  = options$pftId
  idxLat = options$idxLat
  idxLon = options$idxLon
  varNames = options$variables
  experimentName = options$experimentName
  
  # preallocate output arrays [ndays x 1]
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', ndays,',','ncol=',ncirc, ')', sep='') ) )
    eval( parse(text= paste('var',ii,'=','matrix(nrow=', ndays,',','ncol=1', ')', sep='') ) )
  }
  
  # Loop over yearly files and concatenate daily data for each variable
  counter_days=0
  for (yy in 1:nyears){
    infile_yy = paste(options$path2file, options$expFilename, '_', options$years[yy], '0101_', options$years[yy], '1231_1M_stomate_history.nc', sep='');
    ncin      = nc_open(infile_yy)
    
    ndays = length( RCHIDEE.day2decyy.calGregorian(options$years[yy], options$years[yy]) )
    
    # loop over variables
    for (vv in 1:length(varNames)){
      eval( parse(text= paste('var',vv,'[(counter_days+1):(counter_days+ndays)] = ncvar_get(ncin, varNames[vv], collapse_degen=F)[idxLon,idxLat,pftId,]',sep='') ) )
    }# end variable loop
    
    counter_days=  counter_days + ndays
    nc_close(ncin)
  }# end year/file loop
  
  # Make time vector
  # TIME = c(1:counter_days)
   #syear=min(options$years); eyear=max(options$years)
  #TIME = RCHIDEE.day2decyy(syear, eyear)
   
  # Make experiment vector
  expName=rep(experimentName, counter_days)
  
  # prepare output
  out=as.character()
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', ndays,',','ncol=',ncirc, ')', sep='') ) )
    out[ii] = paste(paste(varNames[ii],'=var',ii, sep='',collapse = ",") )
  }
  strOut=paste(rbind(out), collapse=', ')
  return( eval(parse( text=paste('data.frame(', 'TIME=TIME', ',', 'EXPERIMENT=expName', ',', strOut ,')', sep='' ))) )
  
} # end function



####################################################################################  
RCHIDEEsto.catSiteDailyVars.cal365 = function(options){
####################################################################################  
  # RCHIDEEsto.catSiteDailyVars. Concatenate daily time series for a list of 
  #    variables from yearly stomate or sechiba history files with fixed years of 365 days. 
  #  
  # 
  # Call:
  #      RCHIDEEsto.catSiteDailyVars.cal365(options)
  #
  # Input:
  #       options: a list object with the names of the input files and run settings.
  #         options$path2file    : Path to history files
  #         options$expFilename  : Experiment name prefix for history files, e.g., 'run' (see example below)
  #         options$years        : Sequence of years for history files, e.g., c(1901:1955)
  #         options$pftId        : Index for the PFT to extract (look at your run.def used)
  #         options$idxLat       : Index for latitude of the gridbox to extract
  #         options$idxLon       : Index for longitude of the gridbox to extract
  #         options$variables    : List of variable names to process, e.g., c('NPP','GPP','TOTAL_M')
  #         options$experimentName: String or numeric ID for the current experiment
  #
  # Output: The output is a list of the variables specified in options$variables plus a column of TIME and EXPERIMENT .
  #         dimensions(nyears*365, 1) [rows, cols+2]
  #
  # Example call:
  #      options=list() 
  #      options$path2file    = '/home/scratch01/jbari/IGCM_OUT/OL2/PROD/secsto/tstrec03/SBG/Output/'
  #      options$expFilename  = 'run';  
  #      options$years        = c(1901:1955)
  #      options$pftId        = 3
  #      options$idxLat       = 1
  #      options$idxLon       = 1
  #      options$variables    = c('NPP','GPP','TOTAL_M')
  #      options$experimentName = 'EXP1'
  #      BA = RCHIDEEsto.catSiteDailyVars.cal365(options)
  #
  # Dependencies:
  #      require(ncdf)
  #      require(lubridate)
  #      RCHIDEE.month2decyy
  #
  # @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
  # @   14 February 2015, Paris
  #  
  
  # set options 
  nyears = length(options$years)
  ndays  = nyears * 365
  pftId  = options$pftId
  idxLat = options$idxLat
  idxLon = options$idxLon
  varNames = options$variables
  experimentName = options$experimentName
  
  # preallocate output arrays [ndays x 1]
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', ndays,',','ncol=',ncirc, ')', sep='') ) )
    eval( parse(text= paste('var',ii,'=','matrix(nrow=', ndays,',','ncol=1', ')', sep='') ) )
  }
  
  # Loop over yearly files and concatenate daily data for each variable
  counter_days=0
  for (yy in 1:nyears){
    infile_yy = paste(options$path2file, options$expFilename, '_', options$years[yy], '0101_', options$years[yy], '1231_1M_stomate_history.nc', sep='');
    ncin      = nc_open(infile_yy)
    
    # loop over variables
    for (vv in 1:length(varNames)){
      eval( parse(text= paste('var',vv,'[(counter_days+1):(counter_days+365)] = ncvar_get(ncin, varNames[vv], collapse_degen=F)[idxLon,idxLat,pftId,]',sep='') ) )
    }# end variable loop
    
    counter_days=  counter_days + 365
    nc_close(ncin)
  }# end year/file loop
  
  # Make time vector
  # TIME = c(1:counter_days)
  syear=min(options$years); eyear=max(options$years)
  TIME = RCHIDEE.day2decyy(syear, eyear)
   
  # Make experiment vector
  expName=rep(experimentName, counter_days)
  
  # prepare output
  out=as.character()
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', ndays,',','ncol=',ncirc, ')', sep='') ) )
    out[ii] = paste(paste(varNames[ii],'=var',ii, sep='',collapse = ",") )
  }
  strOut=paste(rbind(out), collapse=', ')
  return( eval(parse( text=paste('data.frame(', 'TIME=TIME', ',', 'EXPERIMENT=expName', ',', strOut ,')', sep='' ))) )
  
} # end function


####################################################################################  
RCHIDEEsto.readSiteDailyVars.cal365 = function(options){
####################################################################################  
  # RCHIDEEsto.catSiteDailyVars. Concatenate daily time series for a list of 
  #    variables from an ALREADY CONCATENATED stomate or sechiba history file with fixed 365-day calendar. 
  #  
  # 
  # Call:
  #      RCHIDEEsto.readSiteDailyVars.cal365(options)
  #
  # Input:
  #       options: a list object with the names of the input files and run settings.
  #         options$path2file    : Path to history files
  #         options$years        : Sequence of years for history files, e.g., c(1901:1955)
  #         options$pftId        : Index for the PFT to extract (look at your run.def used)
  #         options$idxLat       : Index for latitude of the gridbox to extract
  #         options$idxLon       : Index for longitude of the gridbox to extract
  #         options$variables    : List of variable names to process, e.g., c('NPP','GPP','TOTAL_M')
  #         options$experimentName: String or numeric ID for the current experiment
  #
  # Output: The output is a list of the variables specified in options$variables plus a column of TIME and EXPERIMENT .
  #         dimensions(nyears*365, 1) [rows, cols+2]
  #
  # Example call:
  #      options=list() 
  #      options$path2file    = '/home/users/jjeong/ferret/sigma+.thin.o.power5_all.nc'
  #      options$years        = c(1901:2000)
  #      options$pftId        = 14
  #      options$idxLat       = 1
  #      options$idxLon       = 1
  #      options$variables    = c('NPP','GPP')
  #      options$experimentName = 'EXP1'
  #      OUT = RCHIDEEsto.readSiteDailyVars.cal365(options)
  #
  # Dependencies:
  #      require(ncdf)
  #      require(lubridate)
  #      RCHIDEE.month2decyy
  #
  # @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
  # @   18 February 2015, Paris
  #  
  
  # set options 
  nyears = length(options$years)
  ndays  = nyears * 365
  pftId  = options$pftId
  idxLat = options$idxLat
  idxLon = options$idxLon
  varNames = options$variables
  experimentName = options$experimentName
  
  # preallocate output arrays [ndays x 1]
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', ndays,',','ncol=',ncirc, ')', sep='') ) )
    eval( parse(text= paste('var',ii,'=','matrix(nrow=', ndays,',','ncol=1', ')', sep='') ) )
  }
  
  # Loop over yearly files and concatenate daily data for each variable
    infile = options$path2file;
    ncin   = nc_open(infile)
    
    # loop over variables
    for (vv in 1:length(varNames)){
      eval( parse(text= paste('var',vv,' = ncvar_get(ncin, varNames[vv], collapse_degen=F)[idxLon,idxLat,pftId,]',sep='') ) )
    }# end variable loop
    
     nc_close(ncin)

  
  # Make time vector
  # TIME = c(1:counter_days)
  syear=min(options$years); eyear=max(options$years)
  TIME = RCHIDEE.day2decyy(syear, eyear)
   
  # Make experiment vector
  expName=rep(experimentName, ndays)
  
  # prepare output
  out=as.character()
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', ndays,',','ncol=',ncirc, ')', sep='') ) )
    out[ii] = paste(paste(varNames[ii],'=var',ii, sep='',collapse = ",") )
  }
  strOut=paste(rbind(out), collapse=', ')
  return( eval(parse( text=paste('data.frame(', 'TIME=TIME', ',', 'EXPERIMENT=expName', ',', strOut ,')', sep='' ))) )
  
} # end function




####################################################################################  
RCHIDEEsto.catSiteDailyVars.cal366 = function(options){
####################################################################################  
  # RCHIDEEsto.catSiteDailyVars.cal366. Concatenate daily time series for a list of 
  #    variables from yearly stomate or sechiba history files with constant years of 366 days. 
  #  
  # 
  # Call:
  #      RCHIDEEsto.catSiteDailyVars(options)
  #
  # Input:
  #       options: a list object with the names of the input files and run settings.
  #         options$path2file    : Path to history files
  #         options$expFilename  : Experiment name prefix for history files, e.g., 'run' (see example below)
  #         options$years        : Sequence of years for history files, e.g., c(1901:1955)
  #         options$pftId        : Index for the PFT to extract (look at your run.def used)
  #         options$idxLat       : Index for latitude of the gridbox to extract
  #         options$idxLon       : Index for longitude of the gridbox to extract
  #         options$variables    : List of variable names to process, e.g., c('NPP','GPP','TOTAL_M')
  #         options$experimentName: String or numeric ID for the current experiment
  #
  # Output: The output is a list of the variables specified in options$variables plus a column of TIME and EXPERIMENT .
  #         dimensions(nyears*365, 1) [rows, cols+2]
  #
  # Example call:
  #      options=list() 
  #      options$path2file    = '/home/scratch01/jbari/IGCM_OUT/OL2/PROD/secsto/tstrec03/SBG/Output/'
  #      options$expFilename  = 'run';  
  #      options$years        = c(1901:1955)
  #      options$pftId        = 3
  #      options$idxLat       = 1
  #      options$idxLon       = 1
  #      options$variables    = c('NPP','GPP','TOTAL_M')
  #      options$experimentName = 'EXP1'
  #      BA = RCHIDEEsto.catSiteDailyVars(options)
  #
  # Dependencies:
  #      require(ncdf)
  #      require(lubridate)
  #      RCHIDEE.month2decyy
  #
  # @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
  # @   14 February 2015, Paris
  #  
  
  # set options 
  nyears = length(options$years)
  ndays  = nyears * 366
  pftId  = options$pftId
  idxLat = options$idxLat
  idxLon = options$idxLon
  varNames = options$variables
  experimentName = options$experimentName
  
  # preallocate output arrays [ndays x 1]
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', ndays,',','ncol=',ncirc, ')', sep='') ) )
    eval( parse(text= paste('var',ii,'=','matrix(nrow=', ndays,',','ncol=1', ')', sep='') ) )
  }
  
  # Loop over yearly files and concatenate daily data for each variable
  counter_days=0
  for (yy in 1:nyears){
    infile_yy = paste(options$path2file, options$expFilename, '_', options$years[yy], '0101_', options$years[yy], '1231_1M_stomate_history.nc', sep='');
    ncin      = nc_open(infile_yy)
    
    # loop over variables
    for (vv in 1:length(varNames)){
      eval( parse(text= paste('var',vv,'[(counter_days+1):(counter_days+366)] = ncvar_get(ncin, varNames[vv], collapse_degen=F)[idxLon,idxLat,pftId,]',sep='') ) )
    }# end variable loop
    
    counter_days=  counter_days + 366
    nc_close(ncin)
  }# end year/file loop
  
  # Make time vector
  # TIME = c(1:counter_days)
  syear=min(options$years); eyear=max(options$years)
  #TIME = RCHIDEE.day2decyy(syear, eyear)
  TIME = RCHIDEE.day2decyy.cal366d(syear, eyear)
  
  # Make experiment vector
  expName=rep(experimentName, counter_days)
  
  # prepare output
  out=as.character()
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', ndays,',','ncol=',ncirc, ')', sep='') ) )
    out[ii] = paste(paste(varNames[ii],'=var',ii, sep='',collapse = ",") )
  }
  strOut=paste(rbind(out), collapse=', ')
  return( eval(parse( text=paste('data.frame(', 'TIME=TIME', ',', 'EXPERIMENT=expName', ',', strOut ,')', sep='' ))) )
  
} # end function





####################################################################################  
RCHIDEEsto.catSiteMonthlyVars = function(options){
  ####################################################################################  
  # RCHIDEEsto.catSiteMonthlyVars. Concatenate monthly time series for a list of 
  #    variables from yearly stomate or sechiba history files. 
  #  
  # 
  # Call:
  #      RCHIDEEsto.catSiteMonthlyVars(options)
  #
  # Input:
  #       options: a list object with the names of the input files and run settings.
  #         options$path2file    : Path to history files
  #         options$expFilename  : Experiment name prefix for history files, e.g., 'run' (see example below)
  #         options$years        : Sequence of years for history files, e.g., c(1901:1955)
  #         options$pftId        : Index for the PFT to extract (look at your run.def used)
  #         options$idxLat       : Index for latitude of the gridbox to extract
  #         options$idxLon       : Index for longitude of the gridbox to extract
  #         options$variables    : List of variable names to process, e.g., c('NPP','GPP','TOTAL_M')
  #         options$experimentName: String or numeric ID for the current experiment
  #
  # Output: The output is a list of the variables specified in options$variables.
  #         dimensions(nyears*12, 1) [rows, cols]
  #
  # Example call:
  #      options=list() 
  #      options$path2file    = '/home/scratch01/jbari/IGCM_OUT/OL2/PROD/secsto/tstrec03/SBG/Output/'
  #      options$expFilename  = 'run';  
  #      options$years        = c(1901:1955)
  #      options$pftId        = 3
  #      options$idxLat       = 1
  #      options$idxLon       = 1
  #      options$variables    = c('NPP','GPP','TOTAL_M')
  #      options$experimentName = 'EXP1'
  #      BA = RCHIDEEsto.catSiteMontlyVars(options)
  #
  # Dependencies:
  #      require(ncdf)
  #      RCHIDEE.month2decyy
  #
  # @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
  # @   14 February 2015, Paris
  #  
  
  # set options 
  nyears = length(options$years)
  nmonths  = nyears * 12
  pftId  = options$pftId
  idxLat = options$idxLat
  idxLon = options$idxLon
  varNames = options$variables
  experimentName = options$experimentName
  
  # preallocate output arrays [nmonths x 1]
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', nmonths,',','ncol=',ncirc, ')', sep='') ) )
    eval( parse(text= paste('var',ii,'=','matrix(nrow=', nmonths,',','ncol=1', ')', sep='') ) )
  }
  
  # Loop over yearly files and concatenate daily data for each variable
  counter_months=0
  for (yy in 1:nyears){
    infile_yy = paste(options$path2file, options$expFilename, '_', options$years[yy], '0101_', options$years[yy], '1231_1M_stomate_history.nc', sep='');
    ncin      = nc_open(infile_yy)
    
    # loop over variables
    for (vv in 1:length(varNames)){
      eval( parse(text= paste('var',vv,'[(counter_months+1):(counter_months+12)] = ncvar_get(ncin, varNames[vv], collapse_degen=F )[idxLon,idxLat,pftId,]',sep='') ) )
    }# end variable loop
    
    counter_months=  counter_months + 12
    nc_close(ncin)
  }# end year/file loop

  # Make time vector
  syear=min(options$years); eyear=max(options$years)
  TIME = RCHIDEE.month2decyy(syear, eyear)
  
  # Make experiment vector
  expName=rep(experimentName, counter_months)  
  
  # prepare output
  out=as.character()
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', nmonths,',','ncol=',ncirc, ')', sep='') ) )
    out[ii] = paste(paste(varNames[ii],'=var',ii, sep='',collapse = ",") )
  }
  strOut=paste(rbind(out), collapse=', ')
  #return( eval(parse( text=paste('list(', strOut ,')', sep='' ))) )
  return( eval(parse( text=paste('data.frame(', 'TIME=TIME', ',', 'EXPERIMENT=expName', ',', strOut ,')', sep='' ))) )
  
} # end function




####################################################################################  
RCHIDEEsto.catSiteAnnualVars = function(options){
####################################################################################  
  # RCHIDEEsto.catSiteAnnualVars Concatenate annual time series for a list of 
  #    variables from yearly stomate or sechiba history files. 
  #  
  # 
  # Call:
  #      RCHIDEEsto.catSiteAnnualVars(options)
  #
  # Input:
  #       options: a list object with the names of the input files and run settings.
  #         options$path2file    : Path to history files
  #         options$expFilename  : Experiment name prefix for history files, e.g., 'run' (see example below)
  #         options$years        : Sequence of years for history files, e.g., c(1901:1955)
  #         options$pftId        : Index for the PFT to extract (look at your run.def used)
  #         options$idxLat       : Index for latitude of the gridbox to extract
  #         options$idxLon       : Index for longitude of the gridbox to extract
  #         options$variables    : List of variable names to process, e.g., c('NPP','GPP','TOTAL_M')
  #         options$experimentName: String or numeric ID for the current experiment
  #
  # Output: The output is a list of the variables specified in options$variables.
  #         dimensions(nyears, 1) [rows, cols]
  #
  # Example call:
  #      options=list() 
  #      options$path2file    = '/home/scratch01/jbari/IGCM_OUT/OL2/PROD/secsto/tstrec03/SBG/Output/'
  #      options$expFilename  = 'run';  
  #      options$years        = c(1901:1955)
  #      options$pftId        = 3
  #      options$idxLat       = 1
  #      options$idxLon       = 1
  #      options$variables    = c('NPP','GPP','TOTAL_M')
  #      options$experimentName = 'EXP1'
  #      BA = RCHIDEEsto.catSiteAnnualVars(options)
  #      
  # Dependencies:
  #      require(ncdf)
  #
  # @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
  # @   14 February 2015, Paris
  #  
  
  # set options 
  nyears = length(options$years)
  nmonths  = nyears
  pftId  = options$pftId
  idxLat = options$idxLat
  idxLon = options$idxLon
  varNames = options$variables
  experimentName = options$experimentName
  
  # preallocate output arrays [nmonths x 1]
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', nmonths,',','ncol=',ncirc, ')', sep='') ) )
    eval( parse(text= paste('var',ii,'=','matrix(nrow=', nmonths,',','ncol=1', ')', sep='') ) )
  }
  
  # Loop over yearly files and concatenate daily data for each variable
  counter_years=0
  for (yy in 1:nyears){
    infile_yy = paste(options$path2file, options$expFilename, '_', options$years[yy], '0101_', options$years[yy], '1231_1Y_stomate_history.nc', sep='');
    ncin      = nc_open(infile_yy)
    
    # loop over variables
    for (vv in 1:length(varNames)){
      eval( parse(text= paste('var',vv,'[(counter_years+1):(counter_years+1)] = ncvar_get(ncin, varNames[vv], collapse_degen=F )[idxLon,idxLat,pftId,]',sep='') ) )
    }# end variable loop
    
    counter_years=  counter_years + 1
    nc_close(ncin)
  }# end year/file loop
  
  # Make time vector
  TIME = options$years
  
  # Make experiment vector
  expName=rep(experimentName, counter_years)  
    
  # prepare output
  out=as.character()
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', nyears,',','ncol=',ncirc, ')', sep='') ) )
    out[ii] = paste(paste(varNames[ii],'=var',ii, sep='',collapse = ",") )
  }
  strOut=paste(rbind(out), collapse=', ')
  #return( eval(parse( text=paste('list(', strOut ,')', sep='' ))) )
  return( eval(parse( text=paste('data.frame(', 'TIME=TIME', ',', 'EXPERIMENT=expName', ',', strOut ,')', sep='' ))) )
    
} # end function


####################################################################################  
RCHIDEEsto4dim.catSiteAnnualNCIRVars = function(options){
####################################################################################  
  # RCHIDEEsto4dim.catSiteAnnualNCIRVars Concatenate annual time series for a list of 
  #    variables across NCIR classes from yearly 4dim stomate or sechiba history files. 
  #  
  # 
  # Call:
  #      RCHIDEEsto4dim.catSiteAnnualNCIRVars(options)
  #
  # Input:
  #       options: a list object with the names of the input files and run settings.
  #         options$path2file    : Path to history files
  #         options$expFilename  : Experiment name prefix for history files, e.g., 'run' (see example below)
  #         options$years        : Sequence of years for history files, e.g., c(1901:1955)
  #         options$pftId        : Index for the PFT to extract (look at your run.def used)
  #         options$idxLat       : Index for latitude of the gridbox to extract
  #         options$idxLon       : Index for longitude of the gridbox to extract
  #         options$ncir         : Number of circ classes
  #         options$variables    : List of variable names to process, e.g., c('CCDELTABA','CCN','CCD','CCH')
  #         options$experimentName: String or numeric ID for the current experiment
  #
  # Output: The output is a list of the variables specified in options$variables.
  #         dimensions(nyears, 1) [rows, cols]
  #
  # Example call:
  #      options=list() 
  #      options$path2file    = '/home/scratch01/jbari/IGCM_OUT/OL2/PROD/secsto/tstrec03/SBG/Output/'
  #      options$expFilename  = 'run';  
  #      options$years        = c(1901:1955)
  #      options$pftId        = 3
  #      options$idxLat       = 1
  #      options$idxLon       = 1
  #      options$ncir         = 5
  #      options$variables    = c('CCDELTABA','CCN','CCD','CCH')
  #      options$experimentName = 'EXP1'
  #      TMP = RCHIDEEsto4dim.catSiteAnnualNCIRVars(options)
  #      
  # Dependencies:
  #      require(ncdf)
  #
  # @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
  # @   8 December 2024, Paris
  #  
  
  # set options 
  nyears = length(options$years)
  nmonths  = nyears
  pftId  = options$pftId
  idxLat = options$idxLat
  idxLon = options$idxLon
  ncir   = options$ncir
  varNames = options$variables
  experimentName = options$experimentName
  
  # preallocate output arrays [nmonths x 1]
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', nmonths,',','ncol=',ncirc, ')', sep='') ) )
    eval( parse(text= paste('var',ii,'=','data.frame(matrix(nrow=', nmonths,',','ncol=', ncir,'))', sep='') ) )
  }
  
  # Loop over yearly files and concatenate daily data for each variable
  counter_years=0
  for (yy in 1:nyears){
    infile_yy = paste(options$path2file, options$expFilename, '_', options$years[yy], '0101_', options$years[yy], '1231_1M_stomate_history_4dim.nc', sep='');
    ncin      = nc_open(infile_yy)
    
    # loop over variables and circ classes
    # print(ncin); float CCTRW[lon_domain_landpoints,lat_domain_landpoints,veget,ncirc,time_counter]   (Chunking: [1,1,15,3,1])  (Compression: shuffle,level 2)

    for (vv in 1:length(varNames)){
      for(cc in 1:ncir){
          # set variable name for a given CIRC class
          if(cc < 10){
          #varNames_tmp=paste(varNames[vv],'_00',cc, sep='')
          varNames_tmp=paste(varNames[vv], sep='')
          }else if(cc >=10){
          #varNames_tmp=paste(varNames[vv],'_0',cc, sep='')
          varNames_tmp=paste(varNames[vv], sep='')
          }
          # read variable into the preallocated data frame
          eval( parse(text= paste('var',vv,'[(counter_years+1):(counter_years+1),cc] = ncvar_get(ncin, varNames_tmp,collapse_degen=F )[idxLon,idxLat,pftId,cc,]',sep='') ) )
          # add name to the corresponding column
          if(cc < 10){
             eval( parse(text= paste('colnames(','var',vv,')[',cc,']=("',varNames[vv],'_00',cc,'")' ,sep='')))
          }else if(cc >=10){
             eval( parse(text= paste('colnames(','var',vv,')[',cc,']=("',varNames[vv],'_0',cc,'")' ,sep='')))
          }
      }
    }# end variable loop
    
    counter_years=  counter_years + 1
    nc_close(ncin)
  }# end year/file loop
  
  # Make time vector
  TIME = options$years
  
  # Make experiment vector
  expName=rep(experimentName, counter_years)  
    
  # prepare output
  out=as.character()
  for (ii in 1:length(varNames)){
    ##eval( parse(text= paste('var',ii,'=','matrix(nrow=', nyears,',','ncol=',ncirc, ')', sep='') ) )
    #out[ii] = paste(paste(varNames[ii],'=var',ii, sep='',collapse = ",") )
    out[ii] = paste('var',ii, sep='',collapse = ",")
  }
  strOut=paste(rbind(out), collapse=', ')

  #return( eval(parse( text=paste('list(', strOut ,')', sep='' ))) )
  return( eval(parse( text=paste('data.frame(', 'TIME=TIME', ',', 'EXPERIMENT=expName', ',', strOut ,')', sep='' ))) )
    
} # end function



####################################################################################  
RCHIDEEsto4dim.catSiteMonthlyVars = function(options){
  ####################################################################################  
  # RCHIDEEsto4dim.catSiteMonthlyVars. Concatenate monthly time series for a list of 
  #    variables from yearly 4dim stomate or sechiba history files. 
  #  
  # 
  # Call:
  #      RCHIDEEsto4dim.catSiteMonthlyVars(options)
  #
  # Input:
  #       options: a list object with the names of the input files and run settings.
  #         options$path2file    : Path to history files
  #         options$expFilename  : Experiment name prefix for history files, e.g., 'run' (see example below)
  #         options$years        : Sequence of years for history files, e.g., c(1901:1955)
  #         options$pftId        : Index for the PFT to extract (look at your run.def used)
  #         options$idxLat       : Index for latitude of the gridbox to extract
  #         options$idxLon       : Index for longitude of the gridbox to extract
  #         options$ncir         : Number of circ classes (look at your run.def used)
  #         options$variables    : List of variable names to process, e.g., c('CCTRW','GPP','TOTAL_M')
  #         options$experimentName: String or numeric ID for the current experiment
  #
  # Output: The output is a list of the variables specified in options$variables.
  #         dimensions(nyears*12, 1) [rows, cols]
  #
  # Example call:
  #      options=list() 
  #      options$path2file    = '/home/scratch01/jbari/IGCM_OUT/OL2/PROD/secsto/tstrec03/SBG/Output/'
  #      options$expFilename  = 'run';  
  #      options$years        = c(1901:1955)
  #      options$pftId        = 5
  #      options$idxLat       = 1
  #      options$idxLon       = 1
  #      options$idxLon       = 3
  #      options$variables    = c('CCTRW')
  #      options$experimentName = 'EXP1'
  #      BA = RCHIDEEsto4dim.catSiteMonthlyVars(options)
  #
  # Dependencies:
  #      require(ncdf)
  #      RCHIDEE.month2decyy
  #
  # @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
  # @   8 December 2024, Paris
  #  
  
  # set options 
  nyears = length(options$years)
  nmonths  = nyears * 12
  pftId  = options$pftId
  idxLat = options$idxLat
  idxLon = options$idxLon
  ncir   = options$ncir
  varNames = options$variables
  experimentName = options$experimentName
  
  # preallocate output arrays [nmonths x 1]
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', nmonths,',','ncol=1', ')', sep='') ) )
    eval( parse(text= paste('var',ii,'=','data.frame(matrix(nrow=', nmonths,',','ncol=', ncir,'))', sep='') ) )
  }
  
  # Loop over yearly files and concatenate daily data for each variable
  counter_months=0
  for (yy in 1:nyears){
    infile_yy = paste(options$path2file, options$expFilename, '_', options$years[yy], '0101_', options$years[yy], '1231_1M_stomate_history_4dim.nc', sep='');
    ncin      = nc_open(infile_yy)
    
    for (vv in 1:length(varNames)){
      for(cc in 1:ncir){
          # set variable name for a given CIRC class
          if(cc < 10){
          #varNames_tmp=paste(varNames[vv],'_00',cc, sep='')
          varNames_tmp=paste(varNames[vv], sep='')
          }else if(cc >=10){
          #varNames_tmp=paste(varNames[vv],'_0',cc, sep='')
          varNames_tmp=paste(varNames[vv], sep='')
          }
          # read variable into the preallocated data frame
          eval( parse(text= paste('var',vv,'[(counter_months+1):(counter_months+12),cc] = ncvar_get(ncin, varNames_tmp,collapse_degen=F )[idxLon,idxLat,pftId,cc,]',sep='') ) )
          # add name to the corresponding column
          if(cc < 10){
             #eval( parse(text= paste('colnames(','var',vv,')[',cc,']=("',varNames[vv],'_00',cc,'")' ,sep='')))
             eval( parse(text= paste('colnames(','var',vv,')[',cc,']=("','ncir','_00',cc,'")' ,sep='')))
          }else if(cc >=10){
             #eval( parse(text= paste('colnames(','var',vv,')[',cc,']=("',varNames[vv],'_0',cc,'")' ,sep='')))
             eval( parse(text= paste('colnames(','var',vv,')[',cc,']=("','ncir','_0',cc,'")' ,sep='')))
          }
      }
    }# end variable loop

    counter_months=  counter_months + 12
    nc_close(ncin)
  }# end year/file loop

  # Make time vector
  syear=min(options$years); eyear=max(options$years)
  TIME = RCHIDEE.month2decyy(syear, eyear)
  
  # Make experiment vector
  expName=rep(experimentName, counter_months)  
  
  # prepare output
  out=as.character()
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', nmonths,',','ncol=',ncirc, ')', sep='') ) )
    out[ii] = paste(paste(varNames[ii],'=var',ii, sep='',collapse = ",") )
  }
  strOut=paste(rbind(out), collapse=', ')
  #return( eval(parse( text=paste('list(', strOut ,')', sep='' ))) )
  return( eval(parse( text=paste('data.frame(', 'TIME=TIME', ',', 'EXPERIMENT=expName', ',', strOut ,')', sep='' ))) )
  
} # end function




####################################################################################  
RCHIDEEsto.catSiteAnnualVerticalCanopyVars = function(options){
####################################################################################  
  # RCHIDEEsto.catSiteAnnualVerticalCanopyVars Concatenate annual time series for a list of 
  #    variables describing canopy vertical profiles from yearly stomate history files. 
  #  
  # 
  # Call:
  #      RCHIDEEsto.catSiteAnnualVerticalCanopyVars(options)
  #
  # Input:
  #       options: a list object with the names of the input files and run settings.
  #         options$path2file    : Path to history files
  #         options$expFilename  : Experiment name prefix for history files, e.g., 'run' (see example below)
  #         options$years        : Sequence of years for history files, e.g., c(1901:1955)
  #         options$pftId        : Index for the PFT to extract (look at your run.def used)
  #         options$idxLat       : Index for latitude of the gridbox to extract
  #         options$idxLon       : Index for longitude of the gridbox to extract
  #         options$canopyLevels : Number of vertical canopy levels (usually 10)
  #         options$historyType  : 'sechiba' or 'stomate' history file
  #         options$variables    : List of variable names to process, e.g., c('LAI_PER_LEVEL','CLEVEL_HEIGHT','ABS_CANOPY_LIGHT')
  #         options$experimentName: String or numeric ID for the current experiment
  #
  # Output: The output is a list of the variables specified in options$variables. Each variable/column has the level ID (e.g.,
  #         LAI_PER_LEVEL_01, LAI_PER_LEVEL_02... LAI_PER_LEVEL_010).  
  #         dimensions(nyears, variables) [rows, cols]
  #         Other variables of light per level such as ABS_CANOPY_LIGHT, TRANS_CANOPY_LIGHT and CUM_CANOPY_LIGHT are written
  #         at daily time step. This is accounted for here but the fucntion need to be called separately for daily and annual variables.
  #
  # Example call:
  #      options=list() 
  #      options$path2file    = '/home/scratch01/jbari/IGCM_OUT/OL2/PROD/secsto/tstrec03/SBG/Output/'
  #      options$expFilename  = 'run';  
  #      options$years        = c(1901:1955)
  #      options$pftId        = 2
  #      options$idxLat       = 1
  #      options$idxLon       = 1
  #      options$canopyLevels = 10
  #      options$historyType  = 'sechiba'
  #      options$variables    = c('LAI_PER_LEVEL','CLEVEL_HEIGHT')
  #      options$experimentName = 'EXP1'
  #      X = RCHIDEEsto.catSiteAnnualVerticalCanopyVars(options)
  #      
  # Dependencies:
  #      require(ncdf)
  #
  # @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
  # @   29 March 2015, Paris
  #  
  
  # set options 
  nyears = length(options$years)
  nmonths  = nyears
  pftId  = options$pftId
  idxLat = options$idxLat
  idxLon = options$idxLon
  canopyLevels = options$canopyLevels
  varNames = options$variables
  historyType = options$historyType
  experimentName = options$experimentName
  
  # preallocate output arrays [nmonths x 1]
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', nmonths,',','ncol=',ncirc, ')', sep='') ) )
    if(varNames[ii] %in% c('ABS_CANOPY_LIGHT', 'TRA_CANOPY_LIGHT', 'CUM_CANOPY_LIGHT')){
      eval( parse(text= paste('var',ii,'=','data.frame(matrix(nrow=', nmonths*365,',','ncol=', canopyLevels, '))', sep='') ) )
    }else{
      eval( parse(text= paste('var',ii,'=','data.frame(matrix(nrow=', nmonths,',','ncol=', canopyLevels, '))', sep='') ) )
    }
  }
  
  # Loop over yearly files and concatenate daily data for each variable
  counter_years=0
  counter_days=0
  for (yy in 1:nyears){
    infile_yy = paste(options$path2file, options$expFilename, '_', options$years[yy], '0101_', options$years[yy], '1231_1M_', historyType,'_history.nc', sep='');
    ncin      = nc_open(infile_yy)
    
    # loop over variables and canopy levels
    for (vv in 1:length(varNames)){
       for(ll in 1:canopyLevels){
          varNames_tmp=paste(varNames[vv],'_00',pftId, sep='')
          if(varNames[vv] %in% c('ABS_CANOPY_LIGHT', 'TRA_CANOPY_LIGHT', 'CUM_CANOPY_LIGHT')){
            eval( parse(text= paste('var',vv,'[(counter_days+1):(counter_days+365),ll] = ncvar_get(ncin, varNames_tmp,collapse_degen=T )[idxLon,idxLat,ll,]',sep='') ) )
          }else{
            eval( parse(text= paste('var',vv,'[(counter_years+1):(counter_years+1),ll] = ncvar_get(ncin, varNames_tmp,collapse_degen=T )[idxLon,idxLat,ll]',sep='') ) )
          }
          eval( parse(text= paste('colnames(','var',vv,')[',ll,']=("',varNames[vv],'_0',ll,'")' ,sep='')))
       }
    }# end variable loop
    
    counter_years=  counter_years + 1
    counter_days=  counter_days + 365
    nc_close(ncin)
  }# end year/file loop
  
  # Make time vector
  if(varNames %in% c('ABS_CANOPY_LIGHT', 'TRA_CANOPY_LIGHT', 'CUM_CANOPY_LIGHT')){
    syear=min(options$years); eyear=max(options$years)
    TIME = RCHIDEE.day2decyy(syear, eyear)
  }else{
    TIME = options$years
  }

  # Make experiment vector
  expName=rep(experimentName, counter_years)  
    
  # prepare output
  out=as.character()
  for (ii in 1:length(varNames)){
    ##eval( parse(text= paste('var',ii,'=','matrix(nrow=', nyears,',','ncol=',ncirc, ')', sep='') ) )
    #out[ii] = paste(paste(varNames[ii],'=var',ii, sep='',collapse = ",") )
    out[ii] = paste('var',ii, sep='',collapse = ",")
  }
  strOut=paste(rbind(out), collapse=',')

  #return( eval(parse( text=paste('list(', strOut ,')', sep='' ))) )
  return( eval(parse( text=paste('data.frame(', 'TIME=TIME', ',', 'EXPERIMENT=expName', ',', strOut ,')', sep='' ))) )
    
} # end function





####################################################################################  
RCHIDEE.getCRUNCEPhalfdLandPoint = function(lati, loni){
####################################################################################  
# RCHIDEE.getCRUNCEPhalfdLandPoint. Retrieve corresponding land point ID and its array index  
#        for a given lat/lon from a 0.5-degree CRUN-CEP forcing file. 
#        Function taken from Yiying.
#  
# 
# Call:
#      RCHIDEE.getCRUNCEPhalfdLandPoint(lati, loni)
#
# Input:
#         lati : Latitude for a given site in decimal degrees 
#         loni : Longitude for a given site in decimal degrees
#
# Output: 
#         landpoint_index: Array index for the pixel in 'land' array
#
# Example call:
#      lat_cru=48.4; lon_cru=2.7;
#      lpointINDEX_fon = RCHIDEE.getCRUNCEPhalfdLandPoint(lati, loni)
#      
# Note for postprocessing:
#
#
#
# Dependencies:
#      require(ncdf4)
#
# @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
# @   23 February 2015, Paris
#  
  
# read dimension from simplified forcing file provided by Yiying 
ncin = nc_open("/home/satellites5/jbari/code/cru_ncep/land.nc")

#set dimension & time steps for working
nx = ncin$dim[["x"]]$len
ny = ncin$dim[["y"]]$len
#nt = ncin$dim[["tstep"]]$len  
np = ncin$dim[["land"]]$len  


#======= import dimesions ====================================
# Import  cru-zoom grid data-set from ncin file
#============================================================= 

#d1 = ncvar_get(ncin,"x")
#d2 = ncvar_get(ncin,"y") # <=== COULDNT READ IT IN NCDF4!
d1 = c(1:nx)
d2 = c(1:ny)
nav_lon  = ncvar_get(ncin,"nav_lon",  start=c(1,1), count=c(nx,ny))
nav_lat  = ncvar_get(ncin,"nav_lat",  start=c(1,1), count=c(nx,ny))
land = ncvar_get(ncin,"land")


#==================================================================
#========Start compressing Procedure ==============================	
# start two counter for the land_id  and exact sequency of matrix
#  cont1:
#  for example  a 3 by 3 array the exact secquency of each point is 
#  1 2 3
#  4 5 6
#  7 8 9
#  cont2:
#  for eaxmple with 3 by 3 array only have 3 point at certain exact point 
#  0 0 1
#  0 2 0
#  0 0 3    
#so we need to create a funtion land() for the mapping
#land(1) = 3
#land(2) = 5
#land(3) = 9
#==================================================================

# get points to look for
lat_vector = nav_lat[1,]
lon_vector = nav_lon[,1]
#find_x = as.numeric(lon_cru) # e.g., 45.25 N 
#find_y = as.numeric(lat_cru) # e.g., 4.25 E

# find nearest nodes
idx_y = which(abs(lat_vector-lati)== min(abs(lat_vector-lati)))
idx_x = which(abs(lon_vector-loni)== min(abs(lon_vector-loni)))
find_y = lat_vector[idx_y]
find_x = lon_vector[idx_x]

# throw error if two grids are found at the same distance
if( length(find_x) != 1 ) stop( paste('ERROR: THERE ARE TWO OR MORE LONGITUDE GRID-NODES AT THE SAME DISTANCE FROM YOUR POINT', sep=' '))
if( length(find_y) != 1 ) stop( paste('ERROR: THERE ARE TWO OR MORE LATITUDE GRID-NODES AT THE SAME DISTANCE FROM YOUR POINT', sep=' '))

find_land = 0
cont = 0
ld_find = FALSE
for ( jj in 1:ny )  {
for ( ii in 1:nx )  { 
            cont = cont + 1 
            if (  (nav_lon[ii,jj] == find_x ) && (nav_lat[ii,jj] == find_y) )  {
            ld_find   = TRUE
            find_land = cont
            #print( paste("land_point:", find_land, sep="") )
            #find_land[jj*ii]   # sequencial with the all index  
               for (kk in 1:np) {
                   if (land[[kk]] == find_land ) {
                       print(paste("lon:",nav_lon[ii,jj], "lat:",nav_lat[ii,jj])) 
                       #print(paste("landpoint_id:", kk," <seq_regular_matrix_count>", find_land))
                       print(paste("landpoint_id:", find_land," <index_in_land_array>",kk ))
                       landpoint_index = kk  # this is the array index value for the pixel in the "land" array that I want for extraction with NCKS
                       # check this to use ncks: http://stackoverflow.com/questions/25183867/extracting-variables-from-gridded-netcdf-data-into-ascii 
                     }
               }
            } # for land mass
} # for longitudei
} # for latitiude

if (ld_find != TRUE) {
   print(paste("No land point find ~"))
#   print(paste("lon:",find_x, "lat:",find_y)) 
}

nc_close(ncin)

# check result -- I should get find_y and find_x coordinates!S 
# > dim(nav_lat)
# [1] 720 360
#
# (nav_lat)[find_land] # same as # z=as.vector((nav_lat))
# (nav_lon)[find_land]
# iidx=which(land==189576) # find_land=189576
# land[iidx]

# define output
return(landpoint_index)

} #end function




####################################################################################  
RCHIDEE.extractCRUNCEPhalfdForcing4Point = function(landpoint_id){
####################################################################################  
# RCHIDEE.extractCRUNCEPhalfdForcing4Point. Extract full CRU-NCEP daily forcing for a land point 
#  
# 
# Call:
#      RCHIDEE.extractCRUNCEPhalfdForcing4Point(landpoint_id)
#
# Input:
#         options: a list object with settings.
#         option$landpoint_index  : Land point for a site obtained with RCHIDEE.getCRUNCEPhalfdLandPoint(lati, loni).
#         options$path2forcings   : Path to forcing files
#         options$path4output     : Path to store output files
#         options$siteName        : Short site name without spaces
#         options$years           : Sequence of years to process, e.g., c(1901:1955)
#
# Output: Yearly NetCDF forcing file for the site
#
# Example call:
#         lati=48.4; loni=2.7;
#         lpointIDX_fon = RCHIDEE.getCRUNCEPhalfdLandPoint(lati, loni)
#         options=list()
#         options$landpoint_id  = 25408
#         options$path2forcings = '/home/orchideeshare/igcmg/IGCM/STORAGE/BC/OL2/CRU-NCEP/v5.3/halfdeg/'
#         options$path4output   = '/home/satellites5/jbari/forcing/sites/CRU-NCEP/v5.3/halfdeg/'
#         options$siteName      = 'fon'
#         options$years         = c(1901:2013)
#         RCHIDEE.extractCRUNCEPhalfdForcing4Point(options)
#
#
# Dependencies:
#      ncks in the system
#
# @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
# @   23 February 2015, Paris
#  
  
  nyears=length(options$years)
  setwd(options$path4output)

  # Loop over yearly files and concatenate daily data for each variable
  counter_years=0
  for (yy in 1:nyears){
    infile_yy = paste(options$path2forcings, 'cruncep_halfdeg_', options$years[yy], '.nc', sep='');
    ofile_yy  = paste(options$path4output, 'cruncep_halfdeg_landpoint_', options$siteName , '_', options$landpoint_id, '_' ,options$years[yy], '.nc', sep='');
    print(ofile_yy)
    system( paste('ncks -O -h -a -d land,', options$landpoint_id, ' ',infile_yy, ' -O ', ofile_yy, sep='') );
    
    # clean .tmp files
    system( 'rm -rf *.tmp' )
    counter_years=  counter_years + 1
  }# end year/file loop

} # end function


####################################################################################  
RCHIDEE.getERA5halfdLandPoint = function(lati, loni){
####################################################################################  
# RCHIDEE.getERA5halfdLandPoint. Retrieve corresponding land point ID and its array index  
#        for a given lat/lon from a 0.5-degree CRUN-CEP forcing file. 
#        Function taken from Yiying.
#  
# 
# Call:
#      RCHIDEE.getERA5halfdLandPoint(lati, loni)
#
# Input:
#         lati : Latitude for a given site in decimal degrees 
#         loni : Longitude for a given site in decimal degrees
#
# Output: 
#         landpoint_index: Array index for the pixel in 'land' array
#
# Example call:
#      lat_cru=48.4; lon_cru=2.7;
#      lpointINDEX_fon = RCHIDEE.getERA5halfdLandPoint(lati, loni)
#      
# Note for postprocessing:
#
#
#
# Dependencies:
#      require(ncdf4)
#
# @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
# @   22 November 2017, Paris
#  
  
# read dimension from simplified forcing file provided by Yiying 
#ncin = nc_open("/home/satellites5/jbari/code/cru_ncep/land.nc")
ncin = nc_open("/home/satellites10/jbari/chile/daymean_ERA5_2016.nc")

#set dimension & time steps for working
nx = ncin$dim[["x"]]$len
ny = ncin$dim[["y"]]$len
#nt = ncin$dim[["tstep"]]$len  
np = ncin$dim[["land"]]$len  


#======= import dimesions ====================================
# Import  cru-zoom grid data-set from ncin file
#============================================================= 

#d1 = ncvar_get(ncin,"x")
#d2 = ncvar_get(ncin,"y") # <=== COULDNT READ IT IN NCDF4!
d1 = c(1:nx)
d2 = c(1:ny)
nav_lon  = ncvar_get(ncin,"nav_lon",  start=c(1,1), count=c(nx,ny))+0.25
nav_lat  = ncvar_get(ncin,"nav_lat",  start=c(1,1), count=c(nx,ny))-0.25
land = ncvar_get(ncin,"land")


#==================================================================
#========Start compressing Procedure ==============================	
# start two counter for the land_id  and exact sequency of matrix
#  cont1:
#  for example  a 3 by 3 array the exact secquency of each point is 
#  1 2 3
#  4 5 6
#  7 8 9
#  cont2:
#  for eaxmple with 3 by 3 array only have 3 point at certain exact point 
#  0 0 1
#  0 2 0
#  0 0 3    
#so we need to create a funtion land() for the mapping
#land(1) = 3
#land(2) = 5
#land(3) = 9
#==================================================================

# get points to look for
lat_vector = nav_lat[1,]
lon_vector = nav_lon[,1]
#find_x = as.numeric(lon_cru) # e.g., 45.25 N 
#find_y = as.numeric(lat_cru) # e.g., 4.25 E

# find nearest nodes
idx_y = which(abs(lat_vector-lati)== min(abs(lat_vector-lati)))
idx_x = which(abs(lon_vector-loni)== min(abs(lon_vector-loni)))
find_y = lat_vector[idx_y]
find_x = lon_vector[idx_x]

# throw error if two grids are found at the same distance
if( length(find_x) != 1 ) stop( paste('ERROR: THERE ARE TWO OR MORE LONGITUDE GRID-NODES AT THE SAME DISTANCE FROM YOUR POINT', sep=' '))
if( length(find_y) != 1 ) stop( paste('ERROR: THERE ARE TWO OR MORE LATITUDE GRID-NODES AT THE SAME DISTANCE FROM YOUR POINT', sep=' '))

find_land = 0
cont = 0
ld_find = FALSE
for ( jj in 1:ny )  {
for ( ii in 1:nx )  { 
            cont = cont + 1 
            if (  (nav_lon[ii,jj] == find_x ) && (nav_lat[ii,jj] == find_y) )  {
            ld_find   = TRUE
            find_land = cont
            #print( paste("land_point:", find_land, sep="") )
            #find_land[jj*ii]   # sequencial with the all index  
               for (kk in 1:np) {
                   if (land[[kk]] == find_land ) {
                       print(paste("lon:",nav_lon[ii,jj], "lat:",nav_lat[ii,jj])) 
                       #print(paste("landpoint_id:", kk," <seq_regular_matrix_count>", find_land))
                       print(paste("landpoint_id:", find_land," <index_in_land_array>",kk ))
                       landpoint_index = kk  # this is the array index value for the pixel in the "land" array that I want for extraction with NCKS
                       # check this to use ncks: http://stackoverflow.com/questions/25183867/extracting-variables-from-gridded-netcdf-data-into-ascii 
                     }
               }
            } # for land mass
} # for longitudei
} # for latitiude

if (ld_find != TRUE) {
   print(paste("No land point find ~"))
#   print(paste("lon:",find_x, "lat:",find_y)) 
}

nc_close(ncin)

# check result -- I should get find_y and find_x coordinates!S 
# > dim(nav_lat)
# [1] 720 360
#
# (nav_lat)[find_land] # same as # z=as.vector((nav_lat))
# (nav_lon)[find_land]
# iidx=which(land==189576) # find_land=189576
# land[iidx]



# define output
return(landpoint_index)

} #end function


####################################################################################  
RCHIDEE.extractERA5halfdForcing4Point = function(landpoint_id){
####################################################################################  
# RCHIDEE.extractERA5halfdForcing4Point. Extract full ERA5 daily forcing for a land point 
#  
# 
# Call:
#      RCHIDEE.extractERA5halfdForcing4Point(landpoint_id)
#
# Input:
#         options: a list object with settings.
#         option$landpoint_index  : Land point for a site obtained with RCHIDEE.getERA5halfdLandPoint(lati, loni).
#         options$path2forcings   : Path to forcing files
#         options$path4output     : Path to store output files
#         options$siteName        : Short site name without spaces
#         options$years           : Sequence of years to process, e.g., c(1901:1955)
#
# Output: Yearly NetCDF forcing file for the site
#
# Example call:
#         lati=48.4; loni=2.7;
#         lpointIDX_fon = RCHIDEE.getERA5halfdLandPoint(lati, loni)
#         options=list()
#         options$landpoint_id  = 25408
#         options$path2forcings = '/home/satellites10/jbari/chile/'
#         options$path4output   = '/home/satellites10/jbari/chile/out/'
#         options$siteName      = 'chl'
#         options$years         = c(2010:2016)
#         RCHIDEE.extractERA5halfdForcing4Point(options)
#
#
# Dependencies:
#      ncks in the system
#
# @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
# @   22 November 2017, Paris
#  
  
  nyears=length(options$years)
  setwd(options$path4output)

  # Loop over yearly files and concatenate daily data for each variable
  counter_years=0
  for (yy in 1:nyears){
    #infile_yy = paste(options$path2forcings, 'daymean_ERA5_', options$years[yy], '.nc', sep='');
    #ofile_yy  = paste(options$path4output, 'daymean_ERA5_halfdeg_landpoint_', options$siteName , '_', options$landpoint_id, '_' ,options$years[yy], '.nc', sep='');
    infile_yy = paste(options$path2forcings, 'ERA5_', options$years[yy], '.nc', sep='');
    ofile_yy  = paste(options$path4output, 'hourly_ERA5_halfdeg_landpoint_', options$siteName , '_', options$landpoint_id, '_' ,options$years[yy], '.nc', sep='');
    print(ofile_yy)
    system( paste('ncks -O -h -a -d land,', options$landpoint_id, ' ',infile_yy, ' -O ', ofile_yy, sep='') );
    
    # clean .tmp files
    system( 'rm -rf *.tmp' )
    counter_years=  counter_years + 1
  }# end year/file loop

} # end function




####################################################################################  
RCHIDEEsec.catSiteMonthlyVars = function(options){
  ####################################################################################  
  # RCHIDEEsec.catSiteMonthlyVars. Concatenate monthly time series for a list of 
  #    variables from yearly sechiba history file. 
  #  
  # 
  # Call:
  #      RCHIDEEsec.catSiteMonthlyVars(options)
  #
  # Input:
  #       options: a list object with the names of the input files and run settings.
  #         options$path2file    : Path to history files
  #         options$expFilename  : Experiment name prefix for history files, e.g., 'run' (see example below)
  #         options$years        : Sequence of years for history files, e.g., c(1901:1955)
  #         options$pftId        : Index for the PFT to extract (look at your run.def used)
  #         options$idxLat       : Index for latitude of the gridbox to extract
  #         options$idxLon       : Index for longitude of the gridbox to extract
  #         options$variables    : List of variable names to process, e.g., c('NPP','GPP','TOTAL_M')
  #         options$experimentName: String or numeric ID for the current experiment
  #
  # Output: The output is a list of the variables specified in options$variables.
  #         dimensions(nyears*12, 1) [rows, cols]
  #
  # Example call:
  #      options=list() 
  #      options$path2file    = '/home/scratch01/jbari/IGCM_OUT/OL2/PROD/secsto/tstrec03/SRF/Output/'
  #      options$expFilename  = 'run';  
  #      options$years        = c(1901:1955)
  #      options$pftId        = 3
  #      options$idxLat       = 1
  #      options$idxLon       = 1
  #      options$variables    = c('rainfall','snow')
  #      options$experimentName = 'EXP1'
  #      BA = RCHIDEEsec.catSiteMontlyVars(options)
  #
  # Dependencies:
  #      require(ncdf)
  #      RCHIDEE.month2decyy
  #
  # @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
  # @   18 August 2015, Paris
  #  
  
  # set options 
  nyears = length(options$years)
  nmonths  = nyears * 12
  pftId  = options$pftId
  idxLat = options$idxLat
  idxLon = options$idxLon
  varNames = options$variables
  experimentName = options$experimentName
  
  # preallocate output arrays [nmonths x 1]
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', nmonths,',','ncol=',ncirc, ')', sep='') ) )
    eval( parse(text= paste('var',ii,'=','matrix(nrow=', nmonths,',','ncol=1', ')', sep='') ) )
  }
  
  # Loop over yearly files and concatenate daily data for each variable
  counter_months=0
  for (yy in 1:nyears){
    infile_yy = paste(options$path2file, options$expFilename, '_', options$years[yy], '0101_', options$years[yy], '1231_1M_sechiba_history.nc', sep='');
    ncin      = nc_open(infile_yy)
    
    # loop over variables
    for (vv in 1:length(varNames)){
      eval( parse(text= paste('var',vv,'[(counter_months+1):(counter_months+12)] = ncvar_get(ncin, varNames[vv], collapse_degen=F )[idxLon,idxLat,]',sep='') ) )
    }# end variable loop
    
    counter_months=  counter_months + 12
    nc_close(ncin)
  }# end year/file loop

  # Make time vector
  syear=min(options$years); eyear=max(options$years)
  TIME = RCHIDEE.month2decyy(syear, eyear)
  
  # Make experiment vector
  expName=rep(experimentName, counter_months)  
  
  # prepare output
  out=as.character()
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', nmonths,',','ncol=',ncirc, ')', sep='') ) )
    out[ii] = paste(paste(varNames[ii],'=var',ii, sep='',collapse = ",") )
  }
  strOut=paste(rbind(out), collapse=', ')
  #return( eval(parse( text=paste('list(', strOut ,')', sep='' ))) )
  return( eval(parse( text=paste('data.frame(', 'TIME=TIME', ',', 'EXPERIMENT=expName', ',', strOut ,')', sep='' ))) )
  
} # end function




####################################################################################  
RCHIDEEsec.catSiteMonthlySoilLayerVars = function(options){
  ####################################################################################  
  # RCHIDEEsec.catSiteMonthlySoilLayerVars. Concatenate monthly time series for a list of 
  #    11 soil layer variables from yearly sechiba history file with dimension ( 1 1 11 / lon lat solay). 
  #  
  # 
  # Call:
  #      RCHIDEEsec.catSiteMonthlySoilLayerVars(options)
  #
  # Input:
  #       options: a list object with the names of the input files and run settings.
  #         options$path2file    : Path to history files
  #         options$expFilename  : Experiment name prefix for history files, e.g., 'run' (see example below)
  #         options$years        : Sequence of years for history files, e.g., c(1901:1955)
  #         options$idxLat       : Index for latitude of the gridbox to extract
  #         options$idxLon       : Index for longitude of the gridbox to extract
  #         options$variables    : List of 11-layer variable names to process , e.g., c('psi_moy','RootDist','SoilSat','SoilMoist')
  #         options$experimentName: String or numeric ID for the current experiment
  #
  # Output: The output is a list of the variables specified in options$variables.
  #         dimensions(nyears*12, 1) [rows, cols]
  #
  # Example call:
  #      options=list() 
  #      options$path2file    = '/home/scratch01/jbari/IGCM_OUT/OL2/PROD/secsto/tstrec03/SRF/Output/'
  #      options$expFilename  = 'run';  
  #      options$years        = c(1901:1955)
  #      options$idxLat       = 1
  #      options$idxLon       = 1
  #      options$variables    = c('SoilMoist','psi_moy')
  #      options$experimentName = 'EXP1'
  #      BA = RCHIDEEsec.catSiteMonthlySoilLayerVars(options)
  #
  # Dependencies:
  #      require(ncdf)
  #      RCHIDEE.month2decyy
  #
  # @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
  # @   8 December 2024, Paris
  #  
  
  # set options 
  nyears = length(options$years)
  nmonths  = nyears * 12
  idxLat = options$idxLat
  idxLon = options$idxLon
  varNames = options$variables
  experimentName = options$experimentName
  nsolay = 11 # number of soil layers
  
  # preallocate output arrays [nmonths x 1]
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', nmonths,',','ncol=1', ')', sep='') ) )
    eval( parse(text= paste('var',ii,'=','data.frame(matrix(nrow=', nmonths,',','ncol=', nsolay,'))', sep='') ) )
  }
  
  # Loop over yearly files and concatenate daily data for each variable
  counter_months=0
  for (yy in 1:nyears){
    infile_yy = paste(options$path2file, options$expFilename, '_', options$years[yy], '0101_', options$years[yy], '1231_1M_sechiba_history.nc', sep='');
    ncin      = nc_open(infile_yy)
    
    for (vv in 1:length(varNames)){
      for(ll in 1:nsolay){
          # set variable name for a given soil layer
          varNames_tmp=paste(varNames[vv], sep='')
          # read variable into the preallocated data frame
          eval( parse(text= paste('var',vv,'[(counter_months+1):(counter_months+12),ll] = ncvar_get(ncin, varNames_tmp,collapse_degen=F )[idxLon,idxLat,ll,]',sep='') ) )
          # add name to the corresponding column
          if (ll <10){
             eval( parse(text= paste('colnames(','var',vv,')[',ll,']=("','solay','_0',ll,'")' ,sep='')))
          }else{
          eval( parse(text= paste('colnames(','var',vv,')[',ll,']=("','solay_',ll,'")' ,sep='')))
          } # end if
        } # end layer loop
    }# end variable loop

    counter_months=  counter_months + 12
    nc_close(ncin)
  }# end year/file loop

  # Make time vector
  syear=min(options$years); eyear=max(options$years)
  TIME = RCHIDEE.month2decyy(syear, eyear)
  
  # Make experiment vector
  expName=rep(experimentName, counter_months)  
  
  # prepare output
  out=as.character()
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', nmonths,',','ncol=',ncirc, ')', sep='') ) )
    out[ii] = paste(paste(varNames[ii],'=var',ii, sep='',collapse = ",") )
  }
  strOut=paste(rbind(out), collapse=', ')
  #return( eval(parse( text=paste('list(', strOut ,')', sep='' ))) )
  return( eval(parse( text=paste('data.frame(', 'TIME=TIME', ',', 'EXPERIMENT=expName', ',', strOut ,')', sep='' ))) )
  
} # end function




####################################################################################  
RCHIDEEsec.catSiteMonthlyVars_decade = function(options){
  ####################################################################################  
  # RCHIDEEsec.catSiteMonthlyVars_decade. Concatenate monthly time series for a list of 
  #    variables from yearly sechiba history file. 
  #  
  # 
  # Call:
  #      RCHIDEEsec.catSiteMonthlyVars_decade(options)
  #
  # Input:
  #       options: a list object with the names of the input files and run settings.
  #         options$path2file    : Path to history files
  #         options$expFilename  : Experiment name prefix for history files, e.g., 'run' (see example below)
  #         options$years        : First and last year of history file sequence, e.g., c(2001,2010)
  #         options$pftId        : Index for the PFT to extract (look at your run.def used)
  #         options$idxLat       : Index for latitude of the gridbox to extract
  #         options$idxLon       : Index for longitude of the gridbox to extract
  #         options$variables    : List of variable names to process, e.g., c('NPP','GPP','TOTAL_M')
  #         options$experimentName: String or numeric ID for the current experiment
  #
  # Output: The output is a list of the variables specified in options$variables.
  #         dimensions(nyears*12, 1) [rows, cols]
  #
  # Example call:
  #      options=list() 
  #      options$path2file    = '/home/scratch01/jbari/IGCM_OUT/OL2/PROD/secsto/tstrec03/SRF/Output/'
  #      options$expFilename  = 'run';  
  #      options$years        = c(2001,2010)
  #      options$pftId        = 3
  #      options$idxLat       = 1
  #      options$idxLon       = 1
  #      options$variables    = c('rainfall','snow')
  #      options$experimentName = 'EXP1'
  #      BA = RCHIDEEsec.catSiteMontlyVars(options)
  #
  # Dependencies:
  #      require(ncdf)
  #      RCHIDEE.month2decyy
  #
  # @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
  # @   30 November 2023, Schiermonnikoog, The Netherlands
  #  
  
  # set options 
  nyears = length(options$years[1]:options$years[2])
  nfiles = 1 # only one file for the range of years
  nmonths  = nyears * 12
  pftId  = options$pftId
  idxLat = options$idxLat
  idxLon = options$idxLon
  varNames = options$variables
  experimentName = options$experimentName
  
  # preallocate output arrays [nmonths x 1]
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', nmonths,',','ncol=',ncirc, ')', sep='') ) )
    eval( parse(text= paste('var',ii,'=','matrix(nrow=', nmonths,',','ncol=1', ')', sep='') ) )
  }
  
  # Loop over files and concatenate monthly data for each variable
  counter_months=0
  for (yy in 1: nfiles){
    infile_yy = paste(options$path2file, options$expFilename, '_', options$years[1], '0101_', options$years[2], '1231_1M_sechiba_history.nc', sep='');
    ncin      = nc_open(infile_yy)
    
    # loop over variables
    for (vv in 1:length(varNames)){
      if(varNames[vv] == 'Ca'){ # Ca is not by PFT
        eval( parse(text= paste('var',vv,'[(counter_months+1):(counter_months+nmonths)] = ncvar_get(ncin, varNames[vv], collapse_degen=F )[idxLon,idxLat,]',sep='') ) )
      }else{
        eval( parse(text= paste('var',vv,'[(counter_months+1):(counter_months+nmonths)] = ncvar_get(ncin, varNames[vv], collapse_degen=F )[idxLon,idxLat,pftId,]',sep='') ) )
      }
    }# end variable loop
    
    #counter_months=  counter_months + nmonths
    nc_close(ncin)
  }# end year/file loop

  # Make time vector
  syear=min(options$years); eyear=max(options$years)
  TIME = RCHIDEE.month2decyy(syear, eyear)
  
  # Make experiment vector
  expName=rep(experimentName, counter_months)  
  
  # prepare output
  out=as.character()
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', nmonths,',','ncol=',ncirc, ')', sep='') ) )
    out[ii] = paste(paste(varNames[ii],'=var',ii, sep='',collapse = ",") )
  }
  strOut=paste(rbind(out), collapse=', ')
  #return( eval(parse( text=paste('list(', strOut ,')', sep='' ))) )
  return( eval(parse( text=paste('data.frame(', 'TIME=TIME', ',', 'EXPERIMENT= experimentName', ',', strOut ,')', sep='' ))) )
  
} # end function



####################################################################################  
RCHIDEEsec.catSiteAnnualVars = function(options){
####################################################################################  
  # RCHIDEEsec.catSiteAnnualVars Concatenate annual time series for a list of 
  #    variables from yearly stomate or sechiba history files. 
  #  
  # 
  # Call:
  #      RCHIDEEsto.catSiteAnnualVars(options)
  #
  # Input:
  #       options: a list object with the names of the input files and run settings.
  #         options$path2file    : Path to history files
  #         options$expFilename  : Experiment name prefix for history files, e.g., 'run' (see example below)
  #         options$years        : Sequence of years for history files, e.g., c(1901:1955)
  #         options$pftId        : Index for the PFT to extract (look at your run.def used)
  #         options$idxLat       : Index for latitude of the gridbox to extract
  #         options$idxLon       : Index for longitude of the gridbox to extract
  #         options$variables    : List of variable names to process, e.g., c('NPP','GPP','TOTAL_M')
  #         options$experimentName: String or numeric ID for the current experiment
  #
  # Output: The output is a list of the variables specified in options$variables.
  #         dimensions(nyears, 1) [rows, cols]
  #
  # Example call:
  #      options=list() 
  #      options$path2file    = '/home/scratch01/jbari/IGCM_OUT/OL2/PROD/secsto/tstrec03/SBG/Output/'
  #      options$expFilename  = 'run';  
  #      options$years        = c(1901:1955)
  #      options$pftId        = 3
  #      options$idxLat       = 1
  #      options$idxLon       = 1
  #      options$variables    = c('NPP','GPP','TOTAL_M')
  #      options$experimentName = 'EXP1'
  #      BA = RCHIDEEsec.catSiteAnnualVars(options)
  #      
  # Dependencies:
  #      require(ncdf)
  #
  # @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
  # @   20 October 2023, Paris
  #  
  
  # set options 
  nyears = length(options$years)
  nmonths  = nyears
  pftId  = options$pftId
  idxLat = options$idxLat
  idxLon = options$idxLon
  varNames = options$variables
  experimentName = options$experimentName
  
  # preallocate output arrays [nmonths x 1]
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', nmonths,',','ncol=',ncirc, ')', sep='') ) )
    eval( parse(text= paste('var',ii,'=','matrix(nrow=', nmonths,',','ncol=1', ')', sep='') ) )
  }
  
  # Loop over yearly files and concatenate daily data for each variable
  counter_years=0
  for (yy in 1:nyears){
    infile_yy = paste(options$path2file, options$expFilename, '_', options$years[yy], '0101_', options$years[yy], '1231_1Y_sechiba_history.nc', sep='');
    ncin      = nc_open(infile_yy)
    
    # loop over variables
    for (vv in 1:length(varNames)){
      eval( parse(text= paste('var',vv,'[(counter_years+1):(counter_years+1)] = ncvar_get(ncin, varNames[vv], collapse_degen=F )[idxLon,idxLat,pftId,]',sep='') ) )
    }# end variable loop
    
    counter_years=  counter_years + 1
    nc_close(ncin)
  }# end year/file loop
  
  # Make time vector
  TIME = options$years
  
  # Make experiment vector
  expName=rep(experimentName, counter_years)  
    
  # prepare output
  out=as.character()
  for (ii in 1:length(varNames)){
    #eval( parse(text= paste('var',ii,'=','matrix(nrow=', nyears,',','ncol=',ncirc, ')', sep='') ) )
    out[ii] = paste(paste(varNames[ii],'=var',ii, sep='',collapse = ",") )
  }
  strOut=paste(rbind(out), collapse=', ')
  #return( eval(parse( text=paste('list(', strOut ,')', sep='' ))) )
  return( eval(parse( text=paste('data.frame(', 'TIME=TIME', ',', 'EXPERIMENT=expName', ',', strOut ,')', sep='' ))) )
    
} # end function


