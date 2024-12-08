# SCRIPT. Use of RCHIDEE library to visualise the output of pixel-level ORCHIDEE4 runs
#
############################################
# CHECK RCHIDEE
############################################
# USAGE IN OBELIX (7 Dec 2023):
# * module load R/4.0.3-gcc8
# * R 
#
# Description of PFTs
# PFT_NAME__01=SoilBareGlobal
# PFT_NAME__02=BroadLeavedEvergreenTropical
# PFT_NAME__03=BroadLeavedRaingreenTropical
# PFT_NAME__04=NeedleleafEvergreenTemperate
# PFT_NAME__05=BroadLeavedEvergreenTemperate
# PFT_NAME__06=BroadLeavedSummergreenTemperate
# PFT_NAME__07=NeedleleafEvergreenBoreal
# PFT_NAME__08=BroadLeavedSummergreenBoreal
# PFT_NAME__09=LarixSpBoreal
# PFT_NAME__10=C3GrassTemperate
# PFT_NAME__11=C4GrassTemperate
# PFT_NAME__12=C3AgricultureTemperate
# PFT_NAME__13=C4AgricultureTemperate

# load library
rm(list = ls())
source('/home/satellites5/USER/code/RCHIDEE/RchideeLibrary_v1p3.R')

########################################################################
# 1. READ AND PLOT VARIABLES OF NORMAL STOMATE HISTORY FILE
########################################################################
## First, see list of variables and dimensions in normal STOMATE history file
infile = '/home/satellites5/USER/scratch.output/IGCM_OUT/OL2/PROD/ACF/hist0.trw/SBG/Output/MO/hist0.trw_20200101_20201231_1M_stomate_history.nc';
RCHIDEE.listVariables(infile)

# now set required options for function and variables to read
options=list() # do not touch this line!
options$path2file    = '/home/satellites5/USER/scratch.output/IGCM_OUT/OL2/PROD/ACF/hist0.trw/SBG/Output/MO/';
options$expFilename  = 'hist0.trw';  # it must match the start of the output file
options$years        = c(1901:2020)
options$pftId        = 5
options$idxLat       = 1
options$idxLon       = 1
options$ncir         = 3
options$experimentName = 'ACF.trw'
options$variables    = c('VEGET_MAX','LAI_MEAN','NPP','GPP','HET_RESP','BA','BA_INV','AGE','AGE_STAND', # only variables with dimensions [lon lat veget] (i.e., 1 1 15)
                         'HEIGHT','HEIGHT_DOM','HEIGHT_INV','DIAMETER','DIA_DOM','DIA_INV','DIAMETER_DOM','IND_INV','TREE_MORTALITY','IND_ESTAB','RECRUITS_IND',
                         'RESERVE_M_c', 'LEAF_M_c','ROOT_M_c','TOTAL_M_c', 
                         'TOTAL_AB_M_c','IND','LIGHT_TRAN_SEASON',
                         'KF','C0_ALLOC','GTEMP_ALLOC','K_LATOSA_ADAPT',
                         'VCMAX','VCMAX_NEW','LEAF_AGE')
# then read orchidee4 output files using the options provided above
STO = RCHIDEEsto.catSiteMonthlyVars(options)

# now do some plots
plot(STO$TIME,STO$GPP,t='l', ylim=c(0,5))

plot(STO$TIME,STO$RESERVE_M_c,t='l')

plot(STO$TIME, STO$RECRUITS_IND*10000*365, t='l', xlab='', ylab='# recruits/ha/year')



########################################################################
# 2. READ AND PLOT VARIABLES PER CIRCUMFERENCE CLASS (STOMATE 4DIM)
########################################################################
# first, list all variables in the orchidee 4dim stomate output files
infile = '/home/satellites5/USER/scratch.output/IGCM_OUT/OL2/PROD/ACF/hist0.trw/SBG/Output/MO/hist0.trw_20200101_20201231_1M_stomate_history_4dim.nc';
RCHIDEE.listVariables(infile)

# now set required options for function and variables to read
options=list() 
options$path2file    = '/home/satellites5/USER/scratch.output/IGCM_OUT/OL2/PROD/ACF/hist0.trw/SBG/Output/MO/';
options$expFilename  = 'hist0.trw';  
options$years        = c(1901:2020)
options$pftId        = 5
options$idxLat       = 1
options$idxLon       = 1
options$ncir         = 3
options$experimentName = 'hist0.trw'
options$variables    = c('CCTRW','CCDIAMETER','CCHEIGHT','CCIND','CCBA','CCDELTABA','CCDIAHOR','CCDIAVER','CCCROWN',
                         'CCCROWNVOL','CCMORTALITY','CCMORTALITY_M_AB_c','CCMORTALITY_M_BE_c','CCTOTAL_M_c','CLEVEL_HEIGHT',
                         'LAI_PER_LEVEL','LIGHT_ABS_TO_LEVEL','LIGHT_TRAN_TO_LEVEL','ROOT_PROF_STRUC','ROOT_PROF_FUNC')

# then read orchidee4 output files using the options provided above
S4D = RCHIDEEsto4dim.catSiteMonthlyVars(options)


# plot tree size and growth data
plot(S4D$TIME, S4D$CCTRW.ncir_003*100, t='l')
lines(S4D$TIME, S4D$CCTRW.ncir_002*100, col='red')
lines(S4D$TIME, S4D$CCTRW.ncir_001*100, col='blue')

plot(S4D$TIME, cumsum(S4D$CCTRW.ncir_003*100), t='l')

plot(S4D$TIME,S4D$LIGHT_TRAN_TO_LEVEL.LIGHT_TRAN_TO_LEVEL_003, t='l')
plot(S4D$TIME,S4D$CCDIAHOR.ncir_001, t='l')

plot(S4D$TIME, S4D$CCHEIGHT.ncir_003, t='l')

plot(S4D$TIME,S4D$CCHEIGHT.ncir_003, t='l',ylim=c(0,10))
lines(S4D$TIME, S4D$CCHEIGHT.ncir_002, col='red')
lines(S4D$TIME, S4D$CCHEIGHT.ncir_001, col='blue')


plot(S4D$TIME,S4D$CCCROWNVOL.ncir_003, t='l',ylim=c(0,200))
lines(S4D$TIME, S4D$CCCROWNVOL.ncir_002, col='red')
lines(S4D$TIME, S4D$CCCROWNVOL.ncir_001, col='blue')





########################################################################
# 3. READ AND PLOT VARIABLES OF NORMAL SECHIBA HISTORY FILE
########################################################################
## First, see list of variables and dimensions in normal STOMATE history file
infile = '/home/satellites5/USER/scratch.output/IGCM_OUT/OL2/PROD/ACF/hist0.trw/SRF/Output/MO/hist0.trw_20200101_20201231_1M_sechiba_history.nc';
RCHIDEE.listVariables(infile)

# now set required options for function and variables to read
options=list() # do not touch this line!
options$path2file    = '/home/satellites5/USER/scratch.output/IGCM_OUT/OL2/PROD/ACF/hist0.trw/SRF/Output/MO/';
options$expFilename  = 'hist0.trw';  # it must match the start of the output file
options$years        = c(1901:2020)
options$pftId        = 5
options$idxLat       = 1
options$idxLon       = 1
options$ncir         = 3
options$experimentName = 'ACF'
options$variables    = c('evap','riverflow','runoff','drainage','rain','TWS','tran','humtot_top','mrso', # only variables with dimensions [lon lat time] (i.e., 1 1 12)
                         'snow','snowmelt','humtot','evapot') 
SEC = RCHIDEEsec.catSiteMonthlyVars(options)

plot(SEC$TIME,SEC$humtot_top, t='l')

plot(SEC$TIME,SEC$humtot_top, t='l',xlim=c(1990,2021))
abline(h=22,lty=3)


# Explore soil moistuyre by layers
options=list() # do not touch this line!
options$path2file    = '/home/satellites5/USER/scratch.output/IGCM_OUT/OL2/PROD/ACF/hist0.trw/SRF/Output/MO/';
options$expFilename  = 'hist0.trw';  # it must match the start of the output file
options$years        = c(1901:2020)
options$idxLat       = 1
options$idxLon       = 1
options$ncir         = 3
options$experimentName = 'ACF'
options$variables    = c('SoilMoist') # only variables with dimensions [lon lat solay] (i.e., 1 1 11)
SEC = RCHIDEEsec.catSiteMonthlySoilLayerVars(options)

# plot layer SWC across time
plot(SEC$TIME,SEC$SoilMoist.solay_08, t='l', xlim=c(1990,2020))


# compare soil moisture profiles for Feb 2015 and Feb 2017
monthly_dates = as.POSIXct( seq(ymd("1901-01-01"), ymd("2020-12-01"), by = "month"), tz = "UTC")
idx = which(monthly_dates == as.POSIXct(ymd("2015-02-01"),tz = "UTC"))
idx2 = which(monthly_dates == as.POSIXct(ymd("2017-02-01"),tz = "UTC"))
idx3 = which(monthly_dates == as.POSIXct(ymd("1999-02-01"),tz = "UTC"))
soil_cols = grep("^SoilMoist\\.solay_\\d+$", names(SEC))

plot(as.numeric(SEC[idx,soil_cols]), 1:11, t='l', xlim=c(0,260),ylim=c(11,1))
lines(as.numeric(SEC[idx2,soil_cols]), 1:11, col='blue')
lines(as.numeric(SEC[idx3,soil_cols]), 1:11, col='red')




