# rchidee
RCHIDEE function library is a collection of R functions to easily extract, plot and postprocess 
   variables from stomate and sechiba history files of ORCHIDEE4. Currently, it works only 
   for a SINGLE PIXEL and using a nonleap (365-day) calendar for daily output.

 Version 1.3 -- Added the following functions:
   * RCHIDEEsto4dim.catSiteMonthlyVars        # tree-ring width 
   * RCHIDEEsto4dim.catSiteAnnualNCIRVars     # tree-ring width; replaced deprecated function RCHIDEEsto.catSiteAnnualNCIRVars 

 Version 1.2 -- Added the following functions:
   * RCHIDEE.extractERA5halfdForcing4Point
   * RCHIDEE.getERA5halfdLandPoint

 Functions available:
   GENERAL
   * RCHIDEE.listVariables          -- Print to the console the name and dimensions of all variables in a history file.
   * RCHIDEE.makeWeibull3parDist    -- Write to the console a three-parameter Weibull CIRC_CLASS_DIST for ORCHIDEE-CAN run.def.
   * RCHIDEE.day2decyy
   * RCHIDEE.month2decyy
   * RCHIDEE.day2decyy.cal366d
   * RCHIDEE.extractCRUNCEPhalfdForcing4Point
   * RCHIDEE.getCRUNCEPhalfdLandPoint
   * RCHIDEE.extractERA5halfdForcing4Point
   * RCHIDEE.getERA5halfdLandPoint

   STOMATE
   * RCHIDEEsto.catSiteDailyCCBA    -- CONCATENATE daily basal area, basal area increment and number of trees for each circ class.
   * RCHIDEEsto.readSiteDailyCCBA   -- READ daily basal area, basal area increment and number of trees for each circ class from an already concatenated file.
   * RCHIDEEsto.retrieveTRW         -- Compute annual tree-ring width time series for each circ class from annual cincumference increments.
   * RCHIDEEsto.catSiteDailyVars.cal365  -- CONCATENATE daily time series of a list of variables for a run with a fixed 365-day calendar 
   * RCHIDEEsto.readSiteDailyVars.cal365 -- READ daily time series of a list of variables for a run with a fixed 365-day calendar 
   * RCHIDEEsto.catSiteDailyVars.cal366  -- CONCATENATE daily time series of a list of variables for a run with a fixed 366-day calendar 
   * RCHIDEEsto.catSiteMonthlyVars  -- CONCATENATE monthly time series of a list of variables
   * RCHIDEEsto.catSiteAnnualVars   -- CONCATENATE annual time series of a list of variables
   * RCHIDEEsto.catSiteAnnualNCIRVars -- CONCATENATE annual time series for variables across NCIR classes (otherwise it is tedious to do with RCHIDEEsto.catSiteAnnualVars for large NCIR)
   * RCHIDEEsto4dim.catSiteMonthlyVars -- CONCATENATE monthly time series of a list of variables by ncric classee from the 4dim stomate history file
   * RCHIDEEsto.catSiteAnnualVerticalCanopyVars -- CONCATENATE annual variables for each vertical canopy level from yearly stomate history files 

   SECHIBA
   * RCHIDEEsec.catSiteMonthlyVars  -- CONCATENATE monthly time series of a list of variables from a sechiba file

 To do:
   * improve general input checking and error management
   * function to concatenate any custom variable from the list
   * function to plot a set of variables from a single experiment
   * function to plot a comparison of several experiments for a set of variables

 Changes
  * v1p1: RCHIDEE.listVariables: adding case when dimension of variable is 1
  * v1p3: 8 dec 2024, added functions to read tree-ring width from new 4dim stomate output file in ORCHIDEE4
    
 @   by Jonathan Barichivich (jonathan.barichivich@lsce.ips.fr), LSCE.
 @   February 2015, revised and published 9 Dec 2024, Paris

