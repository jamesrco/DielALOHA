# KM1513_DielPARcorr.R

# Purpose: Examine correlations between PAR data and lipid concentrations over course of the SCOPE "diel" experiment (HOE-Legacy 2); conducted aboard cruise KM1513 aboard R/V Kilo Moana, 24 Jul - 6 Aug 2015

# Created 1/26/2017 by J.R.C.

### load PAR data ####

# # if reading in raw individual data files
# # per correspondence with Robert Arko at rvdata.us, most of these are corrected values in raw, arbitrary (voltage) units
# 
# # doi:10.7284/119285, from file KM1513_119285_metstation.tar.gz downloaded from http://get.rvdata.us/cruise/KM1513/fileset/119285 on 26 Jan 2017
# 
# # set wd to where we extracted the .tar file downloaded from rvdata.us (KM1513_119285_metstation.tar.gz)
# 
# setwd("/Users/jrcollins/Code/DielPAR/data/met/raw/KM1513/119285/data/")
# 
# KM1513_metfiles = list.files(recursive=TRUE, full.names=FALSE, pattern="_raw") # get list of met data files to be parsed
# 
# # put in order
# 
# KM1513_metfiles = KM1513_metfiles[order(KM1513_metfiles)]
# 
# # iterate through list, parse and collect into single dasta frame
# 
# for (i in 1:length(KM1513_metfiles)) {
#   
#   this.metfile = read.delim(KM1513_metfiles[i], header = FALSE, sep="", colClasses = "character")
#   
#   if (i==1) { # it's the first file
#     
#     allKM1513.met = this.metfile
#     
#   } else { # it's not the first file, so rbind
#     
#     allKM1513.met = rbind(allKM1513.met,this.metfile)
#     
#   }
#   
# }

# if loading in post-processed data (has values in actual units)

# doi:10.7284/123652, from file KM1513_123652_metstation.tar.gz downloaded from http://get.rvdata.us/cruise/KM1513/fileset/123652 on 27 Jan 2017

# set wd to where we extracted the .tar file downloaded from rvdata.us (KM1513_123652_metstation.tar.gz)

setwd("/Users/jrcollins/Code/DielPAR/data/met/processed/KM1513/123652/data/")

allKM1513.met = read.delim("km1513_met", header = FALSE, skip = 0, sep = "") # load in

# assuming this is the processed met data, should have 21 columns, per "/Users/jrcollins/Code/DielPAR/data/met/processed/KM_Formats_of_data.pdf", p. 8)

ncol(allKM1513.met)

# label PAR data (col. 19, in W/m2, per above pdf doc)

colnames(allKM1513.met)[19] = c("PAR_W_m2")

# create, format timestamp

allKM1513.met$Timestamp_POSIXct = as.POSIXct(
  strptime(
    paste(as.Date(as.numeric(allKM1513.met$V2), origin=as.Date("2015-01-01")),
          # *** assumes all of these data were from 2015 ***
          allKM1513.met$V3,
          allKM1513.met$V4,
          allKM1513.met$V5,
          sep = "-"),
    "%Y-%m-%d-%H-%M-%S"),format='%Y-%m-%d %T', tz = "HST")

# plot the PAR data to make sure it makes sense

plot(allKM1513.met$Timestamp_POSIXct, allKM1513.met$PAR_W_m2)

### load lipid data ###

# from spreadsheet from KB on 1/26/17; all concentrations are corrected using response factors
# concentrations, by lipid class, in ng per L

KM1513_lipids = read.csv("/Users/jrcollins/Code/DielPAR/data/Lipid concentration_diel (from KB).csv",
                         header = T, skip = 1)

# a version of the data with "outliers" included

KM1513_lipids.w_out = read.csv("/Users/jrcollins/Code/DielPAR/data/Lipid concentration_diel (from KB)_w_outliers.csv",
                         header = T, skip = 1)

# format timestamp

KM1513_lipids$Timestamp_POSIXct = as.POSIXct(strptime(KM1513_lipids$Timestamp,
                                                      "%m/%d/%y %H:%M"),format='%Y-%m-%d %T', tz = "HST")

### simplifying the PAR data ###

# first, can simplify the dataset a bit (don't need it at 1 s intervals; in addition, there is some obviously spurious data, based on inspection of full plot, above)

# define a simple 10 s moving average function

mav.10 <- function(x,n=10){filter(x,rep(1/n,n), sides=2)}

allKM1513.PAR.smooth = mav.10(allKM1513.met$PAR_W_m2) # apply function

# now, subset to every 10th data point; create matching timestamp series

allKM1513.PAR.red = allKM1513.PAR.smooth[seq(from = 1, to = length(allKM1513.PAR.smooth), by = 10)]
allKM1513.time.red = allKM1513.met$Timestamp_POSIXct[seq(from = 1, to = length(allKM1513.PAR.smooth), by = 10)]

# plot the reduced, simplified series

plot(allKM1513.time.red,allKM1513.PAR.red)

### some daily PAR integrals ###

library(caTools)

# generate vector of times in seconds, with t = 0 being timepoint at beginning of simplified series

allKM1513.red.times_s = as.numeric(rev(allKM1513.time.red[length(allKM1513.time.red)]-
                                         allKM1513.time.red))

# preallocate vector for results

KM1513.metdates = unique(strftime(allKM1513.time.red, format = "%D")) # list of dates in the time series
KM1513_PARint_J_d = vector(mode = "numeric", length = length(KM1513.metdates))

# integrals by day

for (i in 1:length(KM1513_PARint_J_d)) {
  
  # get index to data for date i
  ind_todays.data = strftime(allKM1513.time.red, format = "%D")==KM1513.metdates[i]
  
  KM1513_PARint_J_d[i] =
  caTools::trapz(allKM1513.red.times_s[ind_todays.data], # there's at least one other pkg with a trapz function
                 allKM1513.PAR.red[ind_todays.data])
  
}

# take a look at the integrals

plot(as.POSIXct(KM1513.metdates, format = "%m/%d/%y"),KM1513_PARint_J_d)

### some daily net lipid fluxes, by class ###

# first, take a look at TAG data and chl data to see if there's any sort of trend (will have to de-trend before calculating integrals, if so)

plot(KM1513_lipids$Timestamp_POSIXct, KM1513_lipids$TAG, type = "p")
plot(KM1513_lipids$Timestamp_POSIXct, KM1513_lipids$Chl, type = "p")

# no apparent trend, so let's calculate some daily integrals

# vector of times, in seconds, first measurement considered t = 0
KM1513_lipids.times_s = as.numeric(rev(KM1513_lipids$Timestamp_POSIXct[
  length(KM1513_lipids$Timestamp_POSIXct)]-
    KM1513_lipids$Timestamp_POSIXct))

# preallocate data frame for results

KM1513.lipiddates = unique(strftime(KM1513_lipids$Timestamp_POSIXct, format = "%D")) # list of dates in the time series
KM1513_intlipids.ng_L_d = as.data.frame(matrix(nrow = length(KM1513.lipiddates),
                                 ncol = ncol(KM1513_lipids) - 4 + 1))
colnames(KM1513_intlipids.ng_L_d) =
  c("Date",colnames(KM1513_lipids)[4:(ncol(KM1513_lipids)-1)])
KM1513_intlipids.ng_L_d$Date = as.POSIXct(KM1513.lipiddates, format = "%m/%d/%y")

# integrals by day, by lipid

for (i in 1:nrow(KM1513_intlipids.ng_L_d)) {
  
  for (j in 2:ncol(KM1513_intlipids.ng_L_d)) {
    
    # get time index to data
    ind_todays.data = strftime(KM1513_lipids$Timestamp_POSIXct, format = "%D")==KM1513.lipiddates[i]
    
    KM1513_intlipids.ng_L_d[i,j] =
      caTools::trapz(KM1513_lipids.times_s[ind_todays.data], # there's at least one other pkg with a trapz function
                     KM1513_lipids[ind_todays.data,colnames(KM1513_lipids)==colnames(KM1513_intlipids.ng_L_d)[j]])
    
  }
  
}

### exploring some correlations ###

# date for date

# a vector of PAR integrals to match dates in KM1513.lipiddates
KM1513_PARint_J_d.match = KM1513_PARint_J_d[KM1513.metdates %in% KM1513.lipiddates]

# scatterplot matrix
pairs(~KM1513_PARint_J_d.match+KM1513_intlipids.ng_L_d$DGCC+
        KM1513_intlipids.ng_L_d$DGTS.DGTS+
        KM1513_intlipids.ng_L_d$DGDG+
        KM1513_intlipids.ng_L_d$MGDG+
        KM1513_intlipids.ng_L_d$PC+
        KM1513_intlipids.ng_L_d$PG+
        KM1513_intlipids.ng_L_d$TAG+
        KM1513_intlipids.ng_L_d$SQDG+
        KM1513_intlipids.ng_L_d$PQ+
        KM1513_intlipids.ng_L_d$UQ+
        KM1513_intlipids.ng_L_d$Chl)

# lagged by one day

lagtime = 1 # define our lag time

# a vector of PAR integrals to match dates in KM1513.lipiddates
KM1513_PARint_J_d.match.lagged = KM1513_PARint_J_d[KM1513.metdates %in% 
                                                     strftime(
                                                     as.POSIXct(KM1513.lipiddates, format = "%m/%d/%y")-
                                                     24*60*60*lagtime,
                                                     format = "%D")]

# scatterplot matrix
pairs(~KM1513_PARint_J_d.match.lagged+KM1513_intlipids.ng_L_d$DGCC+
        KM1513_intlipids.ng_L_d$DGTS.DGTS+
        KM1513_intlipids.ng_L_d$DGDG+
        KM1513_intlipids.ng_L_d$MGDG+
        KM1513_intlipids.ng_L_d$PC+
        KM1513_intlipids.ng_L_d$PG+
        KM1513_intlipids.ng_L_d$TAG+
        KM1513_intlipids.ng_L_d$SQDG+
        KM1513_intlipids.ng_L_d$PQ+
        KM1513_intlipids.ng_L_d$UQ+
        KM1513_intlipids.ng_L_d$Chl)
