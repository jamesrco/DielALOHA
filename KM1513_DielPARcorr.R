# KM1513_DielPARcorr.R

# Purpose: Examine correlations between PAR data and lipid concentrations over course of the SCOPE "diel" experiment (HOE-Legacy 2); conducted aboard cruise KM1513 aboard R/V Kilo Moana, 24 Jul - 6 Aug 2015

# Created 1/26/2017 by J.R.C.

### load PAR data ####

# # if reading in raw individual data files
# # per correspondence with Robert Arko at rvdata.us, most of these are corrected values in raw, arbitrary (voltage) units
# 
# # doi:10.7284/119285, from file KM1513_119285_metstation.tar.gz downloaded from http://get.rvdata.us/cruise/KM1513/fileset/119285 on 26 Jan 2017
# 
# # set data to where we extracted the .tar file downloaded from rvdata.us (KM1513_119285_metstation.tar.gz)
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

### load lipid data ###

# from spreadsheet from KB on 1/26/17; all concentrations are corrected using response factors
# concentrations, by lipid class, in ng per L

LM1513_lipids = read.csv("/Users/jrcollins/Code/DielPAR/data/Lipid concentration_diel (from KB).csv",
                         header = T, skip = 1)

# format timestamp

LM1513_lipids$Timestamp_POSIXct = as.POSIXct(strptime(LM1513_lipids$Timestamp,
                                                      "%m/%d/%y %H:%M"),format='%Y-%m-%d %T', tz = "HST")

