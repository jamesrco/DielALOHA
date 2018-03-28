# KM1513_extractData.R

# Purpose: Extract desired CTD and underway/met ship data corresponding to exact
# CTD bottle closure times during cruise KM1513 aboard R/V Kilo Moana,
# 24 Jul - 6 Aug 2015; this was the SCOPE "diel" experiment (HOE-Legacy 2)

# Created 3/26/2018 by J.R.C.

# dependencies:
# *** PAR data is pulled from a data object created with KM1513_DielPARcorr.R
# run KM1513_DielPARcorr.R first if KM1513_shipboard_PAR.30smav.W_m2.csv
# doesn't exist yet

#### load libraries ####

library(repmis) # for sourcing files from web ("repmis" is awesome)

#### load, clean up data sources ####

# 1. PAR
# we'll use a 30s-moving average product I produced from the raw met station data
# using script KM1513_DielPARcorr.R, *not* the PAR data available at
# http://hahana.soest.hawaii.edu/hoelegacy/data/hoelegacy2A.lic
KM1513_PAR_W_m2 = read.csv("~/Code/DielPAR/data/met/processed/KM1513_shipboard_PAR.30smav.W_m2.csv",
                           header = TRUE)

# 2. bottle summary data (for position, some other information)
# if source_data() doesn't work because the file has disappeared, a copy was 
# downloaded to data/bottle_data/hoelegacy2A_bottledata.gof
KM1513_bottledata = source_data("ftp://ftp.soest.hawaii.edu/dkarl/scope/water/hoelegacy/hoelegacy2A.gof",
            rdata = FALSE,
            header = FALSE,
            sep = " ")

# assign some column names
colnames(KM1513_bottledata) =
  c("Station_no","Cast_no","Bottle_no","Lat","Lon","CTD_pres","CTD_temp","CTD_salinity",
"CTD_oxy","CTD_chl","Theta","Sigma","Oxygen","DIC","Alk","Phosphate","Nitrite_plus_nitrate",
"Silicate","DOP","DON","DOC","TDP","TDN","LLN","LLP","PC","PN","PP","PSi","Chl_a",
"Pheo","Chlda_a","Chl_c1_c2_c3","

PERID  19 BUT    FUCO  19 HEX PRASINO    VIOL DIADINO   ALLOX  LUTEIN  ZEAXAN   CHL B   A.CAR   B.CAR DV.CHLA MV.CHLA HPLCchl  H.BACT  P.BACT  S.BACT  E.BACT     ATP     CH4     N2O   QUALT1  QUALT2  QUALT3  QUALT4  QUALT5  QUALT6  QUALT7  QUALT8"