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

# convert date to an R time object
KM1513_PAR_W_m2$Timestamp_PAR_HST =
strptime(as.character(KM1513_PAR_W_m2$Timestamp_PAR_HST),
         "%Y-%m-%dT%H:%M:%S%z", tz="UTC")

# 2. The bottle summary data (for position, some other information)

# if source_data() doesn't work because the file has disappeared, a copy was 
# downloaded to data/bottle_data/hoelegacy2A_bottledata.gof
KM1513_bottledata = source_data("ftp://ftp.soest.hawaii.edu/dkarl/scope/water/hoelegacy/hoelegacy2A.gof",
            rdata = FALSE,
            header = FALSE,
            sep = " ")

# assign some column names
colnames(KM1513_bottledata) =
  c("Station_no","Cast_no","Bottle_no","Lat","Lon","CTD_pres_dbar","CTD_temp_degC","CTD_salinity",
"CTD_oxy_umol_kg","CTD_chl_ug_L","Theta","Sigma","Oxygen","DIC","Alk","Phosphate","Nitrite_plus_nitrate",
"Silicate","DOP","DON","DOC","TDP","TDN","LLN","LLP","PC","PN","PP","PSi","Chl_a",
"Pheo","Chlda_a","Chl_plus","Peridino","Nineteen_but","Fuco","Nineteen_hex","Prasino",
"Viol","Diadino","Allox","Lutein","Zeaxan","Chl_b","Alpha_car","Beta_aar","DV_Chla",
"MV_Chla","HPLC_chl","HetBact","ProBact","SynBact","Euks","ATP","CH4","N2O","QUALT1",
"QUALT2","QUALT3","QUALT4","QUALT5","QUALT6","QUALT7","QUALT8")

# 3. File containing the diel RNA sample metadata

# Exact bottle closure times were assembled from CTD cast tearsheets by @sachacoesel

KM1513_RNA_sample_metadata = read.csv("~/Code/DielPAR/metadata/Diel_RNA_sample_metadata.csv",
                                      header = TRUE,
                                      skip = 3)

#### perform matching and data extraction ####

# From 19 March 2018 email from Bryn Durham, bpdurham@uw.edu, we know that all the RNA
# samples were taken from CTD rosette bottles 8-13

RNAbottles = c(8:13)

# will use this information, along with exact bottle closure times and the fields
# Station_no and Cast_no to extract the most closely matched data

# specify the data fields to be matched & retrieved from various sources

fields.getfromKM1513_bottledata =
  c("Lat","Lon","CTD_pres_dbar","CTD_temp_degC","CTD_salinity","CTD_oxy_umol_kg",
    "CTD_chl_ug_L")

fields.getfromKM1513_PAR_W_m2 =
  c("KM1513_shipboard_PAR_30s.mav_W_m2")

numExtraDielfields = length(fields.getfromKM1513_bottledata) + length(fields.getfromKM1513_PAR_W_m2)
  
# preallocate a data frame based on KM1513_RNA_sample_metadata

KM1513_RNA_sample_matchedDielData =
cbind(KM1513_RNA_sample_metadata,
      as.data.frame((matrix(data = NA, nrow = nrow(KM1513_RNA_sample_metadata),
                            ncol = numExtraDielfields))))
colnames(KM1513_RNA_sample_matchedDielData)[(ncol(KM1513_RNA_sample_matchedDielData)-(numExtraDielfields-1)):
                                              ncol(KM1513_RNA_sample_matchedDielData)] =
  c(fields.getfromKM1513_bottledata,fields.getfromKM1513_PAR_W_m2)

# do the actual matching & extraction

for (i in 1:nrow(KM1513_RNA_sample_matchedDielData)) {
  
  this.Station_no = KM1513_RNA_sample_matchedDielData$Station_no[i]
  this.Cast_no = KM1513_RNA_sample_matchedDielData$Cast_no[i]
  this.Bottletime = 
    strptime(as.character(KM1513_RNA_sample_matchedDielData$Timestamp_bottle_closure_GMT_ISO8601[i]),
             "%Y-%m-%dT%H:%M:%S%z", tz="UTC")
  
  # first, match and get data contained KM1513_bottledata based on station and
  # cast number, using an average of the values collected for bottles 8-13 if
  # there is any difference
  
  for (j in 1:length(fields.getfromKM1513_bottledata)) {
    
    matched.CTDdata = 
      KM1513_bottledata[KM1513_bottledata$Station_no == this.Station_no &
                          KM1513_bottledata$Cast_no == this.Cast_no & 
                          KM1513_bottledata$Bottle_no %in% RNAbottles,
                        fields.getfromKM1513_bottledata[j]]
    
    # take the mean value, if there's any variance
    
    if (length(unique(matched.CTDdata)) > 1) {
      
      KM1513_RNA_sample_matchedDielData[i,ncol(KM1513_RNA_sample_metadata)+j] =
        mean(matched.CTDdata, na.rm = T)
      
      # provide some feedback so we can decide whether the variance in any
      # of the parameters should concern us
      
      cat(this.Station_no,this.Cast_no,fields.getfromKM1513_bottledata[j],
            mean(matched.CTDdata, na.rm = T),sd(matched.CTDdata, na.rm = T),"\n")
      
    } else if (length(unique(matched.CTDdata)) == 1) {
      
      KM1513_RNA_sample_matchedDielData[i,ncol(KM1513_RNA_sample_metadata)+j] =
        matched.CTDdata[1]
      
    }
    
  }
  
  # now, match and extract the PAR data
  
  matched.PARdata = KM1513_PAR_W_m2$KM1513_shipboard_PAR_30s.mav_W_m2[which.min(
    abs(KM1513_PAR_W_m2$Timestamp_PAR_HST-this.Bottletime))]
  
  KM1513_RNA_sample_matchedDielData[i,ncol(KM1513_RNA_sample_metadata)+length(fields.getfromKM1513_bottledata)+1] =
    matched.PARdata
  
}

# finally, export the matched data to a CSV file

write.csv(KM1513_RNA_sample_matchedDielData,
          file = "~/Code/DielPAR/data/data_products/KM1513_RNA_sample_matchedDielData.csv",
          row.names = FALSE)