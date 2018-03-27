# KM1513_extractData.R

# Purpose: Extract desired CTD and underway/met ship data corresponding to exact
# CTD bottle closure times during cruise KM1513 aboard R/V Kilo Moana,
# 24 Jul - 6 Aug 2015; this was the SCOPE "diel" experiment (HOE-Legacy 2)

# Created 3/26/2018 by J.R.C.

#### load libraries ####

library(repmis) # for sourcing files from web ("repmis" is awesome)

source_data