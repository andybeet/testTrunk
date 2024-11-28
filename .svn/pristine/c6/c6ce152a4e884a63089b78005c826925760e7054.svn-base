rm(list = ls()) # clear memory

setwd("/Users/ful083/AtlantisRepository/AtlantisTrunk/example/")

library(shinyrAtlantis)

# Step through using shinyrAtlantis
bgm.file <- system.file("extdata", "VMPA_setas.bgm", package = "shinyrAtlantis")
grp.file <- system.file("extdata", "SETasGroupsDem.csv", package = "shinyrAtlantis")
prm.file <- system.file("extdata", "VMPA_setas_biol_fishing_Trunk.prm", package = "shinyrAtlantis")
nc.file <- system.file("extdata", "INIT_VMPA_Jan2015.nc", package = "shinyrAtlantis")

salinity.file    <- "/Users/ful083/AtlantisRepository/AtlantisTrunk/example/inputs/forcisets/SETAS_VMPAsalt.nc"
temperature.file <- "/Users/ful083/AtlantisRepository/AtlantisTrunk/example/inputs/forcisets/SETAS_VMPAtemp.nc"       # this file is not included in the package
cum.depth <- c(0,20,50,100,250,700,2000)  # cumulative water layer depths

#Built in examples
shinyrAtlantis::SpatialDistributionsExample()
#shinyrAtlantis::DisplayParametersExample()
shinyrAtlantis::DisplayInitializationExample()

