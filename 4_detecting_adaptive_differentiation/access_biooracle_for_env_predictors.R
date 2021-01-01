######################################################################
# Accessing env data for my sites from Bio-Oracle
# Rerunning bay env to include same variables as Xuereb et al (and more if interested)
# 20200226 NL

# Amanda's team picked variables at surface and at bottom, to try and get that range
# Specifically, ones that I don't really cover yet include:
# - bottom temp (min, max, mean)
# - surface current velocity
# - bottom dissolved oxygen
# - bottom chlorophyll

# Bottom temp was sig for them, so definitely should include

# http://www.bio-oracle.org/

######################################################################

# # install packages
# install.packages("sdmpredictors")
# install.packages("leaflet")

# call packages
library(sdmpredictors)
library(leaflet)
#library(tidyverse)


######################################################################
### Extracting environmental information for a set of sites

# Learn about layers
layerinfo <- list_layers()

# List layers avaialble in Bio-ORACLE v2
layers.bio2 <- list_layers( datasets="Bio-ORACLE" )

# check out this df to look at possible environmental variables
layers.bio2

# Download environmental data layers of interest
environmental_vars <- load_layers( layercodes = c( # MEAN VAR AT SURFACE + DEPTH
                                                  "BO2_salinitymean_ss", # mean salinity at surface
                                                  "BO2_ppmean_ss", # mean primary production at surface
                                                  "BO2_nitratemean_ss", # mean nitrate at surface
                                                  "BO2_phosphatemean_ss", # mean nitrate at surface
                                                  "BO2_dissoxmean_ss", # mean dissox at surface

                                                  "BO2_salinitymean_bdmean", # mean salinity at mean bottom
                                                  "BO2_ppmean_bdmean", # primary production at mean bottom
                                                  "BO2_nitratemean_bdmean", # mean nitrate at mean depth
                                                  "BO2_phosphatemean_bdmean", # mean phosphate at mean depth
                                                  "BO2_dissoxmean_bdmean", # mean dissox at max depth
                                                  
                                                  "BO_bathymean", # bathymetry
                                                  "BO_ph", # surface pH 
                                                  "BO_calcite", # mean calcite

                                                  # CURRENT VELOCITY AT SURFACE AND MAX DEPTH
                                                  "BO2_curvelmax_ss", # surface max current velocity
                                                  "BO2_curvelmin_ss", # surface min current velocity
                                                  "BO2_curvelmean_ss", # surface mean current velocity
                                                  "BO2_curvelrange_ss", # surface range current velocity
                                                  "BO2_curvelmax_bdmean", # mean bottom max current velocity
                                                  "BO2_curvelmin_bdmean", # mean bottom min current velocity
                                                  "BO2_curvelrange_bdmean", # mean bottom range current velocity
                                                  "BO2_curvelmean_bdmean", # mean bottom mean current velocity
                          
                                                  # TEMP AT SURFACE AND MAX DEPTH
                                                  "BO2_tempmax_ss", # max SST
                                                  "BO2_tempmean_ss", # mean SST
                                                  "BO2_tempmin_ss", # min SST
                                                  "BO2_temprange_ss", # range SST
                                                  "BO2_tempmax_bdmean", # maximum temp at mean bottom
                                                  "BO2_tempmean_bdmean", # mean temp at mean bottom
                                                  "BO2_tempmin_bdmean", # min temp at mean bottom
                                                  "BO2_temprange_bdmean" # range temp at mean bottom
                                                  ), 
                                   equalarea=FALSE, 
                                   rasterstack=TRUE)


# set wd
setwd("E:/dD_cragig_mox_20191205/filtering/analyses/gene_env_assn")

# Read in a csv file with site name, long, and lat, like this:
# Name	         Lon	      Lat
# ChiniakBayAK	-152.35652	57.712427
# SewardAK	    -149.3997222	60.10027778
# YakutatBayAK 	-139.8463889	59.7158333
my.sites <- read.csv("CG_PC_site_locs.csv", stringsAsFactors = FALSE)
my.sites.locs <- my.sites[,2:3] # get just longs and lats

# Visualise sites of interest in google maps, make sure placing where we think we are
check_map <- leaflet()
check_map <- addTiles(check_map)
check_map <- addMarkers(check_map, lng=my.sites$Lon, lat=my.sites$Lat, popup=my.sites$Name)
check_map

# Extract environmental values from layers, stick into a dataframe
my.sites.env.data <- data.frame(Name=my.sites$Name, 
                                   Var=extract(environmental_vars, my.sites.locs))
my.sites.env.data

# write it to file
write.table(my.sites.env.data, file = "my_sites_env_data_xuereb2.tsv", quote=FALSE, sep = " ")

######################################################################

# 3 of my sites weren't in env data sets, so trying points nearby (hopefully okay since resolution is terrible anyway)
my.tricky.sites <- read_csv("tricky_locs.csv")
my.tricky.sites.locs <- my.tricky.sites[,2:3]

# Visualise sites of interest in google maps, make sure placing where we think we are
check_map <- leaflet()
check_map <- addTiles(check_map)
check_map <- addMarkers(check_map, lng=my.tricky.sites$Lon, lat=my.tricky.sites$Lat, popup=my.tricky.sites$Name)
check_map

# Extract environmental values from layers, stick into a dataframe
my.tricky.sites.env.data <- data.frame(Name=my.tricky.sites$Name, 
                                Var=extract(environmental_vars, my.tricky.sites.locs))
my.tricky.sites.env.data

# What I learned from ^ this exercise is that there is data from bio-oracle to places within a reasonable distance
# to all of my sites EXCEPT Eld Inlet, for which the closest I could get was Commencement Bay. Those
# are probably still more similar than Sekiu vs. south sound, etc., but potentially different too
# I'll note it, and still run the gene by env association test. In the case that there are 
# significant results, I'll decide later whether to include this population at all.

######################################################################
### Generate a map of a specific geographic space, of a layer

# Easy download of raster file (Maximum Temperature at the sea bottom)
ph <- load_layers("BO_ph")

# Crop raster to fit the North Atlantic
ne.pacific.ext <- extent(-160,-117,30,60) # get extent of desired range 
ph.crop <- crop(ph, ne.pacific.ext) # crop raster

# Generate a nice color ramp and plot the map
my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))
plot(ph.crop,col=my.colors(1000),axes=FALSE, box=FALSE)
title(cex.sub = 1.25, sub = "Surface pH")
