# R code for EcMF SDM analysis
# Clara Qin, 2024

# E3.2-rgee-range-maps.R: Use rgee package to calculate species range maps.
# This would take too much computation and memory on the RStudio Server.

# Load and initialize rgee

library(rgee)
ee_Initialize()
task_list <- list()

# Parameters to toggle

nMcmcSubsample <- 10 # number of MCMC samples to include in the subsample; each will require its own prediction raster
thinningGridKm <- 500 # feasible: 500
thinningGridKmAMF <- 200
krigingMaxDistKm <- 5000 # feasible: 5000
krigingResKm <- 50 # feasible: 50, 100
rangeMinimumP <- 0.01 # range defined by mean P greater than this
rangeMinimumPConservative <- 0.01 # range defined by 95% confidence that P is greater than this

# Pointers to raster data layers

image <- ee$Image("projects/spun-geospatial/assets/data-archive/full_stack_v15")
# image_n <- ee$Image("projects/spun-geospatial/assets/nitrogen_5-15cm_mean_1000_reproj")
image_soc <- ee$Image("projects/spun-geospatial/assets/SoilGrids/SoilGrids_SOC_0-15cm_corrected")
image_ph <- ee$Image("projects/spun-geospatial/assets/SoilGrids/SoilGrids_pH_0-15cm_corrected")
image_nitro <- ee$Image("projects/spun-geospatial/assets/SoilGrids/SoilGrids_nitrogen_0-15cm_corrected")

i_mat <- image$select("CHELSA-mean-annual-air-temp")$multiply(0.1)$add(-273.15)
i_mat2 <- i_mat$pow(2)
i_trng <- image$select("CHELSA-annual-range-of-air-temp")$multiply(0.1)
i_map <- image$select("CHELSA-annual-precipitation-amount")$multiply(0.1)
i_map2 <- i_map$pow(2)
i_psea <- image$select("CHELSA-precipitation-seasonality")$multiply(0.1)
i_npp <- image$select("MODIS-NPP-2022")
i_npp2 <- i_npp$pow(2)
i_abio <- image$select("harmonized-aboveground-biomass")
i_abio2 <- i_abio$pow(2)
i_soc <- image_soc$multiply(0.1)
i_soc2 <- i_soc$pow(2)
i_phos <- image$select("phosphorus")
i_phos2 <- i_phos$pow(2)
i_ph <- image_ph$multiply(0.1)
i_ph2 <- i_ph$pow(2)
i_nitro <- image_nitro$multiply(0.01)
i_nitro2 <- i_nitro$pow(2)
i_logits <- ee$Image$constant(11.4)
i_logsamples <- ee$Image$constant(4.4)
i_logits <- i_logits$setDefaultProjection(i_mat$projection())
# i_logsamples <- i_logsamples$setDefaultProjection(i_mat$projection())

# Define ROI for image exports 

ee_ROI <- ee$Geometry$BBox(-180, -60, 180, 84)

# From ChatGPT:
# Define the standard normal CDF (probit CDF) function for an image
standard_normal_cdf_image <- function(x) {
  # Constants for the approximation
  a1 <- 0.319381530
  a2 <- -0.356563782
  a3 <- 1.781477937
  a4 <- -1.821255978
  a5 <- 1.330274429
  pi <- pi
  
  # Compute t, a function of x
  t <- ee$Image(1)$divide(ee$Image(1)$add(ee$Image(0.2316419)$multiply(x$abs())))
  
  # Polynomial approximation
  poly <- ee$Image(a1)$multiply(t)$add(ee$Image(a2)$multiply(t$pow(2)))$
    add(ee$Image(a3)$multiply(t$pow(3)))$
    add(ee$Image(a4)$multiply(t$pow(4)))$
    add(ee$Image(a5)$multiply(t$pow(5)))
  
  # Calculate the main CDF approximation
  approx <- ee$Image(1)$subtract(
    ee$Image(1)$divide(ee$Image(2 * pi)$sqrt())$multiply(x$pow(2)$multiply(-0.5)$exp())$multiply(poly)
  )
  
  # Adjust for negative x using symmetry of the normal distribution
  x$gte(0)$multiply(approx)$add(x$lt(0)$multiply(ee$Image(1)$subtract(approx)))
}




# EcMF --------------------------------------------------------------------

# Model files

model_files <- list.files("./models/singletaxon-EcMF-GFv5", pattern="mpost_.*_2025-02-2(4|5).rds", full.names=TRUE)
subset_OTUs <- str_extract(model_files, "OTU[0-9]+")
length(model_files)

# Existing assets
assets <- ee_manage_assetlist("projects/spun-geospatial/assets/EcMF-ranges-GFv5")
image_ids_ecmf <- assets$ID
asset_OTUs <- str_extract(image_ids_ecmf, "OTU[0-9]+")
head(asset_OTUs)
length(asset_OTUs)

# Remove the models for which we already have assets (range maps)
model_files <- model_files[which(!subset_OTUs %in% asset_OTUs)]
subset_OTUs <- subset_OTUs[which(!subset_OTUs %in% asset_OTUs)]
length(model_files)
head(subset_OTUs)


## Spatial thinning to prep for kriging -------------------------------------

# Grid for spatial thinning for kriging

extent_global <- ext(-17367530, 17367530, -7314540, 7314540)  # Approx. global extent in meters for Equal Earth
resolution <- thinningGridKm * 1e+3 # kriging points will be separated by __ km grid
r_grid <- rast(extent_global, resolution = resolution, crs = "EPSG:8857")  # EPSG 8857 is Equal Earth

# Conduct spatial thinning of eta to reduce computation time of kriging
# Algorithm inspired by dismo::gridSample
mpost <- readRDS(model_files[1])
pts_sf <- st_as_sf(mpost$rL[[1]]$s)
pts_sf <- st_transform(pts_sf, "EPSG:8857")
pts <- st_coordinates(pts_sf)
cell <- cellFromXY(r_grid, pts)
keep_spthin <- !duplicated(cell)

## Loop over OTUs ----------------------------------------------------------

# for (i in seq_along(model_files)) {
for (i in 1:200) {
  message("Started: OTU ", i, " out of ", length(model_files), " (", subset_OTUs[i], ")")
  
  ## Extract parameter estimates ---------------------------------------------
  
  mpost <- readRDS(model_files[i])
  
  thin <- floor(mpost$samples * length(mpost$postList) / nMcmcSubsample)
  post <- poolMcmcChains(mpost$postList, thin=thin)
  postBeta <- lapply(post, function(c) c$Beta)
  postEta <- lapply(post, function(c) c$Eta[[1]])
  postLambda <- lapply(post, function(c) c$Lambda[[1]])
  postAlpha <- lapply(post, function(c) c$Alpha[[1]])
  
  ## Run the following for each MCMC sample ----------------------------------
  
  for(j in 1:length(post)) {
    
    beta <- postBeta[[j]]
    eta <- postEta[[j]][,1]
    lambda <- postLambda[[j]][1,1]
    alpha <- mpost$rL[[1]]$alphapw[postAlpha[[j]],1]
    
    ### Calculate fixed effects -------------------------------------------------
    
    calculate_fixed_effects <- function(beta) {
      L_fe <- ee$Image(beta[1])
      L_fe <- i_mat$multiply(beta[2])$add(L_fe)
      L_fe <- i_mat2$multiply(beta[3])$add(L_fe)
      L_fe <- i_trng$multiply(beta[4])$add(L_fe)
      L_fe <- i_map$multiply(beta[5])$add(L_fe)
      L_fe <- i_map2$multiply(beta[6])$add(L_fe)
      L_fe <- i_psea$multiply(beta[7])$add(L_fe)
      L_fe <- i_npp$multiply(beta[8])$add(L_fe)
      L_fe <- i_npp2$multiply(beta[9])$add(L_fe)
      L_fe <- i_abio$multiply(beta[10])$add(L_fe)
      L_fe <- i_abio2$multiply(beta[11])$add(L_fe)
      L_fe <- i_soc$multiply(beta[12])$add(L_fe)
      L_fe <- i_soc2$multiply(beta[13])$add(L_fe)
      L_fe <- i_phos$multiply(beta[14])$add(L_fe)
      L_fe <- i_phos2$multiply(beta[15])$add(L_fe)
      L_fe <- i_ph$multiply(beta[16])$add(L_fe)
      L_fe <- i_ph2$multiply(beta[17])$add(L_fe)
      L_fe <- i_nitro$multiply(beta[18])$add(L_fe)
      L_fe <- i_nitro2$multiply(beta[19])$add(L_fe)
      L_fe <- i_logits$multiply(beta[20])$add(L_fe)
      L_fe
    }
    
    L_fe <- calculate_fixed_effects(beta)
    
    ### Calculate random effects ------------------------------------------------
    
    if(alpha == 0) {
      
      L_re <- ee$Image$constant(0)
      
    } else {
      eta_sf <- st_as_sf(mpost$rL[[1]]$s[keep_spthin])
      eta_sf$eta <- eta[keep_spthin]
      eta_ee <- sf_as_ee(eta_sf)
      
      # data(meuse)
      # coordinates(meuse) <- ~ x+y
      # meuse <- meuse["zinc"]
      # variogram <- autofitVariogram(zinc~1,meuse, model = "Sph")
      
      # Kriging
      i_eta <- eta_ee$kriging(
        shape = "exponential",
        propertyName = "eta",
        range = alpha * 1e+3,
        # sill = 1,
        sill = var(eta),
        # nugget = 0.05,
        nugget = 0,
        maxDistance = krigingMaxDistKm * 1e+3
      )$reproject(
        crs = "EPSG:4326",
        scale = krigingResKm * 1e+3
      )
      
      # Multiply by lambda to get random effects component of the linear predictor
      L_re <- i_eta$multiply(lambda)
    }
    
    ### Calculate linear predictor ----------------------------------------------
    
    # Linear predictor
    L <- L_fe$add(L_re)
    
    # Append to ImageCollection
    new_collection <- ee$ImageCollection(L)
    if(j==1) {
      L_collection <- new_collection
    } else {
      L_collection <- L_collection$merge(new_collection)
    }
    
    # ### Optionally visualize components of prediction ---------------------------
    # 
    # p <- standard_normal_cdf_image(L)
    # p_without_randomeffects <- standard_normal_cdf_image(L_fe)
    # 
    # # Calculate species range, defined as > __% probability of presence
    # sp_range <- p$gt(rangeMinimumP)
    # sp_range_without_randomeffects <- p_without_randomeffects$gt(rangeMinimumP)
    # 
    # # Get points denoting presences
    # presences_sf <- st_as_sf(mpost$rL[[1]]$s[mpost$Y[,1] == 1])
    # presences_ee <- sf_as_ee(presences_sf)
    # 
    # viridisPal <- c(
    #   '0D0887', '5B02A3',
    #   '9A179B', 'CB4678',
    #   'EB7852', 'FBB32F',
    #   'F0F921'
    # )
    # 
    # # Visualize
    # # eta_lims <- quantile(eta, probs=c(0.01, 0.99))
    # # L_re_lims <- sort(eta_lims * lambda)
    # mapcenter00 <- ee$Geometry$Point(c(0.01, 0.01))
    # Map$centerObject(mapcenter00, zoom=2)
    # Map$addLayer(L_re, name="L_re", shown=FALSE, list(palette = viridisPal)) +
    #   Map$addLayer(L_fe, name="L_fe", shown=FALSE, list(palette = viridisPal)) +
    #   Map$addLayer(eta_ee, name="kriging points", shown=FALSE, list(color="FFFFFF")) +
    #   Map$addLayer(p_without_randomeffects, name="P w/out REs", shown=FALSE, list(min=0, max=1, palette = viridisPal)) +
    #   Map$addLayer(p, name="Probability", list(min=0, max=1, palette=viridisPal)) +
    #   Map$addLayer(sp_range_without_randomeffects, name=paste0("Range w/out RE (P > ", round(rangeMinimumP, 2) * 100, "%)"), shown=FALSE, list(palette=c("000000", "FF0000"))) +
    #   Map$addLayer(sp_range, name=paste0("Range (P > ", round(rangeMinimumP, 2) * 100, "%)"), shown=FALSE, list(palette=c("000000", "FF0000"))) +
    #   Map$addLayer(presences_ee, name="Presences", list(color="FF00FF"))
  }
  
  ## Take mean probability across MCMC samples -------------------------------
  
  L_mean <- L_collection$mean()
  # L_sd <- L_collection$reduce(ee$Reducer$stdDev())
  # L_lower <- L_mean$subtract(L_sd$multiply(1.65)) # 5th %ile
  
  ## Apply probit transform to linear predictor -------------------------------
  
  p_mean <- standard_normal_cdf_image(L_mean)
  # p_lower <- standard_normal_cdf_image(L_lower)
  
  ## Calculate species range -------------------------------------------------

  # Calculate species range, defined as > __% probability of presence
  sp_range <- p_mean$gt(rangeMinimumP)
  # sp_range_conservative <- p_lower$gt(rangeMinimumPConservative)
  
  ## Optionally visualize the mean prediction --------------------------------

  # Get points denoting presences
  # presences_sf <- st_as_sf(mpost$rL[[1]]$s[mpost$Y[,1] == 1])
  # presences_ee <- sf_as_ee(presences_sf)
  # 
  # # Visualize
  # mapcenter00 <- ee$Geometry$Point(c(0.01, 0.01))
  # Map$centerObject(mapcenter00, zoom=2)
  # Map$addLayer(p_mean, name="Probability", list(min=0, max=1, palette=c('0D0887', '5B02A3', '9A179B', 'CB4678', 'EB7852', 'FBB32F', 'F0F921'))) +
  #   # Map$addLayer(p_sd, name="Uncertainty (SD of prob.)", list(palette=c('0D0887', '5B02A3', '9A179B', 'CB4678', 'EB7852', 'FBB32F', 'F0F921'))) +
  #   # Map$addLayer(sp_range_conservative, name=paste0("Conservative range (>95% chance that P > ", round(rangeMinimumPConservative, 2) * 100, "%)"), shown=FALSE, list(palette=c("000000", "FF0000"))) +
  #   Map$addLayer(sp_range, name=paste0("Range (P > ", round(rangeMinimumP, 2) * 100, "%)"), shown=FALSE, list(palette=c("000000", "FF0000"))) +
  #   Map$addLayer(eta_ee, name="kriging points", shown=FALSE, list(color="FFFFFF")) +
  #   Map$addLayer(presences_ee, name="Presences", list(color="FF00FF"))
  
  ## Add taxon name as a property --------------------------------------------
  
  sp_range <- sp_range$set('taxon', subset_OTUs[i])
  # sp_range_conservative <- sp_range_conservative$set('taxon', subset_OTUs[i])
  
  ## Export image as GEE asset -----------------------------------------------
  
  # asset_name_distr <- paste0("Non-Antarc distr of ", subset_OTUs[i])
  # asset_id_distr <- paste0("projects/spun-geospatial/assets/EcMF-distributions/prob_nonantarctic_mean_", subset_OTUs[i])
  # asset_name_range <- paste0("Non-Antarc range of ", subset_OTUs[i])
  # asset_id_range <- paste0("projects/spun-geospatial/assets/EcMF-ranges-GFv5/range_nonantarctic_95gt", str_pad(round(rangeMinimumPConservative*100), 2, pad="0"), "_", subset_OTUs[i])
  asset_name_range0 <- paste0("Range map of ", subset_OTUs[i])
  asset_id_range0 <- paste0("projects/spun-geospatial/assets/EcMF-ranges-GFv5/range_nonantarctic_meangt", str_pad(round(rangeMinimumP*100), 2, pad="0"), "_", subset_OTUs[i])
  
  # task_distr <- ee_image_to_asset(
  #   p_mean,
  #   description = asset_name_distr,
  #   assetId = asset_id_distr,
  #   region = ee_ROI,
  #   scale = 1000,
  #   overwrite = TRUE,
  #   maxPixels = 1e+13
  # )
  # task_range <- ee_image_to_asset(
  #   sp_range_conservative,
  #   description = asset_name_range,
  #   assetId = asset_id_range,
  #   region = ee_ROI,
  #   scale = 1000,
  #   overwrite = TRUE,
  #   maxPixels = 1e+13
  # )
  task_range0 <- ee_image_to_asset(
    sp_range,
    description = asset_name_range0,
    assetId = asset_id_range0,
    region = ee_ROI,
    scale = 1000,
    overwrite = TRUE,
    maxPixels = 1e+13
  )
  
  # task_distr$start()
  # task_range$start()
  task_range0$start()
  
  # Append the task to the list
  # task_list <- c(task_list, task_distr)
  # task_list <- c(task_list, task_range)
  task_list <- c(task_list, task_range0)
}
