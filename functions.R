load(file = "~/SIG/Geo_util/Functions.RData")

prj.str <- CRS("+proj=utm +zone=20 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

geo.str <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

hyb.param <- read.csv("~/SIG/Research_Partners/hyb_param.csv", stringsAsFactors = F)
names(hyb.param) <- gsub("\\.", "", names(hyb.param))

#Red-Yellow-Darkgreen
cols <- colorRampPalette(c("red", "yellow", "darkgreen"), space = "Lab")

# Elevation color ramp
elev_cols <- colorRampPalette(c("#AFF0E9", "#FFFFB3", "#008040", "#FCBA03",
                                "#800000", "#69300D", "#ABABAB", "#FFFCFF"),
                              space = "Lab")

# Electrical conductivity color ramp
ec_cols <- colorRampPalette(c("#2892C7", "#FAFA64", "#E81014"), space = "Lab")

# Organic matter color ramp
om_cols <- colorRampPalette(c("#D6D69C", "#734D00", "#000000"), space = "Lab")

# SWI color ramp
swi_cols <- colorRampPalette(c("#C2523C", "#EDA113", "#FFFF00", "#00DB00", 
                               "#20998F", "#0B2C7A"), space = "Lab")

# CEC color ramp
cec_cols <- colorRampPalette(c("#BA1414", "#FFFFBF", "#369121"), space = "Lab")

# Build landsat list of polygons for matching scenes
lndst_01 <- readOGR("~/SIG/Geo_util/raster/arg/", "lndst_scn_g")
for (a in names(lndst_01@data)) {
  lndst_01@data[, a] <- as.character(lndst_01@data[, a])
}
lndst_02 <- list()
for (a in seq_along(lndst_01)) {
  lndst_02[a] <- lndst_01[a,]
}
lndst.pol <- lndst_02
rm(a, lndst_01, lndst_02)

#Function to get Landsat scene Path and Row from spatial object
scn_pr <- function(sp.layer) {
  require(sp)
  # Check projection and change to WGS84
  if (is.projected(sp.layer)) {
    library(rgdal)
    sp.layer <- spTransform(sp.layer, geo.str)
  }
  require(rgeos)
  pth_rw_lst <- character()
  # For each polygon of landsat scenes check which of them covers the layer
  for (a in lndst.pol) {
    if (gCovers(a, sp.layer)) {
      pth_rw_lst <- c(pth_rw_lst, paste0(sub("P", "", a@data["PATH"]),
                                         sub("R", "", a@data["ROW"])))
    }
  }
  #Return vector of landsat path and rows
  return(pth_rw_lst)
}

# Build SRTM list of polygons for matching scenes
hgt.lst <- list.files("~/SIG/Geo_util/raster/arg/srtm_1s",
                      ".hgt$", full.names = T)
srtm.pol <- list()
for (a in seq_along(hgt.lst)) {
  r <- raster(hgt.lst[a])
  proj4string(r) <- geo.str
  e <- r@extent
  m <- matrix(c(e[1], e[1], e[2], e[2],
                e[3], e[4], e[4], e[3]),
              ncol = 2)
  b <- spPolygons(m, crs = proj4string(r),
                  attr = data.frame(LL = sub(".hgt", "", basename(hgt.lst[a])),
                                    stringsAsFactors = F))
  srtm.pol[a] <- b
  rm(a, r, e, m, b)
}
rm(hgt.lst)

#Function to get SRTM scenes from spatial object
srtm_pr <- function(sp.layer) {
  require(sp)
  # Check projection and change to WGS84
  if (is.projected(sp.layer)) {
    library(rgdal)
    sp.layer <- spTransform(sp.layer, geo.str)
  }
  require(rgeos)
  pth_rw_lst <- character()
  # For each polygon of srtm check which of them overlaps with the layer
  for (a in srtm.pol) {
    if (gIntersects(a, sp.layer)) {
      pth_rw_lst <- c(pth_rw_lst, a@data[1, "LL"])
    }
  }
  #Return vector of landsat path and rows
  return(pth_rw_lst)
}

#Function to get and filter landsat images from selected paths
mk_vi_stk <- function(sp.layer, vindx = "EVI", buff = 30,
                      st.year = 1990, vi.thr = 1500, cv.lim = 15) {
  if (substr(class(sp.layer), 1, 15)[1] != "SpatialPolygons") {
    stop("sp.layer isn't a SpatialPolygon* object")
  }
  require(sp)
  require(rgdal)
  require(rgeos)
  require(raster)
  # Save current directory to return later
  curr.wd <- getwd()
  # Set current directory to the one that has the VI images
  setwd(paste0("~/SIG/Geo_util/raster/arg/" , vindx, "_Landsat/"))
  # Create a list of available images 
  img.lst <- list.files(".", ".tif$")
  # Check projection of layer and project to measure distances
  if (is.projected(sp.layer) == F) {
    sp.layer <- spTransform(sp.layer, prj.str)
  }
  # Assign ownership of holes to parent polygons
  sp.comm <- createSPComment(sp.layer)
  # Create buffer of polygon for border effect
  sp.layer <- gBuffer(sp.comm, width = -buff)
  # Reproject buffered layer to WGS84
  sp.layer <- spTransform(sp.layer, geo.str)
  # Get on which landsat path rows the layer intersects
  scn.pr <- scn_pr(sp.layer)
  # Create empty stacks and data.frames to store information
  r.stk <- stack()
  r.stk2 <- stack()
  df1 <- data.frame()
  df2<- data.frame()
  i <- 1
  # Go through the images until finding the one that fully covers the layer
  vi.lst <- grep(paste(scn.pr, collapse = "|"), img.lst, value = T)
  r.base <- raster(vi.lst[i])
  proj4string(r.base) <- geo.str
  r.crp.bs <- crop(r.base, sp.layer, snap = "near")
  while(any(is.na(r.crp.bs[]))) {
    i <- i + 1
    r.base <- raster(vi.lst[i])
    proj4string(r.base) <- geo.str
    r.crp.bs <- crop(r.base, sp.layer, snap = "near")
  }
  # Mask the selected image to use as base for the next ones
  r.crp.bs <- mask(r.crp.bs, sp.layer)
  # For each image in the list of intersecting ones, do...
  for (c in vi.lst) {
    # Get the year of current image
    scn.year <- as.numeric(substr(c, 10, 13))
    # Check if passes the year limit
    if (scn.year >= st.year) {
      r <- raster(c)
      proj4string(r) <- geo.str
      # Crop raster with polygon
      r.crp <- crop(r, sp.layer, snap = "near")
      # Check if the are any NAs in the cropped raster, go on if there aren't
      if (any(is.na(r.crp[])) == F) {
        # Make mask of the crop
        r.crp <- mask(r.crp, sp.layer)
        # Calculate median and cv
        vi.mdn <- cellStats(r.crp, median)
        r.cv <- (cellStats(r.crp, sd) /
                   cellStats(r.crp, mean)) * 100
        # Check if stat meet requirements
        if (vi.mdn > vi.thr & r.cv < cv.lim){
          if (compareRaster(r.crp, r.crp.bs, extent = T, rowcol = T,
                            crs = F, res = F, orig = F, rotation = T,
                            values = F, stopiffalse = F,
                            showwarning = F) == F) {
            # If this mask doesn't match the base it is resampled
            r.crp <- resample(r.crp, r.crp.bs, method = "bilinear")
          }
          # Add the current mask to the stack
          r.stk <- stack(r.stk, r.crp)
          # Add this mask values to a reference data frame
          df1 <- rbind.data.frame(df1, data.frame(SCN = c,
                                                  Year = scn.year,
                                                  VI = vi.mdn,
                                                  stringsAsFactors = F))
        }
      }
    }
  }
  # The following will leave only one image per year
  if (length(unique(df1$Year)) < length(df1$Year)) {
    for (d in unique(df1$Year)) {
      # Leave the one with highest median
      vi.max <- max(df1[df1$Year == d, "VI"])
      df2 <- rbind.data.frame(df2, df1[df1$VI == vi.max,])
    }
    for (e in df2$SCN) {
      r.stk2 <- stack(r.stk2, r.stk[[which(df1$SCN == e)]])
    }
  } else {
    r.stk2 <- r.stk
  }
  # Project stack
  r.stk2 <- projectRaster(r.stk2, crs = prj.str, method = "bilinear")
  # Return to original working directory
  setwd(curr.wd)
  return(r.stk2)
}

#Function to get and filter srtm images from selected lat/long
dem_srtm <- function(sp.layer, buff = 30) {
  if (substr(class(sp.layer), 1, 15)[1] != "SpatialPolygons") {
    stop("sp.layer isn't a SpatialPolygon* object")
  }
  require(sp)
  require(rgdal)
  require(rgeos)
  require(raster)
  if (buff != 0) {
    # Check projection of layer and project to measure distances
    if (is.projected(sp.layer) == F) {
      sp.layer <- spTransform(sp.layer, prj.str)
    }
    # Assign ownership of holes to parent polygons
    sp.comm <- createSPComment(sp.layer)
    # Create buffer of polygon for border effect
    sp.layer <- gBuffer(sp.comm, width = -buff)
  }
  if (is.projected(sp.layer)) {
    # Reproject buffered layer to WGS84
    sp.layer <- spTransform(sp.layer, geo.str)
  }
  # Get on which srtm polygon the layer intersects
  r.pr <- srtm_pr(sp.layer)
  # Load and mask each raster
  for (b in seq_along(r.pr)) {
    r <- raster(paste0("~/SIG/Geo_util/raster/arg/srtm_1s/", r.pr[b], ".hgt"))
    proj4string(r) <- geo.str
    crp <- crop(r, sp.layer)
    msk <- mask(crp, sp.layer)
    assign(paste0("rmsk", b), msk)
    rm(b, r, crp, msk)
  }
  if (length(r.pr) > 1) {
    # Build a first mosaic from the first two
    msc <- mosaic(rmsk1, rmsk2, fun = mean)
    # If there are more than two add htem to the mosaic
    if (length(r.pr) > 2) {
      # Get masks in the environment
      r.lst <- grep("rmsk[0-9]", ls(), value = T)
      # For each mask add it to the mosaic
      for (c in 3:length(r.lst)) {
        msc <- mosaic(msc, get(r.lst[c]), fun = mean)
      }
    }
  } else {
    msc <- rmsk1
  }
  msc.pnt <- rasterToPoints(msc, spatial = T)
  msc.p <- spTransform(msc.pnt, CRSobj = prj.str)
  names(msc.p) <- "elev"
  return(msc.p)
}

#Function to reclassify a raster in n classes by jenks
rstr_rcls <- function(raster.lyr, n.class = 3, val = 1:3, style = "fisher") {
  if (class(raster.lyr) != "RasterLayer"){
    stop("Input object isn't a RasterLayer object")
  }
  if (n.class != length(val)) {
    stop("Number of classes doesn't match number of values")
  }
  require(classInt)
  require(raster)
  nw.vls <- val
  # Cut the data in the selected number of classes
  cut.vals <- classIntervals(raster.lyr[!is.na(raster.lyr)], n = n.class, style = style)$brks
  # Reclassification matrix
  mat <- as.matrix(data.frame(from = cut.vals[1:n.class],
                              to = cut.vals[2:(n.class + 1)],
                              beco = nw.vls))
  # Reclassification according to matrix
  raster.rcls <- reclassify(raster.lyr, mat, include.lowest = T)
  return(raster.rcls)
}

#Rsaga DEM Covariates
dem_cov <- function(DEM.layer, dem.attr = "DEM") {
  if (class(DEM.layer)[1] != "SpatialPointsDataFrame") {
    if (class(DEM.layer)[1] != "RasterLayer") {
      stop("DEM.layer isn't a SpatialPointsDataFrame or RasterLayer object")
    }
  }
  require(sp)
  require(RSAGA)
  require(raster)
  if ("./DEM_derivates" %in% list.dirs()) {
    unlink("./DEM_derivates", recursive = T, force = T)
  }
  # Topography derivates folder creation
  dir.create("DEM_derivates")
  # Save current directory
  curr.wd <- getwd()
  setwd("./DEM_derivates")
  base.lyr <- DEM.layer
  # Store dem layer CRS
  lyr.crs <- CRS(proj4string(base.lyr))
  # If layer is SPDF convert to raster
  if (class(base.lyr)[1] == "SpatialPointsDataFrame") {
    base.rstr <- pnt2rstr(base.lyr, dem.attr)
  } else {
    base.rstr <- base.lyr
  }
  dem.file <- "dem.tif"
  # Save raster as GeoTIFF
  writeRaster(base.rstr, dem.file)
  # Convert raster layer to points to add the layers
  if (class(DEM.layer)[1] == "RasterLayer") {
    base.lyr <- rasterToPoints(base.lyr, spatial = T)
  }
  # Convert GeoTIFF to Saga Grid
  rsaga.import.gdal(dem.file, show.output.on.console = F)
  # Calculate Slope and Aspect
  rsaga.geoprocessor("ta_morphometry", module = 0,
                     param = list(ELEVATION = "dem.sgrd",
                                  SLOPE = "Slope",
                                  ASPECT = "Aspect",
                                  METHOD = 1),
                     show.output.on.console = F)
  # Calculate Curvatures
  rsaga.geoprocessor("ta_morphometry", module = 0,
                     param = list(ELEVATION = "dem.sgrd",
                                  CURV = "Curv",
                                  HCURV = "PrCurv",
                                  VCURV = "PlCurv",
                                  METHOD = 5),
                     show.output.on.console = F)
  # Calculate Catchment Area
  rsaga.geoprocessor("ta_hydrology", module = 18,
                     param = list(DEM = "dem.sgrd",
                                  AREA = "Catch_Area"),
                     show.output.on.console = F)
  # Calculate CTI
  rsaga.geoprocessor("ta_hydrology", module = 20,
                     param = list(SLOPE = "Slope.sgrd",
                                  AREA = "Catch_Area.sgrd",
                                  TWI = "CTI"),
                     show.output.on.console = F)
  # Calculate Convergence Index
  rsaga.geoprocessor("ta_morphometry", module = 1,
                     param = list(ELEVATION = "dem.sgrd",
                                  RESULT = "Conv_Index"),
                     show.output.on.console = F)
  # Calculate LS Factor
  rsaga.geoprocessor("ta_hydrology", module = 22,
                     param = list(SLOPE = "Slope.sgrd",
                                  AREA = "Catch_Area.sgrd",
                                  LS = "LS_Factor",
                                  METHOD = 2),
                     show.output.on.console = F)
  # Calculate Saga Wetness Index
  rsaga.wetness.index("dem.sgrd", "SWI.sgrd",
                      show.output.on.console = F)
  # List all grid files created
  sgrd.lst <- list.files(".", pattern = ".sgrd$")
  sgrd.lst <- sgrd.lst[!(sgrd.lst %in% "dem.sgrd")]
  # Iterate over list of grids
  for (a in sgrd.lst) {
    # Get name for tiff
    tif.name <- paste0(sub(".sgrd", "", a), ".tif")
    # Convert Saga grid to GeoTIFF
    rsaga.geoprocessor("io_gdal", module = 2,
                       param = list(GRIDS = a,
                                    FILE = tif.name),
                       show.output.on.console = F)
    # Open tiff as RasterLayer
    terr.tif <- raster(tif.name)
    # Get data from the tiff that intersects with points
    terr.data <- extract(terr.tif, base.lyr)
    # Get name for column in attribute table
    attr.name <- sub(".sgrd", "", a)
    # Add column to SPDF
    base.lyr@data[attr.name] <- terr.data
  }
  if (any(is.na(base.lyr@data))) {
    require(DMwR)
    base.lyr@data <- knnImputation(base.lyr@data)
  }
  setwd(curr.wd)
  return(base.lyr)
}

# Defining the MBA interpolation function
int_fx <- function(base.pnts, obs.pnts, vrbl, moran = F, dist = 20, clean = T, krig = F) {
  if (substr(class(base.pnts), 1, 13)[1] != "SpatialPoints" |
        class(obs.pnts)[1] != "SpatialPointsDataFrame") {
    stop("at least one of the inputs isn't a SpatialPoints* object")
  }
  require(sp)
  require(MBA)
  base.map <- base.pnts
  # Creation of prediction grid from base points
  pred.grid <- base.map@coords
  # Duplicate deletion
  if (length(zerodist(obs.pnts)) == 0) {
    obs.pnts <- remove.duplicates(obs.pnts)
  }
  obs.map <- obs.pnts
  # For every variable to interpolate do...
  for (a in vrbl) {
    # Leave positive values
    obs.vrbl <- subset(obs.map, eval(parse(text = a)) > 0)
    # Delete NAs
    obs.vrbl <- obs.vrbl[!is.na(obs.vrbl@data[, a]),]
    # Get vector of numbers
    vrbl.data <- obs.vrbl@data[, a]
    # Cleaning by spatial association
    if (moran) {
      obs.vrbl <- moran_cln(obs.vrbl, a, dist = dist)
      vrbl.data <- obs.vrbl@data[, a]
    }
    # cleaning by data distribution
    if (clean) {
      require(robustbase)
      # Outlier detection
      vrbl.out <- adjboxStats(vrbl.data)$out
      # If there are outliers they will be removed
      if (length(vrbl.out)>0) {
        obs.vrbl <- obs.vrbl[!(vrbl.data %in% vrbl.out),]
        vrbl.data <- obs.vrbl@data[, a]
      }
    }
    if (krig) {
      # Try fit fit variogram
      vg.fit <- try(var_fit(obs.vrbl, a, cln = F), silent = T)
      # If it didn't result in error, krige
      if (class(vg.fit)[1] != "try-error") {
        # Krig with fitted lambda
        vrbl.krig <- krige((vrbl.data ^ vg.fit[[1]]) ~ 1, locations = obs.vrbl,
                           newdata = base.map, model = vg.fit[[3]])
        # Reverse conversion with lambda to obtain original units
        vrbl.krig@data[, 1] <- vrbl.krig@data[, 1] ^ (1 / vg.fit[[1]])
        # Add data to layer
        base.map@data[, a] <- vrbl.krig@data[, 1]
        # Create new layer to plot variable and prediction variance
        map.plot <- base.map[a]
        map.plot@data[paste0(a, "_var")] <- vrbl.krig@data[, 2]
        # save a pdf of variogram and kriging
        pdf(paste0(a, ".pdf"), paper = "letter", width = 0, height = 0)
        #layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
        # Plot exp and fitted variogram
        print(plot(vg.fit[[2]], model = vg.fit[[3]], pch = 21, cex = 2, lty = 2, lwd = 2,
                   col = "red", xlab = "Distance", ylab = "Semivariance",
                   main = paste0("Emp. & fitted semivariogram for ", a)))
        # Plot prediction and prediction variance
        par(mfrow = c(2, 1))
        plot(pnt2rstr(map.plot, c(a, paste0(a, "_var"))), y = 1, col = cols(255))
        plot(pnt2rstr(map.plot, c(a, paste0(a, "_var"))), y = 2, col = rev(cols(255)))
        dev.off()
        # If variogram fitting wasn't successful go with MBA
      } else {
        # Create predicion matrix needed by MBA
        obs.mat <- cbind(obs.vrbl@coords[, 1:2], vrbl.data)
        # Interpolate
        obs.mba <- mba.points(obs.mat, pred.grid, verbose = F)
        # transform prediction to data.frame
        obs.mba <- data.frame(obs.mba$xyz.est)
        # Replace possible negative values with zeroes
        obs.mba[obs.mba$z < 0, 3] <- 0
        # Add data to layer
        base.map@data[, a] <- obs.mba$z
      }
    } else {
      obs.mat <- cbind(obs.vrbl@coords[, 1:2], vrbl.data)
      obs.mba <- mba.points(obs.mat, pred.grid, verbose = F)
      obs.mba <- data.frame(obs.mba$xyz.est)
      obs.mba[obs.mba$z < 0, 3] <- 0
      base.map@data[, a] <- obs.mba$z
    }
  }
  # Return original base layer with added columns
  return(base.map)
}

# plant population response by hybrid
hyb_pp <- function(hybrid, exp.yld, step = 1235) {
  num.param <- as.numeric(hyb.param[1:8, hybrid])
  pl.pop <- vector()
  min.kn.yld <- 7
  if (hyb.param[9, hybrid] == "L") {
    seed.rates <- seq(30000, 150000, step)
    for (a in seq_along(exp.yld)) {
      hyb.pp <- ((-num.param[2] + 2 * num.param[3] * num.param[7] -
                    (num.param[5] * exp.yld[a] - num.param[5] * num.param[8]) -
                    2 * num.param[6] * num.param[7] * (num.param[8] - exp.yld[a])) /
                   (2 * num.param[3] + 2 * num.param[6] * (exp.yld[a] - num.param[8]))) * 10000
      hyb.pp <- seed.rates[which(abs(seed.rates - hyb.pp) == min(abs(seed.rates - hyb.pp)))]
      pl.pop[a] <- hyb.pp
    }
    return(pl.pop)
  } else {
    pp.min <- 35000
    pp.kn <- ((-num.param[2] + 2 * num.param[3] * num.param[7] -
                 (num.param[5] * min.kn.yld - num.param[5] * num.param[8]) -
                 2 * num.param[6] * num.param[7] * (num.param[8] - min.kn.yld)) /
                (2 * num.param[3] + 2 * num.param[6] * (min.kn.yld - num.param[8]))) * 10000
    seed.rates <- seq(35000 + step, 150000, step)
    pp.lm <- lm(y ~ x, data = data.frame(x = c(3, min.kn.yld), y = c(pp.min, pp.kn)))
    for (a in seq_along(exp.yld)) {
      if (exp.yld[a] < 3) {
        hyb.pp <- pp.min
      }
      if (exp.yld[a] >= 3 & exp.yld[a] < min.kn.yld) {
        hyb.pp <- predict(pp.lm, newdata = data.frame(x = exp.yld[a]))
        hyb.pp <- seed.rates[which(abs(seed.rates - hyb.pp) == min(abs(seed.rates - hyb.pp)))]
      }
      if (exp.yld[a] >= min.kn.yld) {
        hyb.pp <- ((-num.param[2] + 2 * num.param[3] * num.param[7] -
                      (num.param[5] * exp.yld[a] - num.param[5] * num.param[8]) -
                      2 * num.param[6] * num.param[7] * (num.param[8] - exp.yld[a])) /
                     (2 * num.param[3] + 2 * num.param[6] * (exp.yld[a] - num.param[8]))) * 10000
        hyb.pp <- seed.rates[which(abs(seed.rates - hyb.pp) == min(abs(seed.rates - hyb.pp)))]
      }
      pl.pop[a] <- hyb.pp
    }
    return(pl.pop)
  }
}

presc_grid <- function(sp.layer, pred.model, hybrid, points = T,
                       fill = F, quantile = NULL, step = 1235) {
  if (class(sp.layer)[1] != "SpatialPointsDataFrame") {
    stop("sp.layer isn't a SpatialPointsDataFrame object")
  }
  require(sp)
  # If one wants the NAs can be filled
  if (fill) {
    if (any(is.na(sp.layer@data))) {
      require(DMwR)
      # Impute de missing data
      sp.layer@data <- knnImputation(sp.layer@data)
    }
  }
  # If the results are needed in polygon
  if (points == F) {
    require(raster)
    lyr.crs <- CRS(proj4string(sp.layer))
    sp.spix <- SpatialPixelsDataFrame(sp.layer,
                                      data = data.frame(1:nrow(sp.layer@data)),
                                      proj4string = lyr.crs)
    # Convert to raster
    sp.r <- raster(sp.spix)
    # Convert to polygons
    sp.poly <- rasterToPolygons(sp.r, dissolve = F)
  } else {
    sp.poly <- sp.layer
  }
  # The next steps according to the prediction model selet the used variables
  if (class(pred.model)[1] == "randomForest") {
    require(randomForest)
    usd.var <- dimnames(pred.model$importance)[[1]]
  }
  if (class(pred.model)[1] == "quantregForest") {
    require(quantregForest)
    usd.var <- names(pred.model$forest$ncat)
  }
  if (class(pred.model) == "cubist") {
    require(Cubist)
    usd.var <- pred.model$vars$all
  }
  if (class(pred.model) == "gbm") {
    require(gbm)
    usd.var <- pred.model$var.names
  }
  if (class(datab.extr) == "extraTrees") {
    usd.var <- names(sp.layer)
  }
  # Create column with specific hybrid
  sp.poly@data <- data.frame(Hybrid = rep(hybrid, nrow(sp.poly@data)))
  # Create expected yield column with prediction
  sp.poly@data["Exp_Yld"] <- predict(pred.model, sp.layer@data[usd.var],
                                     quantiles = quantile)
  # If germplasm is an inbred line divide prediction by 2
  if (hyb.param[9, hybrid] == "L") {
    sp.poly@data["Exp_Yld"] <- sp.poly@data["Exp_Yld"] / 2
  }
  # Calculate plant population according to expected yield
  sp.poly@data["PP"] <- hyb_pp(hybrid, sp.poly@data[, "Exp_Yld"], step)
  return(sp.poly)
}

grd_m <- function(sp.layer, dist = 10) {
  if (substr(class(sp.layer), 1, 15)[1] != "SpatialPolygons") {
    stop("sp.layer isn't a SpatialPolygons* object")
  }
  # If in WGS84 project to UTM
  if (is.projected(sp.layer) == F) {
    library(rgdal)
    sp.layer <- spTransform(sp.layer, prj.str)
  }
  require(sp)
  # Get bounding box
  lyr.bb <- sp.layer@bbox
  # Calculate regularly spaced coordinates
  grd.1 <- expand.grid(x = seq(lyr.bb[1, 1], lyr.bb[1, 2], by = dist),
                       y = seq(lyr.bb[2, 1], lyr.bb[2, 2], by = dist))
  # Get layer CRS
  grd.crs <- CRS(proj4string(sp.layer))
  # Create SPDF from regular coordinates
  grd.1 <- SpatialPointsDataFrame(grd.1, data = grd.1, proj4string = grd.crs)
  # Remove the ones outside the boundary
  grd.inp <- !is.na(over(grd.1, SpatialPolygons(sp.layer@polygons,
                                                proj4string = grd.crs)))
  grd.1 <- grd.1[grd.inp,]
  return(grd.1)
}

mz_smth <- function(sp.layer, area = 2500) {
  if (substr(class(sp.layer), 1, 15)[1] == "SpatialPolygons" |
        class(sp.layer)[1] == "RasterLayer") {
    require(rgrass7)
    require(raster)
    # If the input is a raster convert to polygons and dissolve by zone
    if (class(sp.layer) == "RasterLayer") {
      sp.layer <- rasterToPolygons(sp.layer, dissolve = T, na.rm = T)
    }
    # If in WGS84 project to UTM
    if (is.projected(sp.layer) == F) {
      library(rgdal)
      sp.layer <- spTransform(sp.layer, prj.str)
    }
    # Check wether GRASS is running, else initialize
    if (nchar(Sys.getenv("GISRC")) == 0) {
      initGRASS(gisBase = "c:/Program Files (x86)/GRASS GIS 7.0.0",
                override = TRUE)
    }
    # Convert multipart to singlepart
    sp.layer <- disaggregate(sp.layer)
    zm.pol <- paste0(sample(letters, 1), substr(basename(tempfile()), 9, 14))
    # Convert name 'layer' to 'Zone'
    names(sp.layer) <- sub("layer", "Zone", names(sp.layer))
    # Write GRASS vector
    writeVECT(sp.layer, zm.pol, v.in.ogr_flags = "o")
    zm.gnrl <- paste0(sample(letters, 1), substr(basename(tempfile()), 9, 14))
    # Smooth lines of polygons
    execGRASS("v.generalize", flags = c("overwrite", "quiet"), input = zm.pol,
              output = zm.gnrl, method = "snakes", threshold = 1)
    zm.cln <- paste0(sample(letters, 1), substr(basename(tempfile()), 9, 14))
    # Remove small/sliver polygons
    execGRASS("v.clean", flags = c("overwrite", "quiet"), input = zm.gnrl,
              output = zm.cln, tool = "rmarea", threshold = area)
    # Read back cleaned layer
    zm.fnl <- readVECT(zm.cln)
    # If no CRS, define one
    if (is.na(zm.fnl@proj4string)) {
      proj4string(zm.fnl) <- prj.str
    }
    # Remove 'cat' column from data.frame
    zm.fnl@data["cat"] <- NULL
    return(zm.fnl)
  } else {
    stop("sp.layer isn't a SpatialPolygons* or RasterLayer object")
  }
}

pnt2rstr <- function(sp.layer, field = names(sp.layer)){
  if (class(sp.layer)[1] != "SpatialPointsDataFrame") {
    stop("sp.layer isn't a SpatialPointsDataFrame")
  }
  require(sp)
  require(raster)
  # Get layer CRS
  lyr.crs <- CRS(proj4string(sp.layer))
  # Check if there's more than one field to convert to raster
  p2g <- try(points2grid(sp.layer))
  if (class(p2g)[1] == "try-error") {
    stop("points aren't regularly spaced")
  }
  if (length(field) > 1) {
    # Create empty stack
    sp.rstr <- stack()
    for (a in field) {
      if (a %in% names(sp.layer@data) == F) {
        stop("field isn't an attribute in SpatialPointsDataFrame")
      }
      sp.spix <- SpatialPixelsDataFrame(sp.layer,
                                        data = sp.layer@data[a],
                                        proj4string = lyr.crs)
      sp.rstr <- stack(sp.rstr, raster(sp.spix))
    }
  } else {
    if (field %in% names(sp.layer@data) == F) {
      stop("field isn't an attribute in SpatialPointsDataFrame")
    }
    sp.spix <- SpatialPixelsDataFrame(sp.layer,
                                      data = sp.layer@data[field],
                                      proj4string = lyr.crs)
    sp.rstr <- raster(sp.spix)
  }
  return(sp.rstr)
}

geo_centroid <- function(sp.layer){
  if (substr(class(sp.layer), 1, 7)[1] != "Spatial") {
    stop("sp.layer isn't a Spatial* object")
  }
  # If in UTM project to WGS84
  if (is.projected(sp.layer) == T) {
    library(rgdal)
    sp.layer <- spTransform(sp.layer, geo.str)
  }
  require(rgeos)
  coord <- gCentroid(sp.layer)@coords
  names(coord) <- c("Lon", "Lat")
  coord <- coord[c(2, 1)]
  return(coord)
}

moran_cln <- function(sp.layer, vrbl, dist = 20) {
  if (class(sp.layer)[1] != "SpatialPointsDataFrame") {
    stop("sp.layer isn't a SpatialPointsDataFrame object")
  }
  require(spdep)
  # Remove NA's
  sp.layer <- sp.layer[-which(is.na(sp.layer@data[,vrbl])),]
  # Identify neighbours points by Euclidean distance
  nb.lst <- dnearneigh(sp.layer, d1 = 0, d2 = dist)
  # Get number of neighbours in the neighbours list
  nb.crd <- card(nb.lst)
  # Remove points with no neighbors
  spl.noznb <- subset(sp.layer, nb.crd > 0)
  # Also in neighbor list
  nb.noznb <- subset(nb.lst, nb.crd > 0)
  # Convert cleaned neighbor list to spatial weighted list
  w.mat <- nb2listw(nb.noznb, style = "W")
  # Get numerical data of variable
  vrbl.dt <- spl.noznb@data[, vrbl]
  # Compute the lag vector V x
  wx <- lag.listw(w.mat, vrbl.dt)
  # Lineal model of lagged vs observed variable
  xwx.lm <- lm(wx ~ vrbl.dt)
  # Compute regression (leave-one-out deletion) diagnostics for linear model
  # and only get get the logical influence matrix
  infl.xwx <- influence.measures(xwx.lm)[["is.inf"]]
  # Convert to numeric, 6 column matrix for computation of sums
  infl.mat <- matrix(as.numeric(infl.xwx), ncol = 6)
  # Calculate moran scatterplot parameters
  #mp <- moran.plot(vrbl.dt, w.mat, quiet = T)
  # Get those rows where at least one index is TRUE
  mp.out <- which(rowSums(infl.mat) != 0)
  # Calculate local moran
  lmo <- localmoran(vrbl.dt, w.mat, p.adjust.method = "bonferroni",
                    alternative = "less")
  # Convert to data.frame to select data
  lmo <- data.frame(lmo)
  # Get rows wheres indices are significative
  lm.out <- which(lmo[, "Ii"] <= 0 | lmo[, "Pr.z...0."] <= 0.05)
  # Get the unique rows to delete
  all.out <- unique(c(mp.out, lm.out))
  # Remove them from SPDF
  spl.noznb <- spl.noznb[-all.out, ]
  return(spl.noznb)
}

var_fit <- function(sp.layer, vrbl, cln = F, plot = F){
  if (class(sp.layer)[1] != "SpatialPointsDataFrame") {
    stop("sp.layer isn't a SpatialPointsDataFrame object")
  }
  require(automap)
  require(MASS)
  require(gstat)
  vrbl.nm <- vrbl
  # Leave only positive observed values
  vrbl.spdt <- subset(sp.layer, eval(parse(text = vrbl.nm)) > 0)
  # Remove NAs
  vrbl.spdt <- vrbl.spdt[!is.na(vrbl.spdt@data[, vrbl.nm]),]
  # Remove extra coordinate dimensions other than x, y
  if (ncol(vrbl.spdt@coords) > 2) {
    vrbl.spdt@coords <- vrbl.spdt@coords[, 1:2]
  }
  # Remove duplicates if any
  if (length(zerodist(vrbl.spdt)) == 0) {
    vrbl.spdt <- remove.duplicates(vrbl.spdt)
  }
  # Numeric data vector
  vrbl.dt <- vrbl.spdt@data[, vrbl.nm]
  # Cleaning from distribution
  if (cln) {
    require(robustbase)
    # Get outliers
    vrbl.out <- adjboxStats(vrbl.dt)$out
    # If any they are removed
    if (length(vrbl.out)>0) {
      vrbl.spdt <- vrbl.spdt[!(vrbl.dt %in% vrbl.out),]
      vrbl.dt <- vrbl.spdt@data[, vrbl.nm]
    }
  }
  # Set number of samples according to size
  ifelse(length(vrbl.spdt) > 5000, n.samp <- 5000, n.samp <- length(vrbl.spdt))
  # Assing object in Global Environment for lineal model
  assign(".vrbl.dt", vrbl.dt, envir = .GlobalEnv)
  # Lambda calculation for normality
  vrbl.bc <- boxcox(lm(.vrbl.dt ~ 1), lambda = seq(-5, 5, 0.01), plotit = F)
  # Delete GlobalEnv object
  remove(".vrbl.dt", envir = .GlobalEnv)
  # Test for normality and assignment of lambda
  ifelse(shapiro.test(sample(vrbl.dt, n.samp))$p.value < 0.05,
         lmbd <- vrbl.bc[["x"]][which.max(vrbl.bc[["y"]])],
         lmbd <- 1)
  # Calculation of maximum distance
  max.dist <- max(spDistsN1(vrbl.spdt, vrbl.spdt[1,]))
  #use.dist <- max.dist / 3
  # Variogram models
  vg.mods <- data.frame(short = c("Exp", "Sph", "Gau", "Mat", "Cir", "Wav"),
                        long = c("exponential", "spherical", "gaussian", "matern",
                                 "circular", "wave"),
                        stringsAsFactors = F)
  # Autofit to get initial parameters
  auto.vgm <- tryCatch(autofitVariogram((vrbl.dt ^ lmbd) ~ 1,
                                        vrbl.spdt,model = vg.mods[, 1],
                                        kappa = seq(0.01, 3, by = 0.01),
                                        cressie = T),
                       error = function(e) return("error"),
                       warning = function(w) return("warning"))
  if (class(auto.vgm)[1] == "character") {
    stop("a variogram could not be fitted")
  }
  if (auto.vgm$var_model[2, 3] > max.dist) {
    stop("a variogram could not be fitted")
  }
  if (auto.vgm$var_model[1, 2] > auto.vgm$var_model[2, 2]) {
    stop("a variogram could not be fitted")
  }
  if ((max(auto.vgm$var_model[, 2]) - min(auto.vgm$var_model[, 2])) < 0.0001) {
    stop("a variogram could not be fitted")
  }
  # Short and long varigogram models
  sh.mod <- as.character(auto.vgm$var_model$model[2])
  lg.mod <- vg.mods[vg.mods[1] == sh.mod, 2]
  # Experimental variogram creation
  vrbl.vgm <- variogram((vrbl.dt ^ lmbd) ~ 1, locations = vrbl.spdt, cressie = T)
  # Initial nugget, sill and range
  inug <- auto.vgm$var_model[1, 2]
  isill <- auto.vgm$var_model[2, 2]
  irange <- auto.vgm$var_model[2, 3]
  # Theoretical variogram fitting with selected parameters
  fit.res <- tryCatch(fit.variogram(vrbl.vgm, vgm(isill, sh.mod, irange, inug),
                                    fit.sills = T, fit.ranges = T,
                                    warn.if.neg = T),
                      error = function(e) return("error"),
                      warning = function(w) return("warning"))
  if (class(fit.res)[1] == "character" | attr(fit.res, "singular") |
        fit.res$psill[1] >= fit.res$psill[2]) {
    stop("a variogram could not be fitted")
  } else {
    vrbl.fit <- fit.res
  }
  # Wether to plot the variograms (exp & fit)
  if (plot) {
    print(
      plot(vrbl.vgm, model = vrbl.fit, pch = 21, cex = 2, lty = 2, lwd = 2,
           col = "red", xlab = "Distance", ylab = "Semivariance",
           main = paste0("Emp. & fitted semivariogram for ", vrbl.nm))
    )
  }
  #Return a list with: lambda, experimental and fitted variogram
  return(list(lmbd, vrbl.vgm, vrbl.fit))
}

kmz_sv <- function(sp.layer, kmz.name){
  require(plotKML)
  raster_layer <- pnt2rstr(sp.layer, c('DEM', 'Pred_OM', 'EC30', 'EC90'))
  kml_open(file.name = paste0(kmz.name, '.kml'),
           folder.name = kmz.name, 
           overwrite = T)
  kml_layer(raster_layer, 
            raster.name = 'DEM',
            subfolder.name='DEM (m)',
            colour = 'DEM',
            colour_scale = elev_cols(255))
  kml_layer(raster_layer, 
            raster.name = 'MO',
            subfolder.name='MO (%)',
            colour = 'Pred_OM',
            colour_scale = om_cols(255))
  kml_layer(raster_layer, 
            raster.name = 'EC30',
            subfolder.name='ECa 0-30 cm',
            colour = 'EC30',
            colour_scale = ec_cols(255))
  kml_layer(raster_layer, 
            raster.name = 'EC90',
            subfolder.name='ECa 0-90 cm',
            colour = 'EC90',
            colour_scale = ec_cols(255))
  
  kml_close(file.name = paste0(kmz.name, '.kml'))
  
  file = paste(kmz.name, '.zip', sep = '')
  zip(file, list.files(".", pattern = "png|kml"))
  file.rename(paste0(kmz.name, ".zip"), paste0(kmz.name, ".kmz"))
  file.remove(list.files(pattern = '.png|kml'))
}

# Import Veris data
veris_import <- function(vrs.fl = 'VSECOM', vrbl = c('EC30', 'EC90', 'Red', 'IR')){
  dat.fls <- list.files(path = "./Veris/", pattern = paste0(vrs.fl, '.*\\.txt'))
  for (b in dat.fls) {
    veris.n <- read.table(paste0("Veris/", b), 
                          header = TRUE, 
                          sep = "\t", 
                          skipNul = T)
    if (exists("veris")) {
      veris <- rbind(veris, veris.n)
    } else {
      veris <- veris.n
    }
  }
  ec30.nm <- grep("EC_SH|EC.SH|EC30", names(veris), ignore.case = T, value = T)
  ec90.nm <- grep("EC_DP|EC.DP|EC90", names(veris), ignore.case = T, value = T)
  names(veris) <- sub(ec30.nm, "EC30", names(veris))
  names(veris) <- sub(ec90.nm, "EC90", names(veris))
  require(sp)
  require(rgdal)
  veris.pnt <- SpatialPointsDataFrame(coords = veris[,c('Long','Lat')],
                                      proj4string = geo.str, 
                                      data = veris[,vrbl])
  veris.pnt <- spTransform(veris.pnt, prj.str)
  veris.pnt <- remove.duplicates(veris.pnt)
  
  elev.poly <- readOGR("./Elev", sub(".shp", "", 
                                   list.files("./Elev", pattern = ".shp$")))
  if (summary(elev.poly)$is.projected == F) {
    elev.poly <- spTransform(elev.poly, prj.str)
  }
  join <- over(veris.pnt, elev.poly['elevM'])
  
  veris.pnt@data['elevM'] <- join[1]
  
  writeOGR(veris.pnt, dsn = "./Veris", layer = "veris", 
           driver = "ESRI Shapefile", overwrite_layer = T, )
  return(veris.pnt)
}

var_cal <- function(sp.layer, var = 'OM', soil.layer = 'soil', pdf = T){
  require(rgdal)
  require(rgeos)
  require(ggplot2)
  require(reshape2)
  require(gridExtra)
  
  soil <- readOGR("./Soil", soil.layer)
  # soil <- readOGR("./Soil", "soil")
  if (summary(soil)$is.projected == F) {
    soil <- spTransform(soil, prj.str)
  }
  soil <- soil[soil@data[,var]>0,]
  # Create buffer of 10 m
  soil.buf <- gBuffer(soil, byid = TRUE, width = 10)
  join <- over(soil.buf, sp.layer, fn = mean)
  cal.db <- cbind(join, soil@data)
  
  lm1 <- as.formula(paste(var, '~', paste(c('Red', 'IR'), collapse = ' + ')))
  lm2 <- as.formula(paste(var, '~', paste(c('Red', 'EC30'), collapse = ' + ')))
  lm3 <- as.formula(paste(var, '~', paste(c('Red', 'EC90'), collapse = ' + ')))
  lm4 <- as.formula(paste(var, '~', paste(c('Red', 'DEM'), collapse = ' + ')))
  lm5 <- as.formula(paste(var, '~', paste(c('IR', 'EC30'), collapse = ' + ')))
  lm6 <- as.formula(paste(var, '~', paste(c('IR', 'EC90'), collapse = ' + ')))
  lm7 <- as.formula(paste(var, '~', paste(c('IR', 'DEM'), collapse = ' + ')))
  lm8 <- as.formula(paste(var, '~', paste(c('EC30', 'EC90'), collapse = ' + ')))
  lm9 <- as.formula(paste(var, '~', paste(c('EC30', 'DEM'), collapse = ' + ')))
  lm10 <- as.formula(paste(var, '~', paste(c('EC90', 'DEM'), collapse = ' + ')))
  lms <- paste0("lm", 1:10)
  LMs <- vector("list", 10)
  lm.summ <- data.frame('model' = numeric(10), 'min'= numeric(10), 'max' = numeric(10), 
                        'median' = numeric(10), 'mean' = numeric(10), 'r2' = numeric(10),
                        'rmse' = numeric(10), 'slope' = numeric(10), 'AIC' = numeric(10))
  
  for(i in 1:length(lms)){
    assign("LM", get(lms[i]), envir = .GlobalEnv)
    model <- lm(LM, data = cal.db)
    LMs[[i]] <- model
    vrbl.lm <- as.character(attr(terms(model), "term.labels"))
    
    if (!summary(model)$adj.r.squared == "NaN"){
      predicted <- as.vector(predict(model))
      actual <- cal.db[,var]
      cal.fit <- lm(predicted ~ actual)
          
      # Prediction and interpolation on veris.shp
      pred.int <- predict(model, newdata = sp.layer@data[vrbl.lm])
      lm.summ[i, 'model'] <- i
      lm.summ[i, 'min'] <- min(pred.int)
      lm.summ[i, 'max'] <- max(pred.int)
      lm.summ[i, 'median'] <- median(pred.int)
      lm.summ[i, 'mean'] <- mean(pred.int)
      lm.summ[i, 'r2'] <- summary(cal.fit)$adj.r.squared
      lm.summ[i, 'rmse'] <- summary(cal.fit)$sigma
      lm.summ[i, 'slope'] <- summary(cal.fit)$coefficients[2]
      lm.summ[i, 'AIC'] <- AIC(cal.fit)
      }
  }
  
  write.csv(lm.summ, paste0('Veris/', var, '_lms_summary.csv'))
  #lm.summ <- subset(lm.summ, lm.summ$min > 0.1)
  #lm.summ <- subset(lm.summ, lm.summ$r2 > 0)
  
  # PDF write
  if (pdf == T){
  pdf(paste0('Veris/', var, "_calibration.pdf"), 
      paper = "letter", width = 6, height = 0)
  
  for(i in 1:nrow(lm.summ)){
    model <- LMs[[lm.summ$model[i]]]
    vrbl.lm <- as.character(attr(terms(model), "term.labels"))
    label <- paste('Model N°', lm.summ$model[i], ": ", var, 
                   " ~ ", paste(vrbl.lm, collapse = " + "), sep="")
    pred.int <- predict(model, newdata = sp.layer@data[vrbl.lm])
    sp.layer@data['Pred'] <- pred.int
    plot1 <- ggplot() +
      geom_raster(data = sp.layer@data, aes(x, y, fill = Pred)) +
      coord_equal() + labs(x = 'Longitud', y = 'Latitud', fill = 'Pred', title = label) +
      scale_fill_gradientn(colours = cols(255)) + theme_bw() + 
      geom_point(data = data.frame(soil@coords), aes(x = coords.x1, y = coords.x2), shape = 19, size = 2) +
      theme(plot.title = element_text(size = 14, face = 'bold'),
            axis.text = element_text(size = 10),
            axis.title.x = element_text(size = 12, face = 'bold'),
            axis.title.y = element_text(size = 12, face = 'bold'))
    
    join <- over(soil.buf, sp.layer, fn = mean)
    cal <- cbind(join, soil@data)
    var.nm <- match(var, names(cal))
    names(cal)[var.nm] <- 'Lab'
    cal.fit <- lm(Pred ~ Lab, data = cal)
    eqn <- paste("r2:", format(summary(cal.fit)$adj.r.squared, digits=2), "/", "RMSE:", 
                 format(summary(cal.fit)$sigma, digits=2))
    plot2 <- ggplot(data = cal, aes(Lab, Pred)) +
      geom_point(size = 2) +
      theme_bw() +
      stat_smooth(method = lm, se = F, size = 0.5) +
      geom_abline(intercept = 0, slope = 1, size = 0.5) +
      labs(x = paste0(var, " Laboratorio"), y = paste0(var, " Predicho")) +
      theme(axis.text = element_text(size = 10),
            axis.title.x = element_text(size = 12, face = 'bold'),
            axis.title.y = element_text(size = 12, face = 'bold')) +
      annotate("text", label = eqn, parse = TRUE, x = Inf, y = -Inf,
               hjust = 1.1, vjust = -.5)
    
    title <- paste0('Min:', format(min(sp.layer$Pred), digits = 2),
                    ' / Median:', format(median(sp.layer$Pred), digits = 2),
                    ' / Mean:', format(mean(sp.layer$Pred), digits = 2),
                    ' / Max:', format(max(sp.layer$Pred), digits = 2))
    plot3 <- ggplot(sp.layer@data, aes(x = Pred, y = ..density..)) + 
      geom_histogram(fill="cornsilk", colour="grey60", size=.2) +
      geom_density() +
      labs(x = paste(var, "Predicho"), title = title) +
      theme(axis.text = element_text(size = 10),
            axis.title.x = element_text(size = 12, face = 'bold'),
            axis.title.y = element_text(size = 12, face = 'bold')) +
      theme_bw()
    
    grid.arrange(plot1, plot2, plot3, nrow=3)
  }
  dev.off()
  }
  
  return(LMs)  
}

trat_grd <- function(sp.layer, largo = 10, ancho, ang = 0, num.trat, n.pas = 1) {
  require(rgdal)
  require(rgeos)
  require(maptools)
  # Check if sp.layer is a spatialpolygon
  if (substr(class(sp.layer), 1, 15)[1] != "SpatialPolygons") {
    stop("sp.layer isn't a SpatialPolygons* object")
  }
  # Assignment of field boundary
  bound <- sp.layer
  # Extraction of bounding box
  bbox <- bound@bbox
  #Determination of cell size
  cell.size <- c("largo" = ancho, "ancho" = largo)
  lon.rng <- bbox[1,2] - bbox[1,1]
  lat.rng <- bbox[2,2] - bbox[2,1]
  # Number of rows and columns plu 50 to make it bigger
  nc <- round(lon.rng / cell.size[1], 0) + 50
  nr <- round(lat.rng / cell.size[2], 0) + 50
  # Assignment of number of treatments
  ntrat <- num.trat
  # Creation of the vector of treatments
  trt.ordr <- vector()
  for (b in 1:ceiling(nc / ntrat)) {
    rndm <- sample(ntrat, ntrat)
    trt.ordr <- c(trt.ordr, rndm)
  }
  pas.ordr <- vector()
  for (f in trt.ordr) {
    pas.ordr <- c(pas.ordr, rep(f, n.pas))
  }
  # Creation of an empty data.frame of the number of rows and columns
  mat <- data.frame(matrix(0, nrow = nr, ncol = nc))
  # Creatior of the treatment matrix
  mat.trt <- mat
  for (i in 1:ncol(mat)) {
    mat.trt[, i] <- pas.ordr[i]
  }
  mat.trt <- t(mat.trt)
  trt.v <- as.vector(mat.trt)
  # Creation of the vector of replications
  rep.ordr <- vector()
  for (f in 1:ceiling(nc / ntrat)) {
    rep.ordr <- c(rep.ordr, rep(f, ntrat * n.pas))
  }
  # Creation of the replications matrix
  mat.rep <- mat
  for (i in 1:ncol(mat)) {
    mat.rep[, i] <- rep.ordr[i]
  }
  mat.rep <- t(mat.rep)
  rep.v <- as.vector(mat.rep)
  # Lat-Lon point of the lower left corner
  min.bbox <- bbox[, "min"]
  # Offset of the corner
  cell.offset <- c(min.bbox[1] - lon.rng / 4,
                   min.bbox[2] - lat.rng / 4)
  # Creation of a rectangular grid of defined dimensions
  grd <- GridTopology(cell.offset, cell.size, c(nc, nr))
  # Conversion to spatial polygons
  pol.1 <- as.SpatialPolygons.GridTopology(grd, proj4string = prj.str)
  # Extraction of the IDs of all the polygons
  pol.lst <- list()
  for (a in seq_along(pol.1)) {
    pol.lst[a] <- slot(pol.1[a,]@polygons[[1]], "ID") 
  }
  # Creation of the data frame of the polygons
  data <- data.frame(Trat = trt.v, Rep = rep.v)
  row.names(data) <- unlist(pol.lst)
  # Adding the data frame to the polygons
  pol.2 <- SpatialPolygonsDataFrame(pol.1, data = data, match.ID = T)
  proj4string(pol.2) <- prj.str
  if (ang > 0) {
    # Rotation of the polygons by the defined angle
    pol.3 <- elide(pol.2, rotate = ang,
                   center = gCentroid(bound)@coords)
    proj4string(pol.3) <- prj.str
  } else {
    pol.3 <- pol.2
  }
  # Clipping of the polygon with the boundary
  pol.4 <- gIntersection(bound, pol.3, byid = T)
  # Creation of spatialpoints to join the attribute table with the clipped polygon
  pol.3.pnt <- SpatialPointsDataFrame(coordinates(pol.3),
                                      pol.3@data,
                                      proj4string = prj.str)
  pol.3.pnt <- gBuffer(pol.3.pnt,
                       width = (min(cell.size) - 0.1) / 2,
                       byid = T)
  # Extraction of the data frame rows that match with the clipped polygons
  df <- over(pol.4, pol.3.pnt)
  # Final spatialpolygondf with attribute table
  pol.5 <- SpatialPolygonsDataFrame(pol.4, data = pol.3@data[df,], match.ID = F)
  gc()
  return(pol.5)
}

multi_mz <- function(sp.layer, vrbls = c("DEM", "Aspect", "CTI", "Slope",
                                         "SWI", "EC30", "EC90", "OM",
                                         "CEC", "EVI_mean"),
                     n.mz = 3, dist = 20, plot = F) {
  if (class(sp.layer)[1] != "SpatialPointsDataFrame") {
    stop("sp.layer isn't a SpatialPointsDataFrame object")
  }
  library(ade4)
  library(spdep)
  if (any(is.na(sp.layer@data[vrbls]))) {
    library(DMwR)
    sp.layer@data <- knnImputation(sp.layer@data[vrbls])
  } else {
    sp.layer@data <- sp.layer@data[vrbls]
  }
  # Creation of the data.frame for PCA
  df <- sp.layer@data
  # Calculation of nearest neighbors based on selected distance
  n.neigh <- dnearneigh(sp.layer, 0, dist)
  # Creation of the spatial weighted neighbor list
  sp.w <- nb2listw(n.neigh, style = "W")
  # Calculation of the PCA on the selected variables
  data.pca <- dudi.pca(df, center = T, scannf = F)
  # Run of multispati function
  sp.mltspt <- multispati(data.pca, sp.w, scannf = F)
  # Creation of a SPDF with the created Spatial components
  sp.pca <- SpatialPointsDataFrame(sp.layer,
                                   data.frame(sp.mltspt$li),
                                   proj4string = proj4string(sp.layer))
  if (plot) {
    par(mfrow = c(1, 2))
    s.corcircle(data.pca$co)
    barplot(data.pca$eig * 10, names.arg = 1:length(data.pca$eig), 
            main = "Variances",
            xlab = "Principal Components",
            ylab = "Percentage of variances",
            col = "grey30", border = NA)
    par(mfrow = c(1, 1))
  }
  cs1 <- data.frame("Variable" = row.names(data.pca$c1), "CS1" = data.pca$c1)
  row.names(cs1) <- NULL
  names(cs1)[2:ncol(cs1)] <- paste0("CS", 1:(ncol(cs1)-1))
  cs1 <- cs1[order(-abs(cs1[, 2])),]
  row.names(cs1) <- NULL
  # Clusterization of the first component in the selected number of clusters
  pca.rast <- rstr_rcls(pnt2rstr(sp.pca, "CS1"), n.class = n.mz, val = 1:n.mz)
  # Smoothing and cleaning of the managemnent zones
  sp.pol <- mz_smth(pca.rast)
  print(cs1)
  return(sp.pol)
}

# Function to read shapefiles with proj info
read_shp <- function(dsn, layer) {
  require(rgdal)
  require(maptools)
  if (grep(".shp$", layer, ignore.case = T) == 1) {
    layer.ogr <- sub(".shp", "", layer)
  }
  shp.info <- ogrInfo(dsn, layer.ogr)[["p4s"]]
  shp.pth <- paste0(dsn, "/", layer)
  shp.sp <- readShapeSpatial(fn = shp.pth, proj4string = CRS(shp.info),
                             verbose = F, delete_null_obj = T)
  return(shp.sp)
}

save(lndst.pol, prj.str, geo.str, scn_pr, mk_vi_stk, rstr_rcls, int_fx, dem_cov,
     cols, elev_cols, ec_cols, om_cols, swi_cols, presc_grid, hyb.param, hyb_pp, grd_m,
     mz_smth, pnt2rstr, geo_centroid, moran_cln, var_fit, kmz_sv, veris_import,
     var_cal, trat_grd, multi_mz, srtm.pol, srtm_pr, dem_srtm, read_shp,
     file = "~/SIG/Geo_util/Functions.RData")
