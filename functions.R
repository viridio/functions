load(file = "~/SIG/Geo_util/Functions.RData")

prj_str <- function(zone) {
  zn.str <- paste0("+proj=utm +zone=", zone,
                   " +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  return(CRS(zn.str))
}

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
lndst_01 <- read_shp("~/SIG/Geo_util/raster/arg/lndst_scn_g")
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
mk_vi_stk <- function(sp.layer, vindx = "EVI", buff = 30, st.year = 1990, vi.thr = 1500,
                      cv.lim = 100, proj.obj = T, obj.fmt = "raster") {
  if (!inherits(sp.layer, "SpatialPolygons")) {
    stop("sp.layer isn't a SpatialPolygon* object")
  }
  require(sp)
  require(rgdal)
  require(rgeos)
  require(raster)
  # Create a list of available images 
  img.lst <- list.files(paste0("~/SIG/Geo_util/raster/arg/",
                               vindx, "_Landsat/"),
                        ".tif$", full.names = T)
  # Check projection of layer and project to measure distances
  prj.crs <- prj_str(utm_zone(sp.layer))
  if (!is.projected(sp.layer)) {
    sp.layer <- spTransform(sp.layer, prj.crs)
  }
  # Assign ownership of holes to parent polygons
  sp.comm <- createSPComment(sp.layer)
  # Create buffer of polygon for border effect
  sp.layer <- gBuffer(sp.comm, width = -buff)
  # Reproject buffered layer to WGS84
  sp.layer <- spTransform(sp.layer, geo.str)
  # Get on which landsat path rows the layer intersects
  scn.pr <- scn_pr(sp.layer)
  # Go through the images until finding the one that fully covers the layer
  vi.lst <- grep(paste(scn.pr, collapse = "|"), img.lst, value = T)
  i <- 1
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
  # Create empty bricks and data.frames to store information
  r.stk.lst <- list()
  df.lst1 <- list()
  rwn <- 1
  for (c in vi.lst) {
    # Get the year of current image
    scn.year <- as.numeric(substr(basename(c), 10, 13))
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
          # Add the current mask to the brick
          r.stk.lst[[rwn]] <- r.crp
          # Add this mask values to a reference data frame
          df.lst1[[rwn]] <- data.frame("SCN" = sub(".tif", "", basename(c)),
                                       "Year" = scn.year, "VI" = vi.mdn,
                                       stringsAsFactors = F)
          rwn <- rwn + 1
        }
      }
    }
  }
  # Check length of obtained images
  if (length(r.stk.lst) < 1) {
    stop("no images left with selected parameters")
  }
  # Create brick and data.frame from list
  r.stk <- brick(r.stk.lst)
  df1 <- do.call(rbind, df.lst1)
  # Order data.frame by year
  df1 <- df1[order(df1$Year),]
  df.lst2 <- list()
  rwn <- 1
  # The following will leave only one image per year
  if (length(unique(df1$Year)) < nrow(df1)) {
    for (d in unique(df1$Year)) {
      # Leave the one with highest median
      vi.max <- max(df1[df1$Year == d, "VI"])
      df.lst2[[rwn]] <- df1[df1$VI == vi.max,]
      rwn <- rwn + 1
    }
    # Final data.frame and brick
    df2 <- do.call(rbind, df.lst2)
    r.stk2 <- subset(r.stk, df2$SCN)
  } else {
    r.stk2 <- subset(r.stk, df1$SCN)
  }
  # Get brick name and create new ones
  nms <- names(r.stk2)
  nw.nms <- gsub("^", vindx, substr(nms, 10, 13))
  # Convert to SpatialPointsDF
  r.stk2 <- rasterToPoints(r.stk2, spatial = T)
  names(r.stk2) <- nw.nms
  r.stk2@data <- round(r.stk2@data, 2)
  # Remove points with NAs
  r.stk2 <- r.stk2[complete.cases(r.stk2@data),]
  # Convert to raster if desired
  if (obj.fmt == "raster") {
    r.stk2 <- pnt2rstr(r.stk2)
  }
  # Project object
  if (proj.obj) {
    if (inherits(r.stk2, "Raster")) {
      require(gdalUtils)
      # Generate temp file names
      tmp1 <- tempfile(fileext = ".tif")
      tmp2 <- tempfile(fileext = ".tif")
      # Write temporary raster
      writeRaster(r.stk2, filename = tmp1)
      # Project raster with cubic convolution resampling
      r.stk2 <- gdalwarp(srcfile = tmp1, dstfile = tmp2,
                         t_srs = prj.crs,
                         r = "cubic", output_Raster = T)
      names(r.stk2) <- nw.nms
    }
    if (inherits(r.stk2, "Spatial")) {
      # Project points
      r.stk2 <- spTransform(x = r.stk2,
                            CRSobj = prj.crs)
    }
  }
  return(r.stk2)
}

#Function to get and filter srtm images from selected lat/long
dem_srtm <- function(sp.layer, buff = 30, format = "point", proj.obj = T) {
  if (!inherits(sp.layer, "SpatialPolygons")) {
    stop("sp.layer isn't a SpatialPolygon* object")
  }
  require(sp)
  require(rgdal)
  require(rgeos)
  require(raster)
  prj.crs <- prj_str(utm_zone(sp.layer))
  if (buff != 0) {
    # Check projection of layer and project to measure distances
    if (!is.projected(sp.layer)) {
      sp.layer <- spTransform(sp.layer, prj.crs)
    }
    # Assign ownership of holes to parent polygons
    sp.comm <- createSPComment(sp.layer)
    # Create buffer of polygon for border effect
    sp.layer <- gBuffer(sp.comm, width = -buff)
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
  if (proj.obj) {
    if (format != "raster"){
      # Return projected SpatialPoints
      msc.pnt <- rasterToPoints(msc, spatial = T)
      msc.p <- spTransform(msc.pnt, CRSobj = prj.crs)
      names(msc.p) <- "elev"
      return(msc.p)
    } else {
      # Return projected Raster
      require(gdalUtils)
      # Generate temp file names
      tmp1 <- tempfile(fileext = ".tif")
      tmp2 <- tempfile(fileext = ".tif")
      # Write temporary raster
      writeRaster(msc, filename = tmp1)
      # Project raster with cubic convolution resampling
      msc.p <- gdalwarp(srcfile = tmp1, dstfile = tmp2, t_srs = prj.crs,
                        r = "cubic", output_Raster = T)
      names(msc.p) <- "elev"
      return(msc.p)
    }
  }
  if (format != "raster") {
    # Return SpatialPoints in Lat-Lon
    msc.pnt <- rasterToPoints(msc, spatial = T)
    names(msc.pnt) <- "elev"
    return(msc.pnt)
  }
  # Return raster in Lat-Lon
  return(msc)
}

#Function to reclassify a raster in n classes by jenks
rstr_rcls <- function(raster.lyr, n.class = 3, val = 1:n.class,
                      style = "fisher") {
  if (!inherits(raster.lyr, "Raster")){
    stop("Input object isn't a Raster* object")
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
dem_cov <- function(DEM.layer, dem.attr = "DEM", deriv = "all", smth = T, save.rst = T) {
  if (!inherits(DEM.layer, "SpatialPointsDataFrame") &
      !inherits(DEM.layer, "Raster")) {
    stop("DEM.layer isn't a SpatialPointsDataFrame or Raster* object")
  }
  require(sp)
  require(RSAGA)
  require(raster)
  # Save current directory
  curr.wd <- getwd()
  on.exit(setwd(curr.wd))
  if (save.rst) {
    if ("./DEM_derivates" %in% list.dirs()) {
      unlink("./DEM_derivates", recursive = T, force = T)
    }
    # Topography derivates folder creation
    dir.create("DEM_derivates")
    setwd("./DEM_derivates")
  } else {
    tmp.dir <- paste0(tempdir(), "\\DEM_derivates")
    if (dir.exists(tmp.dir)) {
      unlink(tmp.dir, recursive = T, force = T)
      dir.create(tmp.dir)
    } else {
      dir.create(tmp.dir)
    }
    setwd(tmp.dir)
  }
  # Store dem layer CRS
  lyr.crs <- CRS(proj4string(DEM.layer))
  # If layer is SPDF convert to raster
  if (inherits(DEM.layer, "SpatialPointsDataFrame")) {
    base.rstr <- pnt2rstr(DEM.layer, dem.attr)
    base.lyr <- DEM.layer
  } else {
    # Convert raster layer to points to add other layers
    base.rstr <- DEM.layer
    base.lyr <- rasterToPoints(DEM.layer, spatial = T)
  }
  if (smth) {
    base.rstr <- focal(base.rstr, w = matrix(1, 3, 3),
                       fun = mean, na.rm = T, pad = T)
  }
  dem.file <- "dem.tif"
  # Save raster as GeoTIFF
  writeRaster(base.rstr, dem.file, options = c("COMPRESS=NONE"))
  # All variables
  full.vars <- c("Slope", "Aspect", "Curv", "Catch_area", "CTI",
                 "Conv_Index", "LS_Factor", "SWI")
  # Get list of derivatives to calculate
  if (length(deriv) > 1) {
    deriv.lst <- deriv
  } else {
    if (deriv == "all") {
      deriv.lst <- full.vars
    } else {
      deriv.lst <- deriv
    }
  }
  # Convert GeoTIFF to Saga Grid
  rsaga.import.gdal(dem.file, show.output.on.console = F)
  # Pitremove the DEM
  system("mpiexec -n 4 PitRemove -z dem.tif -fel pitrem.tif",
         show.output.on.console = F)
  if (length(grep("slo|cti|catch|fact", deriv.lst, ignore.case = T)) > 0) {
    # DInf flow directions and slope
    system("mpiexec -n 4 DinfFlowdir -fel pitrem.tif -ang ang.tif -slp Slope.tif",
           show.output.on.console = F)
    slp <- raster("Slope.tif")
    slp <- slp + 0.00000001
    writeRaster(slp, filename = "Slope.tif", options = c("COMPRESS=NONE"),
                overwrite = T)
    # Dinf contributing area
    system("mpiexec -n 4 AreaDinf -nc -ang ang.tif -sca Catch_Area.tif",
           show.output.on.console = F)
    # Wetness Index
    system("mpiexec -n 4 SlopeAreaRatio -slp Slope.tif -sca Catch_Area.tif -sar sar.tif",
           show.output.on.console = F)
    sar <- raster("sar.tif")
    wi <- sar
    wi[,] <- -log(sar[,])
    writeRaster(wi, filename = "CTI.tif", options = c("COMPRESS=NONE"),
                overwrite = T)
  }
  if (length(grep("asp", deriv.lst, ignore.case = T)) > 0) {
    # Calculate Aspect
    rsaga.aspect("dem.sgrd", "Aspect", method = "maxtriangleslope",
                 show.output.on.console = F)
  }
  if (length(grep("curv", deriv.lst, ignore.case = T)) > 0) {
    # Calculate Curvatures
    rsaga.curvature("dem.sgrd", "Curv", method = "poly2zevenbergen",
                    show.output.on.console = F)
  }
  if (length(grep("conv", deriv.lst, ignore.case = T)) > 0) {
    # Calculate Convergence Index
    rsaga.geoprocessor("ta_morphometry", module = 1,
                       param = list(ELEVATION = "dem.sgrd",
                                    RESULT = "Conv_Index"),
                       show.output.on.console = F)
  }
  if (length(grep("ls|fac", deriv.lst, ignore.case = T)) > 0) {
    # Calculate LS Factor
    rsaga.geoprocessor("ta_hydrology", module = 25,
                       param = list(DEM = "dem.sgrd",
                                    LS_FACTOR = "LS_Factor"),
                       show.output.on.console = F)
  }
  if (length(grep("swi|saga", deriv.lst, ignore.case = T)) > 0) {
    # Calculate Saga Wetness Index
    rsaga.wetness.index("dem.sgrd", "SWI.sgrd",
                        show.output.on.console = F)
  }
  gen.deriv <- grep(paste(deriv.lst, collapse = "|"), full.vars,
                    ignore.case = T, value = T)
  # Iterate over list of grids
  base.lyr@data[dem.attr] <- extract(base.rstr, base.lyr)
  for (a in gen.deriv) {
    # Get name for tiff
    tif.name <- paste0(a, ".tif")
    if (file.exists(paste0(a, ".sgrd"))) {
      # Convert Saga grid to GeoTIFF
      rsaga.geoprocessor("io_gdal", module = 2,
                         param = list(GRIDS = paste0(a, ".sgrd"),
                                      FILE = tif.name),
                         show.output.on.console = F)
    }
    # Open tiff as RasterLayer
    terr.tif <- raster(tif.name)
    # Get data from the tiff that intersects with points
    terr.data <- extract(terr.tif, base.lyr, method = "bilinear")
    # Add column to SPDF
    base.lyr@data[a] <- terr.data
  }
  if (any(is.na(base.lyr@data))) {
    base.lyr@data <- df_impute(base.lyr@data)
  }
  file.remove(list.files(path = ".", pattern = "sdat|sgrd|mgrd|prj"))
  return(base.lyr)
}

# Defining the MBA interpolation function
int_fx <- function(base.pnts, obs.pnts, vrbl, moran = F,
                   dist = 20, clean = T, krig = F) {
  if (!inherits(base.pnts, "SpatialPoints") |
      !inherits(obs.pnts, "SpatialPointsDataFrame")) {
    stop("at least one of the inputs isn't a SpatialPoints* object")
  }
  require(sp)
  base.crs <- prj_str(utm_zone(base.pnts))
  obs.crs <- prj_str(utm_zone(obs.pnts))
  if (!identical(base.crs, obs.crs)) stop("spatial layers are in different UTM zones")
  if (!is.projected(base.pnts)) {
    base.pnts <- spTransform(base.pnts, base.crs)
  }
  if (!is.projected(obs.pnts)) {
    obs.pnts <- spTransform(obs.pnts, obs.crs)
  }
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
                           nmax = 100, newdata = base.map, model = vg.fit[[3]])
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
hyb_pp <- function(hybrid, exp.yld, step = 1235, biol = F) {
  num.param <- as.numeric(hyb.param[1:8, hybrid])
  num.param2 <- as.numeric(hyb.param[9:11, hybrid])
  pl.pop <- vector()
  min.kn.yld <- 3.5
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
    pp.min <- 30000
    pp.kn <- ((-num.param[2] + 2 * num.param[3] * num.param[7] -
                 (num.param[5] * min.kn.yld - num.param[5] * num.param[8]) -
                 2 * num.param[6] * num.param[7] * (num.param[8] - min.kn.yld)) /
                (2 * num.param[3] + 2 * num.param[6] * (min.kn.yld - num.param[8]))) * 10000
    if (!biol) {
      pp.kn <- pp.kn * ((100 - (num.param2[1] * min.kn.yld * min.kn.yld + num.param2[2] *
                                  min.kn.yld + num.param2[3]))/100)
    }
    seed.rates <- seq(pp.min + step, 150000, step)
    pp.lm <- lm(y ~ x, data = data.frame(x = c(3, min.kn.yld), y = c(pp.min, pp.kn)))
    for (a in seq_along(exp.yld)) {
      if (exp.yld[a] < 3) {
        hyb.pp <- pp.min
        hyb.corr <- hyb.pp
      }
      if (exp.yld[a] >= 3 & exp.yld[a] < min.kn.yld) {
        hyb.pp <- predict(pp.lm, newdata = data.frame(x = exp.yld[a]))
        hyb.pp <- seed.rates[which(abs(seed.rates - hyb.pp) == min(abs(seed.rates - hyb.pp)))]
        hyb.corr <- hyb.pp
      }
      if (exp.yld[a] >= min.kn.yld) {
        hyb.pp <- ((-num.param[2] + 2 * num.param[3] * num.param[7] -
                      (num.param[5] * exp.yld[a] - num.param[5] * num.param[8]) -
                      2 * num.param[6] * num.param[7] * (num.param[8] - exp.yld[a])) /
                     (2 * num.param[3] + 2 * num.param[6] * (exp.yld[a] - num.param[8]))) * 10000
        if (!biol) {
          hyb.corr <- hyb.pp * ((100 - (num.param2[1] * exp.yld[a] * exp.yld[a] + num.param2[2] *
                                          exp.yld[a] + num.param2[3]))/100)
        } else {
          hyb.corr <- hyb.pp
        }
        hyb.corr <- seed.rates[which(abs(seed.rates - hyb.corr) == min(abs(seed.rates - hyb.corr)))]
      }
      pl.pop[a] <- hyb.corr
    }
    return(pl.pop)
  }
}

presc_grid <- function(sp.layer, pred.model, hybrid, points = T,
                       fill = F, quantile = NULL, step = 1235) {
  if (!inherits(sp.layer, "SpatialPointsDataFrame")) {
    stop("sp.layer isn't a SpatialPointsDataFrame object")
  }
  require(sp)
  # If one wants the NAs can be filled
  if (fill) {
    if (any(is.na(sp.layer@data))) {
      # Impute missing data
      sp.layer@data <- df_impute(sp.layer@data)
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
  sp.poly@data <- data.frame(Hybrid = rep(hybrid, nrow(sp.poly@data)),
                             stringsAsFactors = F)
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
  if (!inherits(sp.layer, "SpatialPolygons")) {
    stop("sp.layer isn't a SpatialPolygons* object")
  }
  # If in WGS84 project to UTM
  prj.crs <- prj_str(utm_zone(sp.layer))
  if (!is.projected(sp.layer)) {
    library(rgdal)
    sp.layer <- spTransform(sp.layer, prj.crs)
  }
  require(sp)
  # Get bounding box
  lyr.bb <- sp.layer@bbox
  # Calculate regularly spaced coordinates
  grd.1 <- expand.grid(x = seq(lyr.bb[1, 1], lyr.bb[1, 2], by = dist),
                       y = seq(lyr.bb[2, 1], lyr.bb[2, 2], by = dist))
  # Create SPDF from regular coordinates
  grd.1 <- SpatialPointsDataFrame(grd.1, data = grd.1, proj4string = prj.crs)
  # Remove the ones outside the boundary
  grd.inp <- !is.na(over(grd.1, SpatialPolygons(sp.layer@polygons,
                                                proj4string = prj.crs)))
  grd.1 <- grd.1[grd.inp,]
  return(grd.1)
}

mz_smth <- function(sp.layer, area = 2500) {
  if (!inherits(sp.layer, "SpatialPolygons") & !inherits(sp.layer, "Raster")) {
    stop("sp.layer isn't a SpatialPolygons* or Raster* object")
  }
  require(rgrass7)
  require(raster)
  # If the input is a raster convert to polygons and dissolve by zone
  if (inherits(sp.layer, "Raster")) {
    sp.layer <- rstr2pol(sp.layer)
  }
  # If in WGS84 project to UTM
  prj.crs <- prj_str(utm_zone(sp.layer))
  if (!is.projected(sp.layer)) {
    library(rgdal)
    sp.layer <- spTransform(sp.layer, prj.crs)
  }
  # Check wether GRASS is running, else initialize
  if (nchar(Sys.getenv("GISRC")) == 0) {
    initGRASS(gisBase = "c:/Program Files (x86)/GRASS GIS 7.0.0",
              override = TRUE)
  }
  # Convert multipart to singlepart
  #sp.layer <- disaggregate(sp.layer)
  zm.pol <- paste0(sample(letters, 1), substr(basename(tempfile()), 9, 14))
  # Convert name 'layer' to 'Zone'
  names(sp.layer) <- sub("layer", "Zone", names(sp.layer))
  # Write GRASS vector
  writeVECT(sp.layer, zm.pol, v.in.ogr_flags = "o")
  zm.gnrl <- paste0(sample(letters, 1), substr(basename(tempfile()), 9, 14))
  # Smooth lines of polygons
#   execGRASS("v.generalize", flags = c("overwrite", "quiet"), input = zm.pol,
#             output = zm.gnrl, method = "snakes", threshold = 1)
#   zm.cln <- paste0(sample(letters, 1), substr(basename(tempfile()), 9, 14))
  # Remove small/sliver polygons
  execGRASS("v.clean", flags = c("overwrite", "quiet"), input = zm.pol,
            output = zm.gnrl, tool = "rmarea", threshold = area)
  # Read back cleaned layer
  zm.fnl <- readVECT(zm.gnrl)
  # If no CRS, define one
  if (is.na(zm.fnl@proj4string)) {
    proj4string(zm.fnl) <- prj.crs
  }
  # Remove 'cat' column from data.frame
  zm.fnl@data["cat"] <- NULL
  return(zm.fnl)
}

pnt2rstr <- function(sp.layer, field = names(sp.layer)){
  if (!inherits(sp.layer, "SpatialPointsDataFrame")) {
    stop("sp.layer isn't a SpatialPointsDataFrame")
  }
  require(sp)
  require(raster)
  # Get layer CRS
  lyr.crs <- CRS(proj4string(sp.layer))
  # Check for regular spacing
  p2g <- try(points2grid(sp.layer))
  if (class(p2g)[1] == "try-error") {
    stop("points aren't regularly spaced")
  }
  # Check if there's more than one field to convert to raster
  if (length(field) > 1) {
    # Create empty stack
    rstr.lst <- list()
    rwn <- 1
    for (a in field) {
      if (a %in% names(sp.layer@data) == F) {
        stop("field isn't an attribute in SpatialPointsDataFrame")
      }
      rstr.lst[[rwn]] <- raster(SpatialPixelsDataFrame(sp.layer, 
                                                       data = sp.layer@data[a],
                                                       proj4string = lyr.crs))
      rwn <- rwn + 1
    }
    sp.rstr <- brick(rstr.lst)
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
  if (!inherits(sp.layer, "Spatial")) {
    stop("sp.layer isn't a Spatial* object")
  }
  lyr.crs <- CRS(proj4string(sp.layer))
  # Get bbox corners
  corners <- cbind(c(sp.layer@bbox[1], sp.layer@bbox[1],
                     sp.layer@bbox[3], sp.layer@bbox[3]),
                   c(sp.layer@bbox[2], sp.layer@bbox[4],
                     sp.layer@bbox[4], sp.layer@bbox[2]))
  pol.bb <- Polygon(corners)
  pols.bb <- Polygons(list(pol.bb), "pol1")
  sp.pol <- SpatialPolygons(list(pols.bb),
                            proj4string = lyr.crs)
  # If in UTM project to WGS84
  if (is.projected(sp.pol)) {
    library(rgdal)
    sp.pol <- spTransform(sp.pol, geo.str)
  }
  require(rgeos)
  gcent <- gCentroid(sp.pol)
  coord <- gcent@coords
  names(coord) <- c("Lon", "Lat")
  coord <- coord[c(2, 1)]
  return(coord)
}

moran_cln <- function(sp.layer, vrbl, dist = 20, GM = F, LM = T) {
  if (!inherits(sp.layer, "SpatialPointsDataFrame")) {
    stop("sp.layer isn't a SpatialPointsDataFrame object")
  }
  require(spdep)
  if (!GM & !LM) {
    cat("WARNING: no cleaning performed, select GM, LM or both")
    return(sp.layer)
  }
  sp.orig <- sp.layer
  # If in WGS84 project to UTM
  prj.crs <- prj_str(utm_zone(sp.layer))
  if (!is.projected(sp.layer)) {
    library(rgdal)
    sp.layer <- spTransform(sp.layer, prj.crs)
  }
  # Remove NA's
  if (any(is.na(sp.layer@data[,vrbl]))) {
    sp.layer <- sp.layer[-which(is.na(sp.layer@data[,vrbl])),]
  }
  # Identify neighbour points by Euclidean distance
  nb.lst <- dnearneigh(sp.layer, d1 = 0, d2 = dist)
  # Get number of neighbours in the neighbours list
  nb.crd <- card(nb.lst)
  if (all(nb.crd == 0)) {
    cat("WARNING: no cleaning performed, try increasing the neighbor distance\n")
    return(sp.orig)
  }
  # Remove points with no neighbors
  spl.noznb <- subset(sp.layer, nb.crd > 0)
  # Also in neighbor list
  nb.noznb <- subset(nb.lst, nb.crd > 0)
  # Convert cleaned neighbor list to spatial weighted list
  w.mat <- nb2listw(nb.noznb, style = "W")
  # Get numerical data of variable
  vrbl.dt <- spl.noznb@data[, vrbl]
  if (GM) {
    # Compute the lag vector V x
    wx <- lag.listw(w.mat, vrbl.dt)
    # Lineal model of lagged vs observed variable
    xwx.lm <- lm(wx ~ vrbl.dt)
    # Compute regression (leave-one-out deletion) diagnostics for
    # linear model and only get get the logical influence matrix
    infl.xwx <- influence.measures(xwx.lm)[["is.inf"]]
    # Convert to numeric, 6 column matrix for computation of sums
    infl.mat <- matrix(as.numeric(infl.xwx), ncol = 6)
    # Get those rows where at least one index is TRUE
    gm.out <- which(rowSums(infl.mat) != 0)
  }
  if (LM) {
    # Calculate local moran
    lmo <- localmoran(vrbl.dt, w.mat, p.adjust.method = "bonferroni",
                      alternative = "less")
    # Convert to data.frame to select data
    lmo <- data.frame(lmo)
    # Get rows wheres indices are significative
    lm.out <- which(lmo[, "Ii"] <= 0 | lmo[, "Pr.z...0."] <= 0.05)
  }
  if (GM & LM) {
    # Get the unique rows to delete
    all.out <- unique(c(mp.out, lm.out))
  } else if(GM) {
    all.out <- gm.out
  } else if(LM) {
    all.out <- lm.out
  }
  # Remove them from SPDF
  spl.noznb <- spl.noznb[-all.out, ]
  if (!is.projected(sp.layer)) {
    sp.noznb <- spTransform(sp.layer, geo.str)
  }
  return(spl.noznb)
}

var_fit <- function(sp.layer, vrbl, cln = F, plot = F){
  if (!inherits(sp.layer, "SpatialPointsDataFrame")) {
    stop("sp.layer isn't a SpatialPointsDataFrame object")
  }
  require(automap)
  require(MASS)
  require(gstat)
  prj.crs <- prj_str(utm_zone(sp.layer))
  if (!is.projected(sp.layer)) {
    library(rgdal)
    sp.layer <- spTransform(sp.layer, prj.crs)
  }
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
  if (class(fit.res)[1] == "character") {
    stop("a variogram could not be fitted")
  }
  if (!is.null(attr(fit.res, "singular"))) {
    if (attr(fit.res, "singular")) {
      stop("a variogram could not be fitted")
    }
  }
  if (fit.res$psill[1] >= fit.res$psill[2]) {
    stop("a variogram could not be fitted")
  }
  fit.res <- fit.variogram(vrbl.vgm, vgm(isill, sh.mod, irange, inug), 
                           fit.sills = T, fit.ranges = T, warn.if.neg = T)
  vrbl.fit <- fit.res
  # Wether to plot the variograms (exp & fit)
  if (plot) {
    print(
      plot(vrbl.vgm, model = vrbl.fit, pch = 19, cex = 2, lty = 2, lwd = 2,
           col = "red", xlab = "Distance", ylab = "Semivariance",
           main = paste0("Emp. & fitted semivariogram for ", vrbl.nm))
    )
  }
  #Return a list with: lambda, experimental and fitted variogram
  return(list(lmbd, vrbl.vgm, vrbl.fit))
}

kmz_sv <- function(sp.layer = interp.rp, spz = spz, Rev = F){
  require(plotKML)
  kmz.name <- paste0(basename(getwd()), '_Reporte_TDCA')
  rstr_lyr <- pnt2rstr(sp.layer, c('DEM', 'SWI', 'EC30', 'EC90', 'OM', 'CEC'))
#   for (pol in seq_along(spz@polygons)) {
#     if (length(spz@polygons[[pol]]@plotOrder) > 1) {
#       pol1 <- spz@polygons[[1]]
#       row1 <- spz@data[1,]
#       spz@data[1, 1] <- spz@data[pol, 1]
#       spz@data[1, 2] <- spz@data[pol, 2]
#       spz@data[pol, 1] <- row1[1, 1]
#       spz@data[pol, 2] <- row1[1, 2]
#       spz@polygons[[1]] <- spz@polygons[[pol]]
#       spz@polygons[[pol]] <- pol1
#     }
#   }
  kml_open(file.name = paste0(kmz.name, '.kml'),
           folder.name = kmz.name, 
           overwrite = T)
  if(Rev == F){
    kml_layer(spz, 
              raster.name = 'Calidad de Sitio',
              subfolder.name='Calidad de Sitio (s/u)',
              colour = Calidad,
              colour_scale = cols(3),
              outline = F)
  } else {
    kml_layer(spz, 
              raster.name = 'Calidad de Sitio',
              subfolder.name='Calidad de Sitio (s/u)',
              colour = Calidad,
              colour_scale = rev(cols(3)),
              outline = F)
  }
  kml_layer(rstr_lyr, 
            raster.name = 'DEM',
            subfolder.name='DEM (m)',
            colour = 'DEM',
            colour_scale = elev_cols(255),
            plot.legend = T)
  kml_layer(rstr_lyr, 
            raster.name = 'Indice de Humedad',
            subfolder.name='IH (s/u)',
            colour = 'SWI',
            colour_scale = swi_cols(255))
  kml_layer(rstr_lyr, 
            raster.name = 'EC30',
            subfolder.name='ECa Superficial (mS/m)',
            colour = 'EC30',
            colour_scale = ec_cols(255))
  kml_layer(rstr_lyr, 
            raster.name = 'EC90',
            subfolder.name='ECa Profunda (mS/m)',
            colour = 'EC90',
            colour_scale = ec_cols(255))
  kml_layer(rstr_lyr, 
            raster.name = 'MO',
            subfolder.name='MO (%)',
            colour = 'OM',
            colour_scale = om_cols(255))
  kml_layer(rstr_lyr, 
            raster.name = 'CIC',
            subfolder.name='CIC (meq/100g)',
            colour = 'CEC',
            colour_scale = cec_cols(255))
  
  kml_close(file.name = paste0(kmz.name, '.kml'))
  
  file = paste(kmz.name, '.zip', sep = '')
  zip(file, list.files(".", pattern = "png|kml"))
  file.rename(paste0(kmz.name, ".zip"), paste0(kmz.name, ".kmz"))
  file.remove(list.files(pattern = '.png|kml'))
}

# Import Veris data
veris_import <- function(vrs.fl = 'VSECOM', vrbl = c('EC30', 'EC90', 'Red', 'IR'), 
                         clean = T, vrbl.cl = c('Speed', 'Depth', 'OM_ratio'), fl.nm = 'veris'){
  fls <- list.files(path = "./Veris/", pattern = paste0(vrs.fl, '.*\\.txt'))
  for (b in fls) {
    veris.n <- read.table(paste0("Veris/", b), header = TRUE, sep = "\t", skipNul = T)
    if (exists("veris")) {
      veris <- rbind(veris, veris.n)
    } else {
      veris <- veris.n
    }
  }
  ec30.nm <- grep("EC_SH|EC.SH|EC30|EC SH", names(veris), ignore.case = T, value = T)
  ec90.nm <- grep("EC_DP|EC.DP|EC90|EC DP", names(veris), ignore.case = T, value = T)
  if (length(ec30.nm) == 1) { names(veris) <- sub(ec30.nm, "EC30", names(veris)) }
  if (length(ec90.nm) == 1) { names(veris) <- sub(ec90.nm, "EC90", names(veris)) }
  # cleaning by data distribution
  if (clean) {
    require(robustbase)
    for (a in vrbl.cl){
      vrbl.data <- veris[, a]
      # Outlier detection
      vrbl.out <- adjboxStats(vrbl.data)$out
      # If there are outliers they will be removed
      if (length(vrbl.out)>0) {
        veris <- veris[!(vrbl.data %in% vrbl.out),]
      }
    }
  }
  require(sp)
  require(rgdal)
  veris.pnt <- SpatialPointsDataFrame(coords = veris[,c('Long','Lat')],
                                      proj4string = geo.str, 
                                      data = veris[,vrbl])
  prj.crs <- prj_str(utm_zone(veris.pnt))
  veris.pnt <- spTransform(veris.pnt, prj.crs)
  veris.pnt <- remove.duplicates(veris.pnt)
  write_shp(veris.pnt, paste0('Veris/', fl.nm), overwrite = T)
  rm(veris, veris.n)
  return(veris.pnt)
}

elev_import <- function(path = 'Elev') {
  require(rgeos)
  elev <- read_shp(paste0(path, '/', list.files(path, pattern = '.shp$')))
  if (inherits(elev, "SpatialPolygonsDataFrame")) {
    elev.cnt <- gCentroid(elev, byid = T)
    elev.df <- over(elev.cnt, elev)
    elev <- SpatialPointsDataFrame(coords = elev.cnt, data = elev.df,
                                   proj4string = CRS(proj4string(elev)))
  }
  prj.crs <- prj_str(utm_zone(elev))
  if (!is.projected(elev)) {
    elev <- spTransform(elev, prj.crs)
  }
  return(elev)
}

soil_import <- function(path = 'Soil') {
  soil.df <- read.table(paste0(path, '/', list.files(path, pattern = '.txt$')),
                     header = TRUE, sep = "\t", skipNul = T)
  sp.soil <- SpatialPointsDataFrame(coords = soil.df[,c('Long', 'Lat')], 
                                 data = soil.df, proj4string = geo.str)
  prj.crs <- prj_str(utm_zone(sp.soil))
  sp.soil <- spTransform(sp.soil, CRSobj = prj.crs)
  writeOGR(sp.soil, overwrite_layer = T, dsn = "./Soil", driver = "ESRI Shapefile",
           layer = "soil")
  return(sp.soil)
}
  
var_cal <- function(sp.layer, var = 'OM', pdf = T, width = 10, soil = 'soil'){
  require(rgdal)
  require(rgeos)
  require(ggplot2)
  require(reshape2)
  require(gridExtra)
  soil <- read_shp(paste0('Soil/', soil, '.shp'))
  proj4string(soil) <- proj4string(sp.layer)
  # Create buffer of 10 m
  soil.buf <- gBuffer(soil, byid = TRUE, width = width)
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
    label <- paste('Model NÂ°', lm.summ$model[i], ": ", var, 
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

trat_grd <- function(sp.layer, largo = 10, ancho, ang = 0, n.trat,
                     n.pas = 1, random = T) {
  require(rgdal)
  require(rgeos)
  require(maptools)
  # Check if sp.layer is a spatialpolygon
  if (!inherits(sp.layer, "SpatialPolygons")) {
    stop("sp.layer isn't a SpatialPolygons* object")
  }
  prj.crs <- prj_str(utm_zone(sp.layer))
  if (!is.projected(sp.layer)) {
    library(rgdal)
    sp.layer <- spTransform(sp.layer, prj.crs)
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
  ntrat <- n.trat
  # Creation of the vector of treatments
  trt.ordr <- vector()
  for (b in 1:ceiling(nc / ntrat)) {
    if (random) {
      trt <- sample(ntrat, ntrat)
    } else {
      trt <- 1:ntrat
    }
    trt.ordr <- c(trt.ordr, trt)
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
  # Creation of column identifiers
  col.mat <- mat
  for (i in 1:ncol(mat)) {
    col.mat[, i] <- i
  }
  col.t <- t(col.mat)
  col.v <- as.vector(col.t)
  # Creation of row identifiers
  row.mat <- mat
  for (i in 1:ncol(mat)) {
    row.mat[, i] <- 1:nrow(mat)
  }
  row.t <- t(row.mat)
  row.v <- as.vector(row.t)
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
  cell.offset <- c(gCentroid(bound)@coords[1] - (cell.size[1] * nc) / 2,
                   gCentroid(bound)@coords[2] - (cell.size[2] * nr) / 2)
  # Creation of a rectangular grid of defined dimensions
  grd <- GridTopology(cell.offset, cell.size, c(nc, nr))
  # Conversion to spatial polygons
  pol.1 <- as.SpatialPolygons.GridTopology(grd, proj4string = prj.crs)
  # Extraction of the IDs of all the polygons
  pol.lst <- list()
  for (a in seq_along(pol.1)) {
    pol.lst[a] <- slot(pol.1[a,]@polygons[[1]], "ID") 
  }
  # Creation of the data frame of the polygons
  data <- data.frame(Col = col.v, Row = row.v, Trat = trt.v, Rep = rep.v)
  row.names(data) <- unlist(pol.lst)
  # Adding the data frame to the polygons
  pol.2 <- SpatialPolygonsDataFrame(pol.1, data = data, match.ID = T)
  proj4string(pol.2) <- prj.crs
  if (ang > 0) {
    # Rotation of the polygons by the defined angle
    pol.3 <- elide(pol.2, rotate = ang,
                   center = gCentroid(bound)@coords)
    proj4string(pol.3) <- prj.crs
  } else {
    pol.3 <- pol.2
  }
  # Clipping of the polygon with the boundary
  pol.4 <- gIntersection(bound, pol.3, byid = T)
  # Creation of spatialpoints to join the attribute table with the clipped polygon
  pol.3.pnt <- SpatialPointsDataFrame(coordinates(pol.3),
                                      pol.3@data,
                                      proj4string = prj.crs)
  pol.3.pnt <- gBuffer(pol.3.pnt,
                       width = (min(cell.size) - 0.001) / 2,
                       byid = T)
  # Extraction of the data frame rows that match with the clipped polygons
  df <- over(pol.4, pol.3.pnt)
  # Final spatialpolygondf with attribute table
  pol.5 <- SpatialPolygonsDataFrame(pol.4, data = df, match.ID = F)
  pol.6 <- pol.5[!is.na(pol.5@data$Col),]
  gc()
  return(pol.5)
}

multi_mz <- function(sp.layer, vrbls = c("DEM", "Aspect", "CTI", "Slope",
                                         "SWI", "EC30", "EC90", "OM",
                                         "CEC", "EVI_mean"),
                     n.mz = 3, dist = 20, plot = F, sp.layer2 = bound.shp, 
                     area = 3000, style = 'fisher') {
  if (!inherits(sp.layer, "SpatialPointsDataFrame")) {
    stop("sp.layer isn't a SpatialPointsDataFrame object")
  }
  library(ade4)
  library(spdep)
  if (any(is.na(sp.layer@data[vrbls]))) {
    sp.layer@data <- df_impute(sp.layer@data[vrbls])
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
  rast <- disaggregate(pnt2rstr(sp.pca, "CS1"), fact = 5, method = 'bilinear')
  pca.rast <- rstr_rcls(rast, n.class = n.mz, val = 1:n.mz, style)
  pca.rast <- mask(pca.rast, sp.layer2)
  # Cleaning of the managemnent zones
  sp.pol <- mz_smth(pca.rast, area)
  print(cs1)
  return(sp.pol)
}

# Function to read shapefiles with proj info
read_shp <- function(shp.file) {
  require(rgdal)
  require(maptools)
  # Get path and file name from object
  if (dirname(shp.file) == ".") {
    dsn <- "."
    fl.nm <- shp.file
  } else {
    dsn <- dirname(shp.file)
    fl.nm <- basename(shp.file)
  }
  # Check if it has extension  for the two functions
  if (length(grep(".shp$", fl.nm, ignore.case = T)) == 1) {
    layer.ogr <- sub(".shp", "", fl.nm)
    fl.nm2 <- shp.file
  } else {
    layer.ogr <- fl.nm
    fl.nm2 <- paste0(shp.file, ".shp")
  }
  # Get projection information from layer
  shp.info <- ogrInfo(dsn, layer.ogr)[["p4s"]]
  # Load the shapefile with the defined parameters
  shp.sp <- readShapeSpatial(fn = fl.nm2, proj4string = CRS(shp.info),
                             verbose = F, delete_null_obj = T)
  return(shp.sp)
}

# Read kmz/kml into Spatial*DataFrame
read_kmz <- function(kmz.file) {
  require(tools)
  require(rgdal)
  # Define temporary directory
  tmp.dir <- tempdir()
  # Check if it is kmz or kml
  if (file_ext(kmz.file) == "kmz") {
    # If a kmz, extract it
    unzip(kmz.file, exdir = tmp.dir)
  } else {
    file.copy(kmz.file, tmp.dir)
  }
  # Get kml working file
  wrk.fl <- list.files(tmp.dir, ".kml$", full.names = T)
  # Get name/s of layer/s in KML
  lyr <- ogrListLayers(wrk.fl)
  # Read the layer/s in the KML into a Spatial*DataFrame
  sp.lyr <- readOGR(wrk.fl, layer = lyr, verbose = F,
                    stringsAsFactors = F)
  # Delete de KML
  file.remove(wrk.fl)
  return(sp.lyr)
}

# Wrapper around writeOGR for simplification
write_shp <- function(sp.layer, file.name, overwrite = F) {
  if (!inherits(sp.layer, "Spatial")) {
    stop("sp.layer isn't a Spatial* object")
  }
  require(rgdal)
  # Check if it has extension  for the two functions
  if (length(grep(".shp$", file.name, ignore.case = T)) == 1) {
    fl.nm2 <- sub(".shp", "", file.name)
  } else {
    fl.nm2 <- file.name
  }
  writeOGR(obj = sp.layer, dsn = dirname(file.name), layer = basename(fl.nm2),
           overwrite_layer = overwrite, driver = "ESRI Shapefile")
}

# Function to convert raster to polygons
rstr2pol <- function(raster) {
  require(rgdal)
  require(RSAGA)
  if (!inherits(raster, "Raster")) {
    stop("sp.layer isn't a Raster* object")
  }
  # Create temporary files
  tmp.rstr <- tempfile(fileext = ".tif")
  writeRaster(raster, tmp.rstr)
  tmp.sgrd <- tempfile(fileext = ".sgrd")
  tmp.shp <- tempfile(fileext = ".shp")
  # Convert tif to sgrd for SAGA
  rsaga.import.gdal(tmp.rstr, tmp.sgrd, show.output.on.console = F)
  # Convert grid to polygons
  rsaga.geoprocessor("shapes_grid", module = 6,
                     param = list(GRID = tmp.sgrd,
                                  POLYGONS = tmp.shp,
                                  SPLIT = 1),
                     show.output.on.console = F)
  # Read the generated shapefile
  sp.pol <- read_shp(tmp.shp)
  # Leave in the data frame only only zone information
  sp.pol@data <- data.frame(layer = as.numeric(sp.pol$NAME),
                            stringsAsFactors = F)
  return(sp.pol)
}

# Report creation
report_tdec <- function(bound = bound.shp, veris = interp.rp, spz = spz, 
                        vrbl.sl = vrbls, zoom = 16, Rev = F){ 
  require(rgdal)
  require(rgeos)
  require(ggplot2)
  require(gridExtra)
  require(ggmap)
  require(png)
  require(RGraphics)
  require(jpeg)
  require(maptools)
  require(plyr)
  require(tools)
  # Load boundary
  # bound <- read_shp('Boundaries/boundaries')
  bound <- spTransform(bound.shp, geo.str)
  data <- fortify(bound)
  # Load Veris data
  # veris <- read_shp('interp_db_5m')
  # Create Multi-Spati zones
  # spz <- multi_mz(sp.layer = veris, plot = T, n.mz = 3, dist = 20, 
  #                 vrbls = c('EC30', 'EC90', 'DEM', 'SWI', 'OM', 'CEC'))
  spz@data$id <- rownames(spz@data)
  attr <- as.data.frame(spz)
  spz.df <- ldply(spz@polygons, fortify)
  spz.df <- cbind(spz.df,attr[as.character(spz.df$id),])
  # Add lat long coordinates to veris@data data frame
  veris <- spTransform(veris, geo.str)
  veris$lat <- veris@coords[,2]
  veris$long <- veris@coords[,1]
  prj.crs <- prj_str(utm_zone(veris))
  veris <- spTransform(veris, prj.crs)
  # Load soil data
  soil <- read_shp('Soil/soil')
  soil <- spTransform(soil, prj.crs)
  soil$lat <- soil@coords[,2]
  soil$long <- soil@coords[,1]
  soil@data$Muestra <- 1:dim(soil@data)[1]
  # Soil samples data
  # vrbl.sl <- c('Muestra', 'OM', 'pH', 'NO3', 'P', 'K', 'Na', 'Zn','CEC')
  # vrbl.sl <- c('Muestra', 'OM', 'pH', 'N', 'P', 'K', 'Na', 'Zn','CEC')
  col.nm <- c('Muestra', 'MO (%)', 'pH', 'N-NO3 (ppm)', 'P (ppm)', 'K (meq/100g)', 
              'Na (meq/100g)', 'Zn (ppm)', 'CIC (meq/100g)')
#   col.nm <- c('Muestra', 'MO (%)', 'pH', 'N-NO3 (ppm)', 'P (ppm)', 'K (ppm)', 
#               'Na (ppm)', 'Zn (ppm)', 'CIC (meq/100g)')
  var.nm <- match(vrbl.sl, names(soil@data))
  soil.tbl <- tableGrob(format(soil@data[var.nm], digits = 3, scientific = F, big.mark = ","), 
                        cols = col.nm, core.just = "center", col.just = 'center', 
                        gpar.coretext = gpar(fontsize = 11), 
                        gpar.coltext = gpar(fontsize = 10, fontface = 'bold'), 
                        show.rownames = F, h.even.alpha = 0,
                        gpar.rowtext = gpar(col = "black", cex = 0.7, equal.width = TRUE,
                                            show.vlines = TRUE, show.hlines = TRUE, separator = "grey"))
  #grid.draw(soil.tbl) 
  # Text
  txt1 <- "Reporte de Caracterizacion Ambiental"
  t1 <- textGrob(txt1, gp = gpar(fontsize = 20, fontface = 'bold'), just = 'center')
  Prod <- basename(dirname(getwd()))
  Lote <- basename(getwd())
  txt2 <-gsub('_', ' - ', gsub('-', ' ', Lote))
  t2 <- textGrob(txt2, gp = gpar(fontsize = 16, fontface = 'italic'), just = 'center')
  # Logos
  lg1 <- readJPEG("C:/AGG/utils/monsanto/td_logo.jpg")
  lg2 <- readPNG("C:/AGG/utils/monsanto/logo.png")
  df <- data.frame(x=sample(1:64, 1000, replace=T), y=sample(1:64, 1000, replace=T))
  l1 <- ggplot(df, aes(x,y)) + 
    annotation_custom(rasterGrob(lg1), -Inf, Inf, -Inf, Inf) + theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_line(colour = NA),
          panel.grid.minor = element_line(colour = NA),
          plot.margin = unit(c(1, 5, 1, 0), 'lines'))
  l2 <- ggplot(df, aes(x,y)) + 
    annotation_custom(rasterGrob(lg2), 
                      -Inf, Inf, -Inf, Inf) + theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_line(colour = NA),
          panel.grid.minor = element_line(colour = NA),
          plot.margin = unit(c(1, 0, 1, 5), 'lines'))
  l3 <- ggplot(df, aes(x,y)) + 
    annotation_custom(rasterGrob(lg1), 
                      -Inf, Inf, -Inf, Inf) + theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_line(colour = NA),
          panel.grid.minor = element_line(colour = NA))
  l4 <- ggplot(df, aes(x,y)) + 
    annotation_custom(rasterGrob(lg2), 
                      -Inf, Inf, -Inf, Inf) + theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_line(colour = NA),
          panel.grid.minor = element_line(colour = NA))
  # Blank plot
  blankPlot <- ggplot()+geom_blank(aes(1,1))+
    theme(plot.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks = element_blank())
  # Create GE figure
  gmap <- get_map(location = bound@polygons[[1]]@labpt, maptype = "satellite", zoom) 
  gm1 <- ggmap(gmap) +
    geom_polygon(data = data, aes(x = long, y = lat, group = group),
                 colour = 'white', fill = 'black', alpha = .3, size = .5) +
    scale_x_continuous(limits = c(bound@bbox[1,1], bound@bbox[1,2]), expand = c(0.001, 0.001)) +
    scale_y_continuous(limits = c(bound@bbox[2,1], bound@bbox[2,2]), expand = c(0.001, 0.001)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  # Create plots (DEM, SWI, EC30, EC90, MO, CIC, Zones)
  p1 <- ggplot() +
    geom_raster(data = veris@data, aes(x, y, fill = DEM)) +
    coord_equal() + labs(x = 'Longitud', y = 'Latitud', fill = 'DEM', title = 'Altimetria') +
    scale_fill_gradientn(colours = elev_cols(255), name = 'Altura (m)') + theme_bw() + 
    theme(plot.title = element_text(size = 16, face = 'bold'),
          legend.position = 'bottom',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          legend.key.height = unit(0.3, units = 'cm'),
          legend.key.width = unit(1, units = 'cm')) +
    guides(fill=guide_colourbar(title.position = 'top', title.hjust = 0.5)) +
    scale_y_continuous(breaks=seq(min(veris$y), max(veris$y), length = 5),
                       labels=c(round(seq(min(veris$lat), max(veris$lat), length = 5),3))) +
    scale_x_continuous(breaks=seq(min(veris$x), max(veris$x), length = 5),
                       labels=c(round(seq(min(veris$long), max(veris$long), length = 5),3)))
  p2 <- ggplot() +
    geom_raster(data = veris@data, aes(x, y, fill = SWI)) +
    coord_equal() + labs(x = 'Longitud', y = 'Latitud', fill = 'SWI', title = 'Indice de Humedad') +
    scale_fill_gradientn(colours = swi_cols(255), name = 'IH (s/u)') + theme_bw() + 
    theme(plot.title = element_text(size = 16, face = 'bold'),
          legend.position = 'bottom',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          legend.key.height = unit(0.3, units = 'cm'),
          legend.key.width = unit(1, units = 'cm')) +
    guides(fill=guide_colourbar(title.position = 'top', title.hjust = 0.5)) +
    scale_y_continuous(breaks=seq(min(veris$y), max(veris$y), length = 5),
                       labels=c(round(seq(min(veris$lat), max(veris$lat), length = 5),3))) +
    scale_x_continuous(breaks=seq(min(veris$x), max(veris$x), length = 5),
                       labels=c(round(seq(min(veris$long), max(veris$long), length = 5),3)))
  p3 <- ggplot() +
    geom_raster(data = veris@data, aes(x, y, fill = EC30)) +
    coord_equal() + labs(x = 'Longitud', y = 'Latitud', fill = 'EC30', title = 'EC Superficial') +
    scale_fill_gradientn(colours = ec_cols(255), name = 'ECs (mS/m)') + theme_bw() + 
    theme(plot.title = element_text(size = 16, face = 'bold'),
          legend.position = 'bottom',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          legend.key.height = unit(0.3, units = 'cm'),
          legend.key.width = unit(1, units = 'cm')) +
    guides(fill=guide_colourbar(title.position = 'top', title.hjust = 0.5)) +
    scale_y_continuous(breaks=seq(min(veris$y), max(veris$y), length = 5),
                       labels=c(round(seq(min(veris$lat), max(veris$lat), length = 5),3))) +
    scale_x_continuous(breaks=seq(min(veris$x), max(veris$x), length = 5),
                       labels=c(round(seq(min(veris$long), max(veris$long), length = 5),3)))
  p4 <- ggplot() +
    geom_raster(data = veris@data, aes(x, y, fill = EC90)) +
    coord_equal() + labs(x = 'Longitud', y = 'Latitud', fill = 'EC90', title = 'EC Profunda') +
    scale_fill_gradientn(colours = ec_cols(255), name = 'ECp (mS/m)') + theme_bw() + 
    theme(plot.title = element_text(size = 16, face = 'bold'),
          legend.position = 'bottom',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          legend.key.height = unit(0.3, units = 'cm'),
          legend.key.width = unit(1, units = 'cm')) +
    guides(fill=guide_colourbar(title.position = 'top', title.hjust = 0.5)) +
    scale_y_continuous(breaks=seq(min(veris$y), max(veris$y), length = 5),
                       labels=c(round(seq(min(veris$lat), max(veris$lat), length = 5),3))) +
    scale_x_continuous(breaks=seq(min(veris$x), max(veris$x), length = 5),
                       labels=c(round(seq(min(veris$long), max(veris$long), length = 5),3)))
  p5 <- ggplot() +
    geom_raster(data = veris@data, aes(x, y, fill = OM)) +
    coord_equal() + labs(x = 'Longitud', y = 'Latitud', fill = 'OM', title = 'Materia Organica') +
    scale_fill_gradientn(colours = om_cols(255), name = 'MO (%)') + theme_bw() + 
    theme(plot.title = element_text(size = 16, face = 'bold'),
          legend.position = 'bottom',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          legend.key.height = unit(0.3, units = 'cm'),
          legend.key.width = unit(1, units = 'cm')) +
    guides(fill=guide_colourbar(title.position = 'top', title.hjust = 0.5)) +
    scale_y_continuous(breaks=seq(min(veris$y), max(veris$y), length = 5),
                       labels=c(round(seq(min(veris$lat), max(veris$lat), length = 5),3))) +
    scale_x_continuous(breaks=seq(min(veris$x), max(veris$x), length = 5),
                       labels=c(round(seq(min(veris$long), max(veris$long), length = 5),3)))
  p6 <- ggplot() +
    geom_raster(data = veris@data, aes(x, y, fill = CEC)) +
    coord_equal() + labs(x = 'Longitud', y = 'Latitud', fill = 'CEC', title = 'CIC') +
    scale_fill_gradientn(colours = cec_cols(255), name = 'CIC (meq/100g)') + theme_bw() + 
    theme(plot.title = element_text(size = 16, face = 'bold'),
          legend.position = 'bottom',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          legend.key.height = unit(0.3, units = 'cm'),
          legend.key.width = unit(1, units = 'cm')) +
    guides(fill=guide_colourbar(title.position = 'top', title.hjust = 0.5)) +
    scale_y_continuous(breaks=seq(min(veris$y), max(veris$y), length = 5),
                       labels=c(round(seq(min(veris$lat), max(veris$lat), length = 5),3))) +
    scale_x_continuous(breaks=seq(min(veris$x), max(veris$x), length = 5),
                       labels=c(round(seq(min(veris$long), max(veris$long), length = 5),3)))
  y.vals <- seq(min(veris$y), max(veris$y), length = 5) 
  y.labs <- c(round(seq(min(veris$lat), max(veris$lat), length = 5),3))
  x.vals <- seq(min(veris$x), max(veris$x), length = 5) 
  x.labs <- c(round(seq(min(veris$long), max(veris$long), length = 5),3))
  pnts <- list('sp.points', soil, pch=19, cex=.8, col = 'black')
  pnts.lbs <- list('sp.pointLabel', soil, label = soil@data$Muestra, cex=1.2, col = 'black')
  if(Rev == T){
    lut <- rev(cols(3))
    lbls <- c('Alta', 'Media', 'Baja')
  } else {
    lut <- cols(3)
    lbls <- c('Baja', 'Media', 'Alta')
  }
  p7 <- spplot(spz, 'Zone', col.regions = lut, cuts = 2,
               scales = list(y=list(at = y.vals, labels = y.labs), x = list(at = x.vals, labels = x.labs)),
               sp.layout = list(pnts, pnts.lbs),
               colorkey = list(labels=list(at = c(1,2,3), labels = lbls), 
                               height = 0.15))
  # Create histogram plots (DEM, SWI, EC30, EC90, MO, CIC, Zones)
  title <- paste0('Min:', format(min(veris@data$DEM), digits = 1, nsmall = 1),
                  ' / Median:', format(median(veris@data$DEM), digits = 1, nsmall = 1),
                  ' / Mean:', format(mean(veris@data$DEM), digits = 1, nsmall = 1),
                  ' / Max:', format(max(veris@data$DEM), digits = 1, nsmall = 1))
  h1 <- ggplot(veris@data, aes(x = DEM)) + 
    geom_histogram(fill="cornsilk", colour="grey60", size=.2) +
    theme_bw() +
    labs(x = "Altura (m)", y = "N° de observaciones", title = title) +
    theme(title = element_text(size = 8),
          axis.text = element_text(size = 10),
          axis.title.x = element_text(size = 12, face = 'bold'),
          axis.title.y = element_text(size = 12, face = 'bold'))
  title <- paste0('Min:', format(min(veris@data$SWI), digits = 1, nsmall = 1),
                  ' / Median:', format(median(veris@data$SWI), digits = 1, nsmall = 1),
                  ' / Mean:', format(mean(veris@data$SWI), digits = 1, nsmall = 1),
                  ' / Max:', format(max(veris@data$SWI), digits = 1, nsmall = 1))
  h2 <- ggplot(veris@data, aes(x = SWI)) + 
    geom_histogram(fill="cornsilk", colour="grey60", size=.2) +
    theme_bw() +
    labs(x = "Indice de Humedad", y = "N° de observaciones", title = title) +
    theme(title = element_text(size = 8),
          axis.text = element_text(size = 10),
          axis.title.x = element_text(size = 12, face = 'bold'),
          axis.title.y = element_text(size = 12, face = 'bold'))
  title <- paste0('Min:', format(min(veris@data$EC30), digits = 1, nsmall = 1),
                  ' / Median:', format(median(veris@data$EC30), digits = 1, nsmall = 1),
                  ' / Mean:', format(mean(veris@data$EC30), digits = 1, nsmall = 1),
                  ' / Max:', format(max(veris@data$EC30), digits = 1, nsmall = 1))
  h3 <- ggplot(veris@data, aes(x = EC30)) + 
    geom_histogram(fill="cornsilk", colour="grey60", size=.2) +
    theme_bw() +
    labs(x = "ECs (mS/m)", y = "N° de observaciones", title = title) +
    theme(title = element_text(size = 8),
          axis.text = element_text(size = 10),
          axis.title.x = element_text(size = 12, face = 'bold'),
          axis.title.y = element_text(size = 12, face = 'bold'))
  title <- paste0('Min:', format(min(veris@data$EC90), digits = 1, nsmall = 1),
                  ' / Median:', format(median(veris@data$EC90), digits = 1, nsmall = 1),
                  ' / Mean:', format(mean(veris@data$EC90), digits = 1, nsmall = 1),
                  ' / Max:', format(max(veris@data$EC90), digits = 1, nsmall = 1))
  h4 <- ggplot(veris@data, aes(x = EC90)) + 
    geom_histogram(fill="cornsilk", colour="grey60", size=.2) +
    theme_bw() +
    labs(x = "ECp (mS/m)", y = "N° de observaciones", title = title) +
    theme(title = element_text(size = 8),
          axis.text = element_text(size = 10),
          axis.title.x = element_text(size = 12, face = 'bold'),
          axis.title.y = element_text(size = 12, face = 'bold'))
  title <- paste0('Min:', format(min(veris@data$OM), digits = 1, nsmall = 1),
                  ' / Median:', format(median(veris@data$OM), digits = 1, nsmall = 1),
                  ' / Mean:', format(mean(veris@data$OM), digits = 1, nsmall = 1),
                  ' / Max:', format(max(veris@data$OM), digits = 1, nsmall = 1))
  h5 <- ggplot(veris@data, aes(x = OM)) + 
    geom_histogram(fill="cornsilk", colour="grey60", size=.2) +
    theme_bw() +
    labs(x = "MO (%)", y = "N° de observaciones", title = title) +
    theme(title = element_text(size = 8),
          axis.text = element_text(size = 10),
          axis.title.x = element_text(size = 12, face = 'bold'),
          axis.title.y = element_text(size = 12, face = 'bold'))
  title <- paste0('Min:', format(min(veris@data$CEC), digits = 1, nsmall = 1),
                  ' / Median:', format(median(veris@data$CEC), digits = 1, nsmall = 1),
                  ' / Mean:', format(mean(veris@data$CEC), digits = 1, nsmall = 1),
                  ' / Max:', format(max(veris@data$CEC), digits = 1, nsmall = 1))
  h6 <- ggplot(veris@data, aes(x = CEC)) + 
    geom_histogram(fill="cornsilk", colour="grey60", size=.2) +
    theme_bw() +
    labs(x = "CIC (meq/100g)", y = "N° de observaciones", title = title) +
    theme(title = element_text(size = 8),
          axis.text = element_text(size = 10),
          axis.title.x = element_text(size = 12, face = 'bold'),
          axis.title.y = element_text(size = 12, face = 'bold'))
  # PDF creation
  pdf(paste0(Lote, '_Reporte_TDCA.pdf'), paper = "letter", height = 0, width = 0)
  grid.arrange(arrangeGrob(l3, ncol=1),
               arrangeGrob(textGrob(txt2, gp = gpar(fontsize = 26, fontface = 'bold'),
                                    just = 'center'), ncol = 1),
               arrangeGrob(textGrob(txt1, gp = gpar(fontsize = 18, fontface = 'bold'),
                                    just = 'center'), ncol = 1),
               arrangeGrob(gm1, ncol=1),
               arrangeGrob(l4, ncol=1),
               nrow = 5, heights = c(1/10, 1/10, 1.5/10, 5.5/10, 1/10))
  grid.arrange(arrangeGrob(l1, blankPlot, l2, ncol=3),
               arrangeGrob(t1, ncol = 1),
               arrangeGrob(t2, ncol=1),
               arrangeGrob(p1, p2, ncol=2),
               arrangeGrob(h1, h2, ncol = 2),
               arrangeGrob(textGrob('Pagina 1', gp = gpar(fontsize = 10), just = 'center'), ncol = 1),
               nrow = 6, heights = c(1/10, 0.5/10, 0.5/10, 5/10, 2.5/10, 0.5/10))
  grid.arrange(arrangeGrob(l1, blankPlot, l2, ncol=3),
               arrangeGrob(t1, ncol = 1),
               arrangeGrob(t2, ncol=1),
               arrangeGrob(p3, p4, ncol=2),
               arrangeGrob(h3, h4, ncol = 2),
               arrangeGrob(textGrob('Pagina 2', gp = gpar(fontsize = 10), just = 'center'), ncol = 1),
               nrow = 6, heights = c(1/10, 0.5/10, 0.5/10, 5/10, 2.5/10, 0.5/10))
  grid.arrange(arrangeGrob(l1, blankPlot, l2, ncol=3),
               arrangeGrob(t1, ncol = 1),
               arrangeGrob(t2, ncol=1),
               arrangeGrob(p5, p6, ncol=2),
               arrangeGrob(h5, h6, ncol = 2),
               arrangeGrob(textGrob('Pagina 3', gp = gpar(fontsize = 10), just = 'center'), ncol = 1),
               nrow = 6, heights = c(1/10, 0.5/10, 0.5/10, 5/10, 2.5/10, 0.5/10))
  grid.arrange(arrangeGrob(l1, blankPlot, l2, ncol=3),
               arrangeGrob(t1, ncol = 1),
               arrangeGrob(t2, ncol=1),
               arrangeGrob(p7, ncol=1),
               arrangeGrob(soil.tbl, ncol = 1),
               arrangeGrob(textGrob('Pagina 4', gp = gpar(fontsize = 10), just = 'center'), ncol = 1),
               nrow = 6, heights = c(1/10, 0.5/10, 0.5/10, 5.5/10, 2/10, 0.5/10))
  dev.off()
}

# Impute missing values in data.frame
df_impute <- function(dt.frm, n.neig = 2) {
  if (!inherits(dt.frm, "data.frame")) {
    stop("input isn't a data.frame")
  }
  if (!any(is.na(dt.frm))) return(dt.frm)
  # Save default value for nearest neighbors
  def.neigh <- n.neig
  # Get columns with NAs
  na.cols <- colnames(dt.frm)[unlist(lapply(dt.frm, function(x) any(is.na(x))))]
  # For every column with NA do...
  for (cl in na.cols) {
    # Get vector of values
    cl.df <- dt.frm[, cl]
    # Get indices of NAs
    na.lns <- which(is.na(cl.df))
    # For each index do...
    for (ln in na.lns) {
      # Get index of neighbors
      mn.lns <- ((ln - n.neig):(ln + n.neig))[-(n.neig + 1)]
      mn.lns <- mn.lns[mn.lns >= 1 & mn.lns <= length(cl.df)]
      # If NAs in neighbors increase number of neighbors
      if (any(is.na(cl.df[mn.lns]))) {
        mn.lns <- ((ln - (n.neig + 1)):(ln + (n.neig + 1)))[-(n.neig + 1 + 1)]
        mn.lns <- mn.lns[mn.lns >= 1 & mn.lns <= length(cl.df)]
      }
      # Compute mean
      if (is.character(cl.df)) {
        imp.mn <- Mode(cl.df[mn.lns], na.rm = T)
      } else {
        imp.mn <- mean(cl.df[mn.lns], na.rm = T)
        if (is.integer(cl.df)) imp.mn <- round(imp.mn)
      }
      # Replace value
      dt.frm[ln, cl] <- imp.mn
    }
  }
  return(dt.frm)
}

# Raster resampling using GdalUtils
r_rsmp <- function(r.layer, fact = 3) {
  if (!inherits(r.layer, "Raster")) {
    stop("sp.layer isn't a Raster* object")
  }
  require(gdalUtils)
  # Store band names for future use
  bnd.nms <- names(r.layer)
  # Generate temp file names
  tmp1 <- tempfile(fileext = ".tif")
  tmp2 <- tempfile(fileext = ".tif")
  # Write temporary raster
  writeRaster(r.layer, filename = tmp1)
  # Define raster size
  cl.sz <- res(r.layer) / fact
  # Project raster with cubic convolution resampling
  rsmp.rstr <- gdalwarp(srcfile = tmp1, dstfile = tmp2, tr = cl.sz,
                        r = "cubic", output_Raster = T)
  # Assign back band names
  names(rsmp.rstr) <- bnd.nms
  return(rsmp.rstr)
}

Mode <- function(x, na.rm = T) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  ux <- unique(x)
  x.mode <- ux[which.max(tabulate(match(x, ux)))]
  return(x.mode)
}

utm_zone <- function(sp.layer) {
  if (!inherits(sp.layer, "Spatial")){
    stop("sp.layer isn't a Spatial* object")
  }
  sp.cent <- geo_centroid(sp.layer)
  long <- sp.cent[2]
  names(long) <- NULL
  utm.zn <- (floor((long + 180) / 6) %% 60) + 1
  return(utm.zn)
}

save(lndst.pol, prj_str, geo.str, scn_pr, mk_vi_stk, rstr_rcls, int_fx, dem_cov, cols,
     elev_cols, ec_cols, om_cols, swi_cols, cec_cols, presc_grid, hyb.param, hyb_pp, grd_m,
     mz_smth, pnt2rstr, geo_centroid, moran_cln, var_fit, kmz_sv, veris_import, elev_import,
     soil_import, var_cal, trat_grd, multi_mz, srtm.pol, srtm_pr, dem_srtm, read_shp, read_kmz, 
     rstr2pol, report_tdec, write_shp, df_impute, r_rsmp, Mode, utm_zone,
     file = "~/SIG/Geo_util/Functions.RData")
