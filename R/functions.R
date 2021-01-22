

#' Load an annotated SVG file
#'
#' This function loads an annotated SVG file, output from QuPath.
#' It can take a minute or so, depending on the file size.
#'
#' @param svgfile Character string name of the input file (ending in .svg)
#' @return An object of class `Picture` from the `grImport2` package.
#' @export
read_svg <- function(svgfile){
  rsvg::rsvg_svg(svgfile, "tmpout.svg")
  p <- grImport2::readPicture("tmpout.svg")
  file.remove("tmpout.svg")
  return(p)
}


#' Check color codes on annotated SVG file
#'
#' This function gives the colors of various annotated objects in an SVG file
#' as well as the count of objects in each color.
#'
#' @param p Picture file output from `read_svg` function
#' @return An object of class `Picture` from the `grImport2` package.
#' @export
colorID <- function(p){
  colvec <- vector()
  for(i in 1:length(p@content[[1]]@content)){
    colvec[i] <- p@content[[1]]@content[[i]]@gp$col
  }
  return(table(colvec))
}


#' Count feature types
#'
#' Gives the number of different feature types in the annotated SVG file based on color codes.
#'
#' @param p Picture file output from `read_svg` function
#' @param box_col hexcode color used for the bounding box
#' @param cam_col hexcode color used for the cambium line
#' @param ray_col hexcode color used for the ray lines
#' @param ves_col hexcode color used for the vessels
#' @return A matrix of feature counts.
#' @export
feature_numbers <- function(p,
                            box_col=c("#44330BFF","#00FF00FF"),
                            cam_col="#00FFFFFF",
                            ray_col="#FF00FFFF",
                            ves_col="#FFFF00FF"){
  ftab <- data.frame(table(unlist(lapply(p@content[[1]]@content, function(x) x@gp$col))))
  fdf <- rbind(ftab[which(ftab$Var1 %in% box_col),],
               ftab[which(ftab$Var1 %in% cam_col),],
               ftab[which(ftab$Var1 %in% ray_col),],
               ftab[which(ftab$Var1 %in% ves_col),])
  colnames(fdf) <- c("Color","Number")
  rownames(fdf) <- c("Bounding box","Cambium line(s)","Ray line(s)","Vessels")
  print(fdf)
}




#' Convert features to sf objects
#'
#' Convert annotated features from the `Picture` object into `sf` objects
#' for spatial analyses.
#'
#' @param p Picture file output from `read_svg` function
#' @param feature Type of feature to be converted.
#' @param col Hexcode color used for the feature. Typically for
#' box=c("#44330BFF","#00FF00FF"),
#' cam="#00FFFFFF",
#' ray="#FF00FFFF",
#' and ves="#FFFF00FF".
#' @param r Size (in micrometers) of each pixel.
#' @return A `sf` object with the annotated features.
#' @export
annotations_to_sf <- function(p,
                              feature=c("box","cam","ray", "ves"),
                              col=NULL,
                              r=0.4424){
  if(feature=="box" & is.null(col)) col <- c("#44330BFF","#00FF00FF")
  if(feature=="cam" & is.null(col)) col <- "#00FFFFFF"
  if(feature=="ray" & is.null(col)) col <- "#FF00FFFF"
  if(feature=="ves" & is.null(col)) col <- "#FFFF00FF"
  if(is.null(col)) {warning('color not specified')}
  selector <- which(unlist(lapply(p@content[[1]]@content, function(x) x@gp$col)) %in% col)
  out_list <- list()
  # selector <- get(paste0(feature,"_sel"))
  if(feature %in% c("box","ves")){
    for(i in 1:length(selector)){
      poly <- do.call(rbind,
                      lapply(p@content[[1]]@content[[selector[i]]]@d@segments,
                             function(x) cbind(x@x, x@y)))
      P1 <- sp::Polygon(poly)
      out_list[[i]] <- sp::Polygons(list(P1), ID=paste(feature, i))
    }
    out_sf <- sf::st_as_sf(sp::SpatialPolygons(out_list))
    out_sf <- sf::st_buffer(out_sf, 0)
  } else {
    for(i in 1:length(selector)){
      line <- do.call(rbind,
                      lapply(p@content[[1]]@content[[selector[i]]]@d@segments,
                             function(x) cbind(x@x, x@y)))
      R1 <- sp::Line(line)
      out_list[[i]] <- sp::Lines(list(R1), ID=paste(feature, i))
    }
    out_sf <- sf::st_as_sf(sp::SpatialLines(out_list))
  }
  out_sf$geometry <- out_sf$geometry * r
  if(feature=="ves"){
    out_sf$area <- sf::st_area(out_sf$geometry)
  }
  return(out_sf)
}




#' Measure distance along rays to cambium
#'
#' Measure distance along rays to cambium line of features.
#'
#' @param csf `sf` object of cambium line
#' @param rsf `sf` object of ray line(s)
#' @param vsf `sf` object of vessel polygons
#' @param parallel Logical if processing should be done in parallel (not yet implemented).
#' @param cores Number of cores to use for parallel processing.
#' @return A vector of distances (in microns) from each feature to the
#' cambium line, along the closest ray.
#' @export
cam_dist <- function(csf,
                     rsf,
                     vsf,
                     parallel=F,
                     cores=detectCores(logical=FALSE)){
  pts <- st_centroid(vsf)
  nr <- length(rsf$geometry)
  nv <- length(pts[[1]])
  ray_nearest <- st_nearest_points(pts, rsf)
  ### Select which ray is closest for each vessel
  ray_selector <- tapply(sf::st_length(ray_nearest),
                         rep(1:nv, each=nr),
                         function(x) which(x==min(x)))
  ### Select the shortest line to nearest ray
  ray_nearest2 <- ray_nearest[seq(0, nr*nv, nr)[-(nv+1)]+ray_selector,]
  # Get the point along the ray line
  ray_nearest_pts <- st_cast(ray_nearest2, "POINT")[seq(2, (2*nv), 2)]
  # Make into polygon to slice rayline
  ray_nearest_polys <- st_buffer(ray_nearest_pts, 0.25)
  # Find nearest point on cambium to vessel-ray intersections
  cam_nearest <- st_nearest_points(csf, rsf)
  cam_nearest_pts <- st_cast(cam_nearest, "POINT")[seq(1, (2 * nr), 2)]
  cam_nearest_poly <- st_buffer(cam_nearest_pts, 0.25)
  cam_nearest_polys <- cam_nearest_poly[ray_selector,]
  ### Get distances to cambium

  if(parallel==F){
    cdist <- vector()
    message(paste("Computing nearest distance to cambium along raylines..."))
    pb <- progress_bar$new(total = nv)
    for(i in 1:nv){
      pb$tick()
      focray <- rsf[ray_selector[i],]
      raypoly <- ray_nearest_polys[i,]
      campoly <- cam_nearest_polys[i,]
      rs_tmp <- st_collection_extract(lwgeom::st_split(focray, raypoly), "LINESTRING")

      ## New way to do this to avoid problems with little loops in the ray lines
      rs <- st_collection_extract(lwgeom::st_split(rs_tmp, campoly), "LINESTRING")
      sint_r <- lengths(st_intersects(st_buffer(rs,2), raypoly))
      sint_c <- lengths(st_intersects(st_buffer(rs,2),campoly))
      cdist[i] <- sum(st_length(rs[which((sint_r + sint_c)==2),]))
    }
  }

  if(parallel==T){
    message(paste("Computing nearest distance to cambium along raylines...(in parallel!)"))
    warning(paste("Error: Parallel not yet implemented!!!"))
    # require(foreach)
    # require(doFuture)
    # require(progressr)
    # registerDoFuture()
    # cl <- parallel::makeCluster(cores, setup_strategy = "sequential")
    # plan(cluster, workers=cl)
    # with_progress({
    #   p <- progressor(along=1:nv)
    #   cdist <- foreach(i=1:nv,
    #                    .packages=c("sf","lwgeom"),
    #                    .combine='c') %dopar% {
    #     p(i)
    #     focray <- ray_sf[ray_selector[i],]
    #     raypoly <- ray_nearest_polys[i,]
    #     campoly <- cam_nearest_polys[i,]
    #     rs_tmp <- st_collection_extract(lwgeom::st_split(focray, raypoly), "LINESTRING")
    #     rs <- st_collection_extract(lwgeom::st_split(rs_tmp, campoly), "LINESTRING")
    #     sint_r <- st_intersects(rs, raypoly)
    #     sint_c <- st_intersects(rs,campoly)
    #     focrs <- rs[min(which(lengths(sint_r)>0)):(max(which(lengths(sint_c)>0))),]
    #     sum(st_length(focrs))
    #     }
    # })
    # parallel::stopCluster(cl)
  }
  return(cdist)
}



#' Create a thumbnail image if annotations
#'
#' Quickly create a thumbnail image to confirm annotations are being
#' properly processed by `qwar` functions
#'
#' @param outfile Character name of output file (end in .pdf), optionally including the path.
#' @param bsf `sf` object of bounding box
#' @param csf `sf` object of cambium line
#' @param rsf `sf` object of ray line(s)
#' @param vsf `sf` object of vessel polygons
#' @param vstat Attribute from vessel polygons to use in plotting.
#' @return Saves a PDF file showing the output annotations.
#' @export
thumbnail_check <- function(outfile=NULL,
                            bsf,
                            csf,
                            rsf,
                            vsf,
                            vstat='area'){
  pdf(outfile)
  plot(bsf$geometry, axes=T, lwd=2,
       main=outfile)
  plot(csf$geometry, col=3, lwd=3, add=T)
  plot(rsf$geometry, col=2, lwd=3, add=T)
  if(is.null(vstat)){
    plot(vsf$geometry, add=T, lwd=0.25, col='yellow')
  } else {
    plot(vsf[vstat], add=T, lwd=0.25)
  }
  dev.off()
}




#' Maximum vessel diameter
#'
#' Compute the maximum diameter across a vessel based on
#' 2x the longest radius from edge of polygon to centroid.
#'
#' @param vsf `sf` object of vessel polygons
#' @param parallel Logical if processing should be done in parallel
#' (implemented but of questionably utility).
#' @param cores Number of cores to use for parallel processing.
#' @return A vector of maximum diameter (in microns) for each vessel.
#' @export
ves_max_diam <- function(vsf,
                         parallel=TRUE,
                         cores=detectCores(logical=FALSE)){
  vc <- st_centroid(vsf$geometry)
  if(parallel==F){
    message(paste("Computing max vessel diameter..."))
    ves_diam <- vector()
    pb <- progress_bar$new(total = length(vsf$geometry))
    for(i in 1:length(vsf$geometry)){
      pb$tick()
      pts <- st_cast(vsf[i,]$geometry, 'POINT')
      ves_diam[i] <- max(st_distance(pts))
    }
  }
  if(parallel==T){
    message(paste("Computing max vessel diameter... (in parallel!)"))
    vc <- st_centroid(vsf$geometry)
    nv <- length(vsf$geometry)
    require(foreach)
    require(doFuture)
    require(progressr)
    registerDoFuture()
    cl <- parallel::makeCluster(cores, setup_strategy = "sequential")
    plan(cluster, workers=cl)
    with_progress({
      p <- progressor(along=1:nv)
      ves_diam <- foreach(i=1:nv, .packages=c("sf","lwgeom"),
                          .combine='c') %dopar% {
                            p(i)
                            pts <- st_cast(vsf[i,]$geometry, 'POINT')
                            max(st_distance(pts))
                          }
    })
    parallel::stopCluster(cl)
    Sys.sleep(4)
  }
  return(ves_diam)
}



#' Circle vessel diameter
#'
#' Compute the diameter of vessels based on a circle with equal area.
#'
#' @param vsf `sf` object of vessel polygons
#' @return A vector of cirular diameter (in microns) for each vessel.
#' @export
ves_diam_circle <- function(vsf){
  message(paste("Computing idealized vessel diameter as: 2 * sqrt(area)"))
  D <- 2 * (sqrt(vsf$area/pi))
  return(D)
}




#' Vessel diameter based on median radius
#'
#' Compute the diameter of vessels based on the median radius from
#' the centroid to vessel polygon edge.
#'
#' @param vsf `sf` object of vessel polygons
#' @param parallel Logical if processing should be done in parallel
#' (implemented but of questionably utility).
#' @param cores Number of cores to use for parallel processing.
#' @return A vector of diameter (in microns) for each vessel.
#' @export
ves_diam_median_radius <- function(vsf,
                                   parallel=TRUE,
                                   cores=detectCores(logical=FALSE)){
  vc <- st_centroid(vsf$geometry)
  if(parallel==F){
    ves_rad <- vector()
    message(paste("Computing median diameter based on all vessel radii..."))
    pb <- progress_bar$new(total = length(vsf$geometry))
    for(i in 1:length(vsf$geometry)){
      pb$tick()
      pts <- st_cast(vsf[i,]$geometry, 'POINT')
      rad <- st_nearest_points(vc[i,], pts)
      ves_rad[i] <- median(st_length(rad))
    }
  }
  if(parallel==T){
    message(paste("Computing median diameter based on all vessel radii...(in parallel!)"))
    nv <- length(vsf$geometry)
    require(foreach)
    require(doFuture)
    require(progressr)
    registerDoFuture()
    cl <- parallel::makeCluster(cores, setup_strategy = "sequential")
    plan(cluster, workers=cl)
    with_progress({
      p <- progressor(along=1:nv)
      ves_rad <- foreach(i=1:nv, .packages=c("sf","lwgeom"),
                         .combine='c') %dopar% {
                           p(i)
                           pts <- st_cast(vsf[i,]$geometry, 'POINT')
                           rad <- st_nearest_points(vc[i,], pts)
                           median(st_length(rad))
                         }
    })
    parallel::stopCluster(cl)
    Sys.sleep(4)
  }
  return(2 * ves_rad)
}




#' Group size
#'
#' Compute the size of the group that each vessels belongs to.
#'
#' @param vsf `sf` object of vessel polygons
#' @param thresh Distance threshold (in microns) to idenfity grouped vessels.
#' @return A vector of group size for each vessel.
#' @export
group_size <- function(vsf,
                       thresh=33.26){
  ves.nb <- spdep::poly2nb(vsf, snap=thresh)
  gs <- card(spdep::nblag_cumul(spdep::nblag(ves.nb, 30)))+1
  return(gs)
}





#' Vessel Grouping Index
#'
#' Function for getting vessel grouping index (average number of vessels per group)
#' and proportion of vessels in groups given a `vsf` object. Can be applied to
#' subsets of the `vsf` to, for example, examine changes in grouping index in
#' different parts of the sample.
#' @param vsf `sf` object of vessel polygons
#' @param thresh Distance threshold (in microns) to identify grouped vessels.
#' @return A data.frame with grouping indices for vessels included in the sample.
#' @export
ves_group_indices <- function(vsf,
                              thresh=33.26){
  if(length(vsf$geometry) > 1){
    ves.nb <- spdep::poly2nb(vsf$geometry, snap=(thresh))
    if(class(vsf)[1]=="sfc_GEOMETRY"){
      nv <- length(vsf)
    } else {
      nv <- length(vsf$geometry)
    }
    # Vessel multiple fraction
    vmf <- (length(unique(unlist(ves.nb)))-1) / nv
    # Vessel grouping index
    tmpvgi <- vector()
    for(v in 1:max(vsf$group_size)){
      tmpvgi[v] <- sum(vsf$group_size==v)/v
    }
    vgi <- nv/sum(tmpvgi)
    # Solitary vessel index
    svi <- sum(vsf$group_size==1)/sum(tmpvgi)
    # Mean group size
    mgs <- mean(card(nblag_cumul(nblag(ves.nb, 30)))+1)
    # Compile grouping parameters for output
    gi_df <- data.frame(N_ves=nv,
                        ves_mult_frac=vmf,
                        ves_grp_index=vgi,
                        sol_ves_index=svi,
                        mean_group_size=mgs)
    return(gi_df)
  } else {
    gi_df <- data.frame(N_ves=length(vsf$geometry),
                        ves_mult_frac=NA,
                        ves_grp_index=NA,
                        sol_ves_index=ifelse(length(vsf$geometry)==1,1,NA),
                        mean_group_size=NA)
  }
}


#' Distance bins of vessel characteristics
#'
#' Get characteristics of vessels by binned distance from cambium.
#' @param vsf `sf` object of vessel polygons
#' @param binsize Size of distance bins to away from cambium.
#' @return A data.frame with grouping indices for vessels included in each distance bin.
#' Includes vessel grouping indices, mean and median vessel area, kurtosis, skewness,
#' hydraulically-weighted diameter (Dh), total theoretical conductivity (Kh_total),
#' and mean theoretical conductivity (Kh_mean).
#' @export
ves_characteristics_bin <- function(vsf,
                                    binsize=200){
  bins <- seq(0, max(vsf$dist)+binsize, binsize)
  lower <- bins[-length(bins)]
  upper <- bins[-1]
  outmat <- matrix(nrow=length(bins)-1, ncol=13)
  pb <- progress_bar$new(total = length(lower))
  for(i in 1:length(lower)){
    pb$tick()
    focves <- vsf[vsf$dist > lower[i] & vsf$dist < upper[i],]
    outmat[i,1] <- mean(c(lower[i], upper[i]))
    vg <- ves_group_indices(focves)
    outmat[i,2:6] <- unlist(vg)
    outmat[i,7] <- mean(focves$area)
    outmat[i,8] <- median(focves$area)
    outmat[i,9] <- moments::kurtosis(focves$area)
    outmat[i,10] <- moments::skewness(focves$area)
    outmat[i,11] <- (sum(focves$diam^4)/length(focves))^(0.25)
    outmat[i,12] <- sum((pi * ((focves$Dmedrad/1000)^4))/(128 * (1.002e-9)))
    outmat[i,13] <- mean((pi * ((focves$Dmedrad/1000)^4))/(128 * (1.002e-9)))
  }
  colnames(outmat) <- c("dist", names(vg), "mean_area", "median_area",
                        "kurtosis", "skewness", "Dh", "Kh_total", "Kh_mean")
  return(as.data.frame(outmat))
}




#' Fractional area of vessels
#'
#' Get fractional area of vessels in distance bins away from cambium based on point sampling.
#' @param bsf `sf` object of bounding box
#' @param csf `sf` object of cambium line
#' @param rsf `sf` object of ray line(s)
#' @param vsf `sf` object of vessel polygons
#' @param npts Number of random points to use for sampling
#' @param binsize Size of distance bins to away from cambium
#' @return A data.frame with the number of random sample points included in each bin, the number of vessels in each bin, the vessel fraction in each bin, the mean vessel area in each bin, and the vessel density in each bin.
#' @export
ves_fraction_pt_sample <- function(bsf,
                                   csf,
                                   rsf,
                                   vsf,
                                   npts=100000,
                                   binsize=200){
  message(paste("Generating", npts, "random points..."))
  x <- st_as_sf(st_sample(bsf, npts))
  nr <- length(rsf$geometry)
  x$dist <- cam_dist(csf, rsf, vsf)
  bins <- seq(0, max(x$dist)+1, binsize)
  lower <- bins[-length(bins)]
  upper <- bins[-1]
  outmat <- matrix(nrow=length(bins)-1, ncol=6)
  message('Computing vessel fraction in bins...')
  pb <- progress_bar$new(total = length(lower))
  for(b in 1:length(lower)){
    pb$tick()
    focpts <- x[x$dist > lower[b] & x$dist < upper[b],]
    focves <- vsf[vsf$dist > lower[b] & vsf$dist < upper[b],]

    if(length(focpts$x)>0){
      focpts$ves <- rowSums(st_intersects(focpts, focves, sparse=F))
      outmat[b,1] <- mean(c(lower[b], upper[b]))
      outmat[b,2] <- length(focpts$dist)
      outmat[b,3] <- length(focves$dist)
      outmat[b,4] <- mean(focpts$ves)
      outmat[b,5] <- mean(focves$area)
      outmat[b,6] <- outmat[b,4] / (mean(focves$area)/1e6)
    } else {
      outmat[b,1] <- mean(c(lower[b], upper[b]))
      outmat[b,2] <- length(focpts$dist)
      outmat[b,3] <- length(focves$dist)
      outmat[b,4] <- NA
      outmat[b,5] <- mean(focves$area, na.rm=T)
      outmat[b,6] <- NA
    }
  }
  colnames(outmat) <- c("dist", "N_pts", "N_ves","ves_fraction","mean_ves_area","ves_density")
  return(as.data.frame(outmat))
}
