
ah_range <- function(xy, species=NA, alpha=6, isol.points.buff=TRUE, buff=10, proj="+proj=longlat +datum=WGS84", 
  exclude_map = NULL, is.write=FALSE, path=NA){
  # the first and second column of xy is x and y, and the third is name of species (optional) 
  
  range <- NULL
  range.type <- data.frame("type"=NA,"nocc"=NA,"ahull.error"=FALSE)
  nocc <- nrow(xy)
  if(ncol(xy)>=3 & is.na(species)) species = as.character(xy[1,3])
  
  if(nocc >= 3){
      hull <- try(ahull(xy[,1:2], alpha = alpha), silent = TRUE)
      range <- try(ah2sp(hull, isol.points.buff=isol.points.buff, buff=buff, proj4string=CRS(proj)), silent = TRUE)
      range.type[1,1:2] <- c("alpha_hull", nocc)
      if(inherits(range, "try-error") | is.null(range)) range.type[1,c(1,3)] <- c(NA, TRUE)
  }
  
  if(nocc<3 | inherits(range, "try-error") |is.null(range)){
    spps <- SpatialPoints(xy[,1:2],proj4string=CRS(proj))
    spps <- sp::spTransform(spps, CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs"))
    range <- rgeos::gBuffer(spps,width=buff)
    # to assure the buffered areas are within the geographical coordinate range
    bound = SpatialPolygons(list(Polygons(list(Polygon(cbind(c(-180,180,180,-180),c(89.99,89.99,-89.99,-89.99)))),"p1")), 
                            proj4string=CRS("+proj=longlat +datum=WGS84"))
    bound <- sp::spTransform(bound,  CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs"))
    range <- gIntersection(range, bound, byid=FALSE)
    # transform to the raw projection
    range <- sp::spTransform(range, CRS(proj))
    range.type[1,1:2] <- c("buffer", nocc)
    }
  
  #union incorrect defined holes
  if(!cleangeo::clgeo_IsValid(range)){
    range <- suppressWarnings(gBuffer(range, byid=FALSE, width=0))
    slot(range, "polygons") <- lapply(slot(range, "polygons"), checkPolygonsHoles) 
    range <- unionSpatialPolygons(range,ID=1)    
  }

  if(!is.null(range)){
    range <- SpatialPolygonsDataFrame(range,data=data.frame("species"=species),match.ID=FALSE)
	if(!is.null(exclude_map)){
	if(raster::projection(exclude_map) != raster::projection(range))  exclude_map <- sp::spTransform(exclude_map,crs(range))
	
	range2 <- suppressWarnings(rgeos::gIntersection(range, exclude_map, byid = FALSE))
	if(inherits(range2, "SpatialCollections")) range2 <- aggregate(range2@polyobj)
	
	range <- SpatialPolygonsDataFrame(range2, data = range@data, match.ID=FALSE)
	}

    if(is.write){
      if(is.na(path)) path = paste("Ranges",paste0("alpha",alpha),species,sep="/")
      dir.create(path, recursive=TRUE, showWarnings=FALSE)
      writeOGR(range, dsn=path, layer=species,driver="ESRI Shapefile")
    }
  }
    return(list(range=range,range.type=range.type))
}


