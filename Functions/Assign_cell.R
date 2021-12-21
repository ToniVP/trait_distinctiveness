#Function to assign the stations to cells of different sizes, you get a list composed with the
#your dataframe with cells assigned and then another dataframe with cell resolution


assign_cell = function(data,num){ #The data needs to have the first column as the name of locations to be included to cells
  
  data[,1] <- as.factor(data[,1])
  uniq <- data[!duplicated(data[,1]),] #obtain the unique names of locations
  try <- data.frame(Long = uniq$lon, Lat = uniq$lat)
  dim <- dim(uniq)
  r<-raster(ncol=round(dim[1]*num), nrow=round(dim[1]*num))
  #here can be modified the size of the cells, it is a grid 
  #with the same number of rows and 
  
  # note that you should use (lon, lat), in that order!
  r <- rasterize(try, r, fun="count") #obtain the number of locations inside each cell, cell size...
  xlim <- c(min(uniq$lon), max(uniq$lon))
  ylim <- c(min(uniq$lat), max(uniq$lat))
  #plot(r, xlim = xlim, ylim = ylim)
  
  tab <- data.frame(rasterToPoints(r)) #pass the raster data to a table
  z <- data.frame(cbind(cell=cellFromXY(r, tab[,1:2]), value=tab[,3]), centroid_lon = tab[,1], centroid_lat = tab[,2]) 
  head(z)
  
  z$xMin <- z$centroid_lon-xres(r)/2 #obtain the edges of each defined cell to cluster locations inside
  z$xMax <- z$centroid_lon+xres(r)/2
  z$yMin <- z$centroid_lat-yres(r)/2
  z$yMax <- z$centroid_lat+yres(r)/2
  
  my_out <- data.frame()
  
  #assign a cell name for each unique station and also the centroids for the cell
  
  system.time(for (i in 1:nrow(uniq)){#start loop cell
    trial <- uniq[i,]
    trial[c("cell", "centroid_lon", "centroid_lat")]<- z[which(trial$lon >= z$xMin &
                                                                 trial$lon < z$xMax & 
                                                                 trial$lat >= z$yMin &
                                                                 trial$lat < z$yMax), c("cell", "centroid_lon", "centroid_lat")]
    
    limits <- z[which(z$cell == trial$cell),]; limits <- limits %>% select(xMin, xMax, yMin, yMax)
    
    trial <- cbind(trial,limits)
    
    my_out <- rbind(my_out, trial);
    
    my_out
    
  }) #end loop cell
  
  my_out[,1] <- as.factor(my_out[,1])
  my_out<- my_out %>% select(station, cell, centroid_lon, centroid_lat, xMin, xMax, yMin, yMax) #add the cell names to the full data frame
  data <- merge(data, my_out, by = "station"); data$cell <- as.factor(data$cell); 
  cellsize <- data.frame("resolution_X" = xres(r), "resolution_Y" = yres(r), 
                         "resolution_X_km" = round(xres(r)*111, digits = 3), "resolution_Y_km"=round(yres(r)*111,digits = 3),
                         "cell_km2" = round((xres(r)*111)*(yres(r)*111), digits = 3))
  rm(r,tab,uniq,try)
  return(list(data, cellsize))
} 

