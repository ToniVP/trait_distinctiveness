#################################################
## Project 1: case of Marenzelleria #############
#################################################
# Last update: November 3, 2021
# Loading a raw database, adding trait information for several species, extracting a matrix of abundances
# Computing a distance matrix based on trait information and projecting species in the trait space of the community (PCoA)

remove(list=ls())

#1: Libraries and functions needed #########
library(multcomp); library(hrbrthemes); library(viridis); library(ggplot2); library(forcats); library(maptools); library(gridExtra)
library(car); library(grid); library(RColorBrewer); library(colorRamps); library(lme4); library(AER); library(ape); library(automap)
library(raster); library(pscl); library(lmtest); library(stats); library(geoR); library(tidyverse); library(sp); library(gstat)
library(ggmap); library(rgdal); library(rdist); library(MASS); library(effects); library(plyr); library(vegan); library(devtools)
library(ggbiplot); library(readxl); library(stringr); library(Rcpp); library(gawdis); library(FD); library(maps); library(funrar)
library(ggrepel); library(fastDummies); library(Rmisc); library(cluster); library(data.table); library(prob); library(splitstackshape); 
library(stringi); library(DT); library(openxlsx); library(MASS); library(parallel);library(foreach);library(doParallel);library(doSNOW)
library(ggthemes); library(colorspace); library(ggsci); library(rangeBuilder);library(ncdf4); library(RANN); library(exactextractr)
#install_github("vqv/ggbiplot")
setwd('C:/Users/avipo/Desktop/PhD/Project 1/MarenzelleriaSweden')

#Load all the functions
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
} #Function to assign the stations to cells of different sizes, you get a list composed with the
                                      #your dataframe with cells assigned and then another dataframe with cell resolution

get_common_species = function(x, y, num) { # x as the abundance matrix fo the year/station that we want, as there are the selected species
  # y as the raw trait matrix, to select species present in the abundance matrix
  # num as the column number where the first species is
  
  if(class(x) == "data.frame"){colnames(x)<- gsub("\\.", " ", colnames(x)); mat_columns = colnames(x)[num:ncol(x)]}
  else {mat_columns = x} #the abundance matrix or a vector containing names
  
  trait_species = rownames(y) #the trait dataset

  # Extract matching species
  common_species = match(trait_species, mat_columns)
  rows = which(!is.na(common_species))
  
  return(y[rows,]) 
  
} #Function to select species of the trait database

standr = function(x){(x-min(x))/(max(x)-min(x))} #Function to standardize

int_distinct = function(trait){ #Here trait is a species x trait information matrix, previously cleaned 
  #and with the traits of interest selected
  
  sensitivity.gow<-trait; sensitivity.gow[1:ncol(sensitivity.gow)] <- NULL
  sum.standard.gow<-matrix(0,nrow(as.data.frame(trait)),nrow(as.data.frame(trait)))
  rownames(sum.standard.gow)<-rownames(trait)
  colnames(sum.standard.gow)<-rownames(trait)
  sum.standard.gow<-as.data.frame(sum.standard.gow)
  
  # Defining the total number of column : N; it MUST be 
  N<-ncol(trait)
  listcol<-seq(1,N,1)
  # Number of columns corresponding to traits
  listk<-seq(1,N,1)
  
  if (ncol(trait) >= 4){#start condition for more than 4 trait categories
    k=4
    for (k in 4:N){ #start loop combination of traits
      
      #Here are determined all possible 4 trait combinations for obtaining the average distance
      #between species, here are expressed only the position that are occupying the traits in the data frame
      
      combination<-data.frame(t(combn(listk,k))) 
      
      if (is.null(combination$X5==T)){
        combination$X5<-NA
      }
      if (is.null(combination$X6==T)){
        combination$X6<-NA
      }
      if (is.null(combination$X7==T)){
        combination$X7<-NA
      }
      if (is.null(combination$X8==T)){
        combination$X8<-NA
      }
      if (is.null(combination$X9==T)){
        combination$X9<-NA
      }
      if (is.null(combination$X10==T)){
        combination$X10<-NA
      }
      if (is.null(combination$X11==T)){
        combination$X11<-NA
      }
      
      #This list of conditionals should be equal or have more elements than your trait number, NEVER less
      
      for (i in 1:nrow(combination)){ 
        data <- trait
        codex<-as.vector(t(combination[i,]))
        
        #Here, each row (in other words, each combination of 4 traits) is selected for a 
        #total calculation of distances based on all combinations possibles of the traits
        data<-dplyr::select(data,all_of(codex))
        
        #Calculation of distinctiveness using Gower's distance
        gow<-compute_dist_matrix(data, metric = "gower") # a distance matrix is calculated for each combination of traits
        m.gow<-as.matrix(gow)
        d.gow<-as.data.frame(m.gow)
        standard.gow<-standr(d.gow) # distances are standardized to fit between 1 and 0
        standard.Di.gow<-colSums(standard.gow)/(nrow(standard.gow)-1) # Mean of the computed distances for each species, to obtain the
        # average distance to the other species based on all trait combinations
        standard.Di.gow<-as.data.frame(standard.Di.gow)
        sensitivity.gow<-cbind(sensitivity.gow,standard.Di.gow)
      }
      
      print (k)
    }
    
  }#end condition for more than 4 trait categories
  
  else{#start condition for less than 3 categories
    
    #Calculation of distinctiveness using Gower's distance
    gow<-compute_dist_matrix(trait, metric = "gower") # a distance matrix is calculated for each combination of traits
    m.gow<-as.matrix(gow)
    d.gow<-as.data.frame(m.gow)
    standard.gow<-standr(d.gow) # distances are standardized to fit between 1 and 0
    standard.Di.gow<-colSums(standard.gow)/(nrow(standard.gow)-1) # Mean of the computed distances for each species, to obtain the
    # average distance to the other species based on all trait combinations
    standard.Di.gow<-as.data.frame(standard.Di.gow)
    sensitivity.gow<-cbind(sensitivity.gow,standard.Di.gow)} #end condition for less than 3 categories
  
  forIntDi<-t(sensitivity.gow); forIntDi<-na.omit(forIntDi); forIntDi<-as.matrix(t(forIntDi))
  outputmatrix<-matrix(0,nrow(forIntDi),2)
  
  for (j in 1:nrow(forIntDi)){# start loop output matrix
    
    subset_sp<-forIntDi[j,1:ncol(forIntDi)]
    subset_sp<-as.numeric(as.character(subset_sp))
    outputmatrix[j,2]<-mean(subset_sp) #compute the mean for all the distances obtained via the combinations
    
  }#end loop output matrix
  
  Int_Di<-data.frame(outputmatrix); Int_Di[,1]<-rownames(forIntDi); colnames(Int_Di)<-c("taxon","int_Di"); 
  #Int_Di<-Int_Di[order(Int_Di$int_Di,decreasing=T),]
  
  return(list(Int_Di, standard.gow))
} #Function to calculate integrated distinctiveness

convexe_hull<-function(poly){
  s <- PCOA_IDi %>%
    split(poly)  # Tranform data.frame in lists of data.frames depending on column 'z'
  ch = s %>%
    # Compute which points are on the convex hull of each data.frame 
    lapply(., function(el) chull(el$PCOA1, el$PCOA2))  # 'ch' now contains the row numbers of points on convex hull per sub-data.frame
  # Get points for each sub-data.frame using names index
  ch = lapply(names(ch), function(PCOA1) s[[PCOA1]][ch[[PCOA1]],]) %>%   
    do.call(rbind, .)  # Join all convex hull points in a single data.frame
  return(ch)
} #Function to make hulls about the similarity of species

format_trait = function(x) { #edit the function based on what should be cleaned from the selected trait dataset
  # to unify the criteria and be able to add different information manually
  
  # Delete non-biological information and the information inside parentheses
  x <- data.frame(lapply(x, function(x) {gsub(" *\\(.*?\\) *|See additional information|Insufficient information|
                                                      Not relevant|Field unresearched|Not researched|See additional Information|
                                                      No information found|None", "", x)}))
  
  x$Size <- stri_replace_all_regex(x$Size, c("Large|Medium-large|Medium","Small-medium","Small", "Very small"),
                                           c(">100","51-100","11-20","0-10"), vectorize=F)
  
  x$LifeSpan <- stri_replace_all_regex(x$LifeSpan, c("11-20 years","3-5 years","6-10 years", "<1 year", "1 year|1-2 years"),
                                                           c(">10yrs","3-6yrs","6-10yrs","<1yr", "1-3yrs"), vectorize=F)
  
  x$Habit <- stri_replace_all_regex(x$Habit, c("Free living","Burrow dwelling","Tubiculous"),
                                                     c("Free","Burrow dweller","Tube dweller"), vectorize=F)
  
  x$envpos <- stri_replace_all_regex(x$envpos, c("Infaunal","Epibenthic;Pelagic"),
                                                       c("Middle;Top","Bentho-pelagic"), vectorize=F)
  
  x$feedingmethod <- stri_replace_all_regex(x$feedingmethod, 
                                                c("Passive suspension feeder|Active suspension feeder","Surface deposit feeder|Sub-surface deposit feeder"),
                                                c("Suspension/filter feeder","Deposit feeder (inkl. Both)"), vectorize=F)
  
  x$ReprodFreq <- stri_replace_all_regex(x$ReprodFreq, c("Annual episodic|Annual protracted|Biannual protracted"), c("Iteroparous-Polytelic"), vectorize=F)
  
  
  x$devmech <- stri_replace_all_regex(x$devmech, c("Direct Development"), c("Direct"), vectorize=F)
  
  x$Bioturbator <- stri_replace_all_regex(x$Bioturbator, c("Not relevant","Conveyor belt transport/Reverse conveyor belt transport"), c("No transport", "Conveyor belt transport;Reverse conveyor belt transport"), vectorize=F)
  
  return(x)
  } #Function to unify categories of the traits

mean_matrix<- function(x){
  y <- array(unlist(x) , c(dim(x[[1]]),length(x)))
  y <- apply( arr , 1:2 , mean )
  colnames(y)<-colnames(x[[1]])
  rownames(y)<-rownames(x[[1]])
  return(y)
} #Obtain the mean of each element from a list of matrices

mean_diss = function(dist_matrix, data, status, Int_Di, num){
  
  if (all(is.na(match(colnames(data), rownames(status[which(status$status == "Non-indigenous"),]))))){NIS <- c()}
  else {NIS <- c(colnames(data)[which(!is.na(match(colnames(data), rownames(status[which(status$status == "Non-indigenous"),]))))])}
  
  Allsp <- c(colnames(data)[num:ncol(data)]); Allsp <- gsub("\\.", " ", Allsp)
  NAT <- setdiff(Allsp,NIS)
  
  #mean distance of NIS to other native species per each location/cell
  NIS_Di <- foreach(i = 1:length(NIS), .combine = cbind, .errorhandling = "pass") %dopar% {
    dis <- c()
    for(j in 1:length(NAT)) {val <- dist_matrix[NIS[i],NAT[j]]; dis <- c(dis,val)}
    dis <- mean(dis); dist <- as.data.frame(dis); colnames(dist) <- NIS[i]; dist}
  
  val <- Int_Di[which(!is.na(match(Int_Di$taxon,NAT))),]
  Di_NAT <- mean(val$int_Di)
  
  val <- Int_Di[which(!is.na(match(Int_Di$taxon,Allsp))),]
  Di_ovrll <- mean(val$int_Di)
  
  return(list(NIS_Di, Di_NAT, Di_ovrll))
} #Function to calculate distinctiveness between different species groups
#dist_matrix accounts for a pairwise matrix of integrated distinctiveness (functional distances)
#data is the species occurrences/abundance data for a certain location/cell where we have different sampling locations
#status is a list for all the species we have in TOTAL where their status (NIS or native) is defined
#IntDi is the result of the integrated distinctiveness for all the species obtained previously

MM <- function(k, z, x){
  y <- k*x/(z+x)
  return(y)
}

MMprime <- function(k,z,x){
  y <- k*z/((z+x)^2)
  return(y)
}

#useful functions to make cells and assign a name to each generated cell, respectively (useful functions by MARCEL)
createGrid <- function(gSize, xlim, ylim){
  # Create grid
  my_grid <- expand.grid(Lon_centroid = seq(xlim[1], xlim[2], by = gSize),
                         Lat_centroid = seq(ylim[1], ylim[2], by = gSize))
  
  # Create longitude names
  my_x <- data.frame(Lon = unique(my_grid$Lon_centroid))
  my_x$letters <- LETTERS[as.numeric(substr(my_x$Lon, 1,2))-round(min(my_x$Lon-1))]
  
  x_name <- data.frame(decimales = unique(as.numeric(substr(my_x$Lon, 3,5))))
  x_name[is.na(x_name)] <- 0
  
  x_name$numbers <- as.numeric(1:length(unique(x_name$decimales)))
  
  my_x$numbers <- 1
  
  # i <- 4
  for(i in unique(x_name$numbers)){
    my_x$numbers[which(as.numeric(substr(my_x$Lon, 3,4)) == x_name$decimales[which(x_name$numbers == i)] )] <- i
  }
  
  # Create latitude names
  my_y <- data.frame(Lat=unique(my_grid$Lat_centroid))
  my_y$numbers2 <- seq(from = 1, to =length(unique(my_y$Lat)))
  my_y$numbers2 <- seq(from=1, to=length(unique(my_y$Lat)), by=1)
  
  if(max(my_y$numbers2) >=10){
    my_y$numbers2 <- sprintf("%02d", as.numeric(my_y$numbers2))
  }
  if(max(my_y$numbers2) >=100){
    my_y$numbers2 <- sprintf("%03d", as.numeric(my_y$numbers2))
  }
  
  my_names <- merge(my_grid, my_y, by.x="Lat_centroid", by.y = "Lat", all = TRUE)
  my_names <- merge(my_names, my_x, by.x="Lon_centroid", by.y = "Lon", all = TRUE)
  # head(my_names)
  # sum(is.na(my_names))
  dGrid <- data.frame(Lon = my_names$Lon_centroid,
                      Lat = my_names$Lat_centroid,
                      Rectangle = paste0(my_names$numbers2, my_names$letters, my_names$numbers))
  return(dGrid)
}; 
separateToGrid <- function(Lat, Lon, myGrid, gSize){
  
  # extractGrid <- expand.grid(Lon_centroid = seq(xlim[1], xlim[2], by = gSize),
  #                        Lat_centroid = seq(ylim[1], ylim[2], by = gSize))
  myGrid$xMin <- myGrid$Lon-gSize/2
  myGrid$xMax <- myGrid$Lon+gSize/2
  myGrid$yMin <- myGrid$Lat-gSize/2
  myGrid$yMax <- myGrid$Lat+gSize/2
  
  my_data <- data.frame(Lat = Lat, Lon = Lon)
  
  
  my_data[,c("Rectangle", "Centroid_lat", "Centroid_lon")] <- NA
  my_out <- data.frame()
  for(i in 1:nrow(my_data)){
    work <- my_data[i,]
    # work$Rectangle <- myGrid$Rectangle[which(work$Lon >= myGrid$xMin &
    #                                            work$Lon < myGrid$xMax &
    #                                            work$Lat >= myGrid$yMin &
    #                                            work$Lat < myGrid$yMax)]
    
    work[c("Rectangle", "Centroid_lat", "Centroid_lon")] <- myGrid[which(
      work$Lon >= myGrid$xMin &
        work$Lon < myGrid$xMax &
        work$Lat >= myGrid$yMin &
        work$Lat < myGrid$yMax),
      c("Rectangle", "Lat", "Lon")]
    
    my_out <- rbind(my_out, work)
  }
  return(my_out)
} 


#Set the number of cores to do functions in parallel
numcores <- detectCores()
registerDoParallel(numcores)
ctrl <- nls.control(warnOnly = TRUE,minFactor = 1/4096) #control parameters for setting the non-linear MM model

########################################
#2: Loading raw data and explore #######!
########################################
#2.0: Start cleaning the raw occurrences data ##########
rawdata<- read_excel("MarenzelleriaSweden_1990_2020.xlsx")
raw<-rawdata
attach(raw) ## Just for exploration of the variables then use detach(data)
#detach(raw)
head(raw)
raw$year <- substr(raw$date, 1,4); raw$month<-substr(raw$date,6,7)#extract the year and month
newdata <- data.frame(station = as.factor(raw$unique_station_temp_id), lon = longitude,
                    lat = latitude, year = as.factor(raw$year), month = as.factor(raw$month), depth = depth,
                    sample_number = sample_number, taxon = taxon,
                    taxon_id = taxonID_dyntaxa, abundance = abundance, wet_weight = wet_weight) #create a new database with variables of interest
head(newdata)
newdata<-newdata[-which(is.na(newdata$taxon)),] #remove NA species
newdata$station <- as.factor(newdata$station)
levels(newdata$station) <- 1:length(levels(newdata$station))
newdata<-newdata[-which(is.na(newdata$lon)),]#remove NA coordinates
#write.table(newdata, "Marenzelleria_1990_2020_FULL.txt", sep="\t")

# Put the stations in the map
ggplot() +
  geom_point(aes(x=newdata$lon,y=newdata$lat), color = "darkgreen") +
  # scale_fill_gradientn(name = "Clusters", colours=colpal, na.value = 'white') +
  borders(fill="gray44",colour="black") +
  coord_quickmap(xlim=c(10, 30),ylim= c(54,67))+
  labs(x = "Lon",y="Lat")+
  theme(panel.background = element_rect(fill=alpha('light blue', 0.4), colour = 'black'))

#2.0.1: Delete those species that are very rare or have very few occurrences ######
data<-read.csv("Marenzelleria_1990_2020_FULL.txt",header=T,dec=".",sep="\t")

data$abundance[is.na(data$abundance)] <- 0
ocur <- data.frame(count(data, vars = c("taxon", "year")))
abund <- data.frame(abundance = tapply(data$abundance, data$taxon, mean)); #abund$abundance <- log(abund$abundance)
ocur <- data.frame(occurrences_mean = tapply(ocur$freq, ocur$taxon, sum))
#ocur <- data.frame(occurrences = tapply(ocur$freq, ocur$taxon, mean))
oc_ab <- cbind(taxon = rownames(abund), abund, ocur); rownames(oc_ab)<-NULL
oc_ab <- oc_ab[-which(oc_ab$abundance <= 0),]
#par(mfrow = c(1,2));hist(oc_ab$abundance);hist(oc_ab$occurrences)

#Find the abundance level that goes below 5% of all abundances for the FULL data
CIa <- CI(oc_ab$abundance, ci = 0.99) # We use a 99% to see which values fall below the 0.5% rather than the 2.5% that would be the 95 conventional confidence interval
CIb <- CI(oc_ab$occurrences, ci = 0.99)
oc_ab <- oc_ab[which(oc_ab$occurrences >= CIb[3]),]
plot(oc_ab$abundance, oc_ab$occurrences)

#ggplot(data, aes(x = year)) + geom_bar() + facet_wrap(vars(taxon)) #then we visualize

# Filter the data by species that are below 2.5% of occurrences for all years
common_species <-  match(data$taxon, oc_ab$taxon)
rows <- which(!is.na(common_species))
data_filter <- data[rows,]; data_filter$taxon <- as.factor(data_filter$taxon)
data <- data_filter
#write.table(data_filter, "Marenzelleria_1990_2020.txt", sep="\t")

# Measure the effort for each year per station
data<-read.csv("Marenzelleria_1990_2020_FULL.txt",header=T,dec=".",sep="\t")
data$station <- as.factor(data$station)
data$sampleID <- paste (data$station,data$sample_number)
data <- data[!duplicated(data[c("sampleID")]),] 

data <- data.frame(data[which(as.numeric(data$station) >= 1001 & as.numeric(data$station) <= 1670),])

# Count the number of occurrences per station/taxon per year
#?count
eff <- data.frame(count(data, vars = c("taxon", "year")))
hist(eff$freq)

#then we visualize all the occurrences
ggplot(data, aes(x = year)) + geom_bar() + facet_wrap(vars(taxon))

#2.0.2: Determine which are the non-indigenous species on the database, assign which species occur occasionally ######
data<-read.csv("Marenzelleria_1990_2020_FULL.txt",header=T,dec=".",sep="\t") #read both the clean (occurrence filtered) and the full database
dclean<-read.csv("Marenzelleria_1990_2020.txt",header=T,dec=".",sep="\t")
introd<- read_excel("introduced_sp.xlsx")
introd$status <- ifelse(as.numeric(introd$status) == 1, "Non-indigenous", NA); introd$status <- ifelse(is.na(introd$status), "Cryptogenic", "Non-indigenous")
introd <- introd %>% select(taxon,status)

#Non-indigenous species from the whole database
common_species <-  match(data$taxon, introd$taxon)
rows <- which(!is.na(common_species))
data_filter <- data[rows,]; data_filter$taxon <- as.factor(data_filter$taxon)
nonnative<-data.frame(taxon = levels(data_filter$taxon))
nonnative <- merge(nonnative, introd, by = "taxon") #see the species names
#write.table(nonnative, "introducedsp_wholespecies.txt", row.names = FALSE,  sep="\t")
see <- data.frame(count(data_filter, vars = c("taxon", "year")))
see <- data.frame(occurrences_sum = tapply(see$freq, see$taxon, sum))
clus <- pam(see, k = 3, metric = "manhattan")
see$frequency <- NA; see$frequency[see$occurrences_sum >= 1000]<- "Established"
                     see$frequency[see$occurrences_sum <= CIb[3]]<- "Rare" #Species that have not been included into the cleaned database
                     see$frequency[is.na(see$frequency)] <- "Sporadic"
see$taxon <- rownames(see); rownames(see)<-NULL
                     
try <- merge(data_filter, see, by = "taxon")

#Non-indigenous species from the cleaned database
common_species <-  match(dclean$taxon, introd$taxon)
rows <- which(!is.na(common_species))
data_filter <- dclean[rows,]; data_filter$taxon <- as.factor(data_filter$taxon)
nonnativeCL<-data.frame(SpeciesName = levels(data_filter$taxon))
nonnativeCL <- merge(nonnativeCL, introd, by = "taxon")

#Obtain the occurrences of the species defined as cryptogenic and non-indigenous
#select some years
data <- try

#See the occurrences by year
head(data)
xlim <- c(min(data$lon), max(data$lon))
ylim <- c(min(data$lat), max(data$lat))

p<-ggplot(legend=FALSE)+ #make a default plot to include the points in
  coord_equal()+#equla lat and long scale
  xlab("Longitude") + 
  ylab("Latitude")+
  borders(fill="gray44",colour="black") +
  coord_quickmap(xlim= xlim,ylim= ylim)+
  geom_point(data=data,aes(x=lon,y=lat,colour=factor(taxon), shape = factor(frequency)), size = 4)+
  labs(title="Non-indigenous species occurrences across time",colour = "Frequency")+
  facet_wrap(~year, nrow = 4)+
  theme(legend.position = "right",
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90,size=10,color="black",vjust=0.5),
        axis.text.y = element_text(size=10,color="black",vjust=0.5),
        axis.title.x = element_text(size=12,color="black"),
        axis.title.y = element_text(size=12,color="black"),
        plot.title = element_text(size=14),
        panel.border = element_rect(fill=NA,color="black", linetype="solid"),
        strip.text.x = element_text(colour = "black", size = 10),
        strip.background = element_rect(colour = "black", fill = "gray",linetype = "solid"),
        axis.ticks=element_line(colour = "black"))
p

#Group by similar names

# x <- c(levels(data$taxon))
# groups <- list()
# i <- 1
# while(length(x) > 0)
# {
#   id <- agrep(x[1], x, ignore.case = TRUE, max.distance = 0.1)
#   groups[[i]] <- x[id]
#   x <- x[-id]
#   i <- i + 1
# }

#2.0.3: Specific database cleaning according to the case ############
#In this case, we will clean all the species names, and those species that are currently present on the trait database
data<- read.csv("Marenzelleria_1990_2020.txt",header=T,dec=".",sep="\t")
head(data)
data$taxon <- as.factor(data$taxon); taxons <- c(levels(data$taxon))

#Exclude taxonomic groups that are too general, poorly specified
excl <- data[which(data$taxon == "Amphipoda"|data$taxon =="Polychaeta"|data$taxon =="Lepidoptera"|data$taxon == "Aulodrilus pluriseta"|
                     data$taxon == "Nematoda"|data$taxon == "Pallaseopsis quadrispinosa"|data$taxon == "Ceratopogonidae"),]
data <- anti_join(data,excl)

#Substitute the species name that are present in the trait database but do not match due to different names.
oldname = c("Aoridae","Ampharete lindstroemi","Fabulina fabula","Idotea chelipes","Glycera unicornis","Kurtiella bidentata", 
            "Leucothoe lilljeborgii", "Spiophanes kroeyeri","Tellimya ferruginosa","Tellimya tenella",
            "Chironomini|Chironomus aprilinus|Chironomus plumosus|Dicrotendipes|Hydroptilidae|Leptoceridae|Microtendipes pedellus|Monodiamesa|Polypedilum|Procladius|Stempellina|Tanypodinae|Tanytarsus",
            "Clitellata|Potamothrix hammoniensis",
            "Callinera|Cerebratulus|Cyanophthalma obscura|Hubrechtella dubia|Tubulanus polymorphus","Edwardsia|Edwardsia tuberculata",
            "Ophiocten affinis","Philine quadripartita", "Euspira nitida","Gammarus oceanicus","Hermania scabra",
            "Hydrobiidae","Jaera|Jaera albifrons|Jaera albifrons albifrons","Microdeutopus gryllotalpa","Mysis relicta","Philomedes brenda","Polycirrus plumosus",
            "Polydora coeca", "Sphaerodorum flavum", "Tharyx killariensis","Thysanocardia procera","Tritia nitida","Mysiella bidentata")

newname = c("Aoroidae","Ampharete finmarchia (lindstroemi)","Tellina fabula","Idothea chelipes","Glycera unicornis (rouxii)",
            "Mysiella bidentata", "Leucothoe lilljborgii","Spiophanes kroyeri (kroeyeri)","Tellimya (Montacuta) ferruginosa",
            "Tellimya (Montacuta) tenella","Chironomidae","Oligochaeta","Nemertea","Edwardsiidae", "Ophiura","Philine","Euspira","Gammarus",
            "Philine scabra","Hydrobia","Jaera","Microdeutopus","Mysidae","Philomedes","Polycirrus",
            "Dipolydora coeca","Sphaerodorum gracilis (flavum)","Caulleriella (Tharyx) killariensis","Golfingia","Nassarius incrassatus", "Mysella bidentata")

data$taxon<- stri_replace_all_regex(data$taxon,oldname,newname, vectorize_all = F)


#Change only specific names
subset <- data.frame(data[grep("Marenzelleria", data$taxon), ])#select all rows that contain the text that caused the error
data<-anti_join(data,subset) #delete that rows from the database
subset$taxon<- "Marenzelleria" #give them another valid name
data<- rbind(data, subset) #add again the rows with unified names on the database

#write.table(data, "Marenzelleria_1990_2020_CLEAN.txt", row.names = FALSE,  sep="\t")
taxons <- levels(as.factor(data$taxon))

#2.1: Extract information of WORMS and BIOTIC databases ############
#To extract trait information from the MERP trait explorer is needed the AphiaID, is a unique identificator for each species to be found in WORMS
data<-read.csv("Marenzelleria_1990_2020.txt",header=T,dec=".",sep="\t")
head(data)
data$taxon <- as.factor(data$taxon)
data$station <- as.factor(data$station)
#install.packages("worms"); install.packages("worrms")
library(worrms); library(worms)

taxons <- c(levels(data$taxon))#create a vector of the species names to look for
AphiaID <- wm_name2id_(name= non_traits) #obtain multiple AphiaID (this will allow us to download trait information from MERP trait explorer)

# ONLY IF NECESSARY/APPLICABLE
#Change bad names that are giving errors, do it before running Aphia ID function again; check it with the warnings() function if
# it appears the error 400, change names that are giving errors before converting names into levels

# subset <- data.frame(data[grep("Marenzelleria", data$taxon), ])#select all rows that contain the text that caused the error
# data<-anti_join(data,subset) #delete that rows from the database
# subset$taxon<- "Marenzelleria viridis" #give them another valid name
# data<- rbind(data, subset) #add again the rows with unified names on the database

AphiaID <- t(data.frame (AphiaID)) #convert the given list to a dataframe and then transpose it to obtain a column with the values
#AphiaID <- AphiaID[which(AphiaID >=0)]
#write.table(non_traits, "AphiaID.txt", row.names = FALSE,  sep="\t") #Export the IDs to find the traits

#Obtain all the species names to obtain trait data from BIOTIC
taxons_name<- data.frame(taxon = levels(data$taxon))
#write.table(taxons_name, "200sp_names.txt", row.names = FALSE,  sep="\t") #Export the IDs to find the traits

#2.2: Cleaning trait data and working with trait dummy variables #########

#Load trait data, from different databases and also the occurrences data
traitsFULL <- read_excel("traitsFULLdummy.xlsx")
colnames(traitsFULL) <- traitsFULL[1,]; colnames(traitsFULL)[1:3] <- c("phylum", "family", "taxon"); traitsFULL <- traitsFULL[-1,]#select the original names
names<- colnames(traitsFULL); 
#write.table(traitsFULL, "traitsFULLdummy.txt", row.names = FALSE,  sep="\t")

Int_traits<-read.csv("traitsFULLdummy.txt",header=T,dec=".",sep="\t")
head(Int_traits); colnames(Int_traits) <- names
Int_traits[is.na(Int_traits)]<-0
rownames(Int_traits)<-Int_traits$taxon; Int_traits$taxon <- NULL
Int_traits[Int_traits == ""] <- NA; rownames(Int_traits) <- gsub(" sp\\.","", rownames(Int_traits)); rownames(Int_traits) <- gsub(" spp\\.","", rownames(Int_traits))

traitbio<-read.csv("BIOTIC_traitdata_253spp.txt",header=T,dec=".",sep="\t")
head(traitbio)
rownames(traitbio)<-traitbio$SpeciesName; traitbio$SpeciesName <- NULL
incom <- read.csv("Incomplete_traits.txt",header=T,dec=".",sep="\t"); incom <- incom[,1:10];incom <- incom[-c(24:25),]
rownames(incom) <- incom$taxon; incom$taxon <- NULL
nrow <- c(match(rownames(incom),rownames(traitbio))); traitbio <- traitbio[-nrow,]
traitbio <- bind_rows(traitbio, incom)
  
data<- read.csv("Marenzelleria_1990_2020_CLEAN.txt",header=T,dec=".",sep="\t")
head(data)
data$taxon <- as.factor(data$taxon); taxons <- c(levels(data$taxon))


#Add trait information for new species manually for BioticID

# [1] "BioticID"           "taxon""        "ResearchedBy"       "DataSuppliedBy"     "RefereedBy"         "Phylum"            
# [7] "Class"              "Ordr"               "Family"             "Genus"              "Species"            "FoodType"          
# [13] "Size"               "Flexibility"        "Fragility"          "Habit"              "Bioturbator"        "GrowthRate"        
# [19] "DispPotAdult"       "Dependency"         "Sociability"        "Regeneration"       "Toxic"              "LifeSpan"          
# [25] "Maturity"           "GenerationTime"     "ReprodFreq"         "ReprodSeason"       "LarvalSettlePeriod" "Fecundity"         
# [31] "EggSize"            "FertilizationType"  "DispPotLarvae"      "LarvalSettlingTime" "ReprodLocation"     "BiogeographicRange"
# [37] "DepthRange"         "Migratory"          "biozone"            "devmech"            "envpos"             "feedingmethod"     
# [43] "growthform"         "mobility"           "physpref"           "reprodtype"         "salinity"           "substratum"        
# [49] "waterflow"          "waveexp"            "status"

colnames(Int_traits)
# traitraw <- traitraw %>% add_row (SpeciesName = "Monoporeia affinis",Phylum = "Arthropoda", Class = "Malacostraca", Ordr = "Amphipoda", Family = "Pontoporeiidae",
#                                   Genus = "Monoporeia", Species = "affinis",FoodType = "Phytoplankton and detritus", Size ="Small-medium",
#                                   Habit = "Free living", Sociability = "Solitary", Fecundity = "10000 - 46000", devmech = "Oviparous", envpos = "Infaunal",
#                                   feedingmethod = "Passive suspension feeder;Active suspension feeder", mobility = "Burrower", physpref = "Open coast;Estuary",
#                                   reprodtype = "Gonochoristic", salinity = "Low;Variable", substratum = "Mixed", status = "non-indigenous")

#write.table(traitraw, "BIOTIC_traitdata_253spp.txt", row.names = FALSE,  sep="\t") #Write the database with the new species

#unify the species for which we have information in the different trait databases

common_species <- get_common_species(taxons,traitbio); sp1 <- rownames(common_species)
common_species2 <- get_common_species(taxons,Int_traits); sp2 <- rownames(common_species2)

sp <- data.frame(commonsp = c(sp1,sp2)); sp <- unique(sp)
non_common <- common_species[setdiff(rownames(common_species), sp2),]
taxons <- as.data.frame(taxons); non_traits <- setdiff(taxons, sp)


#select traits of interest and those traits that could be classified as dummy
colnames(non_common)
dummy<-dplyr::select(non_common,Size,LifeSpan,Habit,envpos,feedingmethod,ReprodFreq, Bioturbator, devmech,mobility,Habit)
dummy <- format_trait(dummy)
#dummy$taxon <- rownames(dummy); openxlsx::write.xlsx(dummy, 'Incomplete_traits.xlsx')

out <- cSplit(dummy, c(colnames(dummy)),";") #split the columns that are falling in more than one category, to convert them to dummy
temp<-c(out[, which(colSums(is.na(out)) != nrow(out))]); out <- out %>% select(c(names(temp))); out[is.na(out)] <- ""#Select the columns for which we have information
colnames(dummy)

dumm<-dummy_cols(out, select_columns = c(colnames(out)), remove_selected_columns = TRUE)
colnames(dumm) <- sub(".*_","", colnames(dumm)); colnames (dumm) <- ifelse(colnames(dumm) == "", "no.data", colnames(dumm))
dumm[,which(colnames(dumm) == "no.data")] <- NULL

#bind the converted dummy data with the previous data
dumm <- sapply(split.default(dumm, sub("\\.\\d+$", "", names(dumm))), rowSums); names <- colnames(dumm); dumm <- ifelse(dumm > 0, 1, 0)
dumm <- as.data.frame(dumm); colnames(dumm)<- names
dumm <- cbind(phylum = non_common$Phylum, family  = non_common$Family, dumm); rownames(dumm) <- rownames(non_common)

colnames(dumm)[colnames(dumm) == "Lecithotrophic"]<-"Lecitotrophic"

try <- bind_rows(common_species2, dumm); rownames(try) <- c(rownames(common_species2), rownames(non_common)); try[is.na(try)]<-0 #here is really important that the columns names match

#Group the columns by attributes to calculate the disticntiveness for each trait
colnames(try)
TRAITS <- list ("Size" = c("0-10","11-20","21-50","51-100",">100"),
                "Reproductive type"= c("Semelparous-Monotelic","Iteroparous-Polytelic","Semi-continous"),
                "Adult life span" = c("1-3yrs","3-6yrs","6-10yrs",">10yrs"),
                "Developmental mechanism" = c("Fragmentation/Fission","Direct","Lecitotrophic","Planktotrophic","Ovoviviparous"),
                "Environmental position" = c("Deep","Middle","Top","Interface","Epibenthic","Bentho-pelagic","Epilithic"),
                "Living Habits" =  c("Attached","Tube dweller","Burrow dweller","Crevic dweller","Case builder","Free","Parasite/commensal"),
                "Feeding method" = c("Suspension/filter feeder","Deposit feeder (inkl. Both)","Predator","Scavenger","Herbivore","Miner/Borer","Parasite", "Grazer"),
                "Mobility" = c("Sessile","Semi-motile","Motile"),
                "Movement method" = c("No movement","Swimmer","Crawler","Rafter/Drifter/Byssus","Tube-builder","Burrower","Temporary attachment"),
                "Bioturbation" = c("No transport","Diffusive mixing","Surface deposition","Conveyor belt transport","Reverse conveyor belt transport"))

Int_traits <- try; Int_traits[,1:ncol(Int_traits)] <- NULL


Int_traits <- foreach (i = 1:length(TRAITS), .combine = cbind, 
                       .packages = c("funrar","dplyr"), .errorhandling = "remove") %dopar% {
                         
                         trait <- try %>% select(TRAITS[[i]]) #select a subset data for each trait
                         see <- int_distinct(trait)
                         
                         see <- as.data.frame(see[1]); result <- as.data.frame(see[,2]); colnames(result) <- names(TRAITS[i])
                         rownames(result) <- rownames(Int_traits)
                         result
                         
                       } #more optimized loop

# system.time(for (i in 1:length(TRAITS)){ #start loop integrated distinctiveness for each trait
# 
# trait <- try %>% select(TRAITS[[i]]) #select a subset data for each trait
# see <- int_distinct(trait)
# 
# see <- as.data.frame(see[1])
# 
# Int_traits <- cbind(Int_traits, see[,2]); colnames(Int_traits)[i] <- names(TRAITS[i])
# 
# print(paste("Lap",i))
# }) #Non-optimized loop

#Add new trait based on the status of the species (non-indigenous, native)
introd<- read.csv("introduced_sp.txt",header=T,dec=".",sep="\t");
# introd$status <- ifelse(as.numeric(introd$status) == 1, "Non-indigenous", NA);introd$status <- ifelse(is.na(introd$status), "Cryptogenic", "Non-indigenous")
# introd <- introd %>% select(taxon,status); colnames(introd)[1] <- "taxon"
taxons <- as.data.frame(taxons); colnames(taxons)[1] <- "taxon"
status <- merge(taxons, introd, by = "taxon"); status <- unique(status)
status <- status[which(!is.na(status$status)),]; Int_traits$status <- NA
NIS <- Int_traits[which(!is.na(match(rownames(Int_traits), status$taxon))),]; NIS$status <- "Non-indigenous"
nonNIS <- Int_traits[which(is.na(match(rownames(Int_traits), status$taxon))),]; nonNIS$status <- "Native"
Int_traits <- rbind(nonNIS,NIS)
#write.table(Int_traits, "Final_traits.txt", row.names = TRUE,  sep="\t")

#2.3: Obtain the matrix of abundance for all the years and the species #######

data<-read.csv("Marenzelleria_1990_2020_CLEAN.txt",header=T,dec=".",sep="\t")
head(data)
data$taxon <- as.factor(data$taxon);

#Observe the 1st occurrence of all the NIS species on the database, to decide which year we should choose
NISsp <- c("Marenzelleria", "Mya arenaria", "Streblospio benedicti", 
           "Potamopyrgus antipodarum", "Polydora cornuta")
intro_year <- data.frame()

for (i in 1:length(NISsp)){
  data_nis <- data[which(data$taxon == NISsp[i]),]; 
  y_occur <- min(data_nis$year); last <-max(data_nis$year); 
  temp <- data.frame(taxon = NISsp[i], first_yr = y_occur, last_yr = last)
  intro_year <- rbind(intro_year,temp)
}

#Obtain the abundance matrix of species

sample.mat <- dcast(data,station+year+ month +lon+lat~taxon, mean, value.var='wet_weight') #with this function we can create a matrix with the rows on the left side
                                                                                   #of this function and the species as column (right side)
sample.mat[is.na(sample.mat)] <- 0

#write.table(sample.mat, "species_abundance.txt", sep="\t")


#Obtain the relative abundance matrix for a specific year
data<- read.csv("species_abundance.txt",header=T,dec=".",sep="\t", check.names = FALSE)
dyea <- data[which(data$year >= "2008"),]
dyea <- dyea[colSums(!is.na(dyea)) > 0]#delete those species that are not present at any station on that year
dyea[,3:4]<-as.numeric(unlist(dyea[,3:4]))
colnames(dyea) <- gsub(" ", ".", colnames(dyea))

data <- dyea
#write.table(sample.mat, "species_abundance.txt", sep="\t")


#function that assigns cell names based on a raster
result <- assign_cell(data, num = 0.61); data <- result[[1]]

# Plot subsample
ggplot() +  geom_tile(data = data, aes(x=centroid_lon, y=centroid_lat, fill="red"), size = 1)+
  borders(fill="grey",colour="black") + 
  coord_quickmap(xlim=c(xlim[1], xlim[2]), ylim=c(ylim[1], ylim[2])) +
  scale_x_continuous(name='Longitude', breaks=c(0, 4, 8),labels=c('0E', '4E', '8E')) +
  scale_y_continuous(name='Latitude', breaks=c(50, 55),labels=c('50N', '55N'))
# geom_tile(data = fishing.slope, aes(x=Long, y=Lat, fill=slope, width=.8, height= 0.45))+
# geom_text(data=subset(fishing.slope, fishing.slope$p.val.year<0.05), aes(x=Long, y=Lat), label="*", colour = "black", size=5)+

#2.4: Make several plots, maps and ratios ######
data<-read.csv("Marenzelleria_1990_2020_CLEAN.txt",header=T,dec=".",sep="\t", check.names = FALSE)
head(data)

traitraw<-read.csv("Final_traits.txt",header=T,dec=".",sep="\t", check.names = FALSE)
traitraw$taxon <- rownames(traitraw); rownames(traitraw)<- NULL; head(traitraw); 

phylum<-read.csv("phylum.txt",header=T,dec=".",sep="\t"); phylum <- phylum[!duplicated(phylum$taxon),]
traitraw <- merge(traitraw, phylum, by ="taxon")
status <- traitraw %>% select(taxon,phylum, status)



#biomass of all the species
biomass_cum <- as.data.frame(data %>% 
                               group_by(taxon) %>% 
                               summarise_at(.vars = "wet_weight", sum, na.rm = TRUE))

biomass <- merge(biomass_cum, status, by="taxon"); biomass$rel_abun <- (biomass$wet_weight/sum(biomass$wet_weight))*100

#occurrences of species
ocur <- data.frame(count(data, vars = c("taxon", "year"))); ocur <- as.data.frame(ocur %>% 
                                                                            group_by(status) %>% 
                                                                            summarise_at(.vars = "freq", mean, na.rm = TRUE))
ocur <- merge(ocur,status, by = "taxon")

#Donut plot for the relative abundances
#Compute the cumulative percentages (top of each rectangle)- Natives vs NIS
biom <- as.data.frame(biomass %>% 
                               group_by(status) %>% 
                               summarise_at(.vars = c("rel_abun","wet_weight"), sum, na.rm = TRUE))
biom$ymax = cumsum(biom$rel_abun) #NIS and Native

nat <- biomass[which(biomass$status == "Native"),] #Ratio groups of native organisms
nat <- as.data.frame(biomass %>% 
                        group_by(phylum) %>% 
                        summarise_at(.vars = c("wet_weight"), sum, na.rm = TRUE))

nat$rel_abun <- (nat$wet_weight/sum(nat$wet_weight))*100; 
nat_other <- nat[which(nat$rel_abun <=1),]; nat <- nat %>% add_row(phylum = "Other", wet_weight = sum(nat_other$wet_weight),
                                                                   rel_abun = sum(nat_other$rel_abun))
nat <- anti_join(nat, nat_other); nat$ymax = cumsum(nat$rel_abun)

nis <- biomass[which(biomass$status == "Non-indigenous"),] #Ratio groups of native organisms

nis$rel_abun <- (nis$wet_weight/sum(nis$wet_weight))*100; nis$ymax = cumsum(nis$rel_abun)


# Compute the bottom of each rectangle
nis$ymin = c(0, head(nis$ymax, n=-1))
#Comptue label position
nis$labelPosition <- (nis$ymax + nis$ymin) / 2
nis$label1 <- paste0(round(nis$rel_abun, digits = 2), "%")


# Make the plot-NIS VS NATIVE
c <- ggplot(biom, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=status)) +
  geom_rect(size = 1, aes(color = status,
                          fill = after_scale(desaturate(lighten(color, .5), .5)))) +
  geom_text(x=2.2, aes(y = labelPosition,label=status, color=status, fontface = "bold"),size=11, check_overlap = TRUE) + # x here controls label position (inner / outer)
  geom_label(x=3.5, aes(y = labelPosition,label=label1, color = status, 
                        fill = after_scale(desaturate(lighten(color, .9), .9))),size=12) +
  scale_fill_jama() +
  scale_color_jama() +
  coord_polar(theta="y") +
  xlim(c(1, 4)) +
  theme_void() +
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank(),
        panel.border = element_blank()) + theme(legend.position = "none")
c

#ggsave(file="NATNIS.svg", plot=c, width=10, height=8)


#Make the plot biomass Natives and/or NIS
c <- ggplot(nis, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=taxon)) +
  geom_rect(size = 1, aes(color = taxon, fill = after_scale(desaturate(lighten(color, .4), .4)))) +
  scale_fill_uchicago(palette = "light") +
   # x here controls label position (inner / outer)
  geom_text_repel(min.segment.length = 0,  nudge_x = -1.8,
                  aes(x = 3, y = labelPosition,label=taxon, colour = factor(taxon), fontface = "bold",
                      segment.size = 0.8, segment.linetype = 3),
                      size=7, max.overlaps = Inf)+
  geom_label_repel(x=3.7, aes(y = labelPosition, label=label1, color = taxon, 
                              fill = after_scale(desaturate(lighten(color, .9), .9))),
                   size=6, fontface = "bold")+
  scale_color_uchicago(palette = "dark") +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank(),
        panel.border = element_blank()) + theme(legend.position = "none")
c

#ggsave(file="NIS.svg", plot=c, width=10, height=8)

#Make the plot occurrences Natives and/or NIS
c = ggplot(nat, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=phylum)) +
  geom_rect(size = 1, aes(color = phylum, fill = after_scale(desaturate(lighten(color, .4), .4)))) +
  scale_fill_uchicago(palette = "light") +
  # x here controls label position (inner / outer)
  geom_text_repel(min.segment.length = 0,  nudge_x = -2,
                  aes(x = 3, y = labelPosition,label=phylum, colour = factor(phylum), fontface = "bold.italic",
                      segment.size = 0.8, segment.linetype = 3), size=14, max.overlaps = Inf)+
  geom_label_repel(x=3.7, aes(y = labelPosition, label=label1, color = phylum, size  = 12,
                              fill = after_scale(desaturate(lighten(color, .9), .9))),
                   size=5, fontface = "bold")+
  scale_color_uchicago(palette = "dark") +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank(),
        panel.border = element_blank()) + theme(legend.position = "none")

c

#ggsave(file="Plot_name.svg", plot=myPlot, width=10, height=8)

###########################################################################################
#3: Compute functional distances matrix and obtain the trait-space for the community#######!
###########################################################################################
#3.0: Prepare the data for the calculation of the distance matrix #########
traitsFULL<-read.csv("Final_traits.txt",header=T,dec=".",sep="\t", check.names = FALSE)
head(traitsFULL)

data<-read.csv("species_abundance.txt",header=T,dec=".",sep="\t", check.names = FALSE)
head(data)

# Function to get same species in presence-absence matrix as trait matrix########!
traitsFULL<- get_common_species(data,traitsFULL, num = 5)
#select traits of interest
trait <- traitsFULL[,1:9]

#Clean up the trait data.frame by sorting traits and ordering them 
#trait[,1:ncol(trait)] <- lapply(trait[,1:ncol(trait)], as.factor) #apply factor only when its not numeric

#3.1: Compute distance matrix #######
#3.1.1: Using a few traits and Gawids function ########
#THIS IS NOT INTEGRATED DISTINCTIVENESS

data <- try
# standardization function
standr = function(x){(x-min(x))/(max(x)-min(x))}

gaw<-gawdis(data, w.type = "optimized", opti.maxiter = 300)
m.gaw<-as.matrix(gaw)
d.gaw<-as.data.frame(m.gaw)

d.gaw <- d.gaw[-which(is.na(d.gaw)),]
d.gaw[is.na(d.gaw)] <- 0

#The number of species ALWAYS has to be the same between the abundance matrix and the trait data frame
# it has to be reviewed several times trough the process
cols <- match(colnames(d.gaw),rownames(d.gaw)); cols <- which(!is.na(cols))
d.gaw<-d.gaw[,cols] 

standard.gaw<-standr(d.gaw) #standardized distance matrix

Di.gaw<-(apply(d.gaw,1,sum))/(nrow(d.gaw)+1) #mean distance of one species to the other species
mean.gaw<-as.data.frame(Di.gaw)
standard.Di.gaw<-mean.gaw
#rownames(standard.Di.gaw)<-NULL

#standardize mean distances from one species to the others
standard.Di.gaw$Di<-standr(standard.Di.gaw$Di.gaw)

Di.gaw<-data.frame(species = rownames(standard.Di.gaw), Di.gaw = standard.Di.gaw$Di)

spe_index <- as.data.frame(Di.gaw)
spe_index<-spe_index[order(spe_index$Di.gaw,decreasing=T),]
spe_index <- spe_index[-which(spe_index$Di.gaw == 0),]


# obtain the place of the species in percentiles, if in the 4th percentile it means similarity,
# when in the 1st means dissimilarity; the same with the deciles
spe_index <- mutate(spe_index, quartile = as.factor(ntile(spe_index$Di.gaw,4)),decile = as.factor(ntile(spe_index$Di.gaw,10)))
spe_index$species <- as.character(spe_index$sp)

#The number of species has to be the same in the distance dataframe and in the trait database
row <- match(rownames(standard.gaw),spe_index$species);row <- which(!is.na(row));standard.gaw<-standard.gaw[row,] 
cols <- match(colnames(standard.gaw),spe_index$species); cols <- which(!is.na(cols));standard.gaw<-standard.gaw[,cols]

#3.1.2: Using integrated distinctiveness #########

#Here, "TRAIT" is a species x trait information matrix, previously cleaned and with the traits of interest selected
results <- int_distinct(trait); dist_matrix <- as.data.frame(results[2]); Int_Di <- as.data.frame(results[1])
colnames(dist_matrix) <- rownames(dist_matrix); dist_matrix[1:5,1:5]
head(Int_Di)
#write.table(dist_matrix,file="dist_matrix_ALL.txt",sep="\t",row.names = TRUE)
#write.table(Int_Di,file="Int_Di_gaw.txt",sep="\t")
Int_Di_gaw <- Int_Di

#Calculate axis for PCOA
#Int_Di_gaw<- read.csv("Int_Di_gaw.txt",header=T,dec=".",sep="\t")
colnames(Int_Di_gaw)[1] <- "taxon"

spe_index = mutate(Int_Di_gaw, quartile = as.factor(ntile(Int_Di_gaw$int_Di,4)),decile = as.factor(ntile(Int_Di_gaw$int_Di,10)))
#write.table(spe_index,file="spe_index.txt",sep="\t")

#3.2: Calculate the axis of the PCOA to set the trait space + plots (TRAIT SPACE for the whole community) ##########
dist_matrix<- read.csv("dist_matrix_ALL.txt",header=T,dec=".",sep="\t")
gaw.pco<-cmdscale(dist_matrix, eig = T) #change the name here depending if you used gawdis (standard.gaw) 
                                         #or integrated disctinctiveness (standard.gow)

var<-round(gaw.pco$eig/sum(gaw.pco$eig)*100,1)

PCOA1<-data.frame(gaw.pco$points[,1])[,1]
PCOA2<-data.frame(gaw.pco$points[,2])[,1]

sum(gaw.pco$values$Cum_corr_eig)
PCOA<-data.frame(PCOA1,PCOA2)

PCOA_IDi<-cbind(spe_index,PCOA)

row <- match(rownames(trait),spe_index$taxon) #The number of species has to be the same in the distance dataframe and in the trait database
row <- which(!is.na(row));trait<-trait[row,]

gaw.env <- envfit(gaw.pco, trait, permutations = 999, na.rm = T) # With na.rm = TRUE it removes the rows containing NA values

# Plot to see the PCOA biplot of traits influencing the position of species
funct<-spe_index
funct$group[funct$decile==1]<-"common"
funct$group[funct$decile==10]<-"distinct"
funct$group[funct$decile!=10 & funct$decile!=1]<-"intermediate"
grp=as.data.frame(funct$group)
colnames(grp)="group"

df1<-PCOA
df2<-data.frame(gaw.env$vectors$arrows)
df2<-rbind(df2,data.frame(gaw.env$factors$centroids[1:nrow(gaw.env$factors$centroids),]))
#These are the arrows for the biplot of how are influencing the variables on the distinctiveness, here we have a lot of
# categories, with one species falling in multiple levels for a certain trait, 
#this can be corrected or putting only one category or using dummy variables

var 

b<-ggplot() +
  geom_segment(data = df2,aes(x = 0, y = 0, xend = Dim1, yend = Dim2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),type = "closed"),
               linetype=1, size=0.6,colour = "black"
               )+
  geom_text_repel(data = df2,aes(Dim1,Dim2,label=row.names(df2)),max.overlaps=27,size=4)+
  labs(x=paste(var[1],"%"),y=paste(var[2],"%"))+
  geom_hline(yintercept=0,linetype=3,size=1) + 
  geom_vline(xintercept=0,linetype=3,size=1)+
  theme_bw(15)+
  theme(panel.grid=element_blank())
b

#ggsave(file="biplot.svg", plot=b, width=10, height=8)

# Compute the convex hulls to cluster the species
# decile
ch_decile<-convexe_hull(PCOA_IDi$decile)

# PCOA with convex hulls - decile
funct_ch<-ch_decile
funct_ch$group[funct_ch$decile==1]<-"common"
funct_ch$group[funct_ch$decile==10]<-"distinct"
funct_ch$group[funct_ch$decile!=10 & funct_ch$decile!=1]<-"intermediate"
grp_ch=as.data.frame(funct_ch$group)
colnames(grp_ch)="group"

#Here are classified the species on based on if they are distinct or common within the trait-space, by deciles
#In decile 1 there are the common species, and in decile 10 are the most distant

# Plot with the species in the trait-space with convex hulls of most common and distinct species
c <- ggplot(aes(PCOA1,PCOA2,color = funct$group), data=df1)+
  geom_point(size=3)+
  scale_color_manual(values=c("#FA4848","grey","#00AFBB"), limits=c("common","intermediate","distinct"))+
  geom_text(data = PCOA_IDi[which(!is.na(match(PCOA_IDi$taxon, rownames(traitsFULL[which(traitsFULL$status == "Non-indigenous"),])))),],
  aes(label=taxon),
  nudge_x = 0.03, nudge_y = 0.03,
  color = "darkgreen", size = 5)+
  #geom_point(data = PCOA_IDi[which(PCOA_IDi$species == "Mya arenaria"),], colour = "green")+
  geom_point(data = PCOA_IDi[which(!is.na(match(PCOA_IDi$taxon, rownames(traitsFULL[which(traitsFULL$status == "Non-indigenous"),])))),], colour = "green", size = 3)+
  labs(x=paste(var[1],"%"),y=paste(var[2],"%"),color="Functional distinctiveness group")+
  geom_hline(yintercept=0,linetype=3,size=1) + 
  geom_vline(xintercept=0,linetype=3,size=1)+
  theme_bw(15)+
  theme(panel.grid=element_blank())+
  geom_polygon(data = ch_decile[which(ch_decile$decile==1|ch_decile$decile==10),],
                            aes(colour=(grp_ch$group[grp_ch$group=="common"|grp_ch$group=="distinct"])), alpha = 0.1)
c

ggsave(file="Trait_space.svg", plot=c, width=4200, height= 3000, units = "px")

#3.3: Calculate the effect that each trait has on the distance between species #######
taxon<-rownames(trait)
traits.effects<- trait; traits.effects[,1:ncol(traits.effects)] <- NULL

# traits to be removed
trait.removed<-matrix(0,ncol(trait)+1,ncol(trait))
colnames<-c(paste("X",sep ="", seq(1:ncol(trait))))

trait.removed<-as.data.frame(trait.removed)

for (a in (1: ncol(trait.removed))){
  trait.removed[a]<-a
  trait.removed[a,a]<-NA
}
#Here we set a cell as NA, it will be the column that will be removed when computing the distance repeatedly
#Then, the columns of the trait database will be selected based on these "rows" values; thus one trait will be removed each time

## Calculating all the distances removing a trait in each case

traits.effects <- foreach (m= 1:nrow(trait.removed), .combine = cbind,.packages = c("funrar","dplyr"), .errorhandling = "remove") %dopar% {
  traits.selected<-as.vector(t(trait.removed[m,]))
  traits.selected<-na.omit(traits.selected)
  traits.selected<-as.vector(traits.selected)
  # data<-select(trait,habitat,feeding.mode,tl,offspring.size,spawning.type,age.maturity,fecundity,length.infinity,length.max,length.maturity,growth.coefficient) 
  data<-trait
  codex<-traits.selected
  data<-select(data,all_of(codex))  
  # calculating Di with gawdis
  results <- int_distinct(data); teff <- as.data.frame(results[1]); teff <- as.data.frame(teff[,2]); rownames(teff) <- rownames(data)
  teff
}

colnames(traits.effects)[1:(ncol(traits.effects)-1)]<-colnames(trait)#Here we obtain the mean distance of one species to the others when deleting one trait
colnames(traits.effects)[ncol(traits.effects)] <- "None"             #Each column represents the trait that was removed, and the distances without that trait
#write.table(traits.effects,file="traits.effects.txt",sep="\t", row.names = TRUE)


#3.4: Metrics for each single station ########
#Load all the data sets
data<-read.csv("species_abundance.txt",header=T,dec=".",sep="\t", check.names = FALSE)
head(data)

traitraw<-read.csv("Final_traits.txt",header=T,dec=".",sep="\t", check.names = FALSE)
head(traitraw)
trait <- traitraw[,1:9]

dist_matrix <- read.csv("dist_matrix_ALL.txt",header=T,dec=".",sep="\t", check.names = FALSE)

Int_Di<- read.csv("Int_Di_gaw.txt",header=T,dec=".",sep="\t", check.names = FALSE)

common_sp <- get_common_species(data, traitraw, 5); rownames(common_sp) <- gsub("\\.", " ", rownames(common_sp))
cols <- c(rownames(common_sp)); data <- cbind(data[1:4],data[,c(cols)])

#create a column of unique names: stations_year
data$station_year <- paste(data$station,data$year, sep = "_"); data$station_year <- as.factor(data$station_year)
#st_yr <- c(levels(data$station_year))
stationdi <- data.frame(matrix(ncol = 6)); 
colnames(stationdi) <- c("station", "year", "lon","lat","overall_Di", "non_NIS_Di")

#calculate all the metrics for each individual station per year
system.time(for(j in 1:nrow(data)){#start loop inside each station
  
  tryCatch({
    sample.mat<- subset(data, station_year=="509_1995")
    sample.mat <- data[j,]
    cols <- sample.mat[,6:ncol(sample.mat)-1]; num = 1
    
    #calculate NIS mean distinctiveness, distinctiveness x cell and distinctiveness x cell excluding NIS
    if (length(cols[,which(colSums(cols)> 0)]) > 1){cols <- cols[,which(colSums(cols)> 0)]} else {name = colnames(cols)[which(colSums(cols)> 0)]; 
    cols <- as.data.frame(cols[,which(colSums(cols)> 0)]); colnames(cols) <- name}
    
    results <- mean_diss(dist_matrix, cols, traitraw, Int_Di, num)
    Di <- as.data.frame(results[1]); Di_NAT <- data.frame(non_NIS_Di =results[[2]]); Di_ovrll <- data.frame(overall_Di = results[[3]])
    distances <- cbind(Di,Di_NAT,Di_ovrll)
    
    #calculate diversity metrics, species richness, FRichness, FEveness, FDivergence, Rao's Q...
    traitcomm <- trait[c(colnames(cols)),]
    if (ncol(cols) > 1){Fmetrics <- dbFD(traitcomm,cols,messages = FALSE,m = "max",
                                         calc.FRic = TRUE, calc.FDiv = TRUE, stand.FRic = T); 
    Richness <- ncol(cols); 
    FD <- lapply(Fmetrics[c(1:2,5:8)], mean, na.rm = TRUE);
    FDiv <- as.data.frame(do.call(cbind,FD))} else {FDiv <- NA}
    
    metrics <- cbind(distances, Richness, FDiv)
    metrics <- metrics %>% select_if(~!all(is.na(.)))
    
    temp <- data.frame(cbind (station = as.factor(sample.mat$station[1]), year = as.factor(sample.mat$year[1]),
                              lon = sample.mat$lon[1],lat = sample.mat$lat[1], metrics));
    stationdi <- dplyr::bind_rows(stationdi, temp)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  print(paste(round((j/length(st_yr))*100,digits = 3),"%")); rm(sample.mat, traitcomm,results,FD,FDiv,Di,Di_NAT,
                                                                Di_ovrll, distances, metrics, Fmetrics, cols)
  

})#end loop inside each station

system.time(stationDi <- foreach(j = 1:nrow(data), .packages = c("dplyr", "vegan",
                                                                    "data.table","raster",
                                                                    "doParallel","foreach", "FD"), .errorhandling = "pass")%dopar%{#start loop inside each station
  
  tryCatch({
    #sample.mat<- subset(data, station_year==st_yr[j])
    sample.mat <- data[j,]
    cols <- sample.mat[,6:ncol(sample.mat)-1]; num = 1
    
    #calculate NIS mean distinctiveness, distinctiveness x cell and distinctiveness x cell excluding NIS
    if (length(cols[,which(colSums(cols)> 0)]) > 1){cols <- cols[,which(colSums(cols)> 0)]} else {name = colnames(cols)[which(colSums(cols)> 0)]; 
    cols <- as.data.frame(cols[,which(colSums(cols)> 0)]); colnames(cols) <- name}
    
    results <- mean_diss(dist_matrix, cols, traitraw, Int_Di, num)
    Di <- as.data.frame(results[1]); Di_NAT <- data.frame(non_NIS_Di =results[[2]]); Di_ovrll <- data.frame(overall_Di = results[[3]])
    distances <- cbind(Di,Di_NAT,Di_ovrll)
    
    #calculate diversity metrics, species richness, FRichness, FEveness, FDivergence, Rao's Q...
    traitcomm <- trait[c(colnames(cols)),]
    if (ncol(cols) > 1){Fmetrics <- dbFD(traitcomm,cols,messages = FALSE,m = "max",
                                         calc.FRic = TRUE, calc.FDiv = TRUE, stand.FRic = T); 
    #Richness <- ncol(cols); 
    FD <- lapply(Fmetrics[c(1:3,5:8)], mean, na.rm = TRUE);
    FDiv <- as.data.frame(do.call(cbind,FD))} else {FDiv <- NA}
    
    Richness <- ncol(cols)
    
    metrics <- cbind(distances, FDiv)
    metrics <- metrics %>% select_if(~!all(is.na(.)))
    
    temp <- data.frame(cbind (station = as.factor(sample.mat$station[1]), year = as.factor(sample.mat$year[1]),
                              lon = sample.mat$lon[1],lat = sample.mat$lat[1], Richness = as.numeric(Richness), metrics))}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  print(paste(round((j/length(st_yr))*100,digits = 3),"%")); rm(sample.mat, traitcomm,results,FD,FDiv,Di,Di_NAT,
                                                                Di_ovrll, distances, metrics, Fmetrics, cols)
  temp
})
stationdi <- do.call(bind_rows, stationDi)

#write.table(stationdi,file="Di_metrics_station.txt",sep="\t", row.names = TRUE)

#3.5: Distinctiveness of non-indigenous species based on cell size; local and regional scale ################
data<-read.csv("species_abundance.txt",header=T,dec=".",sep="\t", check.names = FALSE)
head(data)
dyea <- data[which(data$year >= "2008"),]
dyea <- dyea[colSums(!is.na(dyea)) > 0]#delete those species that are not present at any station on that year
dyea[,3:4]<-as.numeric(unlist(dyea[,3:4]))
data <- dyea

traitraw<-read.csv("Final_traits.txt",header=T,dec=".",sep="\t", check.names = FALSE)
head(traitraw)
trait <- traitraw[,1:9]

dist_matrix <- read.csv("dist_matrix_ALL.txt",header=T,dec=".",sep="\t", check.names = FALSE)

Int_Di<- read.csv("Int_Di_gaw.txt",header=T,dec=".",sep="\t", check.names = FALSE)

common_sp <- get_common_species(data, traitraw, 5); rownames(common_sp) <- gsub("\\.", " ", rownames(common_sp))
cols <- c(rownames(common_sp)); data <- cbind(data[1:4],data[,c(cols)])

sizes <- c(seq(from = 0, to = 1, by = 0.05)); sizes <- sizes[-1]#set the sequence of values by which increment cell size
sizes1.1 <- c(seq(from = 1, to= 5, by = 0.5)); sizes <- c(sizes, 1.5)
# uniq <- data[!duplicated(data[,1]),]; dim <- dim(uniq) #obtain the unique names of locations
# sizes <- append(sizes, 1/dim[1], after = 0) # Add a last value corresponding to the number of unique locations in order 
# that at the end you only have one cell containing all the stations
ncell <- data.frame()

#With this loop we will obtain a table that summarizes the % of cell size reductions and the number of cells that there are in the database
ncell <- foreach (i = 1:length(sizes), .combine = rbind, 
                  .packages = c("raster","dplyr"), .errorhandling = "pass") %dopar% {#start loop cell
  result <- assign_cell(data, num = sizes[i]); trial<- result[[1]]; cellsize <- result[[2]]
  number <- data.frame(ncell = length(levels(trial$cell)), size_reduction = sizes[i])
  number <- cbind(number, cellsize)
} #end loop cell

ncell <- ncell[!duplicated(ncell[,1]),]; sizes_end <- c(ncell$size_reduction)

#write.table(ncell,file="ncell.txt",sep="\t", row.names = TRUE)

#Dixcell <- data.frame(matrix(ncol = 4)); colnames(Dixcell) <- c("size_reduction", "cell", "overall_Di", "non_NIS_Di")

#This loop is composed of certain steps, on the first one it assigns the stations/sampling locations from the data to a certain cell
#Then, we determine which would be the minimum number of samples needed to obtain a certain % of existing species in the database by adjusting a
#Species accumulation curve (SAC). Finally, we select a specific cell from the whole database and in the last step, we resample a minimum number of 
#locations based on the results from SAC and calculate several metrics such as Di, Fev... inside each cell. We iterate the resampling process around 100
#times and then we obtain the mean for each metric. 

system.time(Dixcell <- foreach (t = 1:length(sizes_end),
                                .packages = c("dplyr", "vegan","data.table","raster","doParallel","foreach", "FD"), 
                                .errorhandling = "pass")%dopar%{#start loop metrics x cell
  
  result <- assign_cell(data, num = sizes_end[t]); trial <- result[[1]]
  trial <- trial %>% select(station, cell, year, centroid_lon, centroid_lat, xMin, xMax, yMin, yMax, 5:(ncol(trial)-7))
  cells <- c(levels(trial$cell))
  
  mat.deriv <- foreach (i = 1:length(cells),.combine = rbind, 
                        .packages = c("vegan", "data.table"), .errorhandling = "pass") %dopar% {#start loop SAC's
    
    identity <- as.character(cells[i])
    #mat.deriv$year[i] <- celly$year[i]
    
    tryCatch({
      # Transform into a matrix for SAC function
      sample.mat<- trial[which(trial$cell == cells[i]),]
      
      #sample.mat <- dcast(sub.surv, cell+station+lon+lat+year~taxon, mean, value.var='abundance') #with this function it is transformed into an abundance matrix
      #rownames(sample.mat) <- sample.mat[,1]
      #sample.mat[,1] <- NULL; 
      sample.mat[is.na(sample.mat)] <- 0; sample.mat <- sample.mat[,10:ncol(sample.mat)]
      
      # Get the unobserved total number of species
      unobs <- specpool(sample.mat)
      
      # Arrondir les cpue
      #sample.mat <- round(sample.mat)
      
      sac <- specaccum(sample.mat, method = "random", permutations = 100) #Here it fits 100 models of species accummulation curves
      #rare <- specaccum(sample.mat, method='rarefaction', xvar='individuals')
      model9 <- fitspecaccum(sac,  model='michaelis-menten', method='random', control = ctrl)
      #model9r <- fitspecaccum(rare, model='michaelis-menten')
      fitted9 <- as.data.frame(fitted(model9)) #All the models fitted as a dataframe
      aic9 <- sapply(model9$models, AIC) #See the AIC of all the fitted SAC
      num9 <- which(aic9==min(aic9)) #Choose the best model by minimum AIC
      result9 <- fitted9[,num9]
      means9 <- apply(fitted9, 1, FUN=mean)
      
      #Here Vm accounts for number of species found, only reaching the maximum number of species on the cell
      #Then, we calculate how many species do we obtain as we increase the number of samples
      w <- as.data.frame(coef(model9))
      v <- apply(w, 1, FUN=mean)
      x <- seq(from=1, to=nrow(fitted9), by=1)
      
      y <- MM(v[1], v[2], x)
      
      
      ### Get first derivative at max nhauls
      y2 <- MMprime(v[1],v[2],x)
      
      # plot(rare, xlab='', ylab='')
      # plot(model9r, add=TRUE, col='grey')
      # plot(rare, xlab='', ylab='', add=TRUE)
      # points(y~x, type='l', lwd=5, col='orange')
      
      ### Get all values
      if(length(y2)>0){
        yprim <- y2[length(y2)]
        nb.hauls <- nrow(sample.mat) #Total number of locations/inventories/samples
        chao <- round(unobs$chao)
        sr <- unobs$Species
        asympt <- v[1]
        
        ### Get nb of hauls necessary to reach % of chao
        limit <- 0.8*round(v[1])
        obj <- which(y>limit) #From which number of samples do we obtain a 80% of the species
        if(length(obj)>0){limit0.8 <- obj[1]}
        else{limit0.8 <- NA}
        
        limit <- 0.65*round(v[1])
        obj <- which(y>limit)
        if(length(obj)>0){limit0.65 <- obj[1]}
        else{limit0.65 <- NA}
        
        limit <- 0.5*round(v[1])
        obj <- which(y>limit)
        if(length(obj)>0){limit0.5 <- obj[1]}
        else{limit0.5 <- NA}
        
        #Keep or not if y'<0.5 meaning decreasing trend
        x2 <- seq(from=1, to=100, by=1)
        y0 <- MM(v[1], v[2], x2)
        yprim2 <- MMprime(v[1],v[2],x2)
        lets <- which(yprim2<0.5)
        if(length(lets)>0){
          if(lets[1] <= length(x)){keep <- 'YES'}
          if(lets[1] > length(x)){keep <- 'NO'}
          #Min number of hauls necessary
          nb.hauls.min  <-lets[1]
        }
      }
      if(length(lets)==0){keep <- 'NO'
      nb.hauls.min <- NA}
      
      temp <- data.frame(identity=identity, yprim=yprim, nb.hauls=nb.hauls, nb.hauls.min=nb.hauls.min, keep=keep, sr=sr, chao=chao, limit0.8=limit0.8,
                         limit0.65=limit0.65,limit0.5=limit0.5, asympt=asympt)
      rm(fitted9, model9, aic9, lets, ids, i, num9, result9, yprim2, y0, y, x, x2, v, sample.mat, samples, sac, sub.surv, 
         unobs, w, means9, limit, obj)
      temp
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    
  } #end loop SACs
  
  mat.deriv <- mat.deriv[!is.na(mat.deriv$limit0.65),] %>% select(-yprim)
  
  mat.deriv <- mat.deriv[which(mat.deriv$keep == "YES"),]
  
  celldi <- data.frame(matrix(ncol = 12)); colnames(celldi) <- c("size_reduction", "cell", "centroid_lon","centroid_lat","xMin", "xMax",
                                                                 "yMin", "yMax", "chao", "asympt", "overall_Di", "non_NIS_Di")
  
  for(j in 1:nrow(mat.deriv)){#start loop inside each individual cell
    
    sample.mat<- subset(trial, cell==mat.deriv$identity[j])
    sample.mat$cell <- as.numeric(as.character(sample.mat$cell)); sample.mat$station <- as.numeric(as.character(sample.mat$station))
    
    q <- 1; NIS_dist <- data.frame()
    
    repeat {
      #print(q)
      #calculate NIS mean distinctiveness, distinctiveness x cell and distinctiveness x cell excluding NIS
      random.h<- sample_n(sample.mat, mat.deriv$limit0.65[j]); cols <- random.h[,which(colSums(random.h)> 0)]
      results <- mean_diss(dist_matrix, cols, traitraw, Int_Di)
      Di <- as.data.frame(results[1]); Di_NAT <- data.frame(non_NIS_Di =results[[2]]); Di_ovrll <- data.frame(overall_Di = results[[3]])
      distances <- cbind(Di,Di_NAT,Di_ovrll)
      
      #calculate diversity metrics, species richness, FRichness, FEveness, FDivergence, Rao's Q...
      spss <- as.data.frame(cols[,10:ncol(cols)]); traitcomm <- trait[c(colnames(spss)),]
      spss <- as.data.frame(spss[which(rowSums(spss)> 0),])
      if (ncol(spss) > 1){Fmetrics <- dbFD(traitcomm,spss,messages = FALSE,m = "max",calc.FRic = FALSE, calc.FDiv = FALSE); #Richness <- ncol(spss); 
      FD <- lapply(Fmetrics[c(3:5)], mean, na.rm = TRUE);
      FDiv <- as.data.frame(do.call(cbind,FD))} 
      else {FDiv <- NA}

      metrics <- cbind(distances, FDiv)
      
      NIS_dist <- bind_rows(NIS_dist,metrics)
      #NIS_dist <- NIS_dist[,which(!all(is.na(NIS_dist)))]
      NIS_dist <- NIS_dist %>% select_if(~!all(is.na(.)))
      
      q <- q+1
      
      if (q > 50){break}
    }
    
    #calculate diversity metrics, Richness, FRichness and FDivergence
    tryCatch({spss <- as.data.frame(sample.mat[,10:ncol(sample.mat)]); spss <- spss[,which(colSums(spss)> 0)]
    traitcomm <- trait[c(colnames(spss)),]; spss <- as.data.frame(spss[which(rowSums(spss)> 0),])
    Richness <- ncol(spss)
    Fmetrics <- dbFD(traitcomm,spss,messages = FALSE,m = "max", stand.FRic = TRUE);
    FD <- lapply(Fmetrics[c(3,6)], mean, na.rm = TRUE);
    FRic <- as.data.frame(do.call(cbind,FD))}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    FRic <- cbind(FRic, Richness)
    
    NIS_dist <-  NIS_dist %>% summarise_all(mean, na.rm = TRUE); NIS_dist <- cbind(NIS_dist,FRic)
    
    temp <- data.frame(cbind (size_reduction = as.character(sizes_end[t]), cell = as.character(cells[j]),
                              centroid_lon = sample.mat$centroid_lon[1],centroid_lat = sample.mat$centroid_lat[1],
                              xMin = sample.mat$xMin[1],xMax = sample.mat$xMax[1],yMin = sample.mat$yMin[1],
                              yMax = sample.mat$yMax[1], chao = mat.deriv$chao[j],asympt = mat.deriv$asympt[j], NIS_dist));
    celldi <- dplyr::bind_rows(celldi, temp)
    print(paste ("lap",j))
    
  } #end loop inside each individual cell
  
  rm(results, Di, Di_NAT, Di_ovrll, distances, Fmetrics, FD, FDiv, metrics, Richness, FRic)
  
  print(paste(t/length(sizes_end)*100,"%"))
  rm(mat.deriv)
  
  celldi
  
}) #end loop metrics x cell 

#Here we obtain several indexes x cell x different cell sizes

#write.table(try,file="dixcell.txt",sep="\t", row.names = TRUE)

Dixcell <- Dixcell[rowSums(is.na(Dixcell)) != ncol(Dixcell), ]

meanDi <- as.data.frame(Dixcell %>%
                          group_by(size_reduction) %>% 
                          summarise_all(mean, na.rm = TRUE))

#meanDi <- melt(meanDi, id = c("size_reduction")); colnames(meanDi) <- c("cell_reduction", "taxon","meanDi"); meanDi$cell_reduction <- as.factor(meanDi$cell_reduction)

#3.6: Taxon scarcity overall ################

spe_index<-read.csv("spe_index.txt",header=T,dec=".",sep="\t", check.names = FALSE)

data<-read.csv("species_abundance.txt",header=T,dec=".",sep="\t", check.names = FALSE)
head(data)
common_sp <- get_common_species(data, traitsFULL); cols <- c(rownames(common_sp)); data <- cbind(data[1:4],data[,cols])

N<-nrow(spe_index)
ab_tot_sp <- apply(data[,5:(ncol(data))], 2, sum,na.rm=T)#total abundance of species 
ab_tot<-sum(ab_tot_sp)
ab_rel_sp<-ab_tot_sp/ab_tot

Abi<-data.frame(taxon=spe_index$taxon,decile=spe_index$decile,abi=ab_rel_sp); rownames(Abi) <- NULL
#Abi<-Abi[which(Abi$decile==1|Abi$decile==10),]

med_Abi<- ddply(Abi, "decile", summarise, grp.med=median(abi))#Split data frame, apply function, and return results in a data frame
6.995030e-08/1.068977e-06
mean_Abi<- ddply(Abi, "decile", summarise, grp.mean=mean(abi))
0.0000108675/0.0012587719
sum_Abi<- ddply(Abi, "decile", summarise, grp.sum=sum(abi))
0.0004564351/0.0541271919

scarcity<-data.frame(taxon=spe_index$taxon,ab_rel_sp)
scarcity$Sci<-NA
for (k in 1:nrow(scarcity)){
  scarcity$Sci[k]<-exp(-N*log(2)*scarcity$ab_rel_sp[k])
}

spe_index<-full_join(spe_index,scarcity)

ploty<-spe_index[which(spe_index$decile==1|spe_index$decile==10),]

med<- ddply(ploty, "decile", summarise, grp.med=median(Sci))

b<-ggplot(aes(x=Sci,col=as.factor(quartile)),data=spe_index)+
  theme_classic()+
  geom_density(aes(y=..scaled..))+
  labs(x="Scarcity",y="Normalized\number of species",color="Functional distinctiveness group",labs=c("Common (D1)","Distinct (D10)"))+
  theme_bw(15)+
  theme(panel.grid=element_blank(),legend.position = "bottom")+
  scale_color_discrete(labels=c(levels(as.factor(spe_index$quartile))))

b
#3.7: Taxon restrictedness by different cell sizes #################

data<-read.csv("species_abundance.txt",header=T,dec=".",sep="\t", check.names = FALSE)
head(data)
dat <- read.csv("ncell.txt",header=T,dec=".",sep="\t", check.names = FALSE)
spe_index<-read.csv("spe_index.txt",header=T,dec=".",sep="\t", check.names = FALSE)

system.time(Ri.mean_end <- foreach (t = 1:length(sizes_end), .combine = cbind, 
                                    .packages = c("dplyr","tidyr", "vegan","data.table","raster","doParallel","foreach"), 
                                    .errorhandling = "pass")%dopar%{#start loop cell
  
  result <- assign_cell(data, num = sizes_end[t]); trial <- result[[1]]
  trial <- trial %>% select(station, cell, year, centroid_lon, centroid_lat, 5:(ncol(trial)-7)) #here select the species that you have on the matrix; these MUST be
  
  #species for which we have trait information
  cells <- c(levels(trial$cell))
  
  mat <- foreach (i = 1:length(cells),.combine = rbind, .packages = c("vegan", "data.table"), .errorhandling = "pass") %do% {
    
    identity <- as.character(cells[i])
    #mat.deriv$year[i] <- celly$year[i]
    
    tryCatch({
      # Transform into a matrix for SAC function
      sample.mat<- trial[which(trial$cell == cells[i]),]
      
      #sample.mat <- dcast(sub.surv, cell+station+lon+lat+year~taxon, mean, value.var='abundance') #with this function it is transformed into an abundance matrix
      #rownames(sample.mat) <- sample.mat[,1]
      #sample.mat[,1] <- NULL; 
      sample.mat[is.na(sample.mat)] <- 0; sample.mat <- sample.mat[,8:ncol(sample.mat)]
      
      # Get the unobserved total number of species
      unobs <- specpool(sample.mat)
      
      # Arrondir les cpue
      sample.mat <- round(sample.mat)
      
      sac <- specaccum(sample.mat, method = "random", permutations = 100) #Here it fits 100 models of species accummulation curves
      #rare <- specaccum(sample.mat, method='rarefaction', xvar='individuals')
      model9 <- fitspecaccum(sac,  model='michaelis-menten', method='random', control = ctrl)
      #model9r <- fitspecaccum(rare, model='michaelis-menten')
      fitted9 <- as.data.frame(fitted(model9)) #All the models fitted as a dataframe
      aic9 <- sapply(model9$models, AIC) #See the AIC of all the fitted SAC
      num9 <- which(aic9==min(aic9)) #Choose the best model by minimum AIC
      result9 <- fitted9[,num9]
      means9 <- apply(fitted9, 1, FUN=mean)
      
      #Here Vm accounts for number of species found, only reaching the maximum number of species on the cell
      #Then, we calculate how many species do we obtain as we increase the number of samples
      w <- as.data.frame(coef(model9))
      v <- apply(w, 1, FUN=mean)
      x <- seq(from=1, to=nrow(fitted9), by=1)
      
      y <- MM(v[1], v[2], x)
      
      
      ### Get first derivative at max nhauls
      y2 <- MMprime(v[1],v[2],x)
      
      # plot(rare, xlab='', ylab='')
      # plot(model9r, add=TRUE, col='grey')
      # plot(rare, xlab='', ylab='', add=TRUE)
      # points(y~x, type='l', lwd=5, col='orange')
      
      ### Get all values
      if(length(y2)>0){
        yprim <- y2[length(y2)]
        nb.hauls <- nrow(sample.mat) #Total number of locations/inventories/samples
        chao <- round(unobs$chao)
        sr <- unobs$Species
        asympt <- v[1]
        
        ### Get nb of hauls necessary to reach % of chao
        limit <- 0.8*round(v[1])
        obj <- which(y>limit) #From which number of samples do we obtain a 80% of the species
        if(length(obj)>0){limit0.8 <- obj[1]}
        else{limit0.8 <- NA}
        
        limit <- 0.65*round(v[1])
        obj <- which(y>limit)
        if(length(obj)>0){limit0.65 <- obj[1]}
        else{limit0.65 <- NA}
        
        limit <- 0.5*round(v[1])
        obj <- which(y>limit)
        if(length(obj)>0){limit0.5 <- obj[1]}
        else{limit0.5 <- NA}
        
        #Keep or not if y'<0.5 meaning decreasing trend
        x2 <- seq(from=1, to=100, by=1)
        y0 <- MM(v[1], v[2], x2)
        yprim2 <- MMprime(v[1],v[2],x2)
        lets <- which(yprim2<0.5)
        if(length(lets)>0){
          if(lets[1] <= length(x)){keep <- 'YES'}
          if(lets[1] > length(x)){keep <- 'NO'}
          #Min number of hauls necessary
          nb.hauls.min  <-lets[1]
        }
      }
      if(length(lets)==0){keep <- 'NO'
      nb.hauls.min <- NA}
      
      temp <- data.frame(identity=identity, yprim=yprim, nb.hauls=nb.hauls, nb.hauls.min=nb.hauls.min, keep=keep, sr=sr, chao=chao, limit0.8=limit0.8,
                         limit0.65=limit0.65,limit0.5=limit0.5, asympt=asympt)
      temp
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
    
    rm(fitted9, model9, aic9, lets, ids, i, num9, result9, yprim2, y0, y, x, x2, v, sample.mat, samples, sac, sub.surv, 
       unobs, w, means9, limit, obj)
    
  }
  
  mat <- mat[!is.na(mat$limit0.65),] %>% select(-yprim, -asympt)
  
  data.k<-dplyr::select(trial,cell,6:(ncol(trial))); sel.cell <- c(levels(as.factor(mat$identity))); Ktot <- length(sel.cell); k = 1
  
  Ri.sensitivity<-data.frame(taxon=colnames(trial[6:ncol(trial)]))
  Ri.mean<-data.frame(taxon=colnames(trial[6:ncol(trial)]))
  
  system.time(repeat{
    print (k)
    l<-split(data.k,data.k$cell) #divide the database into a list of different databases each one for one location/cell
    for (i in 1:length(sel.cell)){
      l[[i]]<-l[[i]] %>% slice_sample(n=mat$limit0.65[i]) # randomly selecting sampling needed per rectangle
    }
    dat<-do.call(rbind,l) # creating a dataframe based on the list with one haul per rectangle
    dat<-pivot_longer(dat,2:ncol(dat)) #puts in one column all the species in a certain cell and their abundance in other column
    ab.tot.rect.sp<-aggregate(value~cell+name,dat,sum,na.rm=T)#sums all the abundances for each species in each cell
    presence.rect.sp<-cbind(ab.tot.rect.sp[,1:2],pa=decostand(ab.tot.rect.sp$value, method = "pa"))#transforms abundance data into presence/absence
    Ki<-aggregate(pa~name,presence.rect.sp,sum,na.rm=T) #in how many cells/locations do each species appear
    Ri<-Ki
    for (j in (1:nrow(Ki))){
      Ri$value[j]<-1-(Ri$pa[j]/Ktot)
    }
    Ri.sensitivity<-cbind(Ri.sensitivity,Ri$value)
    k <- k+1; if (k >= 100){break}
  })
  rm(mat)
  Ri.mean<-cbind(Ri.mean,apply(Ri.sensitivity[,-1],1,mean,na.rm=T)); value <- data.frame(Ri.mean[,2]); colnames(value)[1] <- sizes_end[t]
  rownames(value) <- Ri.sensitivity$taxon; 
  value
  #colnames(Ri.mean)<-c("taxon","Ri_mean_100")
  
})


#write.table(Ri.mean_end,file="Restrictedness.txt",sep="\t", row.names = TRUE)
Ri.mean <- Ri.mean_end
#dat<-as.data.frame(do.call(cbind,Ri.mean))
Ri.mean$taxon <- rownames(Ri.mean); rownames(Ri.mean) <- NULL
Ri.mean <- Ri.mean[,c(7,1:6)]; #Ri.mean <- Ri.mean %>% select(taxon, '0.81')
Ri.mean <- pivot_longer(Ri.mean, 1:ncol(Ri.mean))

spe_index<-full_join(spe_index,Ri.mean); spe_index <- spe_index[1:177,]

ploty<-spe_index[which(spe_index$decile==1|spe_index$decile==10),]
med<- ddply(ploty, "quartile", summarise, grp.med=median(`0.21`))

a<-ggplot(aes(x=`0.21`,col=as.factor(quartile)),data=spe_index)+
  theme_classic()+
  geom_density(aes(y=..scaled..))+
  labs(x="Integrated taxonomic restrictedness",y="Normalized\nnumber of species",color="Functional distinctiveness group",labs=c("Common (D1)","Distinct (D10)"))+
  theme_bw(15)+
  theme(panel.grid=element_blank(),legend.position = "bottom")
  #scale_color_discrete(labels=c("Common (1)","Distinct (4)"))
a

#3.8: Plots from the several results #########
data<- read.csv("dixcell.txt",header=T,dec=".",sep="\t", check.names = FALSE)


######################################################################
#4: Use environmental data to describe distinctiveness patterns#######!
######################################################################
#4.1: Determine only one cell size to extract the environmental variables for the sea bottom ##########
#area coordinates: lon = 10.56; 23.94 lat = 65.89, 54.99
#Load the database with the coordinates
setwd('C:/Users/avipo/Desktop/PhD/Project 1/MarenzelleriaSweden')
data <- read.csv("Di_metrics_station.txt",header=T,dec=".",sep="\t", check.names = FALSE)

season<-read.csv("Marenzelleria_1990_2020_CLEAN.txt",header=T,dec=".",sep="\t")
sample.mat <- dcast(season,station+year+ month +lon+lat~taxon, mean, value.var='wet_weight') #with this function we can create a matrix with the rows on the left side

season <- sample.mat %>% select(station, year, month,lon, lat); season$stat_year <- paste(season$station, season$year, sep = "_")
season <- season[which(!duplicated(season$stat_year)),]; season <- season %>% select(station, year, month)
env <- data %>% select(station, year, lon, lat)
env <- merge(env, season, by = c("station", "year")); env$month <- ifelse(env$month > 9, env$month, sprintf("%02d", env$month))
env$year_month <- paste(env$year,env$month, sep = ".")

#Prepare the dataset and extract the coordinates for which we want information
env <- env[which(env$year >= "2005"),] #select the years of interest
datasp <- env[,3:4]; coordinates(datasp)=~lon+lat #Select the coordinates to extract
d2020 <- env[which(env$year == 2020),]; env <- env[which(!env$year == 2020),] #the raster for 2020 data has a different extent, so it has to be separated

#Dissolved oxygen at sea bottom - 2005-2019
setwd("./env_data")
oxy05 <- 'o2_o2bot_2005-2010.nc'; oxy05 <- brick(oxy05, var = "o2b") #open the environmental variables in a raster
oxy10 <- 'o2_o2bot_2011-2015.nc'; oxy10 <- brick(oxy10, var = "o2b")
oxy15 <- 'o2_o2bot_2016-2019.nc'; oxy15 <- brick(oxy15, var = "o2b")
oxy <- raster::stack(oxy05,oxy10,oxy15); #join all of them
rna <- is.na(oxy); plot(rna[[1]])
system.time(var<- raster::extract(oxy,datasp, method = 'bilinear', small = T, fun= mean, na.rm = TRUE, df = T, exact = T, cellnumbers = F));rm(oxy, oxy05,
                                                                                                                                                              oxy10,oxy15)
#var<- exactextractr::exact_extract(oxy,datasp, method = 'bilinear', small = T, na.rm = TRUE, df = T, exact = T, cellnumbers = F); 
colnames(var) <- substr(colnames(var), 2,8); 
for (i in 1:nrow(env)){
a <- as.numeric(var[i,which(colnames(var) == env$year_month[i])])
env$bot_oxy[i] <- ifelse(length(a) == 0,
                           0, var[i,which(colnames(var) == env$year_month[i])])
}

#Dissolved oxygen at sea bottom - 2020
oxy20 <- 'o2_o2bot_2020.nc'; oxy20 <- brick(oxy20, var = "o2b")
var<- raster::extract(oxy20,datasp, method = 'bilinear', small = T, na.rm = T, df = T, exact = T, cellnumbers = F);rm(oxy20)
colnames(var) <- substr(colnames(var), 2,8);
for (i in 1:nrow(d2020)){
a <- as.numeric(var[i,which(colnames(var) == d2020$year_month[i])])
d2020$bot_oxy[i] <- ifelse(length(a) == 0,
                             0, var[i,which(colnames(var) == d2020$year_month[i])])
}

#Salinity at sea bottom- 2005-2019
sal05 <- 'temp_sal_salbot_2005-2010.nc'; sal05 <- brick(sal05, var = "sob")
sal10 <- 'temp_sal_salbot_2011-2015.nc'; sal10 <- brick(sal10, var = "sob")
sal15 <- 'temp_sal_salbot_2016-2019.nc'; sal15 <- brick(sal15, var = "sob")
sal <- raster::stack(sal05,sal10,sal15); 
#coordinates(datasp)=~lon+lat
var<- raster::extract(sal,datasp,method = 'bilinear', small = T, na.rm = T, df = T, exact = T, cellnumbers = F); rm(sal, sal05, sal10,
                                                                                                                  sal15)
colnames(var) <- substr(colnames(var), 2,8); 
for (i in 1:nrow(env)){
  a <- as.numeric(var[i,which(colnames(var) == env$year_month[i])])
  env$bot_sal[i] <- ifelse(length(a) == 0,
                             0, var[i,which(colnames(var) == env$year_month[i])])
}

#Salinity at sea bottom- 2020
sal20 <- 'temp_sal_salbot_2020.nc'; sal20 <- brick(sal20, var = "sob")
var<- raster::extract(sal20,datasp, method = 'bilinear', small = T, na.rm = T, df = T, exact = T, cellnumbers = F);rm(sal20)
colnames(var) <- substr(colnames(var), 2,8);
for (i in 1:nrow(d2020)){
  a <- as.numeric(var[i,which(colnames(var) == d2020$year_month[i])])
  d2020$bot_sal[i] <- ifelse(length(a) == 0,
                               0, var[i,which(colnames(var) == d2020$year_month[i])])
}

#Temperature at sea bottom - 2005-2019
t05 <- 'temp_sal_salbot_2005-2010.nc'; t05 <- brick(t05, var = "bottomT")
t10 <- 'temp_sal_salbot_2011-2015.nc'; t10 <- brick(t10, var = "bottomT")
t15 <- 'temp_sal_salbot_2016-2019.nc'; t15 <- brick(t15, var = "bottomT")
tem <- raster::stack(t05,t10,t15); 
var<- raster::extract(tem,datasp, method = 'bilinear', small = T, na.rm = T, df = T, exact = T, cellnumbers = F); rm(t05, tem, t10, t15)
colnames(var) <- substr(colnames(var), 2,8); 
for (i in 1:nrow(env)){
  a <- as.numeric(var[i,which(colnames(var) == env$year_month[i])])
  env$bot_T[i] <- ifelse(length(a) == 0,
                             0, var[i,which(colnames(var) == env$year_month[i])])
}

#Temperature at sea bottom- 2020
t20 <- 'temp_sal_salbot_2020.nc'; t20 <- brick(t20, var = "bottomT")
var<- raster::extract(t20,datasp, method = 'bilinear', small = T, na.rm = T, df = T, exact = T, cellnumbers = F);rm(t20)
colnames(var) <- substr(colnames(var), 2,8); 
for (i in 1:nrow(d2020)){
  a <- as.numeric(var[i,which(colnames(var) == d2020$year_month[i])])
  d2020$bot_T[i] <- ifelse(length(a) == 0,
                               0, var[i,which(colnames(var) == d2020$year_month[i])])
}

#Chlorophyll - 2005-2019
chl05 <- 'chl_2005-2010.nc'; chl05 <- brick(chl05, var = "chl", level = 3)
chl10 <- 'chl_2011-2015.nc'; chl10 <- brick(chl10, var = "chl", level = 3)
chl15 <- 'chl_2016-2019.nc'; chl15 <- brick(chl15, var = "chl", level = 3)
chl <- raster::stack(chl05,chl10,chl15); 
var<- raster::extract(chl,datasp, method = 'bilinear', small = T, na.rm = T, df = T, exact = T, cellnumbers = F); rm(chl05, chl, chl10, chl15)
colnames(var) <- substr(colnames(var), 2,8); 
for (i in 1:nrow(env)){
  a <- as.numeric(var[i,which(colnames(var) == env$year_month[i])])
  env$chl[i] <- ifelse(length(a) == 0,
                         0, var[i,which(colnames(var) == env$year_month[i])])
}

#Temperature at sea bottom- 2020
chl20 <- 'chl_2020.nc'; chl20 <- brick(chl20, var = "chl", level = 3)
var<- raster::extract(chl20,datasp, method = 'bilinear', small = T, na.rm = T, df = T, exact = T, cellnumbers = F);rm(chl20)
colnames(var) <- substr(colnames(var), 2,8); 
for (i in 1:nrow(d2020)){
  a <- as.numeric(var[i,which(colnames(var) == d2020$year_month[i])])
  d2020$chl[i] <- ifelse(length(a) == 0,
                           0, var[i,which(colnames(var) == d2020$year_month[i])])
}

env <- rbind(env,d2020) #join again the databases

#Bottom sediment type
path = 'C:/Users/avipo/Desktop/PhD/Project 1/MarenzelleriaSweden/env_data/sed_klas'
sed <- raster("C:/Users/avipo/Desktop/PhD/Project 1/MarenzelleriaSweden/env_data/sed_klas/w001001x.adf")
plot(sed)
proj4string(datasp) <- CRS("+proj=longlat +datum=WGS84")
dataspUTM<-spTransform(datasp,CRS("+proj=utm +zone=34 +datum=WGS84 +units=m +no_defs")) #we change the type of projection in the case that would be needed
var<- raster::extract(sed,dataspUTM, method = 'simple', cellnumbers = F)
env<-data.frame(env,Sediment=var); rm(sed)

#Depth
setwd('C:/Users/avipo/Desktop/PhD/Project 1/MarenzelleriaSweden')
depth<-read.csv("Marenzelleria_1990_2020_CLEAN.txt",header=T,dec=".",sep="\t")
head(depth); depth <- depth %>% select(station, depth); depth <- depth[which(!duplicated(depth$station)),]
env <- merge(env, depth, by="station")

#write.table(env,file="env_data.txt",sep="\t", row.names = TRUE)

ggplot() +
  geom_point(aes(x=envNA$lon,y=envNA$lat), color = "darkgreen") +
  # scale_fill_gradientn(name = "Clusters", colours=colpal, na.value = 'white') +
  borders(fill="gray44",colour="black") +
  coord_quickmap(xlim=c(10, 30),ylim= c(54,67))+
  labs(x = "Lon",y="Lat")+
  theme(panel.background = element_rect(fill=alpha('light blue', 0.4), colour = 'black'))

#Deal with NAs, use NEAREST NEIGHBOUR
#Extract all the rows that have NA values for some environmental variables
envNA <- env[which(is.na(env$bot_oxy)),]; envNA <- envNA[,-c(7:10)]; try <- env[which(!is.na(env$bot_oxy)),]
coordinates <- unique(try[,c("lon", "lat")]);rownames(coordinates) <- NULL; coordinates$nR <- rownames(coordinates)
closest <- nn2(coordinates[,c("lon", "lat")], envNA[,c("lon", "lat")],
               k = 4, searchtype = "standard")# fast nearest neighbour search
closest <- sapply(closest, cbind) %>% as_tibble
closest$group <- rep(1:nrow(envNA), each = 4) # group rows in 4 as we have the 4 nearest points to our NAs points

NA4nearest <- data.frame() # obtain the coordinates of the 4 nearest points in order to 
for (i in 1:nrow(closest)){
value <- coordinates[closest$nn.idx[i],]; value <- cbind(value, closest$group[i]) 
NA4nearest <- rbind(NA4nearest, value); rownames(NA4nearest) <- NULL
}
colnames(NA4nearest)[4] <- "group"
prova <- merge(NA4nearest, try, by = c("lon","lat")); prova$group <- as.factor(prova$group); group <- c(levels(prova$group))

for (i in 1:length(group)){
  
  a <- prova[which(prova$group == group[i]),]; b <- envNA[i,]
  
  
}
# nn.idx refers to the rownumber on coordinates dataframe
# nn.dist refers to the distance to that location
neighbour <- cbind(nn.idx = closest$nn.idx, envNA)
# merge the dataframe with the nearest coordinates to the dataframe of interest matching the rows
neighbour <- merge(coordinates, neighbour, by.x = "nR", by.y = "nn.idx", 
                   all.x = FALSE, all.y = TRUE); colnames(neighbour)[2:3] <- c("lon","lat")

# get the coordinates of the nearest points to the NA dataframe and obtain the environmental data from the whole dataframe
coords <- neighbour[,2:3]
prova <- merge(try, coords, by = c("lon", "lat"), all.x = FALSE, all.y = TRUE)
#prova <- prova %>% select(lon, lat, bot_oxy, bot_sal, bot_T)
#neighbour <- neighbour[,-c(1:3)]; colnames(neighbour)[3:4] <- c("lon","lat")

newENV <-  merge(prova, neighbour, by = c("lon", "lat"), all.x = FALSE, all.y = TRUE)
newENV <- newENV[, -c(1:2,6)] # delete the coordinates for the nearest stations and the number of rows names





#4.2: Extract the environmental variables for the water column #######

#Chlorophyll - 2005-2019
nc = nc_open("chl_2005-2010.nc")
# Pull out individual layers and create vectors.
chl <-  ncvar_get(nc,"chl")[,,1:15,] #Depth, getting the levels that account for 50 metres
xlon <- ncvar_get(nc,"longitude"); xlat <- ncvar_get(nc,"latitude"); time <- ncvar_get(nc,"time")
df <- rbind(data.frame(lon = as.vector(xlon),lat = as.vector(xlat),
                       time = as.vector(time), var = as.vector(var)))
df <- df %>% group_by(lon, lat, time) %>%  summarise_at(.vars = "var", mean, na.rm = T)
df <- dcast(df,lon+lat~time, value.var='var')

# Select unique locations with environmental data
coordinates <- unique(df[,c("lon", "lat")]); coordinates$nR <- rownames(coordinates) # Add row number as variable to match the datasets further on (see below)
closest <- nn2(coordinates[,c("lon", "lat")], env[,c("lon", "lat")],
               k = 1, searchtype = "standard")# fast nearest neighbour search
closest <- sapply(closest, cbind) %>% as_tibble
# nn.idx refers to the rownumber on coordinates dataframe
# nn.dist refers to the distance to that location
neighbour <- cbind(nn.idx = closest$nn.idx, env)
neighbour <- merge(coordinates, neighbour, by.x = "nR", by.y = "nn.idx", all.x = FALSE, all.y = TRUE)
coords <- 


#Chlorophyll- 2020
t20 <- 'chl_2020.nc'; t20 <- brick(t20, var = "bottomT")
var<- raster::extract(t20,datasp, method = 'simple', small = T, na.rm = T, df = T, exact = T, cellnumbers = F);rm(t20)
colnames(var) <- substr(colnames(var), 2,8); 
for (i in 1:nrow(d2020)){
  a <- as.numeric(var[i,which(colnames(var) == d2020$year_month[i])])
  d2020$bot_T[[i]] <- ifelse(length(a) == 0,
                             0, var[i,which(colnames(var) == d2020$year_month[i])])
}

#Dissolved oxygen
setwd("./env_data")
oxy05 <- 'o2_o2bot_2005-2010.nc'; oxy05 <- brick(oxy05, var = "o2", lvar = 4)
oxy10 <- 'o2_o2bot_2011-2015.nc'; oxy10 <- brick(oxy10, var = "o2",lvar = 4)
oxy15 <- 'o2_o2bot_2016-2019.nc'; oxy15 <- brick(oxy15, var = "o2",lvar = 4)
oxy <- raster::stack(oxy05,oxy10,oxy15); 
coordinates(datasp)=~lon+lat
var<- raster::extract(oxy,datasp, method = 'simple', small = T, na.rm = T, df = T, exact = T, cellnumbers = F);rm(oxy, oxy05,
                                                                                                                  oxy10,oxy15)
colnames(var) <- substr(colnames(var), 2,8); 
for (i in 1:nrow(env)){
  a <- as.numeric(var[i,which(colnames(var) == env$year_month[i])])
  env$bot_oxy[[i]] <- ifelse(length(a) == 0,
                             0, var[i,which(colnames(var) == env$year_month[i])])
}


#4.3: Set the model to see how environment influences Di ##########

###########################################

int_distinct_opt = function(trait){ #Here trait is a species x trait information matrix, previously cleaned 
  #and with the traits of interest selected
  
  sensitivity.gow<-trait; sensitivity.gow[1:ncol(sensitivity.gow)] <- NULL
  sum.standard.gow<-matrix(0,nrow(as.data.frame(trait)),nrow(as.data.frame(trait)))
  rownames(sum.standard.gow)<-rownames(trait)
  colnames(sum.standard.gow)<-rownames(trait)
  sum.standard.gow<-as.data.frame(sum.standard.gow)
  
  # Defining the total number of column : N; it MUST be the number of traits that we have
  N<-ncol(trait)
  listcol<-seq(1,N,1)
  # Number of columns corresponding to traits
  listk<-seq(1,N,1)
  
  if (ncol(trait) >= 4){#start condition for more than 4 trait categories
      
      k=4; sensitivity.gow <- foreach (k = 4:N, .combine = cbind, .packages = c("funrar","dplyr")) %do% {
      
      #Here are determined all possible 4 trait combinations for obtaining the average distance
      #between species, here are expressed only the position that are occupying the traits in the data frame
      
      combination<-data.frame(t(combn(listk,k))) 
      
      if (is.null(combination$X5==T)){
        combination$X5<-NA
      }
      if (is.null(combination$X6==T)){
        combination$X6<-NA
      }
      if (is.null(combination$X7==T)){
        combination$X7<-NA
      }
      if (is.null(combination$X8==T)){
        combination$X8<-NA
      }
      if (is.null(combination$X9==T)){
        combination$X9<-NA
      }
      if (is.null(combination$X10==T)){
        combination$X10<-NA
      }
      if (is.null(combination$X11==T)){
        combination$X11<-NA
      }
      
      result <- foreach (i = 1:nrow(combination), .combine = cbind, .packages = c("funrar","dplyr"), .errorhandling = "remove") %dopar% {
        
        data <- trait
        codex<-as.vector(t(combination[i,]))
        
        #Here, each row (in other words, each combination of 4 traits) is selected for a 
        #total calculation of distances based on all combinations possibles of the traits
        data<-dplyr::select(data,all_of(codex))
        
        #Calculation of distinctiveness using Gower's distance
        gow<-compute_dist_matrix(data, metric = "gower") # a distance matrix is calculated for each combination of traits
        m.gow<-as.matrix(gow)
        d.gow<-as.data.frame(m.gow)
        standard.gow<-standr(d.gow) # distances are standardized to fit between 1 and 0
        standard.Di.gow<-colSums(standard.gow)/(nrow(standard.gow)-1) # Mean of the computed distances for each species, to obtain the
        # average distance to the other species based on all trait combinations
        standard.Di.gow<-as.data.frame(standard.Di.gow)
        standard.Di.gow
      }
      
    }
  }#end condition for more than 4 trait categories
  
  else{#start condition for less than 3 categories
    
    #Calculation of distinctiveness using Gower's distance
    gow<-compute_dist_matrix(trait, metric = "gower") # a distance matrix is calculated for each combination of traits
    m.gow<-as.matrix(gow)
    d.gow<-as.data.frame(m.gow)
    standard.gow<-standr(d.gow) # distances are standardized to fit between 1 and 0
    standard.Di.gow<-colSums(standard.gow)/(nrow(standard.gow)-1) # Mean of the computed distances for each species, to obtain the
    # average distance to the other species based on all trait combinations
    standard.Di.gow<-as.data.frame(standard.Di.gow)
    sensitivity.gow<-cbind(sensitivity.gow,standard.Di.gow)
    } #end condition for less than 3 categories
  
  forIntDi<-t(sensitivity.gow); forIntDi<-na.omit(forIntDi); forIntDi<-as.matrix(t(forIntDi))
  outputmatrix<-matrix(0,nrow(forIntDi),2)
  
  for (j in 1:nrow(forIntDi)){# start loop output matrix
    
    subset_sp<-forIntDi[j,1:ncol(forIntDi)]
    subset_sp<-as.numeric(as.character(subset_sp))
    outputmatrix[j,2]<-mean(subset_sp) #compute the mean for all the distances obtained via the combinations
    
  }#end loop output matrix
  
  Int_Di<-data.frame(outputmatrix); Int_Di[,1]<-rownames(forIntDi); colnames(Int_Di)<-c("taxon","int_Di"); 
  #Int_Di<-Int_Di[order(Int_Di$int_Di,decreasing=T),]
  
  return(list(Int_Di, standard.gow))
}

# Function to get the distance between two species
status <- traitraw %>% select(status,BioticID)

try <- trial[which(trial$cell == cells[j]),]; 
try$cell <- as.numeric(as.character(try$cell)) 
try$station <- as.numeric(as.character(try$station))
cols <- try[,which(colSums(try)> 0)]

colnames(dist_matrix) <- gsub("\\.", " ", colnames(dist_matrix))


mean_diss = function(dist_matrix, data, status){
  
  if (all(is.na(match(colnames(cols), rownames(status[which(status$status == "Non-indigenous"),]))))){NIS <- c()}
  else {NIS <- c(colnames(data)[which(!is.na(match(colnames(data), rownames(status[which(status$status == "Non-indigenous"),]))))])}
  
  NAT <- c(colnames(cols)[6:ncol(cols)]); NAT <- gsub("\\.", " ", NAT)
  NAT <- setdiff(NAT,NIS)
  #mean distance of NIS to other native species per each location/cell
  NIS_Di <- foreach(i = 1:length(NIS), .combine = cbind, .errorhandling = "pass") %dopar% {
            dis <- c()
            for(j in 1:length(NAT)) {val <- dist_matrix[NIS[i],NAT[j]]; dis <- c(dis,val)}
            dis <- mean(dis); dist <- as.data.frame(dis); colnames(dist) <- NIS[i]; dist}
  return(NIS_Di)
}

#dist_matrix accounts for a pairwise matrix of integrated distinctiveness (functional distances)
#data is the species occurrences/abundance data for a certain location/cell where we have different sampling locations
#status is a list for all the species we have in TOTAL where their status (NIS or native) is defined

Dixcell <- data.frame(matrix(ncol = 2))
#colnames(Dixcell) <- c("size_reduction", "cell",rownames(traitraw[which(traitraw$status == "Non-indigenous"),]))
colnames(Dixcell) <- c("size_reduction", "cell")

dbFD(dummy$trait, dummy$abun)
ex6 <- dbFD(tussock$trait, tussock$abun, corr = "lingoes")

dat <- data[,5:ncol(data)]
dat = dat[!apply(dat, 1, function(x) all(x == 0)), ]
species_funrar <- dbFD(traitraw, dat)

results <- assign_cell(data,num = 0.61); data <- results[[1]]

# Nearest neighbor
# Loop over files
nc = nc_open("chl_2005-2010.nc")
# Pull out individual layers and create vectors.
temp <-  ncvar_get(nc,'')
chl <-  ncvar_get(nc,"chl")[,,1:15,] #Depth, getting the levels that account for 50 metres
xlon <- ncvar_get(nc,"longitude"); xlat <- ncvar_get(nc,"latitude"); time <- ncvar_get(nc,"time")
df <- rbind(data.frame(lon = as.vector(xlon),lat = as.vector(xlat),
                       time = as.vector(time), chl = as.vector(surf.temp1)))
df <- df %>% group_by(lon, lat, time) %>%  summarise_at(.vars = "chl", mean, na.rm = T)
df <- dcast(df,lon+lat~time, value.var='chl')

# Select unique locations with environmental data
coordinates <- unique(df[,c("lon", "lat")]); coordinates$nR <- rownames(coordinates) # Add row number as variable to match the datasets further on (see below)
closest <- nn2(coordinates[,c("lon", "lat")], env[,c("lon", "lat")],
               k = 1, searchtype = "standard")# fast nearest neighbour search
closest <- sapply(closest, cbind) %>% as_tibble
# nn.idx refers to the rownumber on coordinates dataframe
# nn.dist refers to the distance to that location
neighbour <- cbind(nn.idx = closest$nn.idx, env)
neighbour <- merge(coordinates, neighbour, by.x = "nR", by.y = "nn.idx", all.x = FALSE, all.y = TRUE)





















