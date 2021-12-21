# Obtain the mean of each element from a list of matrices

# x accounts for the list of matrices
mean_matrix<- function(x){
  y <- array(unlist(x) , c(dim(x[[1]]),length(x))) #convert the list to a 3D array
  y <- apply( arr , 1:2 , mean ) #apply the mean function to the same cell on all the arrays, to obtain a mean matrix
  colnames(y)<-colnames(x[[1]])
  rownames(y)<-rownames(x[[1]])
  return(y)
}

#Function to standardize
standr = function(x){(x-min(x))/(max(x)-min(x))} 

#Get the species that are present both in occurrences data and also in trait database

# x as the abundance matrix fo the year/station that we want, as there are the selected species
# y as the raw trait matrix, to select species present in the abundance matrix
# num as the column number where the first species is

get_common_species = function(x, y, num) {
  if(class(x) == "data.frame"){colnames(x)<- gsub("\\.", " ", colnames(x)); mat_columns = colnames(x)[num:ncol(x)]}
  else {mat_columns = x} #the abundance matrix or a vector containing names
  
  trait_species = rownames(y) #the trait dataset
  
  # Extract matching species
  common_species = match(trait_species, mat_columns)
  rows = which(!is.na(common_species))
  
  return(y[rows,]) 
} 

#Function to make hulls about the similarity of species
convexe_hull<-function(x){
  s <- PCOA_IDi %>%
    split(x)  # Transform data.frame in lists of data.frames depending on column 'z'
  ch = s %>%
    # Compute which points are on the convex hull of each data.frame 
    lapply(., function(el) chull(el$PCOA1, el$PCOA2))  # 'ch' now contains the row numbers of points on convex hull per sub-data.frame
  # Get points for each sub-data.frame using names index
  ch = lapply(names(ch), function(PCOA1) s[[PCOA1]][ch[[PCOA1]],]) %>%   
    do.call(rbind, .)  # Join all convex hull points in a single data.frame
  return(ch)
}


#Function to unify categories of the traits

#edit the function based on what should be cleaned from the selected trait dataset
# to unify the criteria and be able to add different information manually

format_trait = function(x) { 
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
}

#Functions to adjust a Michaelis-menten equations

MM <- function(k, z, x){
  y <- k*x/(z+x)
  return(y)
}

MMprime <- function(k,z,x){
  y <- k*z/((z+x)^2)
  return(y)
}

#Function to calculate functional distances between different species groups using the distances matrix

#dist_matrix accounts for a pairwise matrix of integrated distinctiveness (functional distances)
#data is the species occurrences/abundance data for a certain location/cell where we have different sampling locations
#status is a list for all the species we have in TOTAL where their status (NIS or native) is defined
#IntDi is the result of the integrated distinctiveness for all the species obtained previously

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
} 
