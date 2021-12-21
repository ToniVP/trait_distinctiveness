# Integrated distinctiveness function


int_distinct = function(trait){ #Here trait is a species x trait information matrix, previously cleaned 
  #and with the traits of interest selected
  
  sensitivity.gow<-trait; sensitivity.gow[1:ncol(sensitivity.gow)] <- NULL
  sum.standard.gow<-matrix(0,nrow(as.data.frame(trait)),nrow(as.data.frame(trait)))
  rownames(sum.standard.gow)<-rownames(trait)
  colnames(sum.standard.gow)<-rownames(trait)
  sum.standard.gow<-as.data.frame(sum.standard.gow)
  dist.matr <- list()
  
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
        dist.matr <- c(dist.matr,list(standard.gow))
      }
      
      print (k)
    }
    
    dist_matrix <- as.data.frame(mean_matrix(dist.matr))
    
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
  
  return(list(Int_Di,dist_matrix))
}