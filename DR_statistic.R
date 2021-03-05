# A script to calculate DR statistics from Jetz et al. (2012)
rm(list=ls())
library("ape")
library("phytools")

#defined function from Jetz et al. (2012), aslo see Harvey et al. (2016)
DR_statistic <- function(tree, return.mean = FALSE){
  rootnode <- length(tree$tip.label) + 1
  sprates <- numeric(length(tree$tip.label))
  for (i in 1:length(sprates)){
    node <- i
    index <- 1
    qx <- 0
    while (node != rootnode){
      el <- tree$edge.length[tree$edge[,2] == node]
      node <- tree$edge[,1][tree$edge[,2] == node]			
      qx <- qx + el* (1 / 2^(index-1))			
      index <- index + 1
    }
    sprates[i] <- 1/qx
  }
  if (return.mean){
    return(mean(sprates))		
  }else{
    names(sprates) <- tree$tip.label
    return(sprates)
  }
}

#read tree
  tree <- read.tree("data/Saussurea.tre")
  if(is.ultrametric(tree)){
    print("TRUE")
  }else{
    tree <- force.ultrametric(tree, method="extend")
  }
  DR <- DR_statistic(tree)
  write.csv(paste0(names(DR), ",", DR, sep=""), "rosids_5g_whole_tree_tip_DR.csv", row.names = FALSE, quote=FALSE)

  
  #run for orders
  
  order.list <- c("Brassicales", "Celastrales", "Crossosomatales", "Cucurbitales", "Fabales", "Fagales", 
                  "Geraniales", "Huerteales", "Malpighiales", "Malvales", "Myrtales", "Oxalidales", "Picramniales", "Rosales", 
                  "Sapindales", "Zygophyllales", "Vitales")
  
  for (i in 1:length(order.list)){
    
    tryCatch({
      
      Order <- order.list[i]
      
      #read tree
      tree <- read.tree(paste("../../data/Rosid_Ultrametric_Trees/", Order, "_5g.tre", sep=""))
      if(is.ultrametric(tree)){
        cat(paste0(Order, "\t", "TRUE", sep=""), sep = "\n")
      }else{
        tree <- force.ultrametric(tree, method="extend") #just make sure the tree is ultrametric
      }
      
      print(Order)
      DR <- DR_statistic(tree)
      write.csv(paste0(names(DR), ",", DR, sep=""), paste0("rosids_5g", Order, "_tree_tip_DR.csv", row.names = FALSE, quote=FALSE))
    })
  }
  