###############################################################
######## Summarizing tables that are stored in a list #########
###############################################################
### Take apart tables in list, put them into vectors, get means
### and standard errors from the elements across columns (which
### represent table replicates) and then put back into a table
### of the original formal
### INPUT: list of tables, assuming that ALL ARE OF THE SAME
### DIMENSIONS
### OUTPUT: listed above in description of function


###For TreePar

### extracting table information (function)
tables.summary <- function(final){
	tmp <- final[[1]][,2:ncol(final[[1]])] #just numbers here
	tmp2 <- final[[1]] # whole thing
	dt<-nrow(tmp) * ncol(tmp) #number of numeric cells to summarize
	dt2 <- nrow(tmp2) * ncol(tmp2) #number of total cells
	# strip irrelevant information in first two columns of input tables and make numeric and vector (at same time with 'as.numeric')
	output.tmp <- sapply(final,function(x){as.numeric(x[,2:ncol(tmp2)])})
	means <- apply(output.tmp,MARGIN=1,mean,na.rm = TRUE)
	stders <- apply(output.tmp,1,function(x){sqrt(var(x,na.rm=TRUE)/length(x))})
	mf <- matrix(means,nrow=nrow(tmp),byrow=FALSE)
	sf <- matrix(stders,nrow=nrow(tmp),byrow=FALSE)
	#mf <- data.frame(tmp2[,1],tmp2[,2],mf)
	mf <- data.frame(tmp2[,1],mf)
	#sf <- data.frame(tmp2[,1],tmp2[,2],sf)	
	sf <- data.frame(tmp2[,1],sf)	
	rownames(mf)<-rownames(sf)<-rownames(tmp2)
	colnames(mf)<-colnames(sf)<-colnames(tmp2)
	return(list(means=mf,std.err=sf))
}