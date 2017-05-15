##################################################
# list of utility functions
##################################################

read_element <- function(cache_name = "data/RData/data_BRCA_all_lazyReady.R_CACHE", element_name = "counts") {
	require(SOAR)
	oldLC <- Sys.getenv("R_LOCAL_CACHE", unset = ".R_Cache")
	Sys.setenv(R_LOCAL_CACHE=cache_name)
	Objects()
	tmp <- get(element_name)
	Sys.setenv(R_LOCAL_CACHE=oldLC)
	return(tmp)
}

read_methylation_matrix_gene <- function(cache_name = "data/RData/data_BRCA_all_lazyReady.R_CACHE", genename = "A1BG") {
	require(SOAR)
	require(dplyr)
	oldLC <- Sys.getenv("R_LOCAL_CACHE", unset = ".R_Cache")
	Sys.setenv(R_LOCAL_CACHE=cache_name)
	Objects()
	tmp <- as.data.frame(get(genename)) %>% select(starts_with("pr_"), starts_with("gb_")) %>% as.matrix()
	Sys.setenv(R_LOCAL_CACHE=oldLC)
	return(tmp)
}

read_grouping <- function(cache_name = "data/RData/data_BRCA_all_lazyReady.R_CACHE", group1 = "ANs", group2 = "Ts") {
	require(SOAR)
	oldLC <- Sys.getenv("R_LOCAL_CACHE", unset = ".R_Cache")
	Sys.setenv(R_LOCAL_CACHE=cache_name)
	Objects()
	tmp <- list(get(group1), get(group2))
	Sys.setenv(R_LOCAL_CACHE=oldLC)
	return(tmp)
}

read_genedata <- function(cache_name = "data/RData/data_BRCA_all_lazyReady.R_CACHE", genename = "A1BG") {
	require(SOAR)
	oldLC <- Sys.getenv("R_LOCAL_CACHE", unset = ".R_Cache")
	Sys.setenv(R_LOCAL_CACHE=cache_name)
	Objects()
	tmp <- get(genename)
	Sys.setenv(R_LOCAL_CACHE=oldLC)
	return(tmp)
}

fishersMethod <- function(x, df=2*length(x[!is.na(x)]), log.p = FALSE, logs = FALSE) {
	x <- x[!is.na(x)]
	if (!logs) pchisq(-2*sum(log(x)), df, lower.tail = FALSE, log.p = log.p) else pchisq(-2*sum(x), df, lower.tail = FALSE, log.p = log.p)
}

integrand_e <- function(x, k) {dpois(k,x)}
integrand_m <- function(x ,mean) {dnorm(x=mean,mean=x,sd=0.14)}

tensor_product <- function(matrix1, matrix2, smooth_h=0, normalize=c("row","column","no"), kernel=c("gauss","cauchy","minvar")) {
	require(smoothie)
	if (is.matrix(matrix1) && is.matrix(matrix2)) result <- matrix(ncol=ncol(matrix1),nrow=ncol(matrix2),data=rep(0,ncol(matrix1)*ncol(matrix2))) else result <- matrix(ncol=length(matrix1),nrow=length(matrix1),data=rep(0,length(matrix1)*length(matrix2)))
	
	if (is.matrix(matrix1) && is.matrix(matrix2)) for (i in 1:nrow(matrix1)) {
		result <- result + matrix(ncol=ncol(matrix1),nrow=ncol(matrix2),byrow=TRUE,data=apply(expand.grid(matrix1[i,],matrix2[i,]), 1, prod))
	} else result <- result + matrix(nrow=length(matrix1),ncol=length(matrix2),byrow=TRUE,data=apply(expand.grid(matrix1,matrix2), 1, prod))
	
	if (is.matrix(matrix1) && is.matrix(matrix2)) result <- result/nrow(matrix1)
	if (!is.null(kernel) && smooth_h > 0) result <- kernel2dsmooth(result,kernel.type = kernel[1],sigma=smooth_h,nx=ncol(matrix2),ny=ncol(matrix1))
	if (normalize[1] == "row") for (i in 1:nrow(result)) result[i,] <- result[i,]/sum(result[i,]) else if (normalize[1] == "column") for (i in 1:nrow(result)) result[,i] <- result[,i]/sum(result[,i])
	return(result)
}

geo_mean <- function(data) {
	log_data <- log(data)
	gm <- exp(mean(log_data[is.finite(log_data)]))
	return(gm)
}
