#!/usr/bin/env Rscript
library(optparse)
# Command Line Arguments
option_list <- list(
  make_option(c("-d","--data"), type='character', help = 'cached data to use', default = "data/RData/data_LUSC_all_lazyReady.R_CACHE"),
  make_option(c("-r","--results"), type='character', help = 'model results to use', default = "data/latestResults/top_100_each_seed.RData"),
  make_option(c("-f","--folds"), type='integer', help = 'number of folds to split the dataset', default = 4),
  make_option(c("-n","--name"), type='character', help = 'name of the analysis performed', default = "MLs_EBADIMEX_50randomSamples_BRCA_progressing"),
  make_option("--n_top", type='integer', help = 'number of top ranks/genes to train for and combine', default = 100),
  make_option(c("-p","--parallel"), type='integer', help = 'if desired, provide a number of cores to parallelize the execution of the analysis', default = 1)
)
#opt <- parse_args(OptionParser(option_list = option_list))
opt <- list(data = "data/RData/data_BRCA_progressing_lazyReady.R_CACHE", results="data/latestResults/top_100_each_seed.RData",  name = "PINCAGE_EBADIMEX_50randomSamples_BRCA_progressing", folds = 50, n_top = 100, parallel=8)

library(caret)
library(pROC)
library(dplyr)
library(parallel)
source("R/utilities_new.R")
load(opt$results)
load("/home/michal/iCancerGenomics/irksome-fibula/essentials_EBADIMEXgenes.RData")
load("/home/michal/iCancerGenomics/irksome-fibula/UQs_BRCA_progressing.RData")
samplenames <- unlist(read_grouping(cache_name = opt$data, group1="progressed", group2="nonProgressed"))
indiv_predictions_PINCAGE_cv <- matrix(ncol=length(samplenames), nrow=opt$n_top)
colnames(indiv_predictions_PINCAGE_cv) <- samplenames

# configure here
res_pr <- 25 # number of bins for promoter methylation
res_gb <- 25 # number of bins for gene body methylation
res_expr <- 25 # number of bins for gene expression

smooth_1d <- 1 # should expression node parameterization be smoothed?
smooth_2d <- 1 # should methylation 2D factors parameterizations be smoothed?
if (smooth_1d > 0) library(aws)
if (smooth_2d > 0) library(smoothie)


# the following objects are expected by this script: workingList_BRCA, mmatrix_pc, counts_plusOne, G1, G2, factors_ls, TSS1500Ind, TSS200Ind, UTR5Ind, EXON1Ind, GENEBODYInd, UTR3Ind

args <- commandArgs(trailingOnly = TRUE)
beg <- as.numeric(args[1]) # first ID to process
end <- as.numeric(args[2]) # last consequitve ID to process

integrand_e <- function(x, k) {dpois(k, x)}
integrand_m <- function(x, mean) {dnorm(x=mean, mean=x, sd=0.14)}

tensor_product <- function(matrix1, matrix2, smooth_h=0, normalize=c("row", "column", "no"),kernel=c("gauss", "cauchy", "minvar")) {
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

sum_asLogs <- function (x) {
	if (x[1] >= x[2]) x[1] + log(1+ exp(x[2]-x[1])) else x[2] + log(1+ exp(x[1]-x[2]))
}


AUCs <- NULL
CV_folds <- NULL
for (cv_run in 1:50) {
	
	current_genes <- top_100_each_seed[[cv_run]]
	system(command=paste('mkdir ', cv_run, sep=""))
	
	# for (gene in 1:length(current_genes)) {
	if (opt$parallel > 1) {
		out <- mclapply(1:length(current_genes), mc.cores = opt$parallel, mc.preschedule = FALSE, FUN = function(gene) {
			cat(paste("doing run ", cv_run,", ", gene, "\n",sep=""))
			#ptm <- proc.time()[3]
			set.seed(cv_run)
			train_set_g1 <- sample(1:38, 10)
			train_set_g2 <- 38+sample(1:43, 15)
			test_set_g1 <- (1:38)[-train_set_g1]
			test_set_g2 <- 38+(1:43)[-(train_set_g2-38)]
			
			genedata <- as.data.frame(read_genedata(cache = opt$data, genename = current_genes[gene])) %>% 
				dplyr::select(starts_with("pr_"), starts_with("gb_"), matches("read_count"))
			genedata_training <- genedata  %>%
				slice(c(train_set_g1, train_set_g2))
			genedata_testing <- genedata %>%
				slice(c(test_set_g1, test_set_g2))
				
			
			G1 <- rownames(genedata)[1:38]
			G2 <- rownames(genedata)[39:81]
			train_set_g1 <- rownames(genedata)[train_set_g1]
			train_set_g2 <- rownames(genedata)[train_set_g2]
			test_set_g1 <- rownames(genedata)[test_set_g1]
			test_set_g2 <- rownames(genedata)[test_set_g2]
			rownames(genedata_training) <- c(train_set_g1, train_set_g2)
			rownames(genedata_testing) <- c(test_set_g1, test_set_g2)
			
			system(command=paste('mkdir ', cv_run, '/', gene, sep=""))
			system(command=paste('mkdir ./', cv_run, '/', gene, '/G1_model',sep=""))
			system(command=paste('mkdir ./', cv_run, '/', gene, '/G2_model',sep=""))
			###########################################################
			############### identify constitutive CpGs ################
			IDs_promoter <- colnames(genedata)[grep("pr_", colnames(genedata))]
			promoterVars <- promoterVars_template[1:length(IDs_promoter)]
			promoter_CpGs <- template_promoter_CpGs[1:length(IDs_promoter)]
			IDs_body <- colnames(genedata)[grep("gb_", colnames(genedata))]
			geneBodyVars <- geneBodyVars_template[1:length(IDs_body)]
			geneBody_CpGs <- template_body_CpGs[1:length(IDs_body)]
			ncol = length(IDs_promoter) + length(IDs_body) + 1
			
			############################################################
			######### calculate epsilons and smoothing params ##########
			epsilon_pr_G2 <- 1/(length(train_set_g2)*length(IDs_promoter))/res_pr
			epsilon_gb_G2 <- 1/(length(train_set_g2)*length(IDs_body))/res_gb
			epsilon_e_G2 <- 1/length(train_set_g2)/res_expr
			epsilon_pr_G1 <- 1/(length(train_set_g1)*length(IDs_promoter))/res_pr
			epsilon_gb_G1 <- 1/(length(train_set_g1)*length(IDs_body))/res_gb
			epsilon_e_G1 <- 1/length(train_set_g1)/res_expr
			
			smooth_e <- 10/(mean(length(train_set_g1),length(train_set_g2))/(res_expr*2))
			smooth_pr <- trunc(10/(mean(length(train_set_g1),length(train_set_g2))*length(IDs_promoter)/(res_pr*res_expr)))
			smooth_gb <- trunc(10/(mean(length(train_set_g1),length(train_set_g2))*length(IDs_body)/(res_gb*res_expr)))
			###########################################################
			############### define the binning scheme  ################
			all_labels_pr <- as.character(seq(1,res_pr,1))
			all_labels_gb <- as.character(seq(1,res_gb,1))
			all_labels_expr <- as.character(seq(1,res_expr,1))
			pseudo_counts_pr <- matrix(ncol=res_pr,nrow=res_expr,data=rep(1,res_expr*res_pr))
			pseudo_counts_gb <- matrix(ncol=res_gb,nrow=res_expr,data=rep(1,res_expr*res_gb))
			
			# gene body
			x <- t(genedata[c(G2,G1), IDs_body])
			
			density <- density(x,bw=0.14,from=-7,to=7,n=2801,na.rm=TRUE)
			density$y <- density$y/sum(density$y)
			density$y <- cumsum(density$y)
			breaks <- NULL
			noBreaks <- res_gb-1
			for (j in 1:noBreaks) { breaks <- c (breaks, density$x[which(density$y >= j*(1/(1+noBreaks)))][1])}
			breaksBODY <- sort(c(-7.01,breaks,7.01))
			
			# promoter
			x <- t(genedata[c(G2,G1), IDs_promoter])
			
			density <- density(x,bw=0.14,from=-7,to=7,n=2801,na.rm=TRUE)
			density$y <- density$y/sum(density$y)
			density$y <- cumsum(density$y)
			breaks <- NULL
			noBreaks <- res_pr-1
			for (j in 1:noBreaks) { breaks <- c (breaks, density$x[which(density$y >= j*(1/(1+noBreaks)))][1])}
			breaksPROMOTER <- sort(c(-7.01,breaks,7.01))
			
			# expression
			temp_UQs <- matrix(ncol=2)
			colnames(temp_UQs) <- c("UQ","density")
			for (j in 1:nrow(genedata)) {
				lambda <- as.numeric(genedata$read_count[j]) + 1
				X <- seq(round(max(lambda - (4 * lambda * lambda^(-1/2)), 1)), round(lambda + (4 * lambda * lambda^(-1/2))))
				current <- UQs[c(G2, G1)[j]]
				temp_UQs <- rbind(temp_UQs, cbind(X / current, dpois(X, lambda = lambda) * current))
			}
			temp_UQs <- as.data.frame(temp_UQs[-1,],)
			temp_UQs <- temp_UQs[order(temp_UQs$UQ),]
			temp_UQs[,3] <- cumsum(temp_UQs[,2])
			temp_UQs[,3] <- temp_UQs[,3]/max(temp_UQs[,3])
			breaks <- NULL
			noBreaks <- res_expr-1
			for (j in 1:noBreaks) { breaks <- c (breaks, temp_UQs[which(temp_UQs[,3] >= j*(1/(1+noBreaks))),1][1])}
			breaksEXPRESSION <- c(0,breaks,10^6)
			########################################################
			# dynamic generation of model specification files here #
			########################################################
			####################### state map ######################
			stateMaps <- file(paste("./",cv_run, '/', gene, "/stateMaps.txt", sep=""),"w")
			exprMap <- paste("NAME:\texprMap\nSYMBOLS:\t", paste(seq(1, res_expr, 1), collapse=" "), "\nMETA_SYMBOLS:\t.=", paste(seq(1, res_expr, 1), collapse=" "), "; *=", paste(seq(1, res_expr, 1), collapse=" "), ";\n\n", collapse="", sep="")
			prMap <- paste("NAME:\tprMap\nSYMBOLS:\t", paste(seq(1, res_pr, 1), collapse=" "), "\nMETA_SYMBOLS:\t.=" ,paste(seq(1, res_pr, 1), collapse=" "), "; *=", paste(seq(1, res_pr, 1), collapse=" "), ";\n\n", collapse="", sep="")
			gbMap <- paste("NAME:\tgbMap\nSYMBOLS:\t", paste(seq(1, res_gb, 1),collapse=" "), "\nMETA_SYMBOLS:\t.=", paste(seq(1, res_gb, 1), collapse=" "), "; *=", paste(seq(1, res_gb, 1), collapse=" "), ";\n", collapse="", sep="")
			cat(exprMap, prMap, gbMap, file=stateMaps)
			close(stateMaps)
			########################################################
			####################### variables ######################
			variables <- file(paste("./", cv_run, '/', gene, "/variables.txt",sep=""),"w")
			cat("STATE_MAP_NAME:\texprMap\nVAR_NAMES:\tEXPR\n\n",sep="",file=variables)
			cat("STATE_MAP_NAME:\tprMap\nVAR_NAMES:\tM.P\n\n",sep="",file=variables)
			cat("STATE_MAP_NAME:\tgbMap\nVAR_NAMES:\tM.GB\n",sep="",file=variables)
			close(variables)
			########################################################
			##################### factor graph #####################
			factorGraph <- file(paste("./", cv_run, '/', gene,"/factorGraph.txt", sep=""), "w")
			cat("NAME:\tEXPR.likelihood\nNB1:\tEXPR\nPOT:\tpot_EXPR.likelihood\n", file=factorGraph)
			cat("\nNAME:\tEXPR.prior\nNB1:\tEXPR\nPOT:\tpot_EXPR.prior\n", file=factorGraph)
			cat(paste("\nNAME:\tEXPR.M.GB\nNB1:\tEXPR\nNB2:\tM.GB\nPOT:\tpot_EXPR.M.GB\n", sep="", collapse=""), file=factorGraph)
			cat(paste("\nNAME:\tEXPR.M.P\nNB1:\tEXPR\nNB2:\tM.P\nPOT:\tpot_EXPR.M.P\n", sep="", collapse=""), file=factorGraph)
			cat(paste("\nNAME:\t", promoter_CpGs, "\nNB1:\tM.P\nPOT:\tpot_", promoter_CpGs,"\n",sep="", collapse=""), file=factorGraph)
			cat(paste("\nNAME:\t", geneBody_CpGs, "\nNB1:\tM.GB\nPOT:\tpot_", geneBody_CpGs,"\n",sep="", collapse=""), file=factorGraph)
			close(factorGraph)
			
			system(command=paste('cp ./',cv_run, '/', gene,'/*.txt ./',cv_run, '/', gene,'/G1_model/',sep=""))
			system(command=paste('cp ./',cv_run, '/', gene,'/*.txt ./',cv_run, '/', gene,'/G2_model/',sep=""))
			
			
			# initialize parameter precomputation matrices
			promoter_G2 <- matrix(ncol=res_pr,nrow=length(G2))
			promoter_G1 <- matrix(ncol=res_pr,nrow=length(G1))
			body_G2 <- matrix(ncol=res_gb,nrow=length(G2))
			body_G1 <- matrix(ncol=res_gb,nrow=length(G1))
			expr_G2 <- matrix(ncol=res_expr,nrow=length(G2))
			expr_G1 <- matrix(ncol=res_expr,nrow=length(G1))
			rownames(promoter_G2) <- rownames(body_G2) <- rownames(expr_G2) <- G2
			rownames(promoter_G1) <- rownames(body_G1) <- rownames(expr_G1) <- G1
			###########################################################################
			################### G1 model using training data ##########################
			###########################################################################
			
			# generate FacData
			tempS_G1 <- matrix(ncol=ncol, nrow=length(G1))
			rownames(tempS_G1) <- G1
			for (current_sample in 1:length(G1)) {
				# expression
				read_count <- as.numeric(trunc(genedata[G1[current_sample], "read_count"])) # adding a pseudo-count to avoide problems with 0-count data
				lambdas <- breaksEXPRESSION * UQs[G1[current_sample]]
				frequencies_expr <- rep(0,length(breaksEXPRESSION)-1)
				for (freq in 1:res_expr) {
					frequencies_expr[freq] <- integrate(integrand_e, lower = lambdas[freq], upper = lambdas[freq+1], read_count,stop.on.error=FALSE)[1]
				}
				frequencies_expr <- unlist(frequencies_expr)
				if (length(which(frequencies_expr==0))==res_expr) frequencies_expr[length(frequencies_expr)] <- 1
				frequencies_expr <- frequencies_expr + epsilon_e_G1
				frequencies_expr <- frequencies_expr/sum(frequencies_expr)
				
				# gene body
				cpg_list_gb <- NULL
				for (cpg in 1:length(IDs_body)) {
					miu <- genedata[G1[current_sample], IDs_body[cpg]]
					if (!is.na(miu)) {
						frequencies_gb <- rep(0,res_gb)
						for (freq in 1:res_gb) {
							frequencies_gb[freq] <- integrate(integrand_m,lower=breaksBODY[freq],upper=breaksBODY[freq+1],mean=miu)$value
						}
						frequencies_gb <- unlist(frequencies_gb) + epsilon_gb_G1
						frequencies_gb <- frequencies_gb/sum(frequencies_gb)
						cpg_list_gb[[cpg]] <- frequencies_gb
					} else cpg_list_gb[[cpg]] <- rep(1/res_gb,res_gb)
				}
				
				# promoter
				cpg_list_pr <- NULL
				for (cpg in 1:length(IDs_promoter)) {
					miu <- genedata[G1[current_sample], IDs_promoter[cpg]]
					if (!is.na(miu)) {
						frequencies_pr <- rep(0,res_pr)
						for (freq in 1:res_pr) {
							frequencies_pr[freq] <- integrate(integrand_m,lower=breaksPROMOTER[freq],upper=breaksPROMOTER[freq+1],mean=miu)$value
						}
						frequencies_pr <- unlist(frequencies_pr) + epsilon_pr_G1
						frequencies_pr <- frequencies_pr/sum(frequencies_pr)
						cpg_list_pr[[cpg]] <- frequencies_pr
					} else cpg_list_pr[[cpg]] <- rep(1/res_pr,res_pr)
				}
				
				tempS_formated <- matrix(ncol=ncol, nrow=1)
				tempS_formated[1,1] <- paste('[1,', res_expr,']((', paste(frequencies_expr,sep="", collapse=","), '))', sep="", collapse="")
				cur_ncol <- 1
				for (element in 1:length(cpg_list_pr)) {
					tempS_formated[1, cur_ncol+1] <- paste('[1,', res_pr,']((', paste(cpg_list_pr[[element]], sep="", collapse=","),'))', sep="", collapse="")
					cur_ncol <- cur_ncol + 1
				}
				for (element in 1:length(cpg_list_gb)) {
					tempS_formated[1, cur_ncol+1] <- paste('[1,', res_gb,']((', paste(cpg_list_gb[[element]], sep="", collapse=","), '))', sep="", collapse="")
					cur_ncol <- cur_ncol + 1
				}
				tempS_G1[current_sample,] <- tempS_formated
				#start precomputing correct initialization of parameters
				promoter_G1[current_sample,] <- apply(matrix(unlist(cpg_list_pr), ncol=res_pr, byrow=TRUE), 2,geo_mean)/sum(apply(matrix(unlist(cpg_list_pr), ncol=res_pr,byrow=TRUE),2, geo_mean))
				body_G1[current_sample,] <- apply(matrix(unlist(cpg_list_gb), ncol=res_gb, byrow=TRUE), 2, geo_mean)/sum(apply(matrix(unlist(cpg_list_gb),ncol=res_gb,byrow=TRUE), 2, geo_mean))
				expr_G1[current_sample,] <- frequencies_expr
			}
			
			##########################################################################
			################### G2 model using training data #########################
			##########################################################################
			
			# generate FacData
			tempS_G2 <- matrix(ncol=ncol, nrow=length(G2))
			rownames(tempS_G2) <- G2
			for (current_sample in 1:length(G2)) {
				# expression
				read_count <- as.numeric(trunc(genedata[G2[current_sample], "read_count"])) + 1 # adding a pseudo-count to avoide problems with 0-count data
				lambdas <- breaksEXPRESSION * UQs[G2[current_sample]]
				frequencies_expr <- rep(0,length(breaksEXPRESSION)-1)
				for (freq in 1:res_expr) {
					frequencies_expr[freq] <- integrate(integrand_e, lower = lambdas[freq], upper = lambdas[freq+1], read_count,stop.on.error=FALSE)[1]
				}
				frequencies_expr <- unlist(frequencies_expr)
				if (length(which(frequencies_expr==0))==res_expr) frequencies_expr[length(frequencies_expr)] <- 1
				frequencies_expr <- frequencies_expr + epsilon_e_G2
				frequencies_expr <- frequencies_expr/sum(frequencies_expr)
				
				# gene body
				cpg_list_gb <- NULL
				for (cpg in 1:length(IDs_body)) {
					miu <- genedata[G2[current_sample], IDs_body[cpg]]
					if (!is.na(miu)) {
						frequencies_gb <- rep(0,res_gb)
						for (freq in 1:res_gb) {
							frequencies_gb[freq] <- integrate(integrand_m,lower=breaksBODY[freq],upper=breaksBODY[freq+1],mean=miu)$value
						}
						frequencies_gb <- unlist(frequencies_gb) + epsilon_gb_G2
						frequencies_gb <- frequencies_gb/sum(frequencies_gb)
						cpg_list_gb[[cpg]] <- frequencies_gb
					} else cpg_list_gb[[cpg]] <- rep(1/res_gb,res_gb)
				}
				
				# promoter
				cpg_list_pr <- NULL
				for (cpg in 1:length(IDs_promoter)) {
					miu <- genedata[G2[current_sample], IDs_promoter[cpg]]
					if (!is.na(miu)) {
						frequencies_pr <- rep(0,res_pr)
						for (freq in 1:res_pr) {
							frequencies_pr[freq] <- integrate(integrand_m,lower=breaksPROMOTER[freq],upper=breaksPROMOTER[freq+1],mean=miu)$value
						}
						frequencies_pr <- unlist(frequencies_pr) + epsilon_pr_G2
						frequencies_pr <- frequencies_pr/sum(frequencies_pr)
						cpg_list_pr[[cpg]] <- frequencies_pr
					} else cpg_list_pr[[cpg]] <- rep(1/res_pr,res_pr)
				}
				
				tempS_formated <- matrix(ncol=ncol, nrow=1)
				tempS_formated[1,1] <- paste('[1,',res_expr,']((',paste(frequencies_expr,sep="",collapse=","),'))',sep="",collapse="")
				cur_ncol <- 1
				for (element in 1:length(cpg_list_pr)) {
					tempS_formated[1,cur_ncol+1] <- paste('[1,',res_pr,']((',paste(cpg_list_pr[[element]],sep="",collapse=","),'))',sep="",collapse="")
					cur_ncol <- cur_ncol + 1
				}
				for (element in 1:length(cpg_list_gb)) {
					tempS_formated[1,cur_ncol+1] <- paste('[1,',res_gb,']((',paste(cpg_list_gb[[element]],sep="",collapse=","),'))',sep="",collapse="")
					cur_ncol <- cur_ncol + 1
				}
				tempS_G2[current_sample,] <- tempS_formated
				#start precomputing correct initialization of parameters
				promoter_G2[current_sample,] <- apply(matrix(unlist(cpg_list_pr), ncol=res_pr, byrow=TRUE), 2,geo_mean)/sum(apply(matrix(unlist(cpg_list_pr), ncol=res_pr, byrow=TRUE), 2, geo_mean))
				body_G2[current_sample,] <- apply(matrix(unlist(cpg_list_gb), ncol=res_gb, byrow=TRUE), 2,geo_mean)/sum(apply(matrix(unlist(cpg_list_gb), ncol=res_gb, byrow=TRUE), 2, geo_mean))
				expr_G2[current_sample,] <- frequencies_expr
			}
			
			########################################
			# build the G1 model using training data
			########################################
			prior_pr <- apply(promoter_G1[train_set_g1,], 2, mean)
			prior_gb <- apply(body_G1[train_set_g1,], 2, mean)
			if (smooth_1d > 0) {
				prior_expr <- kernsm(apply(expr_G1[train_set_g1,], 2, mean), h=smooth_e)
				prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
			} else prior_expr <- apply(expr_G1, 2, mean)
			
			string <- paste(prior_pr, collapse=",")
			promoterPots <- paste("\nNAME:\t\tpot_", promoter_CpGs, "\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,", res_pr,"]((", string, "))\nPC_MAT:\t\t[1,", res_pr, "]((",paste(rep(1, res_pr), collapse=","), "))\n", sep="", collapse="")
			string <- paste(prior_gb,collapse=",")
			geneBodyPots <- paste("\nNAME:\t\tpot_", geneBody_CpGs, "\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,", res_gb, "]((", string, "))\nPC_MAT:\t\t[1,", res_gb,"]((",paste(rep(1, res_gb), collapse=","), "))\n", sep="" , collapse="")
			string <- paste(prior_expr, collapse=",")
			expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,", res_expr, "]((", string, "))\nPC_MAT:\t\t[1,", res_expr, "]((", paste(rep(1,res_expr), collapse=","), "))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,", res_expr, "]((", string, "))\nPC_MAT:\t\t[1,", res_expr,"]((", paste(rep(1, res_expr), collapse=","),"))\n\n", sep="", collapse="")
			
			result <- tensor_product(body_G1[train_set_g1,], expr_G1[train_set_g1,], smooth_h=smooth_gb)
			expr.m <- paste("NAME:\t\tpot_EXPR.M.GB\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[", paste(c(res_expr,res_gb), collapse=","), "]((" ,paste(apply(result , 1, paste, collapse=","), collapse="),\n\t\t\t(") ,"))\nPC_MAT:\t\t[" ,paste(c(res_expr, res_gb), collapse=","),"]((", paste(apply(pseudo_counts_gb, 1, paste, collapse=","), collapse="),\n\t\t\t("), "))\n\n", sep="", collapse="")
			
			result <- tensor_product(promoter_G1[train_set_g1,], expr_G1[train_set_g1,], smooth_h=smooth_pr)
			expr.m <- c(expr.m, paste("NAME:\t\tpot_EXPR.M.P\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[", paste(c(res_expr, res_pr), collapse=","), "]((", paste(apply(result, 1, paste, collapse=","), collapse="),\n\t\t\t("), "))\nPC_MAT:\t\t[", paste(c(res_expr,res_pr), collapse=","), "]((", paste(apply(pseudo_counts_pr,1,paste, collapse=","), collapse="),\n\t\t\t("), "))\n\n", sep="", collapse=""))
			
			potentials <- file(paste("./", cv_run, '/', gene, "/G1_model/factorPotentials.txt", sep=""), "w")
			cat(expr.m, expr.pots, promoterPots, geneBodyPots, file=potentials)
			close(potentials)
			
			########################################
			# build the G2 model using training data
			########################################
			prior_pr <- apply(promoter_G2[train_set_g2,], 2, mean)
			prior_gb <- apply(body_G2[train_set_g2,], 2, mean)
			if (smooth_1d > 0) {
				prior_expr <- kernsm(apply(expr_G2[train_set_g2,], 2, mean), h=smooth_e)
				prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
			} else prior_expr <- apply(expr_G2[train_set_g2,], 2, mean)
			
			string <- paste(prior_pr, collapse=",")
			promoterPots <- paste("\nNAME:\t\tpot_", promoter_CpGs, "\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,", res_pr,"]((", string,"))\nPC_MAT:\t\t[1,", res_pr,"]((", paste(rep(1,res_pr), collapse=","),"))\n", sep="", collapse="")
			string <- paste(prior_gb, collapse=",")
			geneBodyPots <- paste("\nNAME:\t\tpot_", geneBody_CpGs, "\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,", res_gb,"]((", string,"))\nPC_MAT:\t\t[1,", res_gb,"]((", paste(rep(1, res_gb),collapse=",") ,"))\n", sep="", collapse="")
			string <- paste(prior_expr, collapse=",")
			expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,", res_expr, "]((", string,"))\nPC_MAT:\t\t[1,", res_expr, "]((", paste(rep(1,res_expr), collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr, "]((", string, "))\nPC_MAT:\t\t[1,", res_expr, "]((", paste(rep(1, res_expr), collapse=","), "))\n\n", sep="", collapse="")
			
			result <- tensor_product(body_G2, expr_G2, smooth_h=smooth_gb)
			expr.m <- paste("NAME:\t\tpot_EXPR.M.GB\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[", paste(c(res_expr, res_gb), collapse=","), "]((", paste(apply(result, 1, paste, collapse=","), collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[", paste(c(res_expr, res_gb), collapse=","), "]((", paste(apply(pseudo_counts_gb, 1, paste, collapse=","), collapse="),\n\t\t\t("), "))\n\n", sep="", collapse="")
			
			result <- tensor_product(promoter_G2, expr_G2, smooth_h=smooth_pr)
			expr.m <- c(expr.m, paste("NAME:\t\tpot_EXPR.M.P\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[", paste(c(res_expr, res_pr), collapse=","),"]((", paste(apply(result, 1, paste, collapse=","), collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[", paste(c(res_expr,res_pr), collapse=","),"]((", paste(apply(pseudo_counts_pr, 1, paste, collapse=","), collapse="),\n\t\t\t("), "))\n\n", sep="", collapse=""))
			
			potentials <- file(paste("./", cv_run, '/', gene, "/G2_model/factorPotentials.txt", sep=""), "w")
			cat(expr.m, expr.pots, promoterPots, geneBodyPots, file=potentials)
			close(potentials)
			
			#################################################################
			###### generate "missing" Var data for test samples #############
			#################################################################
			tempVar_G1 <- matrix(rep(".", length(test_set_g1) * 3), nrow=length(test_set_g1), ncol=3)
			tempVar_G2 <- matrix(rep(".", length(test_set_g2)*3), nrow=length(test_set_g2), ncol=3)
			colnames(tempVar_G2) <- colnames(tempVar_G1) <- c("NAME:\tEXPR", "M.GB", "M.P")
			rownames(tempVar_G1) <- test_set_g1
			rownames(tempVar_G2) <- test_set_g2
			
			eval(parse(text = paste('write.table(', paste('tempVar_G1, file = "./', cv_run, '/', gene, '/G1_model/G1_VarData_test.tab", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t", append=FALSE)', sep = ""))))
			eval(parse(text = paste('write.table(', paste('tempVar_G2,file = "./', cv_run, '/', gene, '/G2_model/G2_VarData_test.tab", row.names=TRUE, col.names=TRUE, quote=FALSE,sep="\t", append=FALSE)', sep = ""))))
			
			#################################################################
			######## Write the G1 and G2 FacData for test samples ###########
			#################################################################
			tempFac <- tempS_G1[test_set_g1,]
			colnames(tempFac) <- c("NAME:\tEXPR.likelihood", promoter_CpGs, geneBody_CpGs)
			eval(parse(text = paste('write.table(', paste('tempFac, file ="./', cv_run, '/', gene,'/G1_model/G1_FacData_test.tab", row.names=TRUE, col.names=TRUE, quote=FALSE,sep="\t", append=FALSE)', sep = ""))))
			tempFac <- tempS_G2[test_set_g2,]
			colnames(tempFac) <- c("NAME:\tEXPR.likelihood", promoter_CpGs, geneBody_CpGs)
			eval(parse(text = paste('write.table(', paste('tempFac,file ="./', cv_run, '/', gene,'/G2_model/G2_FacData_test.tab", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t", append=FALSE)', sep = ""))))
			###### Query the G1 and G2 with test samples ###########
			
			# query the G1 model with G1 test samples
			string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./', cv_run, '/', gene, '/G1_model/ -l -n - ./', cv_run, '/', gene, '/G1_model/G1_VarData_test.tab ./', cv_run, '/', gene, '/G1_model/G1_FacData_test.tab', sep=""))
			G1test_G1model_mlogliks <- as.numeric(substring(string[-1], 17))
			
			# query the G2 model with G2 test samples
			string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./', cv_run, '/', gene, '/G2_model/ -l -n - ./', cv_run, '/', gene,'/G2_model/G2_VarData_test.tab ./', cv_run, '/', gene, '/G2_model/G2_FacData_test.tab', sep=""))
			G2test_G2model_mlogliks <- as.numeric(substring(string[-1], 17))
			
			
			
			# query the G1 model with G2 test samples
			string<-system(intern=TRUE, command=paste('./dfgEval_static --dfgSpecPrefix=./', cv_run, '/', gene, '/G1_model/ -l -n - ./', cv_run, '/', gene,'/G2_model/G2_VarData_test.tab ./', cv_run, '/', gene,'/G2_model/G2_FacData_test.tab', sep=""))
			G2test_G1model_mlogliks <- as.numeric(substring(string[-1], 17))
			
			# query the G2 model with G1 test samples
			string<-system(intern=TRUE, command=paste('./dfgEval_static --dfgSpecPrefix=./', cv_run, '/', gene,'/G2_model/ -l -n - ./', cv_run, '/', gene,'/G1_model/G1_VarData_test.tab ./', cv_run, '/', gene,'/G1_model/G1_FacData_test.tab', sep=""))
			G1test_G2model_mlogliks <- as.numeric(substring(string[-1], 17))
			
			samples_G1_likelihoods <- c(G1test_G1model_mlogliks, G2test_G1model_mlogliks)
			samples_G2_likelihoods <- c(G1test_G2model_mlogliks, G2test_G2model_mlogliks)
			# G1 posterior probability calculation
			out <- data.frame(sample_ID=c(test_set_g1, test_set_g2), posterior_G1=exp(-samples_G1_likelihoods) / (exp(-samples_G1_likelihoods) + exp(-samples_G2_likelihoods)), G1_mloglik=samples_G1_likelihoods, G2_mloglik=samples_G2_likelihoods)
			colnames(out) <- c("sample_ID", "posterior_G1", "G1_mloglik", "G2_mloglik")
			out
		})
		# transfer the results from parallel execution
		currentGenes_out <- NULL
		for (gene in 1:length(out)) currentGenes_out[[gene]] <- out[[gene]] 
		CV_folds[[cv_run]] <- currentGenes_out
	} else stop("Serial execution currently unavailable, use --parallel >1")
}
save(CV_folds, file=paste(opt$name, ".RData", sep=""))
