# the following objects are expected by this script: workingList_BRCA, mmatrix_pc, counts_plusOne, G1, G2, factors_ls, TSS1500Ind, TSS200Ind, UTR5Ind, EXON1Ind, GENEBODYInd, UTR3Ind

args <- commandArgs(trailingOnly = TRUE)
beg <- as.numeric(args[1]) # first ID to process
end <- as.numeric(args[2]) # last consequitve ID to process
#beg <- end <- 1

# data <- args[1]
load("./essentials_aggressive.RData") # data to work on

G1_template <- G1
G2_template <- G2

# configure method here
res_pr <- 25 # number of bins for promoter methylation
res_gb <- 25 # number of bins for gene body methylation
res_expr <- 25 # number of bins for gene expression
nruns <- 100 # number of permutations to perform for null distribution calculation

smooth_1d <- 1 # should expression node parameterization be smoothed?
smooth_2d <- 1 # should methylation 2D factors parameterizations be smoothed?
if (smooth_1d > 0) library(aws)
if (smooth_2d > 0) library(smoothie)


integrand_e <- function(x,k) {dpois(k,x)}
integrand_m <- function(x,mean) {dnorm(x=mean,mean=x,sd=0.14)}

tensor_product <- function(matrix1,matrix2,smooth_h=0,normalize=c("row","column","no"),kernel=c("gauss","cauchy","minvar")) {
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

for (i in beg:end){
	cat(paste("doing ",i,"\n",sep=""))
	ptm <- proc.time()[3]
	system(command=paste('mkdir',i,sep=" "))
	system(command=paste('mkdir ./',i,'/G1_model',sep=""))
	system(command=paste('mkdir ./',i,'/G1_model/all',sep=""))
	system(command=paste('mkdir ./',i,'/G2_model',sep=""))
	system(command=paste('mkdir ./',i,'/G2_model/all',sep=""))
	system(command=paste('mkdir ./',i,'/full_model',sep=""))
	system(command=paste('mkdir ./',i,'/null',sep=""))
	system(command=paste('mkdir ./',i,'/null/G1_model',sep=""))
	system(command=paste('mkdir ./',i,'/null/G2_model',sep=""))
	
	########################################################
	# dynamic generation of model specification files here #
	########################################################
	####################### state map ######################
	stateMaps <- file(paste("./",i,"/stateMaps.txt",sep=""),"w")
	exprMap <- paste("NAME:\texprMap\nSYMBOLS:\t",paste(seq(1,res_expr,1),collapse=" "),"\nMETA_SYMBOLS:\t.=",paste(seq(1,res_expr,1),collapse=" "),"; *=",paste(seq(1,res_expr,1),collapse=" "),";\n\n",collapse="",sep="")
	prMap <- paste("NAME:\tprMap\nSYMBOLS:\t",paste(seq(1,res_pr,1),collapse=" "),"\nMETA_SYMBOLS:\t.=",paste(seq(1,res_pr,1),collapse=" "),"; *=",paste(seq(1,res_pr,1),collapse=" "),";\n\n",collapse="",sep="")
	gbMap <- paste("NAME:\tgbMap\nSYMBOLS:\t",paste(seq(1,res_gb,1),collapse=" "),"\nMETA_SYMBOLS:\t.=",paste(seq(1,res_gb,1),collapse=" "),"; *=",paste(seq(1,res_gb,1),collapse=" "),";\n",collapse="",sep="")
	cat(exprMap,prMap,gbMap,file=stateMaps)
	close(stateMaps)
	########################################################
	####################### variables ######################
	variables <- file(paste("./",i,"/variables.txt",sep=""),"w")
	cat("STATE_MAP_NAME:\texprMap\nVAR_NAMES:\tEXPR\n\n",sep="",file=variables)
	cat("STATE_MAP_NAME:\tprMap\nVAR_NAMES:\tM.P\n\n",sep="",file=variables)
	cat("STATE_MAP_NAME:\tgbMap\nVAR_NAMES:\tM.GB\n",sep="",file=variables)
	close(variables)
	###########################################################
	############### identify constitutive CpGs ################
	IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',workingList_BRCA[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',workingList_BRCA[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',workingList_BRCA[i],'"',sep=""))))
	promoterVars <- promoterVars_template[1:length(IDs_promoter)]
	promoter_CpGs <- template_promoter_CpGs[1:length(IDs_promoter)]
	IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',workingList_BRCA[i],'"',sep="")))))
	geneBodyVars <- geneBodyVars_template[1:length(IDs_body)]
	geneBody_CpGs <- template_body_CpGs[1:length(IDs_body)]
	ncol = length(IDs_promoter) + length(IDs_body) + 1
	
	############################################################
	######### calculate epsilons and smoothing params ##########
	epsilon_pr_G2 <- 1/(length(G2)*length(IDs_promoter))/res_pr
	epsilon_gb_G2 <- 1/(length(G2)*length(IDs_body))/res_gb
	epsilon_e_G2 <- 1/length(G2)/res_expr
	epsilon_pr_G1 <- 1/(length(G1)*length(IDs_promoter))/res_pr
	epsilon_gb_G1 <- 1/(length(G1)*length(IDs_body))/res_gb
	epsilon_e_G1 <- 1/length(G1)/res_expr
	
	smooth_e <- 10/(mean(length(G1),length(G2))/(res_expr*2))
	smooth_pr <- trunc(10/(mean(length(G1),length(G2))*length(IDs_promoter)/(res_pr*res_expr)))
	smooth_gb <- trunc(10/(mean(length(G1),length(G2))*length(IDs_body)/(res_gb*res_expr)))
	########################################################
	##################### factor graph #####################
	factorGraph <- file(paste("./",i,"/factorGraph.txt",sep=""),"w")
	cat("NAME:\tEXPR.likelihood\nNB1:\tEXPR\nPOT:\tpot_EXPR.likelihood\n",file=factorGraph)
	cat("\nNAME:\tEXPR.prior\nNB1:\tEXPR\nPOT:\tpot_EXPR.prior\n",file=factorGraph)
	cat(paste("\nNAME:\tEXPR.M.GB\nNB1:\tEXPR\nNB2:\tM.GB\nPOT:\tpot_EXPR.M.GB\n",sep="",collapse=""),file=factorGraph)
	cat(paste("\nNAME:\tEXPR.M.P\nNB1:\tEXPR\nNB2:\tM.P\nPOT:\tpot_EXPR.M.P\n",sep="",collapse=""),file=factorGraph)
	cat(paste("\nNAME:\t",promoter_CpGs,"\nNB1:\tM.P\nPOT:\tpot_",promoter_CpGs,"\n",sep="",collapse=""),file=factorGraph)
	cat(paste("\nNAME:\t",geneBody_CpGs,"\nNB1:\tM.GB\nPOT:\tpot_",geneBody_CpGs,"\n",sep="",collapse=""),file=factorGraph)
	close(factorGraph)
	
	# distribute model definition files to all neccessary directories
	system(command=paste('cp ./',i,'/*.txt ./',i,'/G1_model/all',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/G2_model/all',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/full_model',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/null/G1_model',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/null/G2_model',sep=""))
	
	# begin the 7-fold x-val
	# initialize
	Zs <- vector(mode="numeric",length=7)
	G1_G1model_mlogliks <- vector(mode="numeric",length=length(G1_template))
	G2_G1model_mlogliks <- vector(mode="numeric",length=length(G2_template))
	G1_G2model_mlogliks <- vector(mode="numeric",length=length(G1_template))
	G2_G2model_mlogliks <- vector(mode="numeric",length=length(G2_template))
	
	# start x-val outer loop
	for (fold in 0:6) {
		cat(paste("doing ",i," fold ",fold+1,"\n",sep=""))
		ptm_f <- proc.time()[3]
		# find indexes for current fold
		index_G1_fold <- (fold*2+1):((fold+1)*2)
		index_G2_fold <- (fold*8+1):((fold+1)*8)
		if (fold == 6) index_G2_fold <- c(index_G2_fold,(1:length(G2_template))[-(1:((fold+1)*8))])
		# subset the groups for current fold
		G1_eval <- G1_template[index_G1_fold]
		G2_eval <- G2_template[index_G2_fold]
		G1 <- G1_template[-index_G1_fold]
		G2 <- G2_template[-index_G2_fold]
		
		###########################################################
		############### generate "missing" Var data ################
		tempVar_G1 <- matrix(rep(".",length(G1)*3),nrow=length(G1),ncol=3)
		colnames(tempVar_G1) <- c("NAME:\tEXPR","M.GB","M.P")
		rownames(tempVar_G1) <- G1
		eval(parse(text = paste('write.table(', paste('tempVar_G1,file = "./',i,'/G1_model/all/G1_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		tempVar_G2 <- matrix(rep(".",length(G2)*3),nrow=length(G2),ncol=3)
		colnames(tempVar_G2) <- c("NAME:\tEXPR","M.GB","M.P")
		rownames(tempVar_G2) <- G2
		eval(parse(text = paste('write.table(', paste('tempVar_G2,file = "./',i,'/G2_model/all/G2_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		tempVar_all <- rbind(tempVar_G2,tempVar_G1)
		eval(parse(text = paste('write.table(', paste('tempVar_all,file = "./',i,'/full_model/full_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		###########################################################
		############### define the binning scheme  ################
		all_labels_pr <- as.character(seq(1,res_pr,1))
		all_labels_gb <- as.character(seq(1,res_gb,1))
		all_labels_expr <- as.character(seq(1,res_expr,1))
		promoter_G2 <- matrix(ncol=res_pr,nrow=length(G2))
		promoter_G1 <- matrix(ncol=res_pr,nrow=length(G1))
		body_G2 <- matrix(ncol=res_gb,nrow=length(G2))
		body_G1 <- matrix(ncol=res_gb,nrow=length(G1))
		expr_G2 <- matrix(ncol=res_expr,nrow=length(G2))
		expr_G1 <- matrix(ncol=res_expr,nrow=length(G1))
		pseudo_counts_pr <- matrix(ncol=res_pr,nrow=res_expr,data=rep(1,res_expr*res_pr))
		pseudo_counts_gb <- matrix(ncol=res_gb,nrow=res_expr,data=rep(1,res_expr*res_gb))
		
		# gene body
		x <- mmatrix_pc[IDs_body,c(G2,G1)]
		
		density <- density(x,bw=0.14,from=-7,to=7,n=2801,na.rm=TRUE)
		density$y <- density$y/sum(density$y)
		density$y <- cumsum(density$y)
		breaks <- NULL
		noBreaks <- res_gb-1
		for (j in 1:noBreaks) { breaks <- c (breaks, density$x[which(density$y >= j*(1/(1+noBreaks)))][1])}
		breaksBODY <- sort(c(-7.01,breaks,7.01))
		
		# promoter
		x <- mmatrix_pc[IDs_promoter,c(G2,G1)]
		
		density <- density(x,bw=0.14,from=-7,to=7,n=2801,na.rm=TRUE)
		density$y <- density$y/sum(density$y)
		density$y <- cumsum(density$y)
		breaks <- NULL
		noBreaks <- res_pr-1
		for (j in 1:noBreaks) { breaks <- c (breaks, density$x[which(density$y >= j*(1/(1+noBreaks)))][1])}
		breaksPROMOTER <- sort(c(-7.01,breaks,7.01))
		
		# expression
		temp <- counts_plusOne[workingList_BRCA[i],c(G2,G1)]
		temp_cpms <- matrix(ncol=2)
		colnames(temp_cpms) <- c("cpm","density")
		for (j in 1:length(temp)) {
			lambda <- as.numeric(temp[j])
			X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
			current <- factors_ls[c(G2,G1)[j]]
			temp_cpms <- rbind(temp_cpms,cbind(X/current,dpois(X,lambda=lambda)*current))
		}
		temp_cpms <- as.data.frame(temp_cpms[-1,],)
		temp_cpms <- temp_cpms[order(temp_cpms$cpm),]
		temp_cpms[,3] <- cumsum(temp_cpms[,2])
		temp_cpms[,3] <- temp_cpms[,3]/max(temp_cpms[,3])
		breaks <- NULL
		noBreaks <- res_expr-1
		for (j in 1:noBreaks) { breaks <- c (breaks, temp_cpms[which(temp_cpms[,3] >= j*(1/(1+noBreaks))),1][1])}
		max <- max(counts_plusOne[workingList_BRCA[i],c(G2,G1)])
		max_boundary <- 20+max+(4*max*max^(-1/2))
		breaksEXPRESSION <- c(0,breaks,max_boundary)
		
		###########################################################################
		############################## G1 model ###################################
		## full G1 model developed from here, to obtain likelihoods of G1 ########
		
		# generate FacData for full set of G1
		tempS_G1 <- matrix(ncol=ncol,nrow=length(G1))
		for (current_sample in 1:length(G1)) {
			# expression
			read_count <- as.numeric(trunc(counts_plusOne[workingList_BRCA[i],G1[current_sample]]))
			lambdas <- breaksEXPRESSION * factors_ls[G1[current_sample]]
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
			for (cpg in 1:length(geneBodyVars)) {
				miu <- mmatrix_pc[IDs_body[cpg],G1[current_sample]]
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
			for (cpg in 1:length(promoterVars)) {
				miu <- mmatrix_pc[IDs_promoter[cpg],G1[current_sample]]
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
			
			tempS_formated <- matrix(ncol=ncol,nrow=1)
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
			tempS_G1[current_sample,] <- tempS_formated
			
			#start precomputing correct initialization of parameters
			promoter_G1[current_sample,] <- apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean)/sum(apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean))
			body_G1[current_sample,] <- apply(matrix(unlist(cpg_list_gb),ncol=res_gb,byrow=TRUE),2,geo_mean)/sum(apply(matrix(unlist(cpg_list_gb),ncol=res_gb,byrow=TRUE),2,geo_mean))
			expr_G1[current_sample,] <- frequencies_expr
		}
		
		# build the full G1 model
		# precompute correct initialization of parameters for G1-only model
		prior_pr <- apply(promoter_G1,2,mean)
		prior_gb <- apply(body_G1,2,mean)
		if (smooth_1d > 0) {
			prior_expr <- kernsm(apply(expr_G1,2,mean),h=smooth_e)
			prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
		} else prior_expr <- apply(expr_G1,2,mean)
		
		string <- paste(prior_pr,collapse=",")
		promoterPots <- paste("\nNAME:\t\tpot_",promoter_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="",collapse="")
		string <- paste(prior_gb,collapse=",")
		geneBodyPots <- paste("\nNAME:\t\tpot_",geneBody_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="",collapse="")
		string <- paste(prior_expr,collapse=",")
		expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
		
		result <- tensor_product(body_G1,expr_G1,smooth_h=smooth_gb)
		expr.m <- paste("NAME:\t\tpot_EXPR.M.GB\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(pseudo_counts_gb,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse="")
		
		result <- tensor_product(promoter_G1,expr_G1,smooth_h=smooth_pr)
		expr.m <- c(expr.m,paste("NAME:\t\tpot_EXPR.M.P\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(pseudo_counts_pr,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse=""))
		
		potentials <- file(paste("./",i,"/G1_model/all/factorPotentials.txt",sep=""),"w")
		cat(expr.m,expr.pots,promoterPots,geneBodyPots,file=potentials)
		close(potentials)
		
		tempFac <- tempS_G1
		colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs,geneBody_CpGs)
		rownames(tempFac) <- G1
		eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/G1_model/all/G1_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		# query the full G1 model with G1 samples
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/G1_model/all/ -l -n - ./',i,'/G1_model/all/G1_VarData.tab ./',i,'/G1_model/all/G1_FacData.tab',sep=""))
		G1_G1model_mlogliks_cur <- as.numeric(substring(string[-1],14))
		
		###########################################################################
		############################## G2 model ###################################
		### full G2 model developed from here, to obtain likelihoods of G2 #######

		# generate FacData for full set of G2
		tempS_G2 <- matrix(ncol=ncol,nrow=length(G2))
		for (current_sample in 1:length(G2)) {
			# expression
			read_count <- as.numeric(trunc(counts_plusOne[workingList_BRCA[i],G2[current_sample]]))
			lambdas <- breaksEXPRESSION * factors_ls[G2[current_sample]]
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
			for (cpg in 1:length(geneBodyVars)) {
				miu <- mmatrix_pc[IDs_body[cpg],G2[current_sample]]
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
			for (cpg in 1:length(promoterVars)) {
				miu <- mmatrix_pc[IDs_promoter[cpg],G2[current_sample]]
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
			
			tempS_formated <- matrix(ncol=ncol,nrow=1)
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
			promoter_G2[current_sample,] <- apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean)/sum(apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean))
			body_G2[current_sample,] <- apply(matrix(unlist(cpg_list_gb),ncol=res_gb,byrow=TRUE),2,geo_mean)/sum(apply(matrix(unlist(cpg_list_gb),ncol=res_gb,byrow=TRUE),2,geo_mean))
			expr_G2[current_sample,] <- frequencies_expr
		}
		
		# build the full G2 model
		# precompute correct initialization of parameters for G2-only model
		prior_pr <- apply(promoter_G2,2,mean)
		prior_gb <- apply(body_G2,2,mean)
		if (smooth_1d > 0) {
			prior_expr <- kernsm(apply(expr_G2,2,mean),h=smooth_e)
			prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
		} else prior_expr <- apply(expr_G2,2,mean)
		
		string <- paste(prior_pr,collapse=",")
		promoterPots <- paste("\nNAME:\t\tpot_",promoter_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="",collapse="")
		string <- paste(prior_gb,collapse=",")
		geneBodyPots <- paste("\nNAME:\t\tpot_",geneBody_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="",collapse="")
		string <- paste(prior_expr,collapse=",")
		expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
		
		result <- tensor_product(body_G2,expr_G2,smooth_h=smooth_gb)
		expr.m <- paste("NAME:\t\tpot_EXPR.M.GB\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(pseudo_counts_gb,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse="")
		
		result <- tensor_product(promoter_G2,expr_G2,smooth_h=smooth_pr)
		expr.m <- c(expr.m,paste("NAME:\t\tpot_EXPR.M.P\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(pseudo_counts_pr,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse=""))
		
		potentials <- file(paste("./",i,"/G2_model/all/factorPotentials.txt",sep=""),"w")
		cat(expr.m,expr.pots,promoterPots,geneBodyPots,file=potentials)
		close(potentials)
		
		tempFac <- tempS_G2
		colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs,geneBody_CpGs)
		rownames(tempFac) <- G2
		eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/G2_model/all/G2_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		# query the full model with G2 samples
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/G2_model/all/ -l -n - ./',i,'/G2_model/all/G2_VarData.tab ./',i,'/G2_model/all/G2_FacData.tab',sep=""))
		G2_G2model_mlogliks_cur <- as.numeric(substring(string[-1],14))
		
		###########################################################################
		######################## All data model ###################################
		## Full model developed from here, to obtain likelihoods of G2 and G1 ####
		# precompute correct initialization of parameters for joint G1-G2 model
		promoter_all <- rbind(promoter_G2,promoter_G1)
		body_all <- rbind(body_G2,body_G1)
		expr_all <- rbind(expr_G2,expr_G1)
		
		prior_pr <- apply(promoter_all,2,mean)
		prior_gb <- apply(body_all,2,mean)
		if (smooth_1d > 0) {
			prior_expr <- kernsm(apply(expr_all,2,mean),h=smooth_e)
			prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
		} else prior_expr <- apply(expr_all,2,mean)
		
		string <- paste(prior_pr,collapse=",")
		promoterPots <- paste("\nNAME:\t\tpot_",promoter_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="",collapse="")
		string <- paste(prior_gb,collapse=",")
		geneBodyPots <- paste("\nNAME:\t\tpot_",geneBody_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="",collapse="")
		string <- paste(prior_expr,collapse=",")
		expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
		
		result <- tensor_product(body_all,expr_all,smooth_h=smooth_gb)
		expr.m <- paste("NAME:\t\tpot_EXPR.M.GB\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(pseudo_counts_gb,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse="")
		
		result <- tensor_product(promoter_all,expr_all,smooth_h=smooth_pr)
		expr.m <- c(expr.m,paste("NAME:\t\tpot_EXPR.M.P\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(pseudo_counts_pr,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse=""))
		
		potentials <- file(paste("./",i,"/full_model/factorPotentials.txt",sep=""),"w")
		cat(expr.m,expr.pots,promoterPots,geneBodyPots,file=potentials)
		close(potentials)
		
		tempFac <- rbind(tempS_G2,tempS_G1)
		colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs,geneBody_CpGs)
		rownames(tempFac) <- c(G2,G1)
		eval(parse(text = paste('write.table(', paste('tempFac,file = "./',i,'/full_model/full_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		# query the full model with T and AN samples
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/full_model/ -l -n - ./',i,'/full_model/full_VarData.tab ./',i,'/full_model/full_FacData.tab',sep=""))
		allData_jointModel_mlogliks_cur <- as.numeric(substring(string[-1],14))
		###########################################################################
		######################### D calculation ###################################
		D <- 2*(sum(allData_jointModel_mlogliks_cur) - (sum(G1_G1model_mlogliks_cur)+sum(G2_G2model_mlogliks_cur)))
		###########################################################################
		################# Z-score calculation using null distr. #####################
		Ds <- vector(length=nruns,mode="numeric")
		for (run in 1:nruns) {
			cur <- sample(x=1:(length(G2)+length(G1)),size=length(G2),replace=FALSE)
			# G2
			tempFac_G2 <- tempFac[cur,]
			rownames(tempFac_G2) <- G2
			eval(parse(text = paste('write.table(', paste('tempFac_G2,file = "./',i,'/null/G2_model/G2_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
			
			prior_pr <- apply(promoter_all[cur,],2,mean)
			prior_gb <- apply(body_all[cur,],2,mean)
			if (smooth_1d > 0) {
				prior_expr <- kernsm(apply(expr_all[cur,],2,mean),h=smooth_e)
				prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
			} else prior_expr <- apply(expr_all[cur,],2,mean)
			
			string <- paste(prior_pr,collapse=",")
			promoterPots <- paste("\nNAME:\t\tpot_",promoter_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="",collapse="")
			string <- paste(prior_gb,collapse=",")
			geneBodyPots <- paste("\nNAME:\t\tpot_",geneBody_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="",collapse="")
			string <- paste(prior_expr,collapse=",")
			expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
			
			result <- tensor_product(body_all[cur,],expr_all[cur,],smooth_h=smooth_gb)
			expr.m <- paste("NAME:\t\tpot_EXPR.M.GB\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(pseudo_counts_gb,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse="")
			
			result <- tensor_product(promoter_all[cur,],expr_all[cur,],smooth_h=smooth_pr)
			expr.m <- c(expr.m,paste("NAME:\t\tpot_EXPR.M.P\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(pseudo_counts_pr,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse=""))
			
			potentials <- file(paste("./",i,"/null/G2_model/factorPotentials.txt",sep=""),"w")
			cat(expr.m,expr.pots,promoterPots,geneBodyPots,file=potentials)
			close(potentials)
			
			# query
			string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/null/G2_model/ -l -n - ./',i,'/G2_model/all/G2_VarData.tab ./',i,'/null/G2_model/G2_FacData.tab',sep=""))
			G2_G2model_mlogliks_run <- as.numeric(substring(string[-1],14))
			
			# G1
			tempFac_G1 <- tempFac[-cur,]
			rownames(tempFac_G1) <- G1
			eval(parse(text = paste('write.table(', paste('tempFac_G1,file = "./',i,'/null/G1_model/G1_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
			
			prior_pr <- apply(promoter_all[-cur,],2,mean)
			prior_gb <- apply(body_all[-cur,],2,mean)
			if (smooth_1d > 0) {
				prior_expr <- kernsm(apply(expr_all[-cur,],2,mean),h=smooth_e)
				prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
			} else prior_expr <- apply(expr_all[-cur,],2,mean)
			
			string <- paste(prior_pr,collapse=",")
			promoterPots <- paste("\nNAME:\t\tpot_",promoter_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="",collapse="")
			string <- paste(prior_gb,collapse=",")
			geneBodyPots <- paste("\nNAME:\t\tpot_",geneBody_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="",collapse="")
			string <- paste(prior_expr,collapse=",")
			expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
			
			result <- tensor_product(body_all[-cur,],expr_all[-cur,],smooth_h=smooth_gb)
			expr.m <- paste("NAME:\t\tpot_EXPR.M.GB\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(pseudo_counts_gb,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse="")
			
			result <- tensor_product(promoter_all[-cur,],expr_all[-cur,],smooth_h=smooth_pr)
			expr.m <- c(expr.m,paste("NAME:\t\tpot_EXPR.M.P\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(pseudo_counts_pr,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse=""))
			
			potentials <- file(paste("./",i,"/null/G1_model/factorPotentials.txt",sep=""),"w")
			cat(expr.m,expr.pots,promoterPots,geneBodyPots,file=potentials)
			close(potentials)
			
			# query
			string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/null/G1_model/ -l -n - ./',i,'/G1_model/all/G1_VarData.tab ./',i,'/null/G1_model/G1_FacData.tab',sep=""))
			G1_G1model_mlogliks_run <- as.numeric(substring(string[-1],14))
			
			Ds[run] <- 2*(sum(allData_jointModel_mlogliks_cur) - (sum(G1_G1model_mlogliks_run)+sum(G2_G2model_mlogliks_run)))
		}
		if (sd(Ds) != 0) Zs[fold+1] <- (D - mean(Ds)) / sd(Ds) else Zs[fold+1] <- -6
		
		###########################################################################
		###################### process G1 eval samples ############################
		
		tempS_G1 <- matrix(ncol=ncol,nrow=length(G1_eval))
		for (current_sample in 1:length(G1_eval)) {
			# expression
			read_count <- as.numeric(trunc(counts_plusOne[workingList_BRCA[i],G1_eval[current_sample]]))
			lambdas <- breaksEXPRESSION * factors_ls[G1_eval[current_sample]]
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
			for (cpg in 1:length(geneBodyVars)) {
				miu <- mmatrix_pc[IDs_body[cpg],G1_eval[current_sample]]
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
			for (cpg in 1:length(promoterVars)) {
				miu <- mmatrix_pc[IDs_promoter[cpg],G1_eval[current_sample]]
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
			
			tempS_formated <- matrix(ncol=ncol,nrow=1)
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
			tempS_G1[current_sample,] <- tempS_formated
		}
		
		tempFac <- tempS_G1
		colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs,geneBody_CpGs)
		rownames(tempFac) <- G1_eval
		eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/G1_model/all/G1_eval_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		tempVar_G1 <- matrix(rep(".",length(G1_eval)*3),nrow=length(G1_eval),ncol=3)
		colnames(tempVar_G1) <- c("NAME:\tEXPR","M.GB","M.P")
		rownames(tempVar_G1) <- G1_eval
		eval(parse(text = paste('write.table(', paste('tempVar_G1,file = "./',i,'/G1_model/all/G1_eval_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		
		# query the full G1 model with G1 eval samples
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/G1_model/all/ -l -n - ./',i,'/G1_model/all/G1_eval_VarData.tab ./',i,'/G1_model/all/G1_eval_FacData.tab',sep=""))
		G1_G1model_mlogliks[index_G1_fold] <- as.numeric(substring(string[-1],14))
		# query the full G2 model with G1 samples
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/G2_model/all/ -l -n - ./',i,'/G1_model/all/G1_eval_VarData.tab ./',i,'/G1_model/all/G1_eval_FacData.tab',sep=""))
		G1_G2model_mlogliks[index_G1_fold] <- as.numeric(substring(string[-1],14))
		###########################################################################
		###################### process G2 eval samples ############################
		
		tempS_G2 <- matrix(ncol=ncol,nrow=length(G2_eval))
		for (current_sample in 1:length(G2_eval)) {
			# expression
			read_count <- as.numeric(trunc(counts_plusOne[workingList_BRCA[i],G2_eval[current_sample]]))
			lambdas <- breaksEXPRESSION * factors_ls[G2_eval[current_sample]]
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
			for (cpg in 1:length(geneBodyVars)) {
				miu <- mmatrix_pc[IDs_body[cpg],G2_eval[current_sample]]
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
			for (cpg in 1:length(promoterVars)) {
				miu <- mmatrix_pc[IDs_promoter[cpg],G2_eval[current_sample]]
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
			
			tempS_formated <- matrix(ncol=ncol,nrow=1)
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
		}
		
		tempFac <- tempS_G2
		colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs,geneBody_CpGs)
		rownames(tempFac) <- G2_eval
		eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/G2_model/all/G2_eval_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		tempVar_G2 <- matrix(rep(".",length(G2_eval)*3),nrow=length(G2_eval),ncol=3)
		colnames(tempVar_G2) <- c("NAME:\tEXPR","M.GB","M.P")
		rownames(tempVar_G2) <- G2_eval
		eval(parse(text = paste('write.table(', paste('tempVar_G2,file = "./',i,'/G2_model/all/G2_eval_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		
		# query the full G2 model with G2 eval samples
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/G2_model/all/ -l -n - ./',i,'/G2_model/all/G2_eval_VarData.tab ./',i,'/G2_model/all/G2_eval_FacData.tab',sep=""))
		G2_G2model_mlogliks[index_G2_fold] <- as.numeric(substring(string[-1],14))
		# query the full G1 model with G2 samples
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/G1_model/all/ -l -n - ./',i,'/G2_model/all/G2_eval_VarData.tab ./',i,'/G2_model/all/G2_eval_FacData.tab',sep=""))
		G2_G1model_mlogliks[index_G2_fold] <- as.numeric(substring(string[-1],14))
		
		cat(paste("done fold no. ", fold+1," in ", sprintf("%.2f", (proc.time()[3]-ptm_f)/60)," minutes\n",sep=""))
	} # done here with the 7-fold x-val
	###########################################################################################
	################# put together the output ##########################################
	samples_G1model_likelihoods <- c(G1_G1model_mlogliks,G2_G1model_mlogliks)
	samples_G2model_likelihoods <- c(G1_G2model_mlogliks,G2_G2model_mlogliks)
	out <- cbind(c(G1_template,G2_template),samples_G1model_likelihoods,samples_G2model_likelihoods)
	colnames(out) <- c("sample_ID","G1_mloglik","G2_mloglik")
	eval(parse(text=paste('write.table(x=t(Zs), col.names=FALSE, row.names=FALSE, append=TRUE, file="./',i,'.result")',sep="")))
	eval(parse(text=paste('write.table(as.data.frame(out), col.names=TRUE, row.names=FALSE, quote=FALSE, append=TRUE, file="./',i,'.result")',sep="")))
	
	cat(paste("done 7-fold x-val for ",i," in ", sprintf("%.2f", (proc.time()[3]-ptm)/60)," minutes\n",sep=""))
}
