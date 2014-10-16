# the following objects are expected by this script: workingList_BRCA, mmatrix_pc, counts_plusOne, G1, G2, factors_ls, TSS1500Ind, TSS200Ind, UTR5Ind, EXON1Ind, GENEBODYInd, UTR3Ind

args <- commandArgs(trailingOnly = TRUE)
beg <- as.numeric(args[1]) # first ID to process
end <- as.numeric(args[2]) # last consequitve ID to process

# data <- args[1]
load("./essentials_aggressive.RData") # data to work on

# configure here
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
	system(command=paste('mkdir ./',i,'/LOU',sep=""))
	system(command=paste('mkdir ./',i,'/LOU/G1_model',sep=""))
	system(command=paste('mkdir ./',i,'/LOU/G2_model',sep=""))
	
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
	
	system(command=paste('cp ./',i,'/*.txt ./',i,'/G1_model/all',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/G2_model/all',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/LOU/G1_model',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/LOU/G2_model',sep=""))
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
	#########################################################################
	###### Query the full G1 and G2 models with G2 and G1 samples ###########
	# query the full G1 model with G2 samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/G1_model/all/ -l -n - ./',i,'/G2_model/all/G2_VarData.tab ./',i,'/G2_model/all/G2_FacData.tab',sep=""))
	G2_G1model_mlogliks <- as.numeric(substring(string[-1],14))
	
	# query the full G2 model with G1 samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/G2_model/all/ -l -n - ./',i,'/G1_model/all/G1_VarData.tab ./',i,'/G1_model/all/G1_FacData.tab',sep=""))
	G1_G2model_mlogliks <- as.numeric(substring(string[-1],14))
	
	#########################################################################
	####################### LOU loglik calculations #########################
	# G2 LOU
	G2_G2model_mlogliks <- vector(mode="numeric",length=length(G2))
	for (cur in 1:length(G2)) {
		# G2
		tempVar <- tempVar_G2[cur,]
		tempVar[1] <- paste(G2[cur],"\t",tempVar[1],sep="")
		tempVar <- rbind(names(tempVar),tempVar)
		tempFac_G2 <- tempS_G2[cur,]
		tempFac_G2[1] <- paste(G2[cur],"\t",tempFac_G2[1],sep="")
		tempFac_G2 <- rbind(c("NAME:\tEXPR.likelihood",promoter_CpGs,geneBody_CpGs),tempFac_G2)
		eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/LOU/G2_model/G2_VarData.tab",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		eval(parse(text = paste('write.table(', paste('tempFac_G2,file = "./',i,'/LOU/G2_model/G2_FacData.tab",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		prior_pr <- apply(promoter_G2[-cur,],2,mean)
		prior_gb <- apply(body_G2[-cur,],2,mean)
		if (smooth_1d > 0) {
			prior_expr <- kernsm(apply(expr_G2[-cur,],2,mean),h=smooth_e)
			prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
		} else prior_expr <- apply(expr_G2[-cur,],2,mean)
		
		string <- paste(prior_pr,collapse=",")
		promoterPots <- paste("\nNAME:\t\tpot_",promoter_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="",collapse="")
		string <- paste(prior_gb,collapse=",")
		geneBodyPots <- paste("\nNAME:\t\tpot_",geneBody_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="",collapse="")
		string <- paste(prior_expr,collapse=",")
		expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
		
		result <- tensor_product(body_G2[-cur,],expr_G2[-cur,],smooth_h=smooth_gb)
		expr.m <- paste("NAME:\t\tpot_EXPR.M.GB\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(pseudo_counts_gb,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse="")
		
		result <- tensor_product(promoter_G2[-cur,],expr_G2[-cur,],smooth_h=smooth_pr)
		expr.m <- c(expr.m,paste("NAME:\t\tpot_EXPR.M.P\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(pseudo_counts_pr,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse=""))
		
		potentials <- file(paste("./",i,"/LOU/G2_model/factorPotentials.txt",sep=""),"w")
		cat(expr.m,expr.pots,promoterPots,geneBodyPots,file=potentials)
		close(potentials)
		
		# query
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/LOU/G2_model/ -l -n - ./',i,'/LOU/G2_model/G2_VarData.tab ./',i,'/LOU/G2_model/G2_FacData.tab',sep=""))
		G2_G2model_mlogliks[cur] <- as.numeric(substring(string[-1],14))
	}
	
	G1_G1model_mlogliks <- vector(mode="numeric",length=length(G1))
	for (cur in 1:length(G1)) {
		# G1
		tempVar <- tempVar_G1[cur,]
		tempVar[1] <- paste(G1[cur],"\t",tempVar[1],sep="")
		tempVar <- rbind(names(tempVar),tempVar)
		tempFac_G1 <- tempS_G1[cur,]
		tempFac_G1[1] <- paste(G1[cur],"\t",tempFac_G1[1],sep="")
		tempFac_G1 <- rbind(c("NAME:\tEXPR.likelihood",promoter_CpGs,geneBody_CpGs),tempFac_G1)
		eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/LOU/G1_model/G1_VarData.tab",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		eval(parse(text = paste('write.table(', paste('tempFac_G1,file = "./',i,'/LOU/G1_model/G1_FacData.tab",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		prior_pr <- apply(promoter_G1[-cur,],2,mean)
		prior_gb <- apply(body_G1[-cur,],2,mean)
		if (smooth_1d > 0) {
			prior_expr <- kernsm(apply(expr_G1[-cur,],2,mean),h=smooth_e)
			prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
		} else prior_expr <- apply(expr_G1[-cur,],2,mean)
		
		string <- paste(prior_pr,collapse=",")
		promoterPots <- paste("\nNAME:\t\tpot_",promoter_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="",collapse="")
		string <- paste(prior_gb,collapse=",")
		geneBodyPots <- paste("\nNAME:\t\tpot_",geneBody_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="",collapse="")
		string <- paste(prior_expr,collapse=",")
		expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
		
		result <- tensor_product(body_G1[-cur,],expr_G1[-cur,],smooth_h=smooth_gb)
		expr.m <- paste("NAME:\t\tpot_EXPR.M.GB\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_gb),collapse=","),"]((",paste(apply(pseudo_counts_gb,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse="")
		
		result <- tensor_product(promoter_G1[-cur,],expr_G1[-cur,],smooth_h=smooth_pr)
		expr.m <- c(expr.m,paste("NAME:\t\tpot_EXPR.M.P\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_pr),collapse=","),"]((",paste(apply(pseudo_counts_pr,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\n\n",sep="",collapse=""))
		
		potentials <- file(paste("./",i,"/LOU/G1_model/factorPotentials.txt",sep=""),"w")
		cat(expr.m,expr.pots,promoterPots,geneBodyPots,file=potentials)
		close(potentials)
		
		# query
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/LOU/G1_model/ -l -n - ./',i,'/LOU/G1_model/G1_VarData.tab ./',i,'/LOU/G1_model/G1_FacData.tab',sep=""))
		G1_G1model_mlogliks[cur] <- as.numeric(substring(string[-1],14))
	}
	###########################################################################################
	samples_G1_likelihoods <- c(G1_G1model_mlogliks,G2_G1model_mlogliks)
	samples_G2_likelihoods <- c(G1_G2model_mlogliks,G2_G2model_mlogliks)
	# T posterior probability calculation
	out <- cbind(c(G1,G2),exp(-samples_G1_likelihoods) / (exp(-samples_G1_likelihoods) + exp(-samples_G2_likelihoods)),samples_G1_likelihoods,samples_G2_likelihoods)
	colnames(out) <- c("sample_ID","posterior_G1","G1_mloglik","G2_mloglik")
	eval(parse(text=paste('write.table(as.data.frame(out), col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE, file="./',i,'.predicted")',sep="")))
	
	cat(paste("done predicting for ",i," in ", sprintf("%.2f", (proc.time()[3]-ptm)/60)," minutes\n",sep=""))
}
