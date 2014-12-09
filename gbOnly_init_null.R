args <- commandArgs(trailingOnly = TRUE)
beg <- as.numeric(args[1])
end <- as.numeric(args[2])

load("./essentials_allBRCA.RData")
IDs_length <- nchar(G1[1])+2
G1 <- G1[1:55]
G2 <- G2[1:486]

res_gb <- 25
smooth <- 1
if (smooth > 0) library(aws)

integrand_m <- function(x,mean) {dnorm(x=mean,mean=x,sd=0.14)}

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
	
	IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',workingList_BRCA[i],'"',sep="")))))
	promoterVars <- geneBodyVars_template[1:length(IDs_body)]
	promoter_CpGs <- template_body_CpGs[1:length(IDs_body)]
	# generate "missing" Var data
	ncol = length(IDs_body)
	tempVar_G1 <- matrix(rep(".",length(G1)*1),nrow=length(G1),ncol=1)
	colnames(tempVar_G1) <- "NAME:\tM.GB"
	rownames(tempVar_G1) <- G1
	eval(parse(text = paste('write.table(', paste('tempVar_G1,file = "./',i,'/G1_model/all/G1_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))  
	tempVar_G2 <- matrix(rep(".",length(G2)*1),nrow=length(G2),ncol=1)
	colnames(tempVar_G2) <- "NAME:\tM.GB"
	rownames(tempVar_G2) <- G2
	eval(parse(text = paste('write.table(', paste('tempVar_G2,file = "./',i,'/G2_model/all/G2_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	tempVar_all <- rbind(tempVar_G2,tempVar_G1)
	rownames(tempVar_all) <- c(G2,G1)
	eval(parse(text = paste('write.table(', paste('tempVar_all,file = "./',i,'/full_model/full_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	############################################################
	######### calculate epsilons and smoothing params ##########
	epsilon_pr_t <- 1/(length(G2)*ncol)/res_gb
	epsilon_pr_an <- 1/(length(G1)*ncol)/res_gb
	
	smooth_gb <- 10/(mean(length(G1),length(G2))*length(IDs_body)/(res_gb*2))
	###########################################################
	############### binning scheme defined here ###############
	# promoter
	x <- as.vector(mmatrix_pc[IDs_body,c(G2,G1)])
	
	density <- density(x,bw=0.14,from=-7,to=7,n=2801)
	density$y <- density$y/sum(density$y)
	density$y <- cumsum(density$y)
	breaks <- NULL
	noBreaks <- res_gb-1
	for (j in 1:noBreaks) { breaks <- c(breaks, density$x[which(density$y >= j*(1/(1+noBreaks)))][1])}
	breaksPROMOTER <- sort(c(-7.01,breaks,7.01))
	
	all_labels_pr <- as.character(seq(1,res_gb,1))
	promoter_t <- matrix(ncol=res_gb,nrow=length(G2))
	promoter_an <- matrix(ncol=res_gb,nrow=length(G1))
	########################################################
	# dynamic generation of model specification files here #
	########################################################
	stateMaps <- file(paste("./",i,"/stateMaps.txt",sep=""),"w")
	prMap <- paste("NAME:\tprMap\nSYMBOLS:\t",paste(seq(1,res_gb,1),collapse=" "),"\nMETA_SYMBOLS:\t.=",paste(seq(1,res_gb,1),collapse=" "),"; *=",paste(seq(1,res_gb,1),collapse=" "),";\n\n",collapse="",sep="")
	cat(prMap,file=stateMaps)
	close(stateMaps)
	########################################################
	####################### variables ######################
	variables <- file(paste("./",i,"/variables.txt",sep=""),"w")
	cat("STATE_MAP_NAME:\tprMap\nVAR_NAMES:\tM.GB\n\n",sep="",file=variables)
	close(variables)
	########################################################
	##################### factor graph #####################
	factorGraph <- file(paste("./",i,"/factorGraph.txt",sep=""),"w")
	cat(paste("\nNAME:\tprior.M.GB\nNB1:\tM.GB\nPOT:\tpot_P.M.prior\n",sep="",collapse=""),file=factorGraph)
	cat(paste("\nNAME:\t",promoter_CpGs,"\nNB1:\tM.GB\nPOT:\tpot_",promoterVars,"\n",sep="",collapse=""),file=factorGraph)
	close(factorGraph)
	
	system(command=paste('cp ./',i,'/*.txt ./',i,'/G1_model/all',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/G2_model/all',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/full_model',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/null/G1_model',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/null/G2_model',sep=""))
	
	promoter_CpGs[1] <- paste("NAME:\t",promoter_CpGs[1],sep="")
	###########################################################################
	############################## G1 model ###################################
	## full G1 model developed from here, to obtain likelihoods of G1 ########
	
	# generate FacData for full set of G1
	tempS_G1 <- matrix(ncol=ncol,nrow=length(G1))
	for (current_sample in 1:length(G1)) {
		# promoter
		cpg_list_pr <- NULL
		for (cpg in 1:length(promoterVars)) {
			miu <- mmatrix_pc[IDs_body[cpg],G1[current_sample]]
			frequencies_pr <- rep(0,res_gb)
			for (freq in 1:res_gb) {
				frequencies_pr[freq] <- integrate(integrand_m,lower=breaksPROMOTER[freq],upper=breaksPROMOTER[freq+1],mean=miu)$value
			}
			frequencies_pr <- unlist(frequencies_pr) + epsilon_pr_an
			frequencies_pr <- frequencies_pr/sum(frequencies_pr)
			cpg_list_pr[[cpg]] <- frequencies_pr
		}
		
		tempS_formated <- matrix(ncol=ncol,nrow=1)
		cur_ncol <- 0
		for (element in 1:length(cpg_list_pr)) {
			tempS_formated[1,cur_ncol+1] <- paste('[1,',res_gb,']((',paste(cpg_list_pr[[element]],sep="",collapse=","),'))',sep="")
			cur_ncol <- cur_ncol + 1
		}
		tempS_G1[current_sample,] <- tempS_formated
		#start precomputing correct initialization of parameters
		promoter_an[current_sample,] <- apply(matrix(unlist(cpg_list_pr),ncol=res_gb,byrow=TRUE),2,geo_mean)/sum(apply(matrix(unlist(cpg_list_pr),ncol=res_gb,byrow=TRUE),2,geo_mean))
	}
	# precompute correct initialization of parameters for AN-only model
	if (smooth > 0) {
		prior_pr <- kernsm(apply(promoter_an,2,mean),h=smooth_gb)
		prior_pr <- prior_pr@yhat/sum(prior_pr@yhat)
	} else prior_pr <- apply(promoter_an,2,mean)
	
	# write potentials file
	string <- paste(prior_pr,collapse=",")
	promoterPots <- paste("\nNAME:\t\tpot_",c("P.M.prior",promoterVars),"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="")
	potentials <- file(paste("./",i,"/G1_model/all/factorPotentials.txt",sep=""),"w")
	cat(promoterPots,file=potentials)
	close(potentials)
	
	tempFac <- tempS_G1
	colnames(tempFac) <- promoter_CpGs[1:length(promoterVars)]
	rownames(tempFac) <- G1
	eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/G1_model/all/G1_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	
	# query the full model with AN samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/G1_model/all/ -l -n - ./',i,'/G1_model/all/G1_VarData.tab ./',i,'/G1_model/all/G1_FacData.tab',sep=""))
	G1_G1_likelihoods <- as.numeric(substring(string[-1],17))
	##########################################################################
	############################## G2 model ###################################
	### full G2 model developed from here, to obtain likelihoods of G2 #######
	tempS_G2 <- matrix(ncol=ncol,nrow=length(G2))
	for (current_sample in 1:length(G2)) {
		# promoter
		cpg_list_pr <- NULL
		for (cpg in 1:length(promoterVars)) {
			miu <- mmatrix_pc[IDs_body[cpg],G2[current_sample]]
			frequencies_pr <- rep(0,res_gb)
			for (freq in 1:res_gb) {
				frequencies_pr[freq] <- integrate(integrand_m,lower=breaksPROMOTER[freq],upper=breaksPROMOTER[freq+1],mean=miu)$value
			}
			frequencies_pr <- unlist(frequencies_pr) + epsilon_pr_t
			frequencies_pr <- frequencies_pr/sum(frequencies_pr)
			cpg_list_pr[[cpg]] <- frequencies_pr
		}
		
		tempS_formated <- matrix(ncol=ncol,nrow=1)
		cur_ncol <- 0
		for (element in 1:length(cpg_list_pr)) {
			tempS_formated[1,cur_ncol+1] <- paste('[1,',res_gb,']((',paste(cpg_list_pr[[element]],sep="",collapse=","),'))',sep="")
			cur_ncol <- cur_ncol + 1
		}
		tempS_G2[current_sample,] <- tempS_formated
		#start precomputing correct initialization of parameters
		promoter_t[current_sample,] <- apply(matrix(unlist(cpg_list_pr),ncol=res_gb,byrow=TRUE),2,geo_mean)/sum(apply(matrix(unlist(cpg_list_pr),ncol=res_gb,byrow=TRUE),2,geo_mean))
	}
	# precompute correct initialization of parameters for T-only model
	if (smooth > 0) {
		prior_pr <- kernsm(apply(promoter_t,2,mean),h=smooth_gb)
		prior_pr <- prior_pr@yhat/sum(prior_pr@yhat)
	} else prior_pr <- apply(promoter_t,2,mean)
	
	# write potentials file
	string <- paste(prior_pr,collapse=",")
	promoterPots <- paste("\nNAME:\t\tpot_",c("P.M.prior",promoterVars),"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="")
	potentials <- file(paste("./",i,"/G2_model/all/factorPotentials.txt",sep=""),"w")
	cat(promoterPots,file=potentials)
	close(potentials)
	
	tempFac <- tempS_G2
	colnames(tempFac) <- promoter_CpGs[1:length(promoterVars)]
	rownames(tempFac) <- G2
	eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/G2_model/all/G2_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	
	# query the full model with T samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/G2_model/all/ -l -n - ./',i,'/G2_model/all/G2_VarData.tab ./',i,'/G2_model/all/G2_FacData.tab',sep=""))
	G2_G2_likelihoods <- as.numeric(substring(string[-1],17))
	###########################################################################
	######################## All data model ###################################
	## Full model developed from here, to obtain likelihoods of G2 and G1 ####
	# precompute correct initialization of parameters for full model
	promoter_all <- rbind(promoter_t,promoter_an)
	
	if (smooth > 0) {
		prior_pr <- kernsm(apply(promoter_all,2,mean),h=smooth_gb)
		prior_pr <- prior_pr@yhat/sum(prior_pr@yhat)
	} else prior_pr <- apply(promoter_all,2,mean)
	
	# write potentials file
	string <- paste(prior_pr,collapse=",")
	promoterPots <- paste("\nNAME:\t\tpot_",c("P.M.prior",promoterVars),"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="")
	potentials <- file(paste("./",i,"/full_model/factorPotentials.txt",sep=""),"w")
	cat(promoterPots,file=potentials)
	close(potentials)

	tempFac <- rbind(tempS_G2,tempS_G1)
	colnames(tempFac) <- promoter_CpGs[1:length(promoterVars)]
	rownames(tempFac) <- c(G2,G1)
	eval(parse(text = paste('write.table(', paste('tempFac,file = "./',i,'/full_model/full_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	
	# build and query the full model with G2 and G1 samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/full_model/ -l -n - ./',i,'/full_model/full_VarData.tab ./',i,'/full_model/full_FacData.tab',sep=""))
	allData_full_likelihoods <- as.numeric(substring(string[-1],17))
	###########################################################################
	######################## D calculation ###################################
	D <- 2*(sum(allData_full_likelihoods) - (sum(G1_G1_likelihoods)+sum(G2_G2_likelihoods)))
	###########################################################################
	################# P val calculation using null distr. #####################
	nruns <- 100
	Ds <- vector(length=nruns,mode="numeric")
	for (run in 1:nruns) {
		cur <- sample(x=1:(length(G2)+length(G1)),size=length(G2),replace=FALSE)
		
		# G2
		tempFac_G2 <- as.data.frame(tempFac[cur,])
		colnames(tempFac_G2) <- promoter_CpGs[1:length(promoterVars)]
		rownames(tempFac_G2) <- G2
		eval(parse(text = paste('write.table(', paste('tempFac_G2,file = "./',i,'/null/G2_model/G2_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		if (smooth > 0) {
			prior_pr <- kernsm(apply(promoter_all[cur,],2,mean),h=smooth_gb)
			prior_pr <- prior_pr@yhat/sum(prior_pr@yhat)
		} else prior_pr <- apply(promoter_all[cur,],2,mean)
		
		string <- paste(prior_pr,collapse=",")
		promoterPots <- paste("\nNAME:\t\tpot_",c("P.M.prior",promoterVars),"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="")
		potentials <- file(paste("./",i,"/null/G2_model/factorPotentials.txt",sep=""),"w")
		cat(promoterPots,file=potentials)
		close(potentials)
		# query
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/null/G2_model/ -l -n - ./',i,'/G2_model/all/G2_VarData.tab ./',i,'/null/G2_model/G2_FacData.tab',sep=""))
		G2_G2_likelihoods <- as.numeric(substring(string[-1],17))
		
		# G1
		tempFac_G1 <- as.data.frame(tempFac[-cur,])
		colnames(tempFac_G1) <- promoter_CpGs[1:length(promoterVars)]
		rownames(tempFac_G1) <- G1
		eval(parse(text = paste('write.table(', paste('tempFac_G1,file = "./',i,'/null/G1_model/G1_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		if (smooth > 0) {
			prior_pr <- kernsm(apply(promoter_all[-cur,],2,mean),h=smooth_gb)
			prior_pr <- prior_pr@yhat/sum(prior_pr@yhat)
		} else prior_pr <- apply(promoter_all[-cur,],2,mean)
		
		string <- paste(prior_pr,collapse=",")
		promoterPots <- paste("\nNAME:\t\tpot_",c("P.M.prior",promoterVars),"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_gb,"]((",string,"))\nPC_MAT:\t\t[1,",res_gb,"]((",paste(rep(1,res_gb),collapse=","),"))\n",sep="")
		potentials <- file(paste("./",i,"/null/G1_model/factorPotentials.txt",sep=""),"w")
		cat(promoterPots,file=potentials)
		close(potentials)
		# query
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/null/G1_model/ -l -n - ./',i,'/G1_model/all/G1_VarData.tab ./',i,'/null/G1_model/G1_FacData.tab',sep=""))
		G1_G1_likelihoods <- as.numeric(substring(string[-1],17))
		# Ds calculation
		Ds[run] <- 2*(sum(allData_full_likelihoods) - (sum(G1_G1_likelihoods)+sum(G2_G2_likelihoods)))
	}
	
	if (sd(Ds) != 0) zscore <- (D - mean(Ds)) / sd(Ds) else zscore <- -6
	pval_zscore <- pnorm(zscore,lower.tail=FALSE)
	###########################################################################################
	eval(parse(text=paste('write.table(x=t(c(pval_zscore,D,mean(Ds),sd(Ds),zscore)), col.names=FALSE, row.names=FALSE, append=TRUE, file="./',i,'.result")',sep="")))
	cat(paste("done ",i," in ", sprintf("%.2f", (proc.time()[3]-ptm)/60)," minutes\n",sep=""))
}
