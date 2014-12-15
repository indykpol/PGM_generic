args <- commandArgs(trailingOnly = TRUE)
beg <- as.numeric(args[1])
end <- as.numeric(args[2])

load("./essentials_allBRCA.RData") # data to work on
IDs_length <- nchar(G1[1])+2
G1 <- G1[c(1,2,5,6,7,8,10,14,18,20,21,22,23,24,25,26,27,28,29,30,32,33,34,35,37,38,39,41,42,45,46,47,48,49,51,53,54,56,57,58,59,61,62,64,65,67,68,70,71,75,76,77,78,80,81)]
G2 <- G2[c(1,2,3,5,7,9,10,12,13,14,15,16,19,20,21,22,23,24,25,26,27,30,31,32,34,35,36,37,38,39,41,42,44,46,47,48,52,53,54,55,58,59,61,62,64,65,66,67,69,72,74,76,77,79,81,82,83,87,88,89,90,91,92,93,94,95,97,98,99,101,102,103,104,106,108,112,113,115,116,117,119,120,121,122,125,126,128,129,130,131,132,133,135,136,137,139,141,142,143,145,146,147,148,149,152,155,156,157,158,159,160,164,165,166,168,169,170,172,173,174,176,178,180,181,182,183,184,187,189,190,191,192,193,194,195,197,198,201,202,203,204,205,208,210,212,214,215,217,219,220,221,223,224,225,228,229,232,237,238,239,240,241,244,245,246,247,248,250,251,252,253,254,255,256,259,260,263,264,265,266,268,269,270,272,277,278,279,280,282,283,284,287,289,290,291,292,295,296,297,298,299,300,301,302,303,305,308,309,310,311,313,314,316,317,318,319,320,321,322,324,326,330,332,333,334,335,336,337,339,340,342,345,349,351,352,353,354,355,356,357,358,360,361,362,363,364,365,367,368,369,370,371,372,373,376,377,379,380,382,383,387,390,391,393,397,399,400,402,403,404,405,406,407,409,411,412,413,414,415,417,418,420,422,423,424,426,427,428,429,431,432,433,435,436,437,438,439,442,445,447,448,450,451,452,453,454,455,456,459,461,463,467,468,471,472,473,474,475,477,479,481,482,484,485,486,488,489,490,491,493,494,495,497,500,502,503,505,507,509,510,511,512,513,514,517,518,519,520,521,523,524,525,527,528,529,530,531,532,534,535,536,537,539,540,541,545,546,551,552,553,554,555,556,557,559,561,567,569,570,575,579,580,582,583,584,585,586,589,591,594,596,598,600,602,603,604,605,606,609,610,611,612,613,614,616,618,619,620,621,622,623,625,628,630,631,632,633,634,635,636,638,639,640,641,642,643,645,646,647,648,649,650,652,653,654,658,659,660,661,662,663,664,666,667,669,671,672,673,674,676,677,679,680,681,682,684,685,686,688,689,690,691,696,698,702,703,705,706,707,708,709,711,712,713,714,715,716,718,720,721,722,723,724,726,728,729,730)]

res_pr <- 25
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
	
	IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',workingList_BRCA[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',workingList_BRCA[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',workingList_BRCA[i],'"',sep=""))))
	promoterVars <- promoterVars_template[1:length(IDs_promoter)]
	promoter_CpGs <- template_promoter_CpGs[1:length(IDs_promoter)]
	ncol = length(IDs_promoter)
	# generate "missing" Var data
	
	tempVar_G1 <- matrix(rep(".",length(G1)*1),nrow=length(G1),ncol=1)
	colnames(tempVar_G1) <- "NAME:\tM.P"
	rownames(tempVar_G1) <- G1
	eval(parse(text = paste('write.table(', paste('tempVar_G1,file = "./',i,'/G1_model/all/AN_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))  
	tempVar_G2 <- matrix(rep(".",length(G2)*1),nrow=length(G2),ncol=1)
	colnames(tempVar_G2) <- "NAME:\tM.P"
	rownames(tempVar_G2) <- G2
	eval(parse(text = paste('write.table(', paste('tempVar_G2,file = "./',i,'/G2_model/all/T_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	tempVar_all <- rbind(tempVar_G2,tempVar_G1)
	rownames(tempVar_all) <- c(G2,G1)
	eval(parse(text = paste('write.table(', paste('tempVar_all,file = "./',i,'/full_model/full_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	############################################################
	######### calculate epsilons and smoothing params ##########
	epsilon_pr_t <- 1/(length(G2)*ncol)/res_pr
	epsilon_pr_an <- 1/(length(G1)*ncol)/res_pr
	
	smooth_pr <- 10/(mean(length(G1),length(G2))*length(IDs_promoter)/(res_pr*2))
	###########################################################
	############### binning scheme defined here ###############
	# promoter
	x <- mmatrix_pc[IDs_promoter,c(G2,G1)]
	
	density <- density(x,bw=0.14,from=-7,to=7,n=2801)
	density$y <- density$y/sum(density$y)
	density$y <- cumsum(density$y)
	breaks <- NULL
	noBreaks <- res_pr-1
	for (j in 1:noBreaks) { breaks <- c(breaks, density$x[which(density$y >= j*(1/(1+noBreaks)))][1])}
	breaksPROMOTER <- sort(c(-7.01,breaks,7.01))
	
	all_labels_pr <- as.character(seq(1,res_pr,1))
	promoter_t <- matrix(ncol=res_pr,nrow=length(G2))
	promoter_an <- matrix(ncol=res_pr,nrow=length(G1))
	########################################################
	# dynamic generation of model specification files here #
	########################################################
	stateMaps <- file(paste("./",i,"/stateMaps.txt",sep=""),"w")
	prMap <- paste("NAME:\tprMap\nSYMBOLS:\t",paste(seq(1,res_pr,1),collapse=" "),"\nMETA_SYMBOLS:\t.=",paste(seq(1,res_pr,1),collapse=" "),"; *=",paste(seq(1,res_pr,1),collapse=" "),";\n\n",collapse="",sep="")
	cat(prMap,file=stateMaps)
	close(stateMaps)
	########################################################
	####################### variables ######################
	variables <- file(paste("./",i,"/variables.txt",sep=""),"w")
	cat("STATE_MAP_NAME:\tprMap\nVAR_NAMES:\tM.P\n\n",sep="",file=variables)
	close(variables)
	########################################################
	##################### factor graph #####################
	factorGraph <- file(paste("./",i,"/factorGraph.txt",sep=""),"w")
	cat(paste("\nNAME:\tprior.M.P\nNB1:\tM.P\nPOT:\tpot_P.M.prior\n",sep="",collapse=""),file=factorGraph)
	cat(paste("\nNAME:\t",promoter_CpGs,"\nNB1:\tM.P\nPOT:\tpot_",promoterVars,"\n",sep="",collapse=""),file=factorGraph)
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
			miu <- mmatrix_pc[IDs_promoter[cpg],G1[current_sample]]
			frequencies_pr <- rep(0,res_pr)
			for (freq in 1:res_pr) {
				frequencies_pr[freq] <- integrate(integrand_m,lower=breaksPROMOTER[freq],upper=breaksPROMOTER[freq+1],mean=miu)$value
			}
			frequencies_pr <- unlist(frequencies_pr) + epsilon_pr_an
			frequencies_pr <- frequencies_pr/sum(frequencies_pr)
			cpg_list_pr[[cpg]] <- frequencies_pr
		}
		
		tempS_formated <- matrix(ncol=ncol,nrow=1)
		cur_ncol <- 0
		for (element in 1:length(cpg_list_pr)) {
			tempS_formated[1,cur_ncol+1] <- paste('[1,',res_pr,']((',paste(cpg_list_pr[[element]],sep="",collapse=","),'))',sep="")
			cur_ncol <- cur_ncol + 1
		}
		tempS_G1[current_sample,] <- tempS_formated
		#start precomputing correct initialization of parameters
		promoter_an[current_sample,] <- apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean)/sum(apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean))
	}
	# precompute correct initialization of parameters for AN-only model
	if (smooth > 0) {
		prior_pr <- kernsm(apply(promoter_an,2,mean),h=smooth_pr)
		prior_pr <- prior_pr@yhat/sum(prior_pr@yhat)
	} else prior_pr <- apply(promoter_an,2,mean)
	
	# write potentials file
	string <- paste(prior_pr,collapse=",")
	promoterPots <- paste("\nNAME:\t\tpot_",c("P.M.prior",promoterVars),"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="")
	potentials <- file(paste("./",i,"/G1_model/all/factorPotentials.txt",sep=""),"w")
	cat(promoterPots,file=potentials)
	close(potentials)
	
	tempFac <- tempS_G1
	colnames(tempFac) <- promoter_CpGs[1:length(promoterVars)]
	rownames(tempFac) <- G1
	eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/G1_model/all/G1_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	
	# query the full model with G1 samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/G1_model/all/ -l -n - ./',i,'/G1_model/all/AN_VarData.tab ./',i,'/G1_model/all/G1_FacData.tab',sep=""))
	G1_G1_likelihoods <- as.numeric(substring(string[-1],17))
	##########################################################################
	############################## G2 model ###################################
	### full G2 model developed from here, to obtain likelihoods of G2 #######
	tempS_G2 <- matrix(ncol=ncol,nrow=length(G2))
	for (current_sample in 1:length(G2)) {
		# promoter
		cpg_list_pr <- NULL
		for (cpg in 1:length(promoterVars)) {
			miu <- mmatrix_pc[IDs_promoter[cpg],G2[current_sample]]
			frequencies_pr <- rep(0,res_pr)
			for (freq in 1:res_pr) {
				frequencies_pr[freq] <- integrate(integrand_m,lower=breaksPROMOTER[freq],upper=breaksPROMOTER[freq+1],mean=miu)$value
			}
			frequencies_pr <- unlist(frequencies_pr) + epsilon_pr_t
			frequencies_pr <- frequencies_pr/sum(frequencies_pr)
			cpg_list_pr[[cpg]] <- frequencies_pr
		}
		
		tempS_formated <- matrix(ncol=ncol,nrow=1)
		cur_ncol <- 0
		for (element in 1:length(cpg_list_pr)) {
			tempS_formated[1,cur_ncol+1] <- paste('[1,',res_pr,']((',paste(cpg_list_pr[[element]],sep="",collapse=","),'))',sep="")
			cur_ncol <- cur_ncol + 1
		}
		tempS_G2[current_sample,] <- tempS_formated
		#start precomputing correct initialization of parameters
		promoter_t[current_sample,] <- apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean)/sum(apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean))
	}
	# precompute correct initialization of parameters for T-only model
	if (smooth > 0) {
		prior_pr <- kernsm(apply(promoter_t,2,mean),h=smooth_pr)
		prior_pr <- prior_pr@yhat/sum(prior_pr@yhat)
	} else prior_pr <- apply(promoter_t,2,mean)
	
	# write potentials file
	string <- paste(prior_pr,collapse=",")
	promoterPots <- paste("\nNAME:\t\tpot_",c("P.M.prior",promoterVars),"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="")
	potentials <- file(paste("./",i,"/G2_model/all/factorPotentials.txt",sep=""),"w")
	cat(promoterPots,file=potentials)
	close(potentials)
	
	tempFac <- tempS_G2
	colnames(tempFac) <- promoter_CpGs[1:length(promoterVars)]
	rownames(tempFac) <- G2
	eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/G2_model/all/G2_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	
	# query the full model with G2 samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/G2_model/all/ -l -n - ./',i,'/G2_model/all/T_VarData.tab ./',i,'/G2_model/all/G2_FacData.tab',sep=""))
	G2_G2_likelihoods <- as.numeric(substring(string[-1],17))
	###########################################################################
	######################## All data model ###################################
	## Full model developed from here, to obtain likelihoods of G2 and G1 ####
	# precompute correct initialization of parameters for full model
	promoter_all <- rbind(promoter_t,promoter_an)
	
	if (smooth > 0) {
		prior_pr <- kernsm(apply(promoter_all,2,mean),h=smooth_pr)
		prior_pr <- prior_pr@yhat/sum(prior_pr@yhat)
	} else prior_pr <- apply(promoter_all,2,mean)
	
	# write potentials file
	string <- paste(prior_pr,collapse=",")
	promoterPots <- paste("\nNAME:\t\tpot_",c("P.M.prior",promoterVars),"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="")
	potentials <- file(paste("./",i,"/full_model/factorPotentials.txt",sep=""),"w")
	cat(promoterPots,file=potentials)
	close(potentials)

	tempFac <- rbind(tempS_G2,tempS_G1)
	colnames(tempFac) <- promoter_CpGs[1:length(promoterVars)]
	rownames(tempFac) <- c(G2,G1)
	eval(parse(text = paste('write.table(', paste('tempFac,file = "./',i,'/full_model/full_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	
	# build and query the full model with T and AN samples
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
			prior_pr <- kernsm(apply(promoter_all[cur,],2,mean),h=smooth_pr)
			prior_pr <- prior_pr@yhat/sum(prior_pr@yhat)
		} else prior_pr <- apply(promoter_all[cur,],2,mean)
		
		string <- paste(prior_pr,collapse=",")
		promoterPots <- paste("\nNAME:\t\tpot_",c("P.M.prior",promoterVars),"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="")
		potentials <- file(paste("./",i,"/null/G2_model/factorPotentials.txt",sep=""),"w")
		cat(promoterPots,file=potentials)
		close(potentials)
		# query
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/null/G2_model/ -l -n - ./',i,'/G2_model/all/T_VarData.tab ./',i,'/null/G2_model/G2_FacData.tab',sep=""))
		G2_G2_likelihoods <- as.numeric(substring(string[-1],17))
		
		# G1
		tempFac_G1 <- as.data.frame(tempFac[-cur,])
		colnames(tempFac_G1) <- promoter_CpGs[1:length(promoterVars)]
		rownames(tempFac_G1) <- G1
		eval(parse(text = paste('write.table(', paste('tempFac_G1,file = "./',i,'/null/G1_model/G1_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		if (smooth > 0) {
			prior_pr <- kernsm(apply(promoter_all[-cur,],2,mean),h=smooth_pr)
			prior_pr <- prior_pr@yhat/sum(prior_pr@yhat)
		} else prior_pr <- apply(promoter_all[-cur,],2,mean)
		
		string <- paste(prior_pr,collapse=",")
		promoterPots <- paste("\nNAME:\t\tpot_",c("P.M.prior",promoterVars),"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="")
		potentials <- file(paste("./",i,"/null/G1_model/factorPotentials.txt",sep=""),"w")
		cat(promoterPots,file=potentials)
		close(potentials)
		# query
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/null/G1_model/ -l -n - ./',i,'/G1_model/all/AN_VarData.tab ./',i,'/null/G1_model/G1_FacData.tab',sep=""))
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
args <- commandArgs(trailingOnly = TRUE)
beg <- as.numeric(args[1])
end <- as.numeric(args[2])

load("./essentials_allBRCA.RData") # data to work on
IDs_length <- nchar(G1[1])+2
G1 <- G1[1:55]
G2 <- G2[1:486]

res_pr <- 25
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
	
	IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',workingList_BRCA[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',workingList_BRCA[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',workingList_BRCA[i],'"',sep=""))))
	promoterVars <- promoterVars_template[1:length(IDs_promoter)]
	promoter_CpGs <- template_promoter_CpGs[1:length(IDs_promoter)]
	ncol = length(IDs_promoter)
	# generate "missing" Var data
	
	tempVar_G1 <- matrix(rep(".",length(G1)*1),nrow=length(G1),ncol=1)
	colnames(tempVar_G1) <- "NAME:\tM.P"
	rownames(tempVar_G1) <- G1
	eval(parse(text = paste('write.table(', paste('tempVar_G1,file = "./',i,'/G1_model/all/AN_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))  
	tempVar_G2 <- matrix(rep(".",length(G2)*1),nrow=length(G2),ncol=1)
	colnames(tempVar_G2) <- "NAME:\tM.P"
	rownames(tempVar_G2) <- G2
	eval(parse(text = paste('write.table(', paste('tempVar_G2,file = "./',i,'/G2_model/all/T_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	tempVar_all <- rbind(tempVar_G2,tempVar_G1)
	rownames(tempVar_all) <- c(G2,G1)
	eval(parse(text = paste('write.table(', paste('tempVar_all,file = "./',i,'/full_model/full_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	############################################################
	######### calculate epsilons and smoothing params ##########
	epsilon_pr_t <- 1/(length(G2)*ncol)/res_pr
	epsilon_pr_an <- 1/(length(G1)*ncol)/res_pr
	
	smooth_pr <- 10/(mean(length(G1),length(G2))*length(IDs_promoter)/(res_pr*2))
	###########################################################
	############### binning scheme defined here ###############
	# promoter
	x <- mmatrix_pc[IDs_promoter,c(G2,G1)]
	
	density <- density(x,bw=0.14,from=-7,to=7,n=2801)
	density$y <- density$y/sum(density$y)
	density$y <- cumsum(density$y)
	breaks <- NULL
	noBreaks <- res_pr-1
	for (j in 1:noBreaks) { breaks <- c(breaks, density$x[which(density$y >= j*(1/(1+noBreaks)))][1])}
	breaksPROMOTER <- sort(c(-7.01,breaks,7.01))
	
	all_labels_pr <- as.character(seq(1,res_pr,1))
	promoter_t <- matrix(ncol=res_pr,nrow=length(G2))
	promoter_an <- matrix(ncol=res_pr,nrow=length(G1))
	########################################################
	# dynamic generation of model specification files here #
	########################################################
	stateMaps <- file(paste("./",i,"/stateMaps.txt",sep=""),"w")
	prMap <- paste("NAME:\tprMap\nSYMBOLS:\t",paste(seq(1,res_pr,1),collapse=" "),"\nMETA_SYMBOLS:\t.=",paste(seq(1,res_pr,1),collapse=" "),"; *=",paste(seq(1,res_pr,1),collapse=" "),";\n\n",collapse="",sep="")
	cat(prMap,file=stateMaps)
	close(stateMaps)
	########################################################
	####################### variables ######################
	variables <- file(paste("./",i,"/variables.txt",sep=""),"w")
	cat("STATE_MAP_NAME:\tprMap\nVAR_NAMES:\tM.P\n\n",sep="",file=variables)
	close(variables)
	########################################################
	##################### factor graph #####################
	factorGraph <- file(paste("./",i,"/factorGraph.txt",sep=""),"w")
	cat(paste("\nNAME:\tprior.M.P\nNB1:\tM.P\nPOT:\tpot_P.M.prior\n",sep="",collapse=""),file=factorGraph)
	cat(paste("\nNAME:\t",promoter_CpGs,"\nNB1:\tM.P\nPOT:\tpot_",promoterVars,"\n",sep="",collapse=""),file=factorGraph)
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
			miu <- mmatrix_pc[IDs_promoter[cpg],G1[current_sample]]
			frequencies_pr <- rep(0,res_pr)
			for (freq in 1:res_pr) {
				frequencies_pr[freq] <- integrate(integrand_m,lower=breaksPROMOTER[freq],upper=breaksPROMOTER[freq+1],mean=miu)$value
			}
			frequencies_pr <- unlist(frequencies_pr) + epsilon_pr_an
			frequencies_pr <- frequencies_pr/sum(frequencies_pr)
			cpg_list_pr[[cpg]] <- frequencies_pr
		}
		
		tempS_formated <- matrix(ncol=ncol,nrow=1)
		cur_ncol <- 0
		for (element in 1:length(cpg_list_pr)) {
			tempS_formated[1,cur_ncol+1] <- paste('[1,',res_pr,']((',paste(cpg_list_pr[[element]],sep="",collapse=","),'))',sep="")
			cur_ncol <- cur_ncol + 1
		}
		tempS_G1[current_sample,] <- tempS_formated
		#start precomputing correct initialization of parameters
		promoter_an[current_sample,] <- apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean)/sum(apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean))
	}
	# precompute correct initialization of parameters for AN-only model
	if (smooth > 0) {
		prior_pr <- kernsm(apply(promoter_an,2,mean),h=smooth_pr)
		prior_pr <- prior_pr@yhat/sum(prior_pr@yhat)
	} else prior_pr <- apply(promoter_an,2,mean)
	
	# write potentials file
	string <- paste(prior_pr,collapse=",")
	promoterPots <- paste("\nNAME:\t\tpot_",c("P.M.prior",promoterVars),"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="")
	potentials <- file(paste("./",i,"/G1_model/all/factorPotentials.txt",sep=""),"w")
	cat(promoterPots,file=potentials)
	close(potentials)
	
	tempFac <- tempS_G1
	colnames(tempFac) <- promoter_CpGs[1:length(promoterVars)]
	rownames(tempFac) <- G1
	eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/G1_model/all/G1_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	
	# query the full model with G1 samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/G1_model/all/ -l -n - ./',i,'/G1_model/all/AN_VarData.tab ./',i,'/G1_model/all/G1_FacData.tab',sep=""))
	G1_G1_likelihoods <- as.numeric(substring(string[-1],17))
	##########################################################################
	############################## G2 model ###################################
	### full G2 model developed from here, to obtain likelihoods of G2 #######
	tempS_G2 <- matrix(ncol=ncol,nrow=length(G2))
	for (current_sample in 1:length(G2)) {
		# promoter
		cpg_list_pr <- NULL
		for (cpg in 1:length(promoterVars)) {
			miu <- mmatrix_pc[IDs_promoter[cpg],G2[current_sample]]
			frequencies_pr <- rep(0,res_pr)
			for (freq in 1:res_pr) {
				frequencies_pr[freq] <- integrate(integrand_m,lower=breaksPROMOTER[freq],upper=breaksPROMOTER[freq+1],mean=miu)$value
			}
			frequencies_pr <- unlist(frequencies_pr) + epsilon_pr_t
			frequencies_pr <- frequencies_pr/sum(frequencies_pr)
			cpg_list_pr[[cpg]] <- frequencies_pr
		}
		
		tempS_formated <- matrix(ncol=ncol,nrow=1)
		cur_ncol <- 0
		for (element in 1:length(cpg_list_pr)) {
			tempS_formated[1,cur_ncol+1] <- paste('[1,',res_pr,']((',paste(cpg_list_pr[[element]],sep="",collapse=","),'))',sep="")
			cur_ncol <- cur_ncol + 1
		}
		tempS_G2[current_sample,] <- tempS_formated
		#start precomputing correct initialization of parameters
		promoter_t[current_sample,] <- apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean)/sum(apply(matrix(unlist(cpg_list_pr),ncol=res_pr,byrow=TRUE),2,geo_mean))
	}
	# precompute correct initialization of parameters for T-only model
	if (smooth > 0) {
		prior_pr <- kernsm(apply(promoter_t,2,mean),h=smooth_pr)
		prior_pr <- prior_pr@yhat/sum(prior_pr@yhat)
	} else prior_pr <- apply(promoter_t,2,mean)
	
	# write potentials file
	string <- paste(prior_pr,collapse=",")
	promoterPots <- paste("\nNAME:\t\tpot_",c("P.M.prior",promoterVars),"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="")
	potentials <- file(paste("./",i,"/G2_model/all/factorPotentials.txt",sep=""),"w")
	cat(promoterPots,file=potentials)
	close(potentials)
	
	tempFac <- tempS_G2
	colnames(tempFac) <- promoter_CpGs[1:length(promoterVars)]
	rownames(tempFac) <- G2
	eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/G2_model/all/G2_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	
	# query the full model with G2 samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/G2_model/all/ -l -n - ./',i,'/G2_model/all/T_VarData.tab ./',i,'/G2_model/all/G2_FacData.tab',sep=""))
	G2_G2_likelihoods <- as.numeric(substring(string[-1],17))
	###########################################################################
	######################## All data model ###################################
	## Full model developed from here, to obtain likelihoods of G2 and G1 ####
	# precompute correct initialization of parameters for full model
	promoter_all <- rbind(promoter_t,promoter_an)
	
	if (smooth > 0) {
		prior_pr <- kernsm(apply(promoter_all,2,mean),h=smooth_pr)
		prior_pr <- prior_pr@yhat/sum(prior_pr@yhat)
	} else prior_pr <- apply(promoter_all,2,mean)
	
	# write potentials file
	string <- paste(prior_pr,collapse=",")
	promoterPots <- paste("\nNAME:\t\tpot_",c("P.M.prior",promoterVars),"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="")
	potentials <- file(paste("./",i,"/full_model/factorPotentials.txt",sep=""),"w")
	cat(promoterPots,file=potentials)
	close(potentials)

	tempFac <- rbind(tempS_G2,tempS_G1)
	colnames(tempFac) <- promoter_CpGs[1:length(promoterVars)]
	rownames(tempFac) <- c(G2,G1)
	eval(parse(text = paste('write.table(', paste('tempFac,file = "./',i,'/full_model/full_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	
	# build and query the full model with T and AN samples
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
			prior_pr <- kernsm(apply(promoter_all[cur,],2,mean),h=smooth_pr)
			prior_pr <- prior_pr@yhat/sum(prior_pr@yhat)
		} else prior_pr <- apply(promoter_all[cur,],2,mean)
		
		string <- paste(prior_pr,collapse=",")
		promoterPots <- paste("\nNAME:\t\tpot_",c("P.M.prior",promoterVars),"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="")
		potentials <- file(paste("./",i,"/null/G2_model/factorPotentials.txt",sep=""),"w")
		cat(promoterPots,file=potentials)
		close(potentials)
		# query
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/null/G2_model/ -l -n - ./',i,'/G2_model/all/T_VarData.tab ./',i,'/null/G2_model/G2_FacData.tab',sep=""))
		G2_G2_likelihoods <- as.numeric(substring(string[-1],17))
		
		# G1
		tempFac_G1 <- as.data.frame(tempFac[-cur,])
		colnames(tempFac_G1) <- promoter_CpGs[1:length(promoterVars)]
		rownames(tempFac_G1) <- G1
		eval(parse(text = paste('write.table(', paste('tempFac_G1,file = "./',i,'/null/G1_model/G1_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		if (smooth > 0) {
			prior_pr <- kernsm(apply(promoter_all[-cur,],2,mean),h=smooth_pr)
			prior_pr <- prior_pr@yhat/sum(prior_pr@yhat)
		} else prior_pr <- apply(promoter_all[-cur,],2,mean)
		
		string <- paste(prior_pr,collapse=",")
		promoterPots <- paste("\nNAME:\t\tpot_",c("P.M.prior",promoterVars),"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_pr,"]((",string,"))\nPC_MAT:\t\t[1,",res_pr,"]((",paste(rep(1,res_pr),collapse=","),"))\n",sep="")
		potentials <- file(paste("./",i,"/null/G1_model/factorPotentials.txt",sep=""),"w")
		cat(promoterPots,file=potentials)
		close(potentials)
		# query
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/null/G1_model/ -l -n - ./',i,'/G1_model/all/AN_VarData.tab ./',i,'/null/G1_model/G1_FacData.tab',sep=""))
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
