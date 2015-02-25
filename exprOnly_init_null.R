args <- commandArgs(trailingOnly = TRUE)
beg <- as.numeric(args[1])
end <- as.numeric(args[2])

# data <- args[1]
load("./essentials_allBRCA.RData")
counts_plusOne <- counts_plusOne-1
IDs_length <- nchar(G1[1])+2
G1 <- G1[c(1,2,5,6,7,8,10,14,18,20,21,22,23,24,25,26,27,28,29,30,32,33,34,35,37,38,39,41,42,45,46,47,48,49,51,53,54,56,57,58,59,61,62,64,65,67,68,70,71,75,76,77,78,80,81)]
G2 <- G2[c(1,2,3,5,7,9,10,12,13,14,15,16,19,20,21,22,23,24,25,26,27,30,31,32,34,35,36,37,38,39,41,42,44,46,47,48,52,53,54,55,58,59,61,62,64,65,66,67,69,72,74,76,77,79,81,82,83,87,88,89,90,91,92,93,94,95,97,98,99,101,102,103,104,106,108,112,113,115,116,117,119,120,121,122,125,126,128,129,130,131,132,133,135,136,137,139,141,142,143,145,146,147,148,149,152,155,156,157,158,159,160,164,165,166,168,169,170,172,173,174,176,178,180,181,182,183,184,187,189,190,191,192,193,194,195,197,198,201,202,203,204,205,208,210,212,214,215,217,219,220,221,223,224,225,228,229,232,237,238,239,240,241,244,245,246,247,248,250,251,252,253,254,255,256,259,260,263,264,265,266,268,269,270,272,277,278,279,280,282,283,284,287,289,290,291,292,295,296,297,298,299,300,301,302,303,305,308,309,310,311,313,314,316,317,318,319,320,321,322,324,326,330,332,333,334,335,336,337,339,340,342,345,349,351,352,353,354,355,356,357,358,360,361,362,363,364,365,367,368,369,370,371,372,373,376,377,379,380,382,383,387,390,391,393,397,399,400,402,403,404,405,406,407,409,411,412,413,414,415,417,418,420,422,423,424,426,427,428,429,431,432,433,435,436,437,438,439,442,445,447,448,450,451,452,453,454,455,456,459,461,463,467,468,471,472,473,474,475,477,479,481,482,484,485,486,488,489,490,491,493,494,495,497,500,502,503,505,507,509,510,511,512,513,514,517,518,519,520,521,523,524,525,527,528,529,530,531,532,534,535,536,537,539,540,541,545,546,551,552,553,554,555,556,557,559,561,567,569,570,575,579,580,582,583,584,585,586,589,591,594,596,598,600,602,603,604,605,606,609,610,611,612,613,614,616,618,619,620,621,622,623,625,628,630,631,632,633,634,635,636,638,639,640,641,642,643,645,646,647,648,649,650,652,653,654,658,659,660,661,662,663,664,666,667,669,671,672,673,674,676,677,679,680,681,682,684,685,686,688,689,690,691,696,698,702,703,705,706,707,708,709,711,712,713,714,715,716,718,720,721,722,723,724,726,728,729,730)]

res_expr <- 25
smooth = 1
if (smooth > 0) library(aws)

integrand_e <- function(x,k) {dpois(k,x)}

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
	
	# generate "missing" Var data
	ncol=1
	tempVar_G1 <- matrix(rep(".",length(G1)*ncol),nrow=length(G1),ncol=ncol)
	colnames(tempVar_G1)[1] <- "NAME:\tEXPR"
	rownames(tempVar_G1) <- G1
	eval(parse(text = paste('write.table(', paste('tempVar_G1,file = "./',i,'/G1_model/all/G1_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	tempVar_G2 <- matrix(rep(".",length(G2)*ncol),nrow=length(G2),ncol=ncol)
	colnames(tempVar_G2)[1] <- "NAME:\tEXPR"
	rownames(tempVar_G2) <- G2
	eval(parse(text = paste('write.table(', paste('tempVar_G2,file = "./',i,'/G2_model/all/G2_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	tempVar <- rbind(tempVar_G2,tempVar_G1)
	eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/full_model/full_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	############################################################
	######### calculate epsilons and smoothing params ##########
	epsilon_e_g2 <- 1/length(G2)/res_expr
	epsilon_e_g1 <- 1/length(G1)/res_expr
	
	smooth_e <- 10/(mean(length(G1),length(G2))/(res_expr*2))
	###########################################################
	############### binning scheme defined here ###############
	# expression
	temp <- counts_plusOne[workingList_BRCA[i],c(G2,G1)]
	tempAN <- matrix(ncol=2)
	colnames(tempAN) <- c("cpm","density")
	for (j in 1:length(temp)) {
		lambda <- as.numeric(temp[j])
		if (lambda > 0) X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),0)),round(lambda+(4*lambda*lambda^(-1/2)))) else X <- 0:1
		current <- factors_ls[c(G2,G1)[j]]
		tempAN <- rbind(tempAN,cbind(X/current,dpois(X,lambda=lambda)))
	}
	tempAN <- as.data.frame(tempAN[-1,],)
	tempAN <- tempAN[order(tempAN$cpm),]
	tempAN[,3] <- cumsum(tempAN[,2])
	if (max(tempAN[,3]) != 0) tempAN[,3] <- tempAN[,3]/max(tempAN[,3])
	breaks <- NULL
	noBreaks <- res_expr-1
	for (j in 1:noBreaks) { breaks <- c(breaks, tempAN[which(tempAN[,3] >= j*(1/(1+noBreaks))),1][1])}
	max <- max(temp)
	if (max > 0) max_boundary <- 20+max+(4*max*max^(-1/2)) else max_boundary <- 20
	breaks <- unique(breaks)
	if (length(breaks) == 1) breaks <- seq(max(tempAN[,1]),max_boundary-max(tempAN[,1]),(max_boundary-max(tempAN[,1]))/24) else if (length(breaks) != 24){
		breaks <- breaks[-which(breaks %in% 0)]
		breaks <- c(breaks,seq(max(tempAN[,1])+max(breaks),max_boundary-max(tempAN[,1]),(max_boundary-max(tempAN[,1])-max(breaks))/(24-length(breaks))))
	}
	breaksEXPRESSION <- c(0,breaks,max_boundary)
	
	all_labels_expr <- as.character(seq(1,res_expr,1))
	expr_g2 <- matrix(ncol=res_expr,nrow=length(G2))
	expr_g1 <- matrix(ncol=res_expr,nrow=length(G1))
	########################################################
	# dynamic generation of model specification files here #
	########################################################
	####################### state map ######################
	stateMaps <- file(paste("./",i,"/stateMaps.txt",sep=""),"w")
	exprMap <- paste("NAME:\texprMap\nSYMBOLS:\t",paste(seq(1,res_expr,1),collapse=" "),"\nMETA_SYMBOLS:\t.=",paste(seq(1,res_expr,1),collapse=" "),"; *=",paste(seq(1,res_expr,1),collapse=" "),";\n\n",collapse="",sep="")
	cat(exprMap,file=stateMaps)
	close(stateMaps)
	####################### variables ######################
	variables <- file(paste("./",i,"/variables.txt",sep=""),"w")
	cat("STATE_MAP_NAME:\texprMap\nVAR_NAMES:\tEXPR\n\n",sep="",file=variables)
	close(variables)
	##################### factor graph #####################
	factorGraph <- file(paste("./",i,"/factorGraph.txt",sep=""),"w")
	cat("NAME:\tEXPR.likelihood\nNB1:\tEXPR\nPOT:\tpot_EXPR.likelihood\n",file=factorGraph)
	cat("\nNAME:\tEXPR.prior\nNB1:\tEXPR\nPOT:\tpot_EXPR.prior\n",file=factorGraph)
	close(factorGraph)
	
	system(command=paste('cp ./',i,'/*.txt ./',i,'/G1_model/all',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/G2_model/all',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/full_model/',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/null/G1_model',sep=""))
	system(command=paste('cp ./',i,'/*.txt ./',i,'/null/G2_model',sep=""))
	####################################################################
	############################## G1 model ############################
	## full G1 model developed from here, to obtain likelihoods of G1##
	# generate FacData for full set of G1
	tempFac <- matrix(ncol=ncol,nrow=length(G1))
	for (current_sample in 1:length(G1)) {
		# expression
		read_count <- trunc(as.numeric(counts_plusOne[workingList_BRCA[i],G1[current_sample]]))
		lambdas <- breaksEXPRESSION * factors_ls[G1[current_sample]]
		frequencies_expr <- rep(0,length(breaksEXPRESSION)-1)
		for (freq in 1:res_expr) {
			frequencies_expr[freq] <- integrate(integrand_e, lower = lambdas[freq], upper = lambdas[freq+1], read_count,stop.on.error=FALSE)[1]
		}
		frequencies_expr <- unlist(frequencies_expr)
		if (length(which(frequencies_expr==0))==res_expr) frequencies_expr[length(frequencies_expr)] <- 1
		frequencies_expr <- frequencies_expr + epsilon_e_g1
		frequencies_expr <- frequencies_expr/sum(frequencies_expr)
		
		#start precomputing correct initialization of parameters
		expr_g1[current_sample,] <- frequencies_expr
		tempFac[current_sample,] <- paste('[1,',res_expr,']((',paste(frequencies_expr,sep="",collapse=","),'))',sep="")
	}
	# precompute correct initialization of parameters for AN-only model
	if (smooth > 0) {
		prior_expr <- kernsm(apply(expr_g1,2,mean),h=smooth_e)
		prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
	} else prior_expr <- apply(expr_g1,2,mean)
	
	# write potentials file
	string <- paste(prior_expr,collapse=",")
	expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
	potentials <- file(paste("./",i,"/G1_model/all/factorPotentials.txt",sep=""),"w")
	cat(expr.pots,file=potentials)
	close(potentials)
	
	tempS_G1 <- tempFac
	#tempFac <- cbind(tempFac,paste('[1,',res_expr,']((',paste(potentials_an,sep="",collapse=","),'))',sep=""))
	colnames(tempFac) <- c("NAME:\tEXPR.likelihood")
	rownames(tempFac) <- G1
	eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/G1_model/all/G1_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	# build and query the full model with T and AN samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/G1_model/all/ -l -n - ./',i,'/G1_model/all/G1_VarData.tab ./',i,'/G1_model/all/G1_FacData.tab',sep=""))
	G1_G1_likelihoods <- as.numeric(substring(string[-1],IDs_length))
	####################################################################
	############################# G2 model ##############################
	### full G2 model developed from here, to obtain likelihoods of G2 ##
	# generate FacData for full set of G2
	tempFac <- matrix(ncol=ncol,nrow=length(G2))
	for (current_sample in 1:length(G2)) {
		# expression
		read_count <- trunc(as.numeric(counts_plusOne[workingList_BRCA[i],G2[current_sample]]))
		lambdas <- breaksEXPRESSION * factors_ls[G2[current_sample]]
		frequencies_expr <- rep(0,length(breaksEXPRESSION)-1)
		for (freq in 1:res_expr) {
			frequencies_expr[freq] <- integrate(integrand_e, lower = lambdas[freq], upper = lambdas[freq+1], read_count,stop.on.error=FALSE)[1]
		}
		frequencies_expr <- unlist(frequencies_expr)
		if (length(which(frequencies_expr==0))==res_expr) frequencies_expr[length(frequencies_expr)] <- 1
		frequencies_expr <- frequencies_expr + epsilon_e_g2
		frequencies_expr <- frequencies_expr/sum(frequencies_expr)
		
		#start precomputing correct initialization of parameters
		expr_g2[current_sample,] <- frequencies_expr
		
		tempFac[current_sample,] <- paste('[1,',res_expr,']((',paste(frequencies_expr,sep="",collapse=","),'))',sep="")
	}
	# precompute correct initialization of parameters for T-only model
	if (smooth > 0) {
		prior_expr <- kernsm(apply(expr_g2,2,mean),h=smooth_e)
		prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
	} else prior_expr <- apply(expr_g2,2,mean)
	
	# write potentials file
	string <- paste(prior_expr,collapse=",")
	expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
	potentials <- file(paste("./",i,"/G2_model/all/factorPotentials.txt",sep=""),"w")
	cat(expr.pots,file=potentials)
	close(potentials)
	
	tempS_G2 <- tempFac
	colnames(tempFac) <- c("NAME:\tEXPR.likelihood")
	rownames(tempFac) <- G2
	eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/G2_model/all/G2_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	
	# query the full model with AN samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/G2_model/all/ -l -n - ./',i,'/G2_model/all/G2_VarData.tab ./',i,'/G2_model/all/G2_FacData.tab',sep=""))
	G2_G2_likelihoods <- as.numeric(substring(string[-1],IDs_length))
	###########################################################################
	######################## All data model ###################################
	## Full model developed from here, to obtain likelihoods of G2 and G1 ####
	# precompute correct initialization of parameters for joint model
	expr_all <- rbind(expr_g2,expr_g1)
	
	if (smooth > 0) {
		prior_expr <- kernsm(apply(expr_all,2,mean),h=smooth_e)
		prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
	} else prior_expr <- apply(expr_all,2,mean)
	
	# write potentials file
	string <- paste(prior_expr,collapse=",")
	expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
	potentials <- file(paste("./",i,"/full_model/factorPotentials.txt",sep=""),"w")
	cat(expr.pots,file=potentials)
	close(potentials)
	
	# generate FacData for G2
	tempFac <- rbind(tempS_G2,tempS_G1)
	rownames(tempFac) <- c(G2,G1)
	colnames(tempFac) <- c("NAME:\tEXPR.likelihood")
	eval(parse(text = paste('write.table(', paste('tempFac,file = "./',i,'/full_model/full_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
	
	# build and query the full model with G2 and G1 samples
	string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/full_model/ -l -n - ./',i,'/full_model/full_VarData.tab ./',i,'/full_model/full_FacData.tab',sep=""))
	allData_full_likelihoods <- as.numeric(substring(string[-1],IDs_length))
	###########################################################################
	######################## D calculation ####################################
	D <- 2*(sum(allData_full_likelihoods) - (sum(G1_G1_likelihoods)+sum(G2_G2_likelihoods)))
	###########################################################################
	################# P val calculation using null distr. #####################
	nruns <- 100
	Ds <- vector(length=nruns,mode="numeric")
	for (run in 1:nruns) {
		cur <- sample(x=1:(length(G2)+length(G1)),size=length(G2),replace=FALSE)
		
		# G2
		tempFac_T <- as.data.frame(tempFac[cur,])
		colnames(tempFac_T) <- c("NAME:\tEXPR.likelihood")
		rownames(tempFac_T) <- G2
		eval(parse(text = paste('write.table(', paste('tempFac_T,file = "./',i,'/null/G2_model/G2_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		if (smooth > 0) {
			prior_expr <- kernsm(apply(expr_all[cur,],2,mean),h=smooth_e)
			prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
		} else prior_expr <- apply(expr_all[cur,],2,mean)
		
		string <- paste(prior_expr,collapse=",")
		expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
		potentials <- file(paste("./",i,"/null/G2_model/factorPotentials.txt",sep=""),"w")
		cat(expr.pots,file=potentials)
		close(potentials)
		# query
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/null/G2_model/ -l -n - ./',i,'/G2_model/all/G2_VarData.tab ./',i,'/null/G2_model/G2_FacData.tab',sep=""))
		G2_G2_likelihoods <- as.numeric(substring(string[-1],IDs_length))
		
		# G1
		tempFac_AN <- as.data.frame(tempFac[-cur,])
		colnames(tempFac_AN) <- c("NAME:\tEXPR.likelihood")
		rownames(tempFac_AN) <- G1
		eval(parse(text = paste('write.table(', paste('tempFac_AN,file = "./',i,'/null/G1_model/G1_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
		
		if (smooth > 0) {
			prior_expr <- kernsm(apply(expr_all[-cur,],2,mean),h=smooth_e)
			prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
		} else prior_expr <- apply(expr_all[-cur,],2,mean)
		
		string <- paste(prior_expr,collapse=",")
		expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
		potentials <- file(paste("./",i,"/null/G1_model/factorPotentials.txt",sep=""),"w")
		cat(expr.pots,file=potentials)
		close(potentials)
		# query
		string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/null/G1_model/ -l -n - ./',i,'/G1_model/all/G1_VarData.tab ./',i,'/null/G1_model/G1_FacData.tab',sep=""))
		G1_G1_likelihoods <- as.numeric(substring(string[-1],IDs_length))
		# Ds calculation
		Ds[run] <- 2*(sum(allData_full_likelihoods) - (sum(G1_G1_likelihoods)+sum(G2_G2_likelihoods)))
	}
	if (sd(Ds) != 0 & D > 0.1) zscore <- (D - mean(Ds)) / sd(Ds) else zscore <- -6
	pval_zscore <- pnorm(zscore,lower.tail=FALSE)
	############################################################################
	eval(parse(text=paste('write.table(x=t(c(pval_zscore,D,mean(Ds),sd(Ds),zscore)), col.names=FALSE, row.names=FALSE, append=TRUE, file="./',i,'.result")',sep="")))
	cat(paste("done ",i," in ", sprintf("%.2f", (proc.time()[3]-ptm)/60)," minutes\n",sep=""))
}
