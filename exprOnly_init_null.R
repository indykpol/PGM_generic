args <- commandArgs(trailingOnly = TRUE)
beg <- as.numeric(args[1])
end <- as.numeric(args[2])

# data <- args[1]
load("./essentials_allBRCA.RData")
IDs_length <- nchar(G1[1])+3
G1 <- G1[1:55]
G2 <- G2[1:486]

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
	epsilon_e_t <- 1/length(G2)/res_expr
	epsilon_e_an <- 1/length(G1)/res_expr
	
	smooth_e <- 10/(mean(length(G1),length(G2))/(res_expr*2))
	###########################################################
	############### binning scheme defined here ###############
	# expression
	temp <- counts_plusOne[workingList_BRCA[i],c(G2,G1)]
	tempAN <- matrix(ncol=2)
	colnames(tempAN) <- c("cpm","density")
	for (j in 1:length(temp)) {
		lambda <- as.numeric(temp[j])
		X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
		current <- factors_ls[c(G2,G1)[j]]
		tempAN <- rbind(tempAN,cbind(X/current,dpois(X,lambda=lambda)*current))
	}
	tempAN <- as.data.frame(tempAN[-1,],)
	tempAN <- tempAN[order(tempAN$cpm),]
	tempAN[,3] <- cumsum(tempAN[,2])
	tempAN[,3] <- tempAN[,3]/max(tempAN[,3])
	breaks <- NULL
	noBreaks <- res_expr-1
	for (j in 1:noBreaks) { breaks <- c (breaks, tempAN[which(tempAN[,3] >= j*(1/(1+noBreaks))),1][1])}
	max <- max(counts_plusOne[workingList_BRCA[i],c(G2,G1)])
	max_boundary <- 20+max+(4*max*max^(-1/2))
	breaksEXPRESSION <- c(0,breaks,max_boundary)
	
	all_labels_expr <- as.character(seq(1,res_expr,1))
	expr_t <- matrix(ncol=res_expr,nrow=length(G2))
	expr_an <- matrix(ncol=res_expr,nrow=length(G1))
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
		frequencies_expr <- frequencies_expr + epsilon_e_an
		frequencies_expr <- frequencies_expr/sum(frequencies_expr)
		
		#start precomputing correct initialization of parameters
		expr_an[current_sample,] <- frequencies_expr
		tempFac[current_sample,] <- paste('[1,',res_expr,']((',paste(frequencies_expr,sep="",collapse=","),'))',sep="")
	}
	# precompute correct initialization of parameters for AN-only model
	if (smooth > 0) {
		prior_expr <- kernsm(apply(expr_an,2,mean),h=smooth_e)
		prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
	} else prior_expr <- apply(expr_an,2,mean)
	
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
		frequencies_expr <- frequencies_expr + epsilon_e_t
		frequencies_expr <- frequencies_expr/sum(frequencies_expr)
		
		#start precomputing correct initialization of parameters
		expr_t[current_sample,] <- frequencies_expr
		
		tempFac[current_sample,] <- paste('[1,',res_expr,']((',paste(frequencies_expr,sep="",collapse=","),'))',sep="")
	}
	# precompute correct initialization of parameters for T-only model
	if (smooth > 0) {
		prior_expr <- kernsm(apply(expr_t,2,mean),h=smooth_e)
		prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
	} else prior_expr <- apply(expr_t,2,mean)
	
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
	expr_all <- rbind(expr_t,expr_an)
	
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
	if (sd(Ds) != 0) zscore <- (D - mean(Ds)) / sd(Ds) else zscore <- -6
	pval_zscore <- pnorm(zscore,lower.tail=FALSE)
	############################################################################
	eval(parse(text=paste('write.table(x=t(c(pval_zscore,D,mean(Ds),sd(Ds),zscore)), col.names=FALSE, row.names=FALSE, append=TRUE, file="./',i,'.result")',sep="")))
	cat(paste("done ",i," in ", sprintf("%.2f", (proc.time()[3]-ptm)/60)," minutes\n",sep=""))
}
