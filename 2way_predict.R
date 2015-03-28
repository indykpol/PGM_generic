args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[1])

load("essentials_allBRCA.RData")
IDs_length <- nchar(G1[1])+2
G1 <- G1[-c(1,2,5,6,7,8,10,14,18,20,21,22,23,24,25,26,27,28,29,30,32,33,34,35,37,38,39,41,42,45,46,47,48,49,51,53,54,56,57,58,59,61,62,64,65,67,68,70,71,75,76,77,78,80,81)]
G2 <- G2[-c(1,2,3,5,7,9,10,12,13,14,15,16,19,20,21,22,23,24,25,26,27,30,31,32,34,35,36,37,38,39,41,42,44,46,47,48,52,53,54,55,58,59,61,62,64,65,66,67,69,72,74,76,77,79,81,82,83,87,88,89,90,91,92,93,94,95,97,98,99,101,102,103,104,106,108,112,113,115,116,117,119,120,121,122,125,126,128,129,130,131,132,133,135,136,137,139,141,142,143,145,146,147,148,149,152,155,156,157,158,159,160,164,165,166,168,169,170,172,173,174,176,178,180,181,182,183,184,187,189,190,191,192,193,194,195,197,198,201,202,203,204,205,208,210,212,214,215,217,219,220,221,223,224,225,228,229,232,237,238,239,240,241,244,245,246,247,248,250,251,252,253,254,255,256,259,260,263,264,265,266,268,269,270,272,277,278,279,280,282,283,284,287,289,290,291,292,295,296,297,298,299,300,301,302,303,305,308,309,310,311,313,314,316,317,318,319,320,321,322,324,326,330,332,333,334,335,336,337,339,340,342,345,349,351,352,353,354,355,356,357,358,360,361,362,363,364,365,367,368,369,370,371,372,373,376,377,379,380,382,383,387,390,391,393,397,399,400,402,403,404,405,406,407,409,411,412,413,414,415,417,418,420,422,423,424,426,427,428,429,431,432,433,435,436,437,438,439,442,445,447,448,450,451,452,453,454,455,456,459,461,463,467,468,471,472,473,474,475,477,479,481,482,484,485,486,488,489,490,491,493,494,495,497,500,502,503,505,507,509,510,511,512,513,514,517,518,519,520,521,523,524,525,527,528,529,530,531,532,534,535,536,537,539,540,541,545,546,551,552,553,554,555,556,557,559,561,567,569,570,575,579,580,582,583,584,585,586,589,591,594,596,598,600,602,603,604,605,606,609,610,611,612,613,614,616,618,619,620,621,622,623,625,628,630,631,632,633,634,635,636,638,639,640,641,642,643,645,646,647,648,649,650,652,653,654,658,659,660,661,662,663,664,666,667,669,671,672,673,674,676,677,679,680,681,682,684,685,686,688,689,690,691,696,698,702,703,705,706,707,708,709,711,712,713,714,715,716,718,720,721,722,723,724,726,728,729,730)]
samples <- c(G1,G2)

integrand_e <- function(x,k) {dpois(k,x)}
integrand_m <- function(x,mean) {dnorm(x=mean, mean=x, sd=0.14)}

cat(paste("predicting for ",i,"\n",sep=""))
ptm <- proc.time()[3]
# prepare models and directories
system(intern=TRUE,command=paste('tar xf ',i,'.tar',sep=""))
system(command=paste('mkdir ./',i,'/toPredict',sep=""))

############### generate "missing" Var data ################
tempVar <- matrix(rep(".",length(samples)*3), nrow=length(samples), ncol=3)
colnames(tempVar) <- c("NAME:\tEXPR","M.GB","M.P")
rownames(tempVar) <- samples
eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/toPredict/samples_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
############### identify constitutive CpGs ################
IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',workingList_BRCA[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',workingList_BRCA[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',workingList_BRCA[i],'"',sep=""))))
promoterVars <- promoterVars_template[1:length(IDs_promoter)]
promoter_CpGs <- template_promoter_CpGs[1:length(IDs_promoter)]
IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',workingList_BRCA[i],'"',sep="")))))
geneBodyVars <- geneBodyVars_template[1:length(IDs_body)]
geneBody_CpGs <- template_body_CpGs[1:length(IDs_body)]
ncol = length(IDs_promoter) + length(IDs_body) + 1
############################################################

breaksEXPRESSION <- t(read.table(nrow=1, skip=1, file=paste("./",i,".result", sep="")))[,1]
breaksBODY <- t(read.table(nrow=1, skip=2, file=paste("./",i,".result", sep="")))[,1]
breaksPROMOTER <- t(read.table(nrow=1, skip=3, file=paste("./",i,".result", sep="")))[,1]

res_expr <- length(breaksEXPRESSION)-1
res_pr <- length(breaksPROMOTER)-1
res_gb <- length(breaksBODY)-1

all_labels_pr <- as.character(seq(1, res_pr, 1))
all_labels_gb <- as.character(seq(1, res_gb, 1))
all_labels_expr <- as.character(seq(1, res_expr, 1))
promoter_samples <- matrix(ncol=res_pr, nrow=length(samples))
body_samples <- matrix(ncol=res_gb, nrow=length(samples))
expr_samples <- matrix(ncol=res_expr, nrow=length(samples))

epsilon_pr <- 1/(length(samples)*length(IDs_promoter))/res_pr
epsilon_gb <- 1/(length(samples)*length(IDs_body))/res_gb
epsilon_e <- 1/length(samples)/res_expr

tempS_samples <- matrix(ncol=ncol,nrow=length(samples))
for (current_sample in 1:length(samples)) {
	# expression
	read_count <- trunc(as.numeric(counts[workingList_BRCA[i],samples[current_sample]]))
	lambdas <- breaksEXPRESSION * factors_ls[samples[current_sample]]
	frequencies_expr <- rep(0, length(breaksEXPRESSION)-1)
	for (freq in 1:res_expr) {
		frequencies_expr[freq] <- integrate(integrand_e, lower = lambdas[freq], upper = lambdas[freq+1], read_count, stop.on.error=FALSE)[1]
	}
	frequencies_expr <- unlist(frequencies_expr)
	if (length(which(frequencies_expr==0))==res_expr) frequencies_expr[length(frequencies_expr)] <- 1
	frequencies_expr <- frequencies_expr + epsilon_e
	frequencies_expr <- frequencies_expr/sum(frequencies_expr)
	
	# gene body
	cpg_list_gb <- NULL
	for (cpg in 1:length(geneBodyVars)) {
		miu <- mmatrix_pc[IDs_body[cpg],samples[current_sample]]
		if (!is.na(miu)) {
			frequencies_gb <- rep(0, res_gb)
			for (freq in 1:res_gb) {
			frequencies_gb[freq] <- integrate(integrand_m, lower=breaksBODY[freq], upper=breaksBODY[freq+1], mean=miu)$value
			}
			frequencies_gb <- unlist(frequencies_gb) + epsilon_gb
			frequencies_gb <- frequencies_gb/sum(frequencies_gb)
			cpg_list_gb[[cpg]] <- frequencies_gb
		} else cpg_list_gb[[cpg]] <- rep(1/res_gb, res_gb)
	}
	
	# promoter
	cpg_list_pr <- NULL
	for (cpg in 1:length(promoterVars)) {
		miu <- mmatrix_pc[IDs_promoter[cpg],samples[current_sample]]
		if (!is.na(miu)) {
			frequencies_pr <- rep(0, res_pr)
			for (freq in 1:res_pr) {
			frequencies_pr[freq] <- integrate(integrand_m, lower=breaksPROMOTER[freq], upper=breaksPROMOTER[freq+1], mean=miu)$value
			}
			frequencies_pr <- unlist(frequencies_pr) + epsilon_pr
			frequencies_pr <- frequencies_pr/sum(frequencies_pr)
			cpg_list_pr[[cpg]] <- frequencies_pr
		} else cpg_list_pr[[cpg]] <- rep(1/res_pr, res_pr)
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
	tempS_samples[current_sample,] <- tempS_formated
}

tempFac <- tempS_samples
colnames(tempFac) <- c("NAME:\tEXPR.likelihood", promoter_CpGs, geneBody_CpGs)
rownames(tempFac) <- samples
eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/toPredict/samples_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))

# query the full G1 model
string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/G1_model/all/ -l -n - ./',i,'/toPredict/samples_VarData.tab ./',i,'/toPredict/samples_FacData.tab',sep=""))
samples_G1_likelihoods <- as.numeric(substring(string[-1], IDs_length))

# query the full G2 model
string<-system(intern=TRUE,command=paste('./dfgEval_static --dfgSpecPrefix=./',i,'/G2_model/all/ -l -n - ./',i,'/toPredict/samples_VarData.tab ./',i,'/toPredict/samples_FacData.tab',sep=""))
samples_G2_likelihoods <- as.numeric(substring(string[-1], IDs_length))

# G2 posterior probability calculation
out <- cbind(samples, exp(-samples_G2_likelihoods) / (exp(-samples_G1_likelihoods) + exp(-samples_G2_likelihoods)), samples_G2_likelihoods, samples_G1_likelihoods)
colnames(out) <- c("sample_ID","posterior_G2","G2_mloglik","G1_mloglik")
eval(parse(text=paste('write.table(as.data.frame(out), col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE, file="./',i,'.predicted")', sep="")))

cat(paste("done predicting for ",i," in ", sprintf("%.2f", (proc.time()[3]-ptm)/60)," minutes\n",sep=""))
