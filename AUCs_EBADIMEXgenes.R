library(pROC)
sum_asLogs <- function (x) {
	if (x[1] >= x[2]) x[1] + log(1+ exp(x[2]-x[1])) else x[2] + log(1+ exp(x[1]-x[2]))
}

calculatePosterior <- function(x) {
	1 - exp(x[,1] - apply(x, 1, sum_asLogs))
}

AUCs_indiv <- NULL
AUCs_comb <- NULL

for (cv_run in 1:length(CV_folds)) {
	auc_top_combined_PINCAGE <- auc_top_PINCAGE <- NULL
	
	running_logliks <- CV_folds[[cv_run]][[1]][,3:4]
	for (i in 1:length(CV_folds[[cv_run]])) {
		response <- as.factor(c(rep(1, 28), rep(0, 28)))
		logliks <- CV_folds[[cv_run]][[i]][,3:4]
		auc_top_PINCAGE[[i]] <- tryCatch(auc(response = response, predictor = calculatePosterior(logliks))[1], error = function(e) return(NA))
		auc_top_combined_PINCAGE[[i]] <- tryCatch(auc(response = response, predictor = calculatePosterior(running_logliks))[1], error = function(e) return(NA))
		running_logliks <- running_logliks + logliks
	}
	AUCs_indiv[[cv_run]] <- auc_top_PINCAGE
	AUCs_comb[[cv_run]] <- auc_top_combined_PINCAGE
}

# Analyze AUC tables
auc_top_combined_PINCAGE <- data.frame(model="auc_top_combined_PINCAGE", rank=1:100, AUC_mean =apply(sapply(1:50, FUN = function(x) AUCs_comb[[x]]), 1, mean) , AUC_SD = apply(sapply(1:50, FUN = function(x) AUCs_comb[[x]]), 1, sd))