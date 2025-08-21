library(mgcv)
library(R.matlab)
library(gratia)

MFC <- readMat("E:/dHCP/Figure/364sub/Result_Figure/R/MFC_v.mat")  #MFC_v：Vertex-level gMFC/sMFC/MFC
Age <- readMat("E:/dHCP/Figure/364sub/Result_Figure/R/Age_364.mat") 
GAM_data <- list()
GAM_data$PMA <- Age$PMA364     #PMA at scan
GAM_data$MFC_gs <- MFC$gMFC.v  #gMFC
GAM_data$MFC_ss <- MFC$sMFC.v  #sMFC
GAM_data$MFC_gss <- MFC$MFC.v  #MFC
fit_GAM <- list()              #create an empty list
current_GAM_data <- list()     #create an empty list
current_GAM_data$PMA <- GAM_data$PMA
Pvalue_all <- matrix(1:8589, nrow = 8589, ncol = 1)
#age_term_all <- matrix(1:3126396, nrow = 364, ncol = 8589)
derivs <- matrix(1:8589, nrow = 8589, ncol = 1)


###gMFC
for (i in 1:nrow(GAM_data$MFC_gs)) {
  current_GAM_data$current_MFC <- GAM_data$MFC_gs[i, ]
  fit_GAM <- gam(current_MFC ~ s(PMA, bs="tp", k=3), data=current_GAM_data, method="REML")  #GAM model
  summary_model <- summary(fit_GAM)
  Pvalue_all[i,1] = summary_model$s.table[1,4]
  deriv <- derivatives(fit_GAM, term='s(PMA)', order=1)
  derivs[i,1] <- mean(deriv$.derivative)   #average derivative: growth rate
}
write.table(Pvalue_all,"E:/dHCP/Figure/364sub/Result_Figure/R/Pvalue_gs.csv",row.names=F, col.names=F)
write.table(derivs,"E:/dHCP/Figure/364sub/Result_Figure/R/derivs_gs.csv",row.names=F, col.names=F)


###sMFC
for (i in 1:nrow(GAM_data$MFC_ss)) {
  current_GAM_data$current_MFC <- GAM_data$MFC_ss[i, ]
  fit_GAM <- gam(current_MFC ~ s(PMA, bs="tp", k=3), data=current_GAM_data, method="REML")   #GAM model
  summary_model <- summary(fit_GAM)
  Pvalue_all[i,1] = summary_model$s.table[1,4]
  deriv <- derivatives(fit_GAM, term='s(PMA)', order=1)
  derivs[i,1] <- mean(deriv$.derivative)   #average derivative: growth rate
}

write.table(Pvalue_all,"E:/dHCP/Figure/364sub/Result_Figure/R/Pvalue_ss.csv",row.names=F, col.names=F)
write.table(derivs,"E:/dHCP/Figure/364sub/Result_Figure/R/derivs_ss.csv",row.names=F, col.names=F)


###MFC
for (i in 1:nrow(GAM_data$MFC_gss)) {
  current_GAM_data$current_MFC <- GAM_data$MFC_gss[i, ]
  fit_GAM <- gam(current_MFC ~ s(PMA, bs="tp", k=3), data=current_GAM_data, method="REML")   #GAM model
  summary_model <- summary(fit_GAM)
  Pvalue_all[i,1] = summary_model$s.table[1,4]
  deriv <- derivatives(fit_GAM, term='s(PMA)', order=1)
  derivs[i,1] <- mean(deriv$.derivative)   #average derivative: growth rate
}

write.table(Pvalue_all,"E:/dHCP/Figure/364sub/Result_Figure/R/Pvalue_gss.csv",row.names=F, col.names=F)
write.table(derivs,"E:/dHCP/Figure/364sub/Result_Figure/R/derivs_gss.csv",row.names=F, col.names=F)
