rm(list=ls())

library(splus2R)
library(foreach)
library(doParallel)
library(doSNOW)
library(mice)
library(caret)
library(pROC)
library(rJava)
options(java.parameters = "-Xmx25g")
library(bartMachine)
library(BART)

path<-"/Users/linjungyi/Desktop/"
source("/Users/linjungyi/Google Drive/SWAN Study/Codes and R output/Paper 2/Organized R code/function.R")

# Replications
## 1. Simulate X
## 2. Simulate Y
## 3. Ampute
## 4. Impute and var selection
## 5. Calculate metrics for final model
## 6. Repeate steps 1-4 250 times

# Set parameters
## Time: 3:22, 22 min
n_rep<-250 # Num of replications (250 for n = 1000 and 5000 and 1000 for n = 300 and 650)
n_impute<- 200# Num of imputed data 
maxit_impute<-5 # (5 to be)
ntreex<-50 # Default in bart=50 and bartMachine=50


# Output
comp_dat_list<-list(NA)
var_select_model<-list(NA)
final_var<-list(NA)


# Sim data
n<-1000 # 650, 300, 5000

#setup parallel backend to use many processors
seedx<-123 # seed: 1, 12, 123, 1234, 12345, 6, 67, 678, 6789, 67890
cores=detectCores()
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

# Generate data
# X
x1<-rbinom(n,1,0.5)
x2<-rbinom(n,1,0.5)
x3<-rnorm(n,0,1)
x4<-rnorm(n,0,1)
x5<-rnorm(n,0,1)
x6<-rgamma(n, shape=4, scale=0.6)
x7<-rnorm(n,-0.4*x5+0.4*x6+0.3*x5*x6,1)
x8<-rnorm(n, 0.1*x5*(x6-2)^2-0.1*x7^2, 1)
x9<-rnorm(n, 0.5*x3+0.3*x4-0.3*x5^2+0.2*x3*x4, 1)
x10<-rnorm(n, 0.1*x3^3-0.3*x4-0.4*x5+0.2*x9^2+0.3*x4*x5,1)

# Z: 20 N(0,1) and 20 bin(1,0.5)
z1<-rnorm(n,0,1); z2<-rnorm(n,0,1); z3<-rnorm(n,0,1); z4<-rnorm(n,0,1); z5<-rnorm(n,0,1)
z6<-rnorm(n,0,1); z7<-rnorm(n,0,1); z8<-rnorm(n,0,1); z9<-rnorm(n,0,1); z10<-rnorm(n,0,1)
z11<-rnorm(n,0,1); z12<-rnorm(n,0,1); z13<-rnorm(n,0,1); z14<-rnorm(n,0,1); z15<-rnorm(n,0,1)
z16<-rnorm(n,0,1); z17<-rnorm(n,0,1); z18<-rnorm(n,0,1); z19<-rnorm(n,0,1); z20<-rnorm(n,0,1)

z21<-rbinom(n,1,0.5); z22<-rbinom(n,1,0.5); z23<-rbinom(n,1,0.5); z24<-rbinom(n,1,0.5); z25<-rbinom(n,1,0.5)
z26<-rbinom(n,1,0.5); z27<-rbinom(n,1,0.5); z28<-rbinom(n,1,0.5); z29<-rbinom(n,1,0.5); z30<-rbinom(n,1,0.5)
z31<-rbinom(n,1,0.5); z32<-rbinom(n,1,0.5); z33<-rbinom(n,1,0.5); z34<-rbinom(n,1,0.5); z35<-rbinom(n,1,0.5)
z36<-rbinom(n,1,0.5); z37<-rbinom(n,1,0.5); z38<-rbinom(n,1,0.5); z39<-rbinom(n,1,0.5); z40<-rbinom(n,1,0.5)

p0<-  1.8*x1+0.5*x2+1.1*x3-
  0.4*exp(x5)-0.4*(x6-3.5)^2+0.3*(x7-1)^3+1.1*x8-1.1*x10+
  5*sin(0.1*pi*x4*x9)-
  0.4*x5*x10^2+0.4*x3^2*x8-2.7
p<-exp(p0)/(1+exp(p0))
y<-rbinom(n,1,p)

dat_comp<-data.frame(y,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,
                     z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,
                     z11,z12,z13,z14,z15,z16,z17,z18,z19,z20,
                     z21,z22,z23,z24,z25,z26,z27,z28,z29,z30,
                     z31,z32,z33,z34,z35,z36,z37,z38,z39,z40) # Because ampute can only ampue numeric, create dummy for categorical in data

# 3. Ampute
# Create interaction and polynomial for amputation
dat_comp$x5x6<-dat_comp$x5*dat_comp$x6 # used in x7 x8
dat_comp$x7x7<-dat_comp$x7*dat_comp$x7 # used in x8
dat_comp$x3x4<-dat_comp$x3*dat_comp$x4 # used in x9
dat_comp$x5x5<-dat_comp$x5*dat_comp$x5 # used in x9
dat_comp$x4x5<-dat_comp$x4*dat_comp$x5 # used in x10
dat_comp$x6x6<-dat_comp$x6*(dat_comp$x6) # used in y
dat_comp$x4x9<-dat_comp$x4*(dat_comp$x9) # used in y
dat_comp$x5x10<-dat_comp$x5*(dat_comp$x10) # used in y
dat_comp$x3x8<-(dat_comp$x3)*dat_comp$x8 # used in y

na_pattern<-matrix(c( c(0, rep(1,6), 1, 1, 1, 1, rep(1,40), rep(1,dim(dat_comp)[2]-51)), # y
                      c(1, rep(1,6), 0, 1, 1, 1, rep(1,40), rep(1,dim(dat_comp)[2]-51)), # x7
                      c(1, rep(1,6), 1, 0, 1, 1, rep(1,40), rep(1,dim(dat_comp)[2]-51)), # x8
                      c(1, rep(1,6), 1, 1, 0, 1, rep(1,40), rep(1,dim(dat_comp)[2]-51)), # x9
                      c(1, rep(1,6), 1, 1, 1, 0, rep(1,40), rep(1,dim(dat_comp)[2]-51)), # x10
                      c(0, rep(1,6), 0, 0, 1, 1, rep(1,40), rep(1,dim(dat_comp)[2]-51)), # y x7 x8
                      c(0, rep(1,6), 1, 1, 0, 0, rep(1,40), rep(1,dim(dat_comp)[2]-51)), # y x9 x10
                      c(0, rep(1,6), 1, 0, 1, 0, rep(1,40), rep(1,dim(dat_comp)[2]-51))# y x8 x10
), ncol=dim(dat_comp)[2], byrow=T)

na_wt<-matrix(c(c(0, 5,5,1,0,-1,-1, 1, 1, 0, 1, rep(0,40), rep(0,5), -0.5,1.5,-0.5,0.5), # y
                c(0, 0,0,0,0,1,1, 0, 0, 0, 0, rep(0,40), 1, rep(0,4), rep(0,4)), # x7
                c(5, 0,0,0,0,1,1, 1, 0, 0, 0, rep(0,40), 1, 1, rep(0,3), rep(0,4)), # x8
                c(5, 0,0,1,1,1,0, 0, 0, 0, 0, rep(0,40), 0, 0, 1, 1, 0, rep(0,4)), # x9
                c(5, 0,0,1,1,1,0, 0, 0, 1, 0, rep(0,40), 0, 0, 0, 0, 1, rep(0,4)), # x10
                c(0, 0,0,0,0,1,1, 0, 0, 0, 0, rep(0,40), 0, rep(0,4), rep(0,4)), # y x7 x8
                c(0, 0,0,1,1,1,0, 0, 0, 0, 0, rep(0,40), 0, 0, 0.5, 0.5, 0,rep(0,4)), # y x9 x10
                c(0, 0,0,0,0,1,0, 0, 0, 0, 0, rep(0,40), rep(0,dim(dat_comp)[2]-51))  # y x8 x10
), ncol=dim(dat_comp)[2], byrow=T)
colnames(na_pattern)<-colnames(na_wt)<-names(dat_comp)

dat_na0<-ampute(dat_comp, prop=0.60, mech="MAR", patterns = na_pattern, 
                freq = c(0.30,0.09,0.09,0.08,0.08,0.16,0.10,0.10), weights=na_wt,
                cont = TRUE, 
                type=c("RIGHT","RIGHT", "RIGHT","RIGHT", "RIGHT","TAIL", "TAIL","TAIL"))$amp
#table(dat_na0$y, useNA = "ifany") # 40.2% NA
#nrow(na.omit(dat_na0)) # 58.2% NA in total

dat_na<-dat_na0[,c("y","x1","x2","x3","x4","x5","x6","x7","x8","x9","x10",
                   "z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                   "z11","z12","z13","z14","z15","z16","z17","z18","z19","z20",
                   "z21","z22","z23","z24","z25","z26","z27","z28","z29","z30",
                   "z31","z32","z33","z34","z35","z36","z37","z38","z39","z40")]
dat_na[,c("x1","x2",
          "z21","z22","z23","z24","z25","z26","z27","z28","z29","z30",
          "z31","z32","z33","z34","z35","z36","z37","z38","z39","z40")]<-
  lapply(dat_na[,c("x1","x2",
                   "z21","z22","z23","z24","z25","z26","z27","z28","z29","z30",
                   "z31","z32","z33","z34","z35","z36","z37","z38","z39","z40")],factor)

## Change y into yes/no
dat_na$y<-factor(ifelse(is.na(dat_na$y)==1,NA,
                        ifelse(dat_na$y==1,"Yes",
                               ifelse(dat_na$y==0,"No",NA))))

# Split data
na_traindat_list<-list(NA)
na_validdat_list<-list(NA)

for(i in 1:n_rep) {
  # Split data
  train_id<-sample(c(1:n),n/2)
  na_traindat_list[[i]]<-dat_na[train_id,]
  na_validdat_list[[i]]<-dat_na[-train_id,]
  
}
na_train_valid_dat<-list("na_traindat_list"=na_traindat_list, "na_validdat_list"=na_validdat_list)


# Impute
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)
final_dat<-foreach(i=1:n_rep, .combine="comb", .multicombine = TRUE,  
                   .init=list(list(), list(), list(), list()), 
                   .packages = c("doParallel", "doSNOW"))  %dopar%{
                     library(mice)
                     
                     traindat_na<-na_train_valid_dat[[1]][[i]]
                     validdat_na<-na_train_valid_dat[[2]][[i]]

                     # Impute traindat using mice
                     trainmice_b<-mice(traindat_na, m=n_impute, maxit = maxit_impute)
                     trainimpute_dat0<-complete(trainmice_b, include=FALSE, action = "broad")
                     trainimpute_dat<-list(NA)
                     for(j in 1:n_impute){
                       num_col<-ncol(traindat_na)
                       trainimpute_dat[[j]]<-trainimpute_dat0[,c( (num_col*(j-1)+1):(num_col*j) )]
                       names(trainimpute_dat[[j]])<-names(dat_na)
                     }
                     i#mpute_traindat_list<-trainimpute_dat
                     #trainmice_list<-trainmice_b
                     
                     # Impute validdat
                     validmice_b<-mice(validdat_na, m=1, maxit = maxit_impute)
                     validimpute_dat0<-complete(validmice_b, include=FALSE)
                     
                     list("impute_traindat_list"=trainimpute_dat, "trainmice_list"=trainmice_b,
                          "impute_validdat_list"=validimpute_dat0, "validmice_list"=validmice_b)
                     
                   }

stopCluster(cl)
saveRDS(final_dat, paste(path,"x", seedx, "rrbart_imputedat_102721.rds", sep=""))

impute_traindat_list<-final_dat[[1]]
trainmice_list<-final_dat[[2]]
impute_validdat_list<-final_dat[[3]]
validmice_list<-final_dat[[4]]

# Model
cl <- makeCluster(cores[1]-2) #not to overload your computer
#cl <- makeCluster(25) #not to overload your computer
registerDoParallel(cl)

vs_prop_factor<-list(NA)
vs_se_factor<-list(NA)

finalMatrix<-NA
finalMatrix<-foreach(i=1:n_rep, .combine="comb", .multicombine = TRUE,  .init=list(list(),list(),list()), 
                     .packages = c("doParallel")) %:%
  
  foreach(b=1:n_impute, .combine="comb", .multicombine = TRUE, .init=list(list(),list(),list()),
          .packages = c("doParallel")) %dopar%{
                                    library(caret)
                                    library(pROC)
                                    library(rJava)
                                    options(java.parameters = "-Xmx25g")
                                    library(bartMachine)
                                    library(BART)

                                    mice_dat<-impute_traindat_list[[i]][[b]]
                                    
                                    # Models
                                    cov_mat<-mice_dat[,-which(colnames(mice_dat)%in%c("y"))]
                                    out_vec<-factor(as.character(mice_dat$y), levels=c("Yes","No"))

                                    
                                    ## bartMachine
                                    #Sys.time() 
                                    bm1<-bartMachine(y = out_vec, X=cov_mat, mem_cache_for_speed=FALSE,
                                                     use_missing_data = TRUE, use_missing_data_dummies_as_covars = TRUE, 
                                                     num_iterations_after_burn_in = 2000) # 1min
                                    #Sys.time() 
                                    
                                    ### VIP
                                    p_var_selected1<-get_var_counts_over_chain(bm1)/ntreex
                                    
                                    ### Sum VIP over the levels of categorical variables
                                    library(dplyr)
                                    vs_prop_factor0<-p_var_selected1 %>%  
                                      as_tibble() %>% 
                                      transmute(
                                        x1=x1_0+x1_1, x2=x2_0+x2_1, x3=x3,
                                        x4=x4, x5=x5,
                                        x6=x6, x7=x7, x8=x8,
                                        x9=x9,
                                        x10=x10,
                                        z1=z1, z2=z2, z3=z3, z4=z4, z5=z5,
                                        z6=z6, z7=z7, z8=z8, z9=z9, z10=z10,
                                        z11=z11, z12=z12, z13=z13, z14=z14, z15=z15,
                                        z16=z16, z17=z17, z18=z18, z19=z19, z20=z20,
                                        z21=z21_0+z21_1, z22=z22_0+z22_1, z23=z23_0+z23_1, z24=z24_0+z24_1, z25=z25_0+z25_1,
                                        z26=z26_0+z26_1, z27=z27_0+z27_1, z28=z28_0+z28_1, z29=z29_0+z29_1, z30=z30_0+z30_1,
                                        z31=z31_0+z31_1, z32=z32_0+z32_1, z33=z33_0+z33_1, z34=z34_0+z34_1, z35=z35_0+z35_1,
                                        z36=z36_0+z36_1, z37=z37_0+z37_1, z38=z38_0+z38_1, z39=z39_0+z39_1, z40=z40_0+z40_1)
                                    detach("package:dplyr")
                                    
                                    ### Mean and SD of VIP
                                    vs_prop_factor<-apply(vs_prop_factor0, 2, mean)
                                    vs_se_factor<-apply(vs_prop_factor0, 2, sd)/sqrt(n)
                                    
                                    list("post_prop"=vs_prop_factor0, "prop"=vs_prop_factor, 
                                         "se"=vs_se_factor)
                                  }

stopCluster(cl)
Sys.time() # 1:40
saveRDS(finalMatrix,  paste(path,"x",seedx,"finalMatrix_rrbart_102721.rds",sep=""))


# Diff in VIP
#finalMatrix<-readRDS(paste(path,"x",seedx,"finalMatrix_rrbart_102721.rds",sep=""))
post_prop_factor<-finalMatrix[[1]]
vs_prop_factor<-finalMatrix[[2]]
vs_se_factor<-finalMatrix[[3]]

# Final var & performance metrics
## Output
var_select<-list(NA)
prec_recall_f1_i<-list(NA)
metrics_mat<-data.frame(matrix(NA,nrow=1, ncol=4))
names(metrics_mat)<-c("Precision","Recall","F1", "AUC")


for(i in 1:n_rep){
  vs_prop_factor_i<-vs_prop_factor[[i]]
  post_prop_factor_i<-post_prop_factor[[i]]
  
  ## 1. Identify var w/ smallest VIP
  vs_prop_factor_mat<-matrix(unlist(vs_prop_factor_i), byrow=T, nrow=n_impute)
  colnames(vs_prop_factor_mat)<-names(vs_prop_factor[[1]][[1]])
  ave_vs_prop<-apply(vs_prop_factor_mat, 2, mean)
  var_min_vip<-names(ave_vs_prop)[ave_vs_prop==min(ave_vs_prop)]
  
  ## 2. Calculate diff from var w/ min vip
  post_diff_vip0<-lapply(post_prop_factor_i, function(x){
    apply(x, 2, function(y){
      y-x[,var_min_vip]
    })
  })
  
  post_diff_vip<-lapply(post_diff_vip0, function(x){
    mat<-matrix(unlist(x), byrow=FALSE, ncol=50)
    colnames(mat)<-names(vs_prop_factor[[1]][[1]])
    return(mat)})
  
  
  
  ## 3. Apply RR to diff
  ### Mean within each impute
  ave_diff_vip<-data.frame(matrix(unlist(lapply(post_diff_vip, function(x){apply(x,2,mean)})), nrow=n_impute, byrow=T))
  names(ave_diff_vip)<-names(vs_prop_factor[[1]][[1]])
  
  ### SE within each impute
  se_diff_vip<-data.frame(matrix(unlist(lapply(post_diff_vip, function(x){apply(x,2, function(x1){sd(x1)/sqrt(n)} )})), nrow=n_impute, byrow=T))
  names(se_diff_vip)<-names(vs_prop_factor[[1]][[1]])
  
  
  ### Within var
  Vw<-apply(se_diff_vip, 2, function(x){sum(x^2)})/n_impute
  
  ### Between var
  p_mean<-matrix(rep( apply(ave_diff_vip, 2, mean), n_impute), nrow=n_impute, byrow=T)
  Vb<-apply(ave_diff_vip-p_mean, 2, function(x){sum(x^2)})/(n_impute-1)
  colnames(p_mean)<-c(1:ncol(p_mean))
  
  ## Rubin
  rubin_se<-sqrt(Vw+(1+1/n_impute)*Vb)
  
  # CI 
  ## RR reference t-dist
  lambda<-(Vb+Vb/n_impute)/rubin_se
  df_old<-(n_impute-1)/(lambda^2)
  k<-1
  df_obs<-((n-k)+1)/((n-k)+3)*(n-k)*(1-lambda)
  df_adj<-df_old*df_obs/(df_old+df_obs)
  t_df<-qt(0.975, df=df_adj)
  
  p_up<-p_mean[1,]+t_df*rubin_se
  p_low<-p_mean[1,]-t_df*rubin_se
  
  ## Normal approximation
  #p_up<-p_mean[1,]+1.96*rubin_se
  #p_low<-p_mean[1,]-1.96*rubin_se
  
  var_select[[i]]<-names(vs_prop_factor[[1]][[1]])[p_low>0&is.na(p_low)==0]
  
  ## 4. Performance metrics (precision, recall, F1)
  prec_recall_f1_i[[i]]<-prec_recall_f1(gold_pos=c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"),
                                            gold_neg=c("z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                                                       "z11","z12","z13","z14","z15","z16","z17","z18","z19","z20",
                                                       "z21","z22","z23","z24","z25","z26","z27","z28","z29","z30",
                                                       "z31","z32","z33","z34","z35","z36","z37","z38","z39","z40"),
                                            var_selected=var_select[[i]])

  
  metrics_mat[i,c("Precision","Recall","F1")]<-c(
    prec_recall_f1_i[[i]]$precision, prec_recall_f1_i[[i]]$recall, 
    prec_recall_f1_i[[i]]$f1)
  
  ## 5. AUC
  ### Fit a model on valid data for calculation of AUC
  mice_valid_dat<-impute_validdat_list[[i]]
  cov_mat<-data.frame(mice_valid_dat[,var_select[[i]]])
  out_vec<-factor(as.character(mice_valid_dat$y), levels=c("Yes","No"))
  bm_valid<-bartMachine(y = out_vec, X=cov_mat, mem_cache_for_speed=FALSE,
                        use_missing_data = TRUE, use_missing_data_dummies_as_covars = TRUE) 
  
  mice_valid_dat$phat<-predict(bm_valid, type="prob", new_data = cov_mat)
  
  ### AUC
  metrics_mat[i,"AUC"]<-auc(roc(response=mice_valid_dat$y, predictor=mice_valid_dat$phat, 
                                  levels=c("No", "Yes"), direction="<"))
  
  
}

saveRDS(var_select, paste(path,"x", seedx, "var_select_rrbart_tdist_102721_na07.rds",sep=""))

# Power/type I error
final_var_factor<-factor(unlist(var_select), 
                         levels = c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10",
                                    "z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                                    "z11","z12","z13","z14","z15","z16","z17","z18","z19","z20",
                                    "z21","z22","z23","z24","z25","z26","z27","z28","z29","z30",
                                    "z31","z32","z33","z34","z35","z36","z37","z38","z39","z40"))
metrics_mat[,c("power_x1","power_x2","power_x3","power_x4","power_x5",
               "power_x6","power_x7","power_x8","power_x9","power_x10",
               "type1err_z1","type1err_z2","type1err_z3","type1err_z4","type1err_z5","type1err_z6","type1err_z7","type1err_z8","type1err_z9","type1err_z10",
               "type1err_z11","type1err_z12","type1err_z13","type1err_z14","type1err_z15","type1err_z16","type1err_z17","type1err_z18","type1err_z19","type1err_z20",
               "type1err_z21","type1err_z22","type1err_z23","type1err_z24","type1err_z25","type1err_z26","type1err_z27","type1err_z28","type1err_z29","type1err_z30",
               "type1err_z31","type1err_z32","type1err_z33","type1err_z34","type1err_z35","type1err_z36","type1err_z37","type1err_z38","type1err_z39","type1err_z40")]<-
  matrix(rep(table(final_var_factor)/n_rep, n_rep), ncol=50, byrow = T)
metrics_mat$power_overall<-apply(metrics_mat[,c("power_x1","power_x2","power_x3","power_x4","power_x5",
                                                "power_x6","power_x7","power_x8","power_x9","power_x10")],1,mean)
metrics_mat$type1err_overall<-apply(metrics_mat[,c("type1err_z1","type1err_z2","type1err_z3","type1err_z4","type1err_z5","type1err_z6",
                                                   "type1err_z7","type1err_z8","type1err_z9","type1err_z10",
                                                   "type1err_z11","type1err_z12","type1err_z13","type1err_z14","type1err_z15","type1err_z16","type1err_z17","type1err_z18","type1err_z19","type1err_z20",
                                                   "type1err_z21","type1err_z22","type1err_z23","type1err_z24","type1err_z25","type1err_z26","type1err_z27","type1err_z28","type1err_z29","type1err_z30",
                                                   "type1err_z31","type1err_z32","type1err_z33","type1err_z34","type1err_z35","type1err_z36","type1err_z37","type1err_z38","type1err_z39","type1err_z40")],1,mean)


write.csv(metrics_mat, paste(path,"x", seedx, "metrics_mat_rrbart_tdist_102721.csv",sep=""))
Sys.time() # from 4/23 17:41 (n_rep=200 *40) --> 23:54 6 hr