rm(list=ls())

#library(here)
library(splus2R)
library(foreach)
library(doParallel)
library(mice)
library(xgboost)
library(missForest)
library(pROC)

path<-"/Users/linjungyi/Desktop/"
source("/Users/linjungyi/Google Drive/SWAN Study/Codes and R output/Paper 2/Organized R code/function.R")

# Replications
## 1. Simulate X
## 2. Simulate Y
## 3. Ampute
## 4. Impute and var selection
## 5. Calculate metrics for final model
## 6. Repeate steps 2-4 1000 times

# Set parameters
## Time: 3:22, 22 min
n_rep<-250 # Num of replications (250 for n = 1000 and 5000 and 1000 for n = 300 and 650)
n_boot_impute<-100 # Num of bootstraps to be imputed (100 to be)
n_impute<-1 # Num of imputation (1 to be)
maxit_impute<-5 # (5 to be)
nroundx<-200 # (200 to be)
nthreadx<-1
vars.drop.frac_x<-0.1
rcv_ce<-FALSE
cutpt_boot_varselect<-c(1/n_boot_impute, seq(0.1, 1.0, 0.1))

# Output
var_select_model<-list(NA)
final_var<-list(NA)

# Sim data
n<-1000 # 650, 300, 5000

#setup parallel backend to use many processors
seedx<-123 
set.seed(seedx)
cores=detectCores()

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

# Impute+model
cl <- makeCluster(cores[1]-2) #not to overload your computer

registerDoParallel(cl)
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

impute_boot<-foreach(i=1:n_rep, .combine="comb", .multicombine = TRUE,  .init=list(list(), list()), 
                     .packages = c("doParallel", "doSNOW")) %:%
  
  foreach(b=1:n_boot_impute, .combine="comb", .multicombine = TRUE, .init=list(list(), list()),
          .packages = c("doParallel")) %dopar%{
            library(caret)
            library(pROC)
            library(mice)
            library(missForest)
            library(xgboost)
            
            dat_na_i<-na_traindat_list[[i]]
            # Bootstarp
            id_boot<-sample(c(1:nrow(dat_na_i)), nrow(dat_na_i), replace=TRUE)
            dat_b<-dat_na_i[c(id_boot),]
            
            # Impute
            missfrst_b<-missForest(dat_b)
            mice_b_dat<-missfrst_b$ximp
            
            # Model
            m1<-try(var_select_backward(datax=mice_b_dat, event="y", method="xgb", 
                                        var_drop_rate=vars.drop.frac_x,
                                        nroundx=nroundx, c.sd=1, nthreadx = nthreadx, rcv_ce = rcv_ce))
            if(inherits(m1, "try-error"))
            {
              m1 <- data.frame()
              saveRDS(mice_b_dat, paste(path,"i",i,"b",b,"xgb_mice_b_dat_debug.rds",sep=""))
              selected_var_b <- NA
            } else{
              selected_var_b<-m1$var_in_model[[1]]
            }

            list("selected_var_b"=selected_var_b, m1)
          }
stopCluster(cl)
saveRDS(impute_boot,  paste(path,"x",seedx,"impute_boot_missrfxgb_102721.rds",sep=""))

# Metric
final_var0<-list(NA)
prec_recall_f1_met<-list(NA)
metrics_mat0<-list(NA)

for(i in 1:n_rep) {
  # Select final vars
  selected_var_i<-impute_boot[[1]][[i]]
  select_freq_i<-table(unlist(selected_var_i))
  final_var_i<-boot_finalvar(cutpt_boot_varselect=cutpt_boot_varselect, 
                             varlist=selected_var_i, n_boot=n_boot_impute)
  
  # Metrics
  ## Impute valid data for AUC
  missfrst_valid<-missForest(na_validdat_list[[i]])
  mice_valid_dat<-missfrst_valid$ximp
  
  prec_recall_f1_met_i<-list(NA)
  metrics_mat_i<-data.frame(matrix(NA, nrow=length(cutpt_boot_varselect), ncol=5))
  names(metrics_mat_i)<-c("Cutpt","Precision","Recall","F1", "AUC")
  
  for(j in 1:length(cutpt_boot_varselect)){
    final_var_j<-final_var_i[[j]]
    ## Precision, Recall, F1
    prec_recall_f1_met_i[[j]]<-prec_recall_f1(gold_pos=c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"),
                                              gold_neg=c("z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                                                         "z11","z12","z13","z14","z15","z16","z17","z18","z19","z20",
                                                         "z21","z22","z23","z24","z25","z26","z27","z28","z29","z30",
                                                         "z31","z32","z33","z34","z35","z36","z37","z38","z39","z40"),
                                              var_selected=final_var_j)
    
    ## Precision, recall f1
    metrics_mat_i[j,c(1:4)]<-c(cutpt_boot_varselect[j], prec_recall_f1_met_i[[j]]$precision, 
                         prec_recall_f1_met_i[[j]]$recall, prec_recall_f1_met_i[[j]]$f1)
    
    ## Fit model to calculate AUC
    out_vec<-ifelse(mice_valid_dat$y=="No",0,
                    ifelse(mice_valid_dat$y=="Yes",1,NA))
    cov_mat<-data.frame(mice_valid_dat[,final_var_i[[j]]])
    ff<-~.
    mf<-model.frame(formula = ff, data = cov_mat, na.action = "na.pass")
    cov_mat_v1<-model.matrix(object = ff, data = mf)
    
    xgb_v1 <- xgboost(data = cov_mat_v1, label = out_vec, verbose = F, missing = NA, 
                      objective = "binary:logistic", na.action = "na.pass", nthread=nthreadx, nrounds=nroundx)
    mice_valid_dat$phat<-predict(xgb_v1, newdata = cov_mat_v1, missing = NA)
    
    ## AUC
    metrics_mat_i[j,5]<-auc(roc(response=mice_valid_dat$y, predictor=mice_valid_dat$phat, 
                                    levels=c("No", "Yes"), direction="<"))
    
  }
  final_var0[[i]]<-final_var_i
  metrics_mat0[[i]]<-metrics_mat_i
  prec_recall_f1_met[[i]]<-prec_recall_f1_met_i
}


# Metrics
metrics_mat_list<-list(NA)
for(i in 1:length(cutpt_boot_varselect)){
  # Precision, recall, f1
  metrics_mat_list_i0<- lapply(metrics_mat0, function(x){x[i,]})
  metrics_mat_list_i<-data.frame(matrix(unlist(metrics_mat_list_i0), nrow=n_rep, byrow=T))
  names(metrics_mat_list_i)<-c("Cutpt","Precision","Recall","F1","AUC")
  
  ## Power & type I error
  final_var<-lapply(final_var0, function(x){x[[i]]})
  final_var_factor<-factor(unlist(final_var), 
                           levels = c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10",
                                      "z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                                      "z11","z12","z13","z14","z15","z16","z17","z18","z19","z20",
                                      "z21","z22","z23","z24","z25","z26","z27","z28","z29","z30",
                                      "z31","z32","z33","z34","z35","z36","z37","z38","z39","z40"))
  metrics_mat_list_i[,c("power_x1","power_x2","power_x3","power_x4","power_x5",
                        "power_x6","power_x7","power_x8","power_x9","power_x10",
                        "type1err_z1","type1err_z2","type1err_z3","type1err_z4","type1err_z5","type1err_z6","type1err_z7","type1err_z8","type1err_z9","type1err_z10",
                        "type1err_z11","type1err_z12","type1err_z13","type1err_z14","type1err_z15","type1err_z16","type1err_z17","type1err_z18","type1err_z19","type1err_z20",
                        "type1err_z21","type1err_z22","type1err_z23","type1err_z24","type1err_z25","type1err_z26","type1err_z27","type1err_z28","type1err_z29","type1err_z30",
                        "type1err_z31","type1err_z32","type1err_z33","type1err_z34","type1err_z35","type1err_z36","type1err_z37","type1err_z38","type1err_z39","type1err_z40")]<-
    matrix(rep(table(final_var_factor)/n_rep, n_rep), ncol=50, byrow = T)
  metrics_mat_list_i$power_overall<-apply(metrics_mat_list_i[,c("power_x1","power_x2","power_x3","power_x4","power_x5",
                                                                "power_x6","power_x7","power_x8","power_x9","power_x10")],1,mean)
  metrics_mat_list_i$type1err_overall<-apply(metrics_mat_list_i[,c("type1err_z1","type1err_z2","type1err_z3","type1err_z4","type1err_z5","type1err_z6",
                                                                   "type1err_z7","type1err_z8","type1err_z9","type1err_z10",
                                                                   "type1err_z11","type1err_z12","type1err_z13","type1err_z14","type1err_z15","type1err_z16","type1err_z17","type1err_z18","type1err_z19","type1err_z20",
                                                                   "type1err_z21","type1err_z22","type1err_z23","type1err_z24","type1err_z25","type1err_z26","type1err_z27","type1err_z28","type1err_z29","type1err_z30",
                                                                   "type1err_z31","type1err_z32","type1err_z33","type1err_z34","type1err_z35","type1err_z36","type1err_z37","type1err_z38","type1err_z39","type1err_z40")],1,mean)
  
  
  
  
  metrics_mat_list[[i]]<-metrics_mat_list_i
}
metrics_mat_list_var<-list("metrics"=metrics_mat_list_i, "final_var"=final_var0, 
                       "prec_recall_f1_met"=prec_recall_f1_met)
# Export
saveRDS(metrics_mat_list_var, paste(path,"x",seedx, "metrics_mat_missrfxgb_102721_n1k_na06.rds",sep=""))
Sys.time() # from 1:25, 4 core, nthread=1, n_rep=8, n_impute=10, 10 min