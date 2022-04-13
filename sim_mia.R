rm(list=ls())

library(splus2R)
library(foreach)
library(doParallel)
library(doSNOW)
library(pROC)
library(mice)
library(missForest)
library(rJava)
options(java.parameters = "-Xmx25g")
library(bartMachine)

path<-"/Users/linjungyi/Desktop/"
source("/Users/linjungyi/Google Drive/SWAN Study/Codes and R output/Paper 2/Organized R code/function.R")

# Set parameters
## Time: 3:22, 22 min
n_rep<-250 # Num of replications (250 for n = 1000 and 5000 and 1000 for n = 300 and 650)
n_impute<-1 # Num of imputation (1 to be)
maxit_impute<-5 # (5 to be)


# Sim data
n<-1000 # 650, 300, 5000

#setup parallel backend to use many processors
seedx<-123 
set.seed(seedx)
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

############################################################################################
# Exclude incompleted Y + BART #################################################################################################################################################################
cl <- makeCluster(cores[1]-2) #not to overload your computer
#cl <- makeCluster(25) #not to overload your computer
registerDoParallel(cl)

metrics_mat_local<-data.frame(matrix(NA, nrow=1, ncol=4))
names(metrics_mat_local)<-c("Precision","Recall","F1","AUC")

finalMatrix<-NA
finalMatrix <- foreach(i=1:n_rep, .combine="comb", .multicombine = TRUE,
                       .init=list(list(), list(), list(), list(), list())) %dopar% {
                                    library(caret)
                                    library(pROC)
                                    library(rJava)
                                    options(java.parameters = "-Xmx25g")
                                    library(bartMachine)

                                    dat_na0<-na_train_valid_dat[[1]][[i]]
                                    dat_na<-dat_na0[is.na(dat_na0$y)==0,]
                                    
                                    # Models
                                    cov_mat<-data.frame(dat_na[,-which(colnames(dat_na)%in%c("y"))])
                                    out_vec<-factor(as.character(dat_na$y), levels=c("Yes","No"))
                                    
                                    Sys.time() # c4 13:48
                                    bm1<-bartMachine(y = out_vec, X=cov_mat, mem_cache_for_speed=FALSE,
                                                     use_missing_data = TRUE, use_missing_data_dummies_as_covars = TRUE,
                                                     num_iterations_after_burn_in = 2000) # 1min
                                    var_select1<-var_selection_by_permute(bm1,plot=FALSE) 
                                    
                                    Sys.time()  # 2 min (nr10_nb5_maxit5)
                                    var_select_model<-var_select1
                                    
                                    # H0: combine the probabilities for each factor with different levels, need to be extremely careful here
                                    ## Replace () and space in colnames
                                    colnames(var_select1$permute_mat)<-gsub("\\(","_",colnames(var_select1$permute_mat))
                                    colnames(var_select1$permute_mat)<-gsub("\\)","_",colnames(var_select1$permute_mat))
                                    colnames(var_select1$permute_mat)<-gsub("\\ ","_",colnames(var_select1$permute_mat))
                                    colnames(var_select1$permute_mat)<-gsub("\\$","",colnames(var_select1$permute_mat))
                                    colnames(var_select1$permute_mat)<-gsub("\\,","",colnames(var_select1$permute_mat))
                                    colnames(var_select1$permute_mat)<-gsub("<","",colnames(var_select1$permute_mat))
                                    colnames(var_select1$permute_mat)<-gsub("\\=","",colnames(var_select1$permute_mat))
                                    colnames(var_select1$permute_mat)<-gsub(">","",colnames(var_select1$permute_mat))
                                    colnames(var_select1$permute_mat)<-gsub("\\/","",colnames(var_select1$permute_mat))
                                    colnames(var_select1$permute_mat)<-gsub("-","",colnames(var_select1$permute_mat))
                                    colnames(var_select1$permute_mat)<-gsub("\\,","",colnames(var_select1$permute_mat))
                                    colnames(var_select1$permute_mat)<-gsub("'","",colnames(var_select1$permute_mat))
                                    colnames(var_select1$permute_mat)<-gsub("\\+","",colnames(var_select1$permute_mat))
                                    colnames(var_select1$permute_mat)<-gsub("&","",colnames(var_select1$permute_mat))
                                    
                                    library(dplyr)
                                    vs_permutate_factor<-var_select1$permute_mat %>%  
                                      as_tibble() %>% 
                                      transmute(
                                        x1=x1_0+x1_1, x2=x2_0+x2_1, x3=x3,
                                        x4=x4, x5=x5,
                                        x6=x6, x7=x7+M_x7, x8=x8+M_x8, x9=x9+M_x9, x10=x10+M_x10,
                                        z1=z1, z2=z2, z3=z3, z4=z4, z5=z5,
                                        z6=z6, z7=z7, z8=z8, z9=z9, z10=z10,
                                        z11=z11, z12=z12, z13=z13, z14=z14, z15=z15,
                                        z16=z16, z17=z17, z18=z18, z19=z19, z20=z20,
                                        z21=z21_0+z21_1, z22=z22_0+z22_1, z23=z23_0+z23_1, z24=z24_0+z24_1, z25=z25_0+z25_1,
                                        z26=z26_0+z26_1, z27=z27_0+z27_1, z28=z28_0+z28_1, z29=z29_0+z29_1, z30=z30_0+z30_1,
                                        z31=z31_0+z31_1, z32=z32_0+z32_1, z33=z33_0+z33_1, z34=z34_0+z34_1, z35=z35_0+z35_1,
                                        z36=z36_0+z36_1, z37=z37_0+z37_1, z38=z38_0+z38_1, z39=z39_0+z39_1, z40=z40_0+z40_1)
                                    
                                    # Final check: if probabilities sum to 1 for all 100 permutations, then let's celebrate! Otherwise, we need to check whether we've missed any levels or factors....
                                    ave_check<-mean(vs_permutate_factor %>% apply(1, sum))
                                    #if(ave_check!=1){
                                    #   break
                                    #}
                                    
                                    # Prop included: combine the probabilities for each factor with different levels, need to be extremely careful here
                                    ## Replace () and space in colnames
                                    names(var_select1$var_true_props_avg)<-gsub("\\(","_",names(var_select1$var_true_props_avg))
                                    names(var_select1$var_true_props_avg)<-gsub("\\)","_",names(var_select1$var_true_props_avg))
                                    names(var_select1$var_true_props_avg)<-gsub("\\ ","_",names(var_select1$var_true_props_avg))
                                    names(var_select1$var_true_props_avg)<-gsub("\\$","",names(var_select1$var_true_props_avg))
                                    names(var_select1$var_true_props_avg)<-gsub("\\,","",names(var_select1$var_true_props_avg))
                                    names(var_select1$var_true_props_avg)<-gsub("<","",names(var_select1$var_true_props_avg))
                                    names(var_select1$var_true_props_avg)<-gsub("\\=","",names(var_select1$var_true_props_avg))
                                    names(var_select1$var_true_props_avg)<-gsub(">","",names(var_select1$var_true_props_avg))
                                    names(var_select1$var_true_props_avg)<-gsub("\\/","",names(var_select1$var_true_props_avg))
                                    names(var_select1$var_true_props_avg)<-gsub("-","",names(var_select1$var_true_props_avg))
                                    names(var_select1$var_true_props_avg)<-gsub("\\,","",names(var_select1$var_true_props_avg))
                                    names(var_select1$var_true_props_avg)<-gsub("'","",names(var_select1$var_true_props_avg))
                                    names(var_select1$var_true_props_avg)<-gsub("\\+","",names(var_select1$var_true_props_avg))
                                    names(var_select1$var_true_props_avg)<-gsub("&","",names(var_select1$var_true_props_avg))
                                    library(tibble)
                                    vs_prop_factor<-var_select1$var_true_props_avg %>%  
                                      as_tibble_row() %>% 
                                      transmute(
                                        x1=x1_0+x1_1, x2=x2_0+x2_1, x3=x3,
                                        x4=x4, x5=x5,
                                        x6=x6, x7=x7+M_x7, x8=x8+M_x8, x9=x9+M_x9, x10=x10+M_x10,
                                        z1=z1, z2=z2, z3=z3, z4=z4, z5=z5,
                                        z6=z6, z7=z7, z8=z8, z9=z9, z10=z10,
                                        z11=z11, z12=z12, z13=z13, z14=z14, z15=z15,
                                        z16=z16, z17=z17, z18=z18, z19=z19, z20=z20,
                                        z21=z21_0+z21_1, z22=z22_0+z22_1, z23=z23_0+z23_1, z24=z24_0+z24_1, z25=z25_0+z25_1,
                                        z26=z26_0+z26_1, z27=z27_0+z27_1, z28=z28_0+z28_1, z29=z29_0+z29_1, z30=z30_0+z30_1,
                                        z31=z31_0+z31_1, z32=z32_0+z32_1, z33=z33_0+z33_1, z34=z34_0+z34_1, z35=z35_0+z35_1,
                                        z36=z36_0+z36_1, z37=z37_0+z37_1, z38=z38_0+z38_1, z39=z39_0+z39_1, z40=z40_0+z40_1
                                      )
                                    # Final check: if probabilities sum to 1 for all 100 permutations, then let's celebrate! Otherwise, we need to check whether we've missed any levels or factors....
                                    ave_check2<-mean(vs_prop_factor %>% apply(1, sum))
                                    #if(ave_check2!=1){
                                    #   break
                                    #}
                                    detach("package:dplyr")
                                    # Var seletion (local): if prop > 1-alpha(0.05) quantils in permutation dist
                                    ## Caltulate cutoff
                                    localcutoff<-apply(vs_permutate_factor, 2, quantile, probs=0.90)
                                    
                                    ## Merge cutoff and vs_prop_factor by var names 
                                    cut_mat<-data.frame("localcutoff"=localcutoff, "varnames"=names(localcutoff))
                                    prop_mat<-data.frame("prop"=t(vs_prop_factor), "varnames"=names(vs_prop_factor))
                                    vs_mat<-merge(cut_mat, prop_mat, by="varnames", all=TRUE)
                                    
                                    vs_mat$local_results<-ifelse(vs_mat$prop>vs_mat$localcutoff,1,0)
                                    table(vs_mat$local_results) # 24 selected
                                    vs_mat$varnames[vs_mat$local_results==1]
                                    
                                    # Global max
                                    vs_mat$globalmaxcut<-quantile(apply(vs_permutate_factor,1,max), probs=0.90)
                                    vs_mat$globalmax_results<-ifelse(vs_mat$prop>vs_mat$globalmaxcut,1,0)
                                    
                                    # Global SE
                                    mk<-apply(vs_permutate_factor,2,mean) # k is number of cov
                                    sk<-apply(vs_permutate_factor,2,sd)
                                    cover_constant<-bisectK(tol = 0.01, coverage = 1 - 0.05, 
                                                            permute_mat = vs_permutate_factor, x_left = 1, x_right = 20, 
                                                            countLimit = 100, perm_mean = mk, perm_se = sk)
                                    globalse_cut<-mk+cover_constant*sk
                                    cut_mat2<-data.frame("globalse_cut"=globalse_cut, "varnames"=names(globalse_cut))
                                    vs_mat2<-merge(cut_mat2, vs_mat, by="varnames", all=TRUE)
                                    vs_mat2$globalse_results<-ifelse(vs_mat2$prop>vs_mat2$globalse_cut,1,0)
                                    
                                    
                                    final_var<-list("local_finalvar"=vs_mat2$varnames[vs_mat2$local_results==1],
                                                    "globalmax_finalvar"=vs_mat2$varnames[vs_mat2$globalmax_results==1],
                                                    "globalse_finalvar"=vs_mat2$varnames[vs_mat2$globalse_results==1])
                                    
                                    # 5. Metrics
                                    ## Precision, Recall, F1
                                    local_met<-prec_recall_f1(gold_pos=c("x1","x2","x3","x4","x5"),
                                                              gold_neg=c("z1","z2","z3","z4","z5","z6","z7","z8","z9","z10"),
                                                              var_selected=final_var$local_finalvar)
                                    
                                    prec_recall_f1_i<-list("local_met"=local_met)
                                    
                                    ## Combine metrics
                                    ## Precision, recall f1
                                    metrics_mat_local[,c("Precision","Recall","F1")]<-c(
                                      prec_recall_f1_i$local_met$precision, prec_recall_f1_i$local_met$recall, 
                                      prec_recall_f1_i$local_met$f1)
                                    
                                    ## AUC
                                    validdat_na0<-na_train_valid_dat[[2]][[i]]
                                    missfrst_valid<-missForest(validdat_na0)
                                    mice_valid_dat<-missfrst_valid$ximp
                                    
                                    cov_mat<-data.frame(mice_valid_dat[,final_var$local_finalvar, drop = F])
                                    out_vec<-factor(as.character(mice_valid_dat$y), levels=c("Yes","No"))
                                    bm_valid<-bartMachine(y = out_vec, X=cov_mat, mem_cache_for_speed=FALSE,
                                                                       use_missing_data = TRUE, use_missing_data_dummies_as_covars = TRUE)
                                    
                                    
                                    mice_valid_dat$phat<-predict(bm_valid, type="prob", new_data = cov_mat)
                                    metrics_mat_local[,"AUC"]<-auc(roc(response=mice_valid_dat$y, predictor=mice_valid_dat$phat, 
                                                                    levels=c("No", "Yes"), direction="<"))

                                    
                                    list(final_var, prec_recall_f1_i, metrics_mat_local,
                                         vs_prop_factor, vs_permutate_factor)
                                  }

stopCluster(cl)
Sys.time()
saveRDS(finalMatrix,  paste(path,"x",seedx,"finalMatrix_compy_MIAbart_110121.rds",sep=""))


final_var<-finalMatrix[[1]]
prec_recall_f1_i<-finalMatrix[[2]]
metrics_mat_local0<-finalMatrix[[3]]

metrics_mat_local<-data.frame(matrix(unlist(metrics_mat_local0), nrow=n_rep, byrow=T))
names(metrics_mat_local)<-c("Precision","Recall","F1","AUC")

############################################################################################
# Impute missing Y + BART#########################################################################
############################################################################################
cl <- makeCluster(cores[1]-2) #not to overload your computer
#cl <- makeCluster(25) #not to overload your computer
registerDoParallel(cl)

metrics_mat_local<-data.frame(matrix(NA, nrow=1, ncol=4))
names(metrics_mat_local)<-c("Precision","Recall","F1","AUC")

finalMatrix<-NA
finalMatrix <- foreach(i=1:n_rep, .combine="comb", .multicombine = TRUE,
                       .init=list(list(), list(), list(), list(), list())) %dopar% {
                         library(caret)
                         library(pROC)
                         library(rJava)
                         options(java.parameters = "-Xmx25g")
                         library(bartMachine)
                         library(missForest)
                         
                         dat_na0<-na_train_valid_dat[[1]][[i]]
                         missfrst_train<-missForest(dat_na0)
                         dat_na<-missfrst_train$ximp
                         
                         # Models
                         cov_mat<-data.frame(dat_na[,-which(colnames(dat_na)%in%c("y"))])
                         out_vec<-factor(as.character(dat_na$y), levels=c("Yes","No"))
                         
                         Sys.time() # c4 13:48
                         bm1<-bartMachine(y = out_vec, X=cov_mat, mem_cache_for_speed=FALSE,
                                          use_missing_data = TRUE, use_missing_data_dummies_as_covars = TRUE,
                                          num_iterations_after_burn_in = 2000) # 1min
                         var_select1<-var_selection_by_permute(bm1,plot=FALSE) 
                         
                         Sys.time()  # 2 min (nr10_nb5_maxit5)
                         var_select_model<-var_select1
                         
                         # H0: combine the probabilities for each factor with different levels, need to be extremely careful here
                         ## Replace () and space in colnames
                         colnames(var_select1$permute_mat)<-gsub("\\(","_",colnames(var_select1$permute_mat))
                         colnames(var_select1$permute_mat)<-gsub("\\)","_",colnames(var_select1$permute_mat))
                         colnames(var_select1$permute_mat)<-gsub("\\ ","_",colnames(var_select1$permute_mat))
                         colnames(var_select1$permute_mat)<-gsub("\\$","",colnames(var_select1$permute_mat))
                         colnames(var_select1$permute_mat)<-gsub("\\,","",colnames(var_select1$permute_mat))
                         colnames(var_select1$permute_mat)<-gsub("<","",colnames(var_select1$permute_mat))
                         colnames(var_select1$permute_mat)<-gsub("\\=","",colnames(var_select1$permute_mat))
                         colnames(var_select1$permute_mat)<-gsub(">","",colnames(var_select1$permute_mat))
                         colnames(var_select1$permute_mat)<-gsub("\\/","",colnames(var_select1$permute_mat))
                         colnames(var_select1$permute_mat)<-gsub("-","",colnames(var_select1$permute_mat))
                         colnames(var_select1$permute_mat)<-gsub("\\,","",colnames(var_select1$permute_mat))
                         colnames(var_select1$permute_mat)<-gsub("'","",colnames(var_select1$permute_mat))
                         colnames(var_select1$permute_mat)<-gsub("\\+","",colnames(var_select1$permute_mat))
                         colnames(var_select1$permute_mat)<-gsub("&","",colnames(var_select1$permute_mat))
                         
                         library(dplyr)
                         vs_permutate_factor<-var_select1$permute_mat %>%  
                           as_tibble() %>% 
                           transmute(
                             x1=x1_0+x1_1, x2=x2_0+x2_1, x3=x3,
                             x4=x4, x5=x5,
                             x6=x6, x7=x7, x8=x8, x9=x9, x10=x10,
                             z1=z1, z2=z2, z3=z3, z4=z4, z5=z5,
                             z6=z6, z7=z7, z8=z8, z9=z9, z10=z10,
                             z11=z11, z12=z12, z13=z13, z14=z14, z15=z15,
                             z16=z16, z17=z17, z18=z18, z19=z19, z20=z20,
                             z21=z21_0+z21_1, z22=z22_0+z22_1, z23=z23_0+z23_1, z24=z24_0+z24_1, z25=z25_0+z25_1,
                             z26=z26_0+z26_1, z27=z27_0+z27_1, z28=z28_0+z28_1, z29=z29_0+z29_1, z30=z30_0+z30_1,
                             z31=z31_0+z31_1, z32=z32_0+z32_1, z33=z33_0+z33_1, z34=z34_0+z34_1, z35=z35_0+z35_1,
                             z36=z36_0+z36_1, z37=z37_0+z37_1, z38=z38_0+z38_1, z39=z39_0+z39_1, z40=z40_0+z40_1)
                         
                         # Final check: if probabilities sum to 1 for all 100 permutations, then let's celebrate! Otherwise, we need to check whether we've missed any levels or factors....
                         ave_check<-mean(vs_permutate_factor %>% apply(1, sum))
                         #if(ave_check!=1){
                         #   break
                         #}
                         
                         # Prop included: combine the probabilities for each factor with different levels, need to be extremely careful here
                         ## Replace () and space in colnames
                         names(var_select1$var_true_props_avg)<-gsub("\\(","_",names(var_select1$var_true_props_avg))
                         names(var_select1$var_true_props_avg)<-gsub("\\)","_",names(var_select1$var_true_props_avg))
                         names(var_select1$var_true_props_avg)<-gsub("\\ ","_",names(var_select1$var_true_props_avg))
                         names(var_select1$var_true_props_avg)<-gsub("\\$","",names(var_select1$var_true_props_avg))
                         names(var_select1$var_true_props_avg)<-gsub("\\,","",names(var_select1$var_true_props_avg))
                         names(var_select1$var_true_props_avg)<-gsub("<","",names(var_select1$var_true_props_avg))
                         names(var_select1$var_true_props_avg)<-gsub("\\=","",names(var_select1$var_true_props_avg))
                         names(var_select1$var_true_props_avg)<-gsub(">","",names(var_select1$var_true_props_avg))
                         names(var_select1$var_true_props_avg)<-gsub("\\/","",names(var_select1$var_true_props_avg))
                         names(var_select1$var_true_props_avg)<-gsub("-","",names(var_select1$var_true_props_avg))
                         names(var_select1$var_true_props_avg)<-gsub("\\,","",names(var_select1$var_true_props_avg))
                         names(var_select1$var_true_props_avg)<-gsub("'","",names(var_select1$var_true_props_avg))
                         names(var_select1$var_true_props_avg)<-gsub("\\+","",names(var_select1$var_true_props_avg))
                         names(var_select1$var_true_props_avg)<-gsub("&","",names(var_select1$var_true_props_avg))
                         library(tibble)
                         vs_prop_factor<-var_select1$var_true_props_avg %>%  
                           as_tibble_row() %>% 
                           transmute(
                             x1=x1_0+x1_1, x2=x2_0+x2_1, x3=x3,
                             x4=x4, x5=x5,
                             x6=x6, x7=x7, x8=x8, x9=x9, x10=x10,
                             z1=z1, z2=z2, z3=z3, z4=z4, z5=z5,
                             z6=z6, z7=z7, z8=z8, z9=z9, z10=z10,
                             z11=z11, z12=z12, z13=z13, z14=z14, z15=z15,
                             z16=z16, z17=z17, z18=z18, z19=z19, z20=z20,
                             z21=z21_0+z21_1, z22=z22_0+z22_1, z23=z23_0+z23_1, z24=z24_0+z24_1, z25=z25_0+z25_1,
                             z26=z26_0+z26_1, z27=z27_0+z27_1, z28=z28_0+z28_1, z29=z29_0+z29_1, z30=z30_0+z30_1,
                             z31=z31_0+z31_1, z32=z32_0+z32_1, z33=z33_0+z33_1, z34=z34_0+z34_1, z35=z35_0+z35_1,
                             z36=z36_0+z36_1, z37=z37_0+z37_1, z38=z38_0+z38_1, z39=z39_0+z39_1, z40=z40_0+z40_1
                           )
                         # Final check: if probabilities sum to 1 for all 100 permutations, then let's celebrate! Otherwise, we need to check whether we've missed any levels or factors....
                         ave_check2<-mean(vs_prop_factor %>% apply(1, sum))
                         #if(ave_check2!=1){
                         #   break
                         #}
                         detach("package:dplyr")
                         # Var seletion (local): if prop > 1-alpha(0.05) quantils in permutation dist
                         ## Caltulate cutoff
                         localcutoff<-apply(vs_permutate_factor, 2, quantile, probs=0.90)
                         
                         ## Merge cutoff and vs_prop_factor by var names 
                         cut_mat<-data.frame("localcutoff"=localcutoff, "varnames"=names(localcutoff))
                         prop_mat<-data.frame("prop"=t(vs_prop_factor), "varnames"=names(vs_prop_factor))
                         vs_mat<-merge(cut_mat, prop_mat, by="varnames", all=TRUE)
                         
                         vs_mat$local_results<-ifelse(vs_mat$prop>vs_mat$localcutoff,1,0)
                         table(vs_mat$local_results) # 24 selected
                         vs_mat$varnames[vs_mat$local_results==1]
                         
                         # Global max
                         vs_mat$globalmaxcut<-quantile(apply(vs_permutate_factor,1,max), probs=0.90)
                         vs_mat$globalmax_results<-ifelse(vs_mat$prop>vs_mat$globalmaxcut,1,0)
                         
                         # Global SE
                         mk<-apply(vs_permutate_factor,2,mean) # k is number of cov
                         sk<-apply(vs_permutate_factor,2,sd)
                         cover_constant<-bisectK(tol = 0.01, coverage = 1 - 0.05, 
                                                 permute_mat = vs_permutate_factor, x_left = 1, x_right = 20, 
                                                 countLimit = 100, perm_mean = mk, perm_se = sk)
                         globalse_cut<-mk+cover_constant*sk
                         cut_mat2<-data.frame("globalse_cut"=globalse_cut, "varnames"=names(globalse_cut))
                         vs_mat2<-merge(cut_mat2, vs_mat, by="varnames", all=TRUE)
                         vs_mat2$globalse_results<-ifelse(vs_mat2$prop>vs_mat2$globalse_cut,1,0)
                         
                         
                         final_var<-list("local_finalvar"=vs_mat2$varnames[vs_mat2$local_results==1],
                                         "globalmax_finalvar"=vs_mat2$varnames[vs_mat2$globalmax_results==1],
                                         "globalse_finalvar"=vs_mat2$varnames[vs_mat2$globalse_results==1])
                         
                         # 5. Metrics
                         ## Precision, Recall, F1
                         local_met<-prec_recall_f1(gold_pos=c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"),
                                                   gold_neg=c("z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                                                              "z11","z12","z13","z14","z15","z16","z17","z18","z19","z20",
                                                              "z21","z22","z23","z24","z25","z26","z27","z28","z29","z30",
                                                              "z31","z32","z33","z34","z35","z36","z37","z38","z39","z40"),
                                                   var_selected=final_var$local_finalvar)
                         
                         prec_recall_f1_i<-list("local_met"=local_met)
                         
                         ## Combine metrics
                         ## Precision, recall f1
                         metrics_mat_local[,c("Precision","Recall","F1")]<-c(
                           prec_recall_f1_i$local_met$precision, prec_recall_f1_i$local_met$recall, 
                           prec_recall_f1_i$local_met$f1)
                         
                         ## AUC
                         validdat_na0<-na_train_valid_dat[[2]][[i]]
                         missfrst_valid<-missForest(validdat_na0)
                         mice_valid_dat<-missfrst_valid$ximp
                         
                         cov_mat<-data.frame(mice_valid_dat[,final_var$local_finalvar, drop = F])
                         out_vec<-factor(as.character(mice_valid_dat$y), levels=c("Yes","No"))
                         bm_valid<-bartMachine(y = out_vec, X=cov_mat, mem_cache_for_speed=FALSE,
                                               use_missing_data = TRUE, use_missing_data_dummies_as_covars = TRUE)
                         
                         
                         mice_valid_dat$phat<-predict(bm_valid, type="prob", new_data = cov_mat)
                         metrics_mat_local[,"AUC"]<-auc(roc(response=mice_valid_dat$y, predictor=mice_valid_dat$phat, 
                                                            levels=c("No", "Yes"), direction="<"))
                         
                         
                         list(final_var, prec_recall_f1_i, metrics_mat_local,
                              vs_prop_factor, vs_permutate_factor)
                       }

stopCluster(cl)
Sys.time()
saveRDS(finalMatrix,  paste(path,"x",seedx,"finalMatrix_missrf_MIAbart_110221.rds",sep=""))

############################################################################################
# Impute missing Y + XGB#############################################################
############################################################################################
cl <- makeCluster(cores[1]-2) #not to overload your computer
#cl <- makeCluster(32) #not to overload your computer
registerDoParallel(cl)

metrics_mat<-data.frame(matrix(NA, nrow=1, ncol=4))
names(metrics_mat)<-c("Precision","Recall","F1","AUC")


finalMatrix <- foreach(i=1:n_rep, .combine="comb", .multicombine = TRUE,
                       .init=list(list(), list(), list())) %dopar% {
                         library(caret)
                         library(pROC)
                         library(xgboost)
                         library(missForest)
                         
                         dat_na0<-na_train_valid_dat[[1]][[i]]
                         missfrst_train<-missForest(dat_na0)
                         dat_na<-missfrst_train$ximp
                         
                         
                         # Model
                         m1<-try(var_select_backward(datax=dat_na, event="y", method="xgb", number=5, repeats=1, 
                                                     var_drop_rate=vars.drop.frac_x,
                                                     nroundx=nroundx, c.sd=1, nthreadx = nthreadx, rcv_ce = rcv_ce))
                         if(inherits(m1, "try-error"))
                         {
                           m1 <- data.frame()
                           saveRDS(mice_b_dat, paste(path,"i",i,"dat_na_debug.rds",sep=""))
                           selected_var_b <- NA
                         } else{
                           selected_var_b<-m1$var_in_model[[1]]
                         }
                         
                         
                         # 5. Metrics
                         ## Precision, Recall, F1
                         prec_recall_f1_i<-prec_recall_f1(gold_pos=c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"),
                                                          gold_neg=c("z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                                                                     "z11","z12","z13","z14","z15","z16","z17","z18","z19","z20",
                                                                     "z21","z22","z23","z24","z25","z26","z27","z28","z29","z30",
                                                                     "z31","z32","z33","z34","z35","z36","z37","z38","z39","z40"),
                                                          var_selected=selected_var_b)
                         
                         metrics_mat[,c("Precision","Recall","F1")]<-c(
                           prec_recall_f1_i$precision, prec_recall_f1_i$recall, 
                           prec_recall_f1_i$f1)
                         
                         ## AUC
                         validdat_na0<-na_train_valid_dat[[2]][[i]]
                         missfrst_valid<-missForest(validdat_na0)
                         mice_valid_dat<-missfrst_valid$ximp
                         
                         out_vec<-ifelse(mice_valid_dat$y=="No",0,
                                         ifelse(mice_valid_dat$y=="Yes",1,NA))
                         cov_mat<-data.frame(mice_valid_dat[,selected_var_b])
                         ff<-~.
                         mf<-model.frame(formula = ff, data = cov_mat, na.action = "na.pass")
                         cov_mat_v1<-model.matrix(object = ff, data = mf)
                         
                         xgb_v1 <- xgboost(data = cov_mat_v1, label = out_vec, verbose = F, missing = NA, 
                                           objective = "binary:logistic", na.action = "na.pass", nthread=nthreadx, nrounds=nroundx)
                         mice_valid_dat$phat<-predict(xgb_v1, newdata = cov_mat_v1, missing = NA)
                         
                         metrics_mat[,"AUC"]<-pROC::auc(roc(response=mice_valid_dat$y, predictor=mice_valid_dat$phat, 
                                                            levels=c("No", "Yes"), direction="<"))
                         
                         
                         list("selected_var_b"=selected_var_b, prec_recall_f1_i, metrics_mat)
                       }


stopCluster(cl)
saveRDS(finalMatrix,  paste(path, "x", seedx,"finalMatrix_missrf_MIAxgb_110221y.rds",sep=""))
Sys.time() 

final_var<-finalMatrix[[1]]
prec_recall_f1_i<-finalMatrix[[2]]
metrics_mat<-finalMatrix[[3]]

metrics_mat<-data.frame(matrix(unlist(metrics_mat), nrow=n_rep, byrow=T))
names(metrics_mat)<-c("Precision","Recall","F1", "AUC")

############################################################################################
# Exclude missing Y + XGB#############################################################
############################################################################################
cl <- makeCluster(cores[1]-2) #not to overload your computer
#cl <- makeCluster(32) #not to overload your computer
registerDoParallel(cl)

metrics_mat<-data.frame(matrix(NA, nrow=1, ncol=4))
names(metrics_mat)<-c("Precision","Recall","F1","AUC")


finalMatrix <- foreach(i=1:n_rep, .combine="comb", .multicombine = TRUE,
                       .init=list(list(), list(), list())) %dopar% {
                         library(caret)
                         library(pROC)
                         library(xgboost)
                         library(missForest)
                         
                         dat_na0<-na_train_valid_dat[[1]][[i]]
                         dat_na<-dat_na0[is.na(dat_na0$y)==0,]
                         
                         # Model
                         m1<-try(var_select_backward(datax=dat_na, event="y", method="xgb", number=5, repeats=1, 
                                                     var_drop_rate=vars.drop.frac_x,
                                                     nroundx=nroundx, c.sd=1, nthreadx = nthreadx, rcv_ce = rcv_ce))
                         if(inherits(m1, "try-error"))
                         {
                           m1 <- data.frame()
                           saveRDS(mice_b_dat, paste(path,"i",i,"dat_na_debug.rds",sep=""))
                           selected_var_b <- NA
                         } else{
                           selected_var_b<-m1$var_in_model[[1]]
                         }
                         
                         
                         # 5. Metrics
                         ## Precision, Recall, F1
                         prec_recall_f1_i<-prec_recall_f1(gold_pos=c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"),
                                                          gold_neg=c("z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                                                                     "z11","z12","z13","z14","z15","z16","z17","z18","z19","z20",
                                                                     "z21","z22","z23","z24","z25","z26","z27","z28","z29","z30",
                                                                     "z31","z32","z33","z34","z35","z36","z37","z38","z39","z40"),
                                                          var_selected=selected_var_b)
                         
                         metrics_mat[,c("Precision","Recall","F1")]<-c(
                           prec_recall_f1_i$precision, prec_recall_f1_i$recall, 
                           prec_recall_f1_i$f1)
                         
                         ## AUC
                         validdat_na0<-na_train_valid_dat[[2]][[i]]
                         missfrst_valid<-missForest(validdat_na0)
                         mice_valid_dat<-missfrst_valid$ximp
                         
                         out_vec<-ifelse(mice_valid_dat$y=="No",0,
                                         ifelse(mice_valid_dat$y=="Yes",1,NA))
                         cov_mat<-data.frame(mice_valid_dat[,selected_var_b])
                         ff<-~.
                         mf<-model.frame(formula = ff, data = cov_mat, na.action = "na.pass")
                         cov_mat_v1<-model.matrix(object = ff, data = mf)
                         
                         xgb_v1 <- xgboost(data = cov_mat_v1, label = out_vec, verbose = F, missing = NA, 
                                           objective = "binary:logistic", na.action = "na.pass", nthread=nthreadx, nrounds=nroundx)
                         mice_valid_dat$phat<-predict(xgb_v1, newdata = cov_mat_v1, missing = NA)
                         
                         metrics_mat[,"AUC"]<-pROC::auc(roc(response=mice_valid_dat$y, predictor=mice_valid_dat$phat, 
                                                            levels=c("No", "Yes"), direction="<"))
                         
                         
                         list("selected_var_b"=selected_var_b, prec_recall_f1_i, metrics_mat)
                       }


stopCluster(cl)
saveRDS(finalMatrix,  paste(path, "x", seedx,"finalMatrix_compy_MIAxgb_110221y.rds",sep=""))
Sys.time() 

final_var<-finalMatrix[[1]]
prec_recall_f1_i<-finalMatrix[[2]]
metrics_mat<-finalMatrix[[3]]

metrics_mat<-data.frame(matrix(unlist(metrics_mat), nrow=n_rep, byrow=T))
names(metrics_mat)<-c("Precision","Recall","F1", "AUC")

## Power & type I error
### local
final_var_factor_local<-factor(unlist(lapply(final_var, function(x){x$local_finalvar})), 
                               levels = c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10",
                                          "z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                                          "z11","z12","z13","z14","z15","z16","z17","z18","z19","z20",
                                          "z21","z22","z23","z24","z25","z26","z27","z28","z29","z30",
                                          "z31","z32","z33","z34","z35","z36","z37","z38","z39","z40"))
metrics_mat_local[,c("power_x1","power_x2","power_x3","power_x4","power_x5",
                     "power_x6","power_x7","power_x8","power_x9","power_x10",
                     "type1err_z1","type1err_z2","type1err_z3","type1err_z4","type1err_z5","type1err_z6","type1err_z7","type1err_z8","type1err_z9","type1err_z10",
                     "type1err_z11","type1err_z12","type1err_z13","type1err_z14","type1err_z15","type1err_z16","type1err_z17","type1err_z18","type1err_z19","type1err_z20",
                     "type1err_z21","type1err_z22","type1err_z23","type1err_z24","type1err_z25","type1err_z26","type1err_z27","type1err_z28","type1err_z29","type1err_z30",
                     "type1err_z31","type1err_z32","type1err_z33","type1err_z34","type1err_z35","type1err_z36","type1err_z37","type1err_z38","type1err_z39","type1err_z40")]<-
  matrix(rep(table(final_var_factor_local)/n_rep, n_rep), ncol=50, byrow = T)
metrics_mat_local$power_overall<-apply(metrics_mat_local[,c("power_x1","power_x2","power_x3","power_x4","power_x5",
                                                            "power_x6","power_x7","power_x8","power_x9","power_x10")],1,mean)
metrics_mat_local$type1err_overall<-apply(metrics_mat_local[,c("type1err_z1","type1err_z2","type1err_z3","type1err_z4","type1err_z5","type1err_z6","type1err_z7","type1err_z8","type1err_z9","type1err_z10",
                                                               "type1err_z11","type1err_z12","type1err_z13","type1err_z14","type1err_z15","type1err_z16","type1err_z17","type1err_z18","type1err_z19","type1err_z20",
                                                               "type1err_z21","type1err_z22","type1err_z23","type1err_z24","type1err_z25","type1err_z26","type1err_z27","type1err_z28","type1err_z29","type1err_z30",
                                                               "type1err_z31","type1err_z32","type1err_z33","type1err_z34","type1err_z35","type1err_z36","type1err_z37","type1err_z38","type1err_z39","type1err_z40")],1,mean)

write.csv(metrics_mat_local, paste(path,"x",seedx,"metrics_mat_compy_MIAbart_110221.csv",sep=""))

Sys.time() # 19:32-->19:34 2 MIN (n_rep=5)