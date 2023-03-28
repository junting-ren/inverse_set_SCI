#####
# function for method 2 of simulating data

ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}

simulateData <- function(M, N.ipd, betas.sim,SIGMA) {
  MAF = rep(0.1, M)
  #MAF <- seq(from=.3, to=.1, length.out = M)
  # SIGMA <- createDistAR(M, order=2)*First.Cor
  # diag(SIGMA) <- 1
  c0 <- qnorm((1-MAF)^2)
  c2 <- qnorm(1-MAF^2)
  # individual person data (ipd)
  G.ipd <-mvtnorm::rmvnorm(N.ipd, sigma=SIGMA)
  G.ipd <- lapply(1:M, FUN=function(m) { ifelse(G.ipd[,m] > c2[m], 2, ifelse(G.ipd[,m]<c0[m], 0,1)) })
  G.ipd <- do.call(cbind, G.ipd)
  Y.ipd <- rnorm(N.ipd, mean=G.ipd%*%betas.sim, sd=1)
  betas.ipd.marg <- unlist(lapply(1:M, FUN=function(m) { summary(lm(Y.ipd ~ G.ipd[,m]))$coef[2,1]}))
  se.ipd.marg <- unlist(lapply(1:M, FUN=function(m) { summary(lm(Y.ipd ~ G.ipd[,m]))$coef[2,2]}))
  # maf.ipd.marg <- colMeans(G.ipd)
  # betas <- betas.ipd.marg
  # r.ipd.joint <- summary(lm(Y.ipd ~ G.ipd[,1:5]))$coef
  # betas.ipd.joint <- r.ipd.joint[2:(1+5),1]
  
  # # reference data
  # G.u <-rmvnorm(N.ref, sigma=SIGMA) 
  # W <- lapply(1:M, FUN=function(m) { ifelse(G.u[,m] > c2[m], 2, ifelse(G.u[,m]<c0[m], 0,1)) })
  # W <- do.call(cbind, W)
  return(list(Y=Y.ipd, X =G.ipd, betas_margin = betas.ipd.marg, sd_margin = se.ipd.marg)) # Edit: Return trait variance
}


# library(tidyverse)
# 
# ########################################################
# # Section 0: Specify parameters for the simulation 
# ########################################################
# 
# # Define causal IDs and correlation structure IDs to iterate
# N_SNP = 100    # Total number of SNPs
# block_size = 10    # Block size of LD matrix
# 
# # Initiate the Causal index
# N_Causal = 1
# if(N_Causal != 0){
#   Causal_ID = seq(1,N_SNP,by = block_size)[1:N_Causal]
# }else{
#   Causal_ID = NULL
# }
# 
# # Initiate the number of studies per ethnic group/population
# N_Study_PerPop = 1
# Cor_ID = c(rep(1,N_Study_PerPop),rep(2,N_Study_PerPop),rep(3,N_Study_PerPop))
# 
# # Specify LD correlation 
# LD_corr = 0.9
# 
# # Number of subjects of each population
# Sample_Size_Scenario = "Balanced"
# 
# if(Sample_Size_Scenario == "Balanced"){
#   N_Sample_Per_Study = as.integer(15000/N_Study_PerPop)
#   N_Sample = c(rep(N_Sample_Per_Study,N_Study_PerPop),rep(N_Sample_Per_Study,N_Study_PerPop),rep(N_Sample_Per_Study,N_Study_PerPop))
# }else{
#   if(N_Study_PerPop == 3){
#     N_Sample = c(rep(11000,3), rep(2000,3),  rep(2000,3))
#   }else{
#     stop("N_Study_PerPop != 3. Please change N_Sample so that they sum up to 45,000.")
#   }
# }
# 
# 
# # Number of replications 
# N_Rep = 500
# # Starting point of replication
# Start_Rep = 1
# 
# # Set the effect size 
# Effect_size = 0.01
# 
# # Set the p-value cutoff for index selection 
# ## final cutoff would be Bonp/N_SNP
# Bonp <-  0.05
# 
# # The name of populations 
# N_Study = length(Cor_ID)
# StudyName = paste("S",c(1:N_Study),sep="")
# 
# 
# 
# ##############################################################################
# # Section 1: This script is to simulate dosage and generate marginal results
# ##############################################################################
# # Set desired MAF and compute the cut-off from a standard normal distribution
# MAF = c(rep(0.2,N_Study_PerPop), rep(0.4,N_Study_PerPop), rep(0.6, N_Study_PerPop))
# Cut_L = qnorm((1-MAF)^2)
# Cut_U = qnorm(1-MAF^2)
# # Define the three correlation structures
# # CAUTION: the length of each correlation vector should be N_SNP and the first element should be 1.
# Cor1 = rep(c(0,rep(LD_corr,block_size-1)),N_SNP/block_size)
# Cor2 = rep(c(0,rep(LD_corr,block_size-1)),N_SNP/block_size)
# Cor3 = rep(c(0,rep(LD_corr,block_size-1)),N_SNP/block_size)
# # Cor1 = Cor2 = Cor3 =rep(0,N_SNP)
# # if(N_LD>0){
# #   for(i in 1:length(Causal_ID)){
# #     Cor1[(Causal_ID[i]+1):(Causal_ID[i]+N_LD)] = LD_corr
# #     Cor2[(Causal_ID[i]+1):(Causal_ID[i]+N_LD)] = LD_corr
# #     Cor3[(Causal_ID[i]+1):(Causal_ID[i]+N_LD)] = LD_corr
# #   }
# # }
# 
# 
# # CAUTION!
# # Construct a conversion vector that stores the beta required in the linear regression to generate SNPs with desired correlation
# # Desired correlationis a vector from 0 to 0.9 with 0.1 space between each pair
# CortoBeta_Cor = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
# CortoBeta_Beta = c(0, 0.25 ,0.4 ,0.6 ,0.84 ,1.1 ,1.4 ,1.8 , 2.5, 4, 8)
# set.seed(2022)
# ##################################
# # Iterate through all replications
# ##################################
# # r2_cont_matrix <- matrix(NA, nrow = N_Rep, ncol = N_Study)
# for (k in Start_Rep:(Start_Rep+N_Rep-1)){
#   print(paste("Now running: marginal results replication ",k,sep=""))
#   Dosage_FileName = paste("Dosage_",StudyName,"_",k,".txt",sep="")
#   Y_FileName = paste("Y_",StudyName,"_",k,".txt",sep="")
#   Marg_FileName = paste0("Marg_",k,".txt")
#   MAF_FileName = paste0("MAF_",k,".txt")
#   # Initialize the marginal output with SNP names
#   SNP_Name = paste("rs",c(1:N_SNP),sep="")
#   Marg_Output <-  MAF_Output <-  SNP_Name
#   #################################
#   # Iterate through all populations
#   #################################
#   for (e in 1:N_Study){
#     
#     ##### METHOD 1 of simulating data ######
#     ## Repeatedly generate dosages until the cholesky decomposition is possible
#     itr = 0
#     repeat{
#       itr = itr + 1
#       # print(itr)
#       # Generate SNP_Ref from normal distribution for this specific populations
#       SNP_Ref = rnorm(N_Sample[e],mean=0,sd=1)
#       SNP_Ref = (SNP_Ref >= Cut_U[e]) + (SNP_Ref >= Cut_L[e])
#       # Generate all other SNPs with desired correlation structure
#       Dosage = SNP_Ref
#       for (i in 2:N_SNP){
#         if (eval(parse(text = paste("Cor",Cor_ID[e],sep="")))[i] == 0){
#           # Re-generate reference SNP
#           SNP_Ref = rnorm(N_Sample[e],mean=0,sd=1)
#           SNP_Ref = (SNP_Ref >= Cut_U[e]) + (SNP_Ref >= Cut_L[e])
#           Dosage = cbind(Dosage,SNP_Ref)
#         }
#         if (eval(parse(text = paste("Cor",Cor_ID[e],sep="")))[i]> 0){
#           temp_SNP = rnorm(N_Sample[e],mean=CortoBeta_Beta[which(CortoBeta_Cor==eval(parse(text = paste("Cor",Cor_ID[e],sep="")))[i])]*(SNP_Ref-mean(SNP_Ref)),1)
#           temp_SNP = (temp_SNP >= Cut_U[e]) + (temp_SNP >= Cut_L[e])
#           Dosage = cbind(Dosage,temp_SNP)
#         }
#       }
#       colnames(Dosage) = SNP_Name
#       temp_Chol = NULL
#       temp_Chol = try(chol(t(Dosage)%*%Dosage))
#       # temp_corr = stats::cor(Dosage)
#       # if (!is.null(dim(temp_Chol)) & !is.singular.matrix(temp_corr)){
#       if (!is.null(dim(temp_Chol))){
#         #print("Great! Found Dosage with valid cholesky decomposition!")
#         break
#       }
#     }
#     ##### METHOD 2 of simulating data ######
#     # sigma <- matrix(0,ncol = N_SNP,nrow = N_SNP)
#     # sigma[Causal_ID:(Causal_ID+N_LD),Causal_ID:(Causal_ID+N_LD)] <- matrix(0.75, nrow = N_LD+1, ncol = N_LD+1)
#     # diag(sigma) <- rep(1, N_SNP)
#     # AllData <- simulateData(M = N_SNP, N.ref = N_Sample[e], N.ipd = N_Sample[e],
#     #                         CausalEffects = Effect_size,SIGMA = sigma)
#     # Dosage <- AllData$W
#     #######
#     # Output the dosage to desired location
#     write.table(Dosage,paste(Dosage_Dir,Dosage_FileName[e],sep=""),quote=F,sep="\t",row.name=F,col.name=T,append=F)
# 
#     ######################
#     ## --- Simulate continuous outcome based on selected causal SNP
#     Beta = rep(0, N_SNP)
#     Beta[Causal_ID] = rep(Effect_size, length(Causal_ID))
#     X_cent <- scale(Dosage, center = TRUE, scale = FALSE)
#     Y <- rnorm(N_Sample[e], mean=X_cent%*%Beta, sd = 1)
#     Y_cent <- scale(Y, center = TRUE, scale = FALSE)
# 
#     # Output the dosage to desired location
#     write.table(Y_cent,paste(Dosage_Dir,Y_FileName[e],sep=""),quote=F,sep="\t",row.name=F,col.name=F,append=F)
# 
#     ## --- Get summary statistics
#     temp_Marg_Output <- susieR::univariate_regression(X_cent, Y)
#     temp_MAF_Output <- colMeans(Dosage)/2
#     Marg_Output <- cbind(Marg_Output, temp_Marg_Output$betahat, temp_Marg_Output$sebetahat)
#     MAF_Output <- cbind(MAF_Output, temp_MAF_Output)
# 
#     # Marg_p = 2*pnorm(abs( temp_Marg_Output$betahat/temp_Marg_Output$sebetahat), lower.tail = F)
#   }
#   # Output the Marginal results for all populations
#   write.table(Marg_Output,paste0(Marg_Dir,Marg_FileName),quote=F,sep="\t",row.name=F,col.name=F,append=F)
#   write.table(MAF_Output,paste0(Marg_Dir,MAF_FileName),quote=F,sep="\t",row.name=F,col.name=F,append=F)
# }
# # write.table(r2_cont_matrix,paste0(Analysis_Output_Dir,"r2_whole_matrix.txt"),quote=F,sep="\t",row.name=F,col.name=F,append=F)

