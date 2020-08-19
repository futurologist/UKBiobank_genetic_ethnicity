library(data.table)
library(dplyr)
library(ggplot2)

################# Functions: ##############################################

# UKB_all_AFR_3PCs[order(Norm_stats_Black), ][Black==1, ][1,]

ethn_relabel <- function(tbl_in, lbl, eth){
  relbl <- function(ntry){ return(eth[ paste(ntry) ]) }
  tbl_in <- data.table(tbl_in, sapply(tbl_in[, ..lbl], relbl))
  colnames(tbl_in)[ncol(tbl_in)] <- 'Ethnicity' 
  return(tbl_in)
}


matrix_of_means <- function(X){
  return( matrix(rep(1, length = nrow(X)), ncol = 1) %*% colMeans(X) )
}

covar_matrix <- function(X, M){
  X <- X - M
  return( (t(X) %*% X)/(nrow(X) - 1) )
}

calc_ethn_score <- function(GD, KG, ethn_list, n_PCs){
  
  PC_GD <- which(startsWith(names(GD), "PC"))
  PC_GD <- PC_GD[c(1: min(length(PC_GD), (n_PCs)))]
  PC <- which(startsWith(names(KG), "PC"))
  PC <- PC[c(1: min(length(PC), (n_PCs)))]
  
  Z_GD <- as.matrix(GD[ , ..PC_GD])
  
  for(etn in ethn_list){
    Z <- as.matrix(KG[Ethnicity == etn, ..PC])
    M <- matrix_of_means(Z)
    Sigma_inv <- solve(covar_matrix(Z, M))
    
    M <- t(M[1,])
    
    score <- function(r){ return( (r-M) %*% Sigma_inv %*% t(r-M) ) }
    
    GD <- data.table( GD, Norm_stats = apply(Z_GD, 1, score) )
    colnames(GD)[ncol(GD)] <- paste('Ethn_Score_', etn, sep='') 
  }
  return(GD)  
  
} 


remove_outlier <- function(KGP_, n_PCs){
  KGP_ <- calc_ethn_score(KGP_, KGP_, c('Black'), n_PCs)
  n_b <- nrow(KGP_[Ethnicity == 'Black',])
  id <- KGP_[order(Ethn_Score_Black), ][Ethnicity == 'Black', ][n_b, IID]
  id <- which(KGP_[, IID] == id)
  rows <- c(c(1:(id-1)), c((id+1): nrow(KGP_)))
  cols <- c(1: (ncol(KGP_)-1))
  KGP_ <- KGP_[rows, ..cols]
  return( KGP_ )
}


calc_tresh <- function(demogr_1K, ethn_list, n_PCs){
  
  demogr_1K <- calc_ethn_score(demogr_1K, demogr_1K, ethn_list, n_PCs)
  
  treshs <- c(1:length(ethn_list))
  names(treshs) <- ethn_list
  
  for(et in ethn_list){
    col_ <- paste('Ethn_Score_', et, sep='')
    treshs[et] <- max((demogr_1K[Ethnicity == et, ..col_]))
  }
  
  return(treshs)
}



set_etnicity <- function(demogr_SE, demogr_1, etn_main, etn_sec, n_PCs){
  
  lbl_main <- paste('Ethn_Score_', etn_main, sep="")
  lbl_sec <- paste('Ethn_Score_', etn_sec, sep="")
  
  demogr_SE <- calc_ethn_score(demogr_SE, demogr_1, c(etn_main, etn_sec), n_PCs)
  tr <- calc_tresh(demogr_1, c(etn_main, etn_sec), n_PCs)
  
  colnames(demogr_SE)[c((ncol(demogr_SE)-1), ncol(demogr_SE))] <- c('e_sc_main', 'e_sc_sec')
  t <- demogr_SE[e_sc_sec<=tr[etn_sec], ][order(e_sc_main),][1,e_sc_main]
  
  #return(t)
  
  #return(demogr_SE)
  
  demogr_SE <- data.table(demogr_SE, 
                          V=ifelse(as.matrix(demogr_SE[, .(e_sc_main)]) < t, 1, 0) )
  colnames(demogr_SE)[c((ncol(demogr_SE)-2) : ncol(demogr_SE))] <- c(lbl_main, lbl_sec, etn_main)
  
  return(demogr_SE)

}  
  # last 3 cols are e_sc_main | e_sc_sec | e_sec
  #score_main <- demogr_SE[e_sc_main <= tr[etn_main], .(e_sc_main)]
  #score_main <- score_main[order(e_sc_main),]
  #score_main <- score_main[,e_sc_main]
  #score_main <- score_main[order(score_main)]
  #return(score_main)


determine_ethn <- function(demogr_UK, etn){
  cols_e <- colnames(demogr_UK)[startsWith(names(demogr_UK), 'Ethn_Score_')]
  col_etn <- paste('Ethn_Score_', etn, sep='')
  
  demogr_UK <- 
    data.table(demogr_UK ,
               this_etn = apply(as.matrix(demogr_UK[,..cols_e]), 1, 
                                function(r){if(min(r) == r[col_etn])
                                {return(1)} else {return(0)}}))
  colnames(demogr_UK)[ncol(demogr_UK)] <- etn
  return(demogr_UK)
}



# OLD
determine_ethn_tr <- function(demogr_UK_sc, trsh){
  cols_e <- colnames(demogr_UK_sc)[startsWith(names(demogr_UK), 'Ethn_Score_')]
  col_etn <- paste('Ethn_Score_', etn, sep='')
  
  demogr_UK <- 
    data.table(demogr_UK ,
               this_etn = apply(as.matrix(demogr_UK[,..cols_e]), 1, 
                                function(r){if(min(r) == r[col_etn])
                                {return(1)} else {return(0)}}))
  colnames(demogr_UK)[ncol(demogr_UK)] <- etn
  return(demogr_UK)
}




################# Import bgen sample of UKB: ##############################################

sample <- fread("E:\\Genomics\\Resources\\UKB_sample_file_newest\\ukb_bgen_sample_new.txt")
colnames(sample) <- 'IID'
sample <- sample[order(sample[,IID]),]

#################### PCA data UKB and 1000: ###############################################################

PCA_data <- fread("E:\\Genomics\\Resources\\PCA_Eur\\PCA_results\\ukb_1000_PCA_full.eigenvec")
cols <- c(2:ncol(PCA_data))
PCA_data <- PCA_data[, ..cols]
gc()

################# UPLOAD AND FILTER ETHNIC BACKGROUND  UKB: ###########################################

demogr_UKB <- fread("E:\\Genomics\\Resources\\Output_data_ph2\\demogr_geno.txt")

names(demogr_UKB)[1] <- 'IID'
demogr_UKB <- demogr_UKB[, .(IID, Ethnic_backgr_v0, Gen_ethnic_grp)]
gc()
demogr_UKB <- as.data.table(left_join(sample, demogr_UKB, by="IID"))
gc()

#################### Relable Ethnicity of UKB: ###########################################################

eth <- c('Unknown', 'Unknown', 'Undislosed', 
         'White_any', 'British', 'Irish', 'White_other',
         'Mixed', 'White_Black_Carr', 'White_Black_Afr', 'White_Asian', 'Mixed_other',
         'Asian_any', 'Indian', 'Pakisani', 'Bangladeshi', 'Asian_other',
         'Black_any', 'Carribbean', 'African', 'Black_other',
         'Chinese', 'Other_Group')

names(eth) <-  c('NA', '-1', '-3',
                 '1', '1001', '1002', '1003', 
                 '2', '2001', '2002', '2003', '2004', 
                 '3', '3001', '3002', '3003', '3004', 
                 '4', '4001', '4002', '4003',
                 '5', '6')

demogr_UKB <- ethn_relabel(demogr_UKB, 'Ethnic_backgr_v0', eth)

demogr_UKB <- demogr_UKB[,.(IID, Ethnicity, Gen_ethnic_grp)]

demogr_UKB[, 'IID'] <- as.character(demogr_UKB[, IID])

################ add PCA data to UKB demogr data ###############################

demogr_UKB <- as.data.table(left_join(demogr_UKB, PCA_data, by='IID'))

################# Demographic data of 1KGP: ##############################################

demogr_1KGP <- fread("E:\\Genomics\\Resources\\PCA_Eur\\1000_files\\1000_ethnicity_indep.txt")
demogr_1KGP <- demogr_1KGP[, .(IID, Population, Population.Description)]
gc()

#################### Relable Ethnicity of 1KGP: ####################################################

eth_1K <-  c('White', 'Black_US', 'American', 'Black', 'Asian', 'Asian', 'American', 'White', 'Indian_1K',
             'White', 'Unknown', 'White', 
             'Black', 'Asian', 'Asian', 'White', 'American', 'American', 'Asian', 'Black')

names(eth_1K) <- c("CEU", "ASW", "MXL", "ACB", "CDX", "CHS", "CLM", "GBR", "GIH", "FIN", "NA", "IBS", 
                   "LWK", "CHB", "JPT", "TSI", "PEL", "PUR", "KHV", "YRI")

demogr_1KGP <- ethn_relabel(demogr_1KGP[, .(IID, Population, Population.Description)], 'Population', eth_1K)

####################### Add PCA data to 1KGP demogr data #######################################

demogr_1KGP <- as.data.table(left_join(demogr_1KGP, PCA_data, by='IID'))

####################### Statistical Analysis: #######################################

########## This is the final version of the analysis! 
# AFR + CAR vs all, the most relevant analysis, use this one as final!!!
# USing 1KGP PCA data and ethnicity labeling, 
# compute ethnicity scores (chi-square stats, based on 
# maximum likelihood estimator of multivariate nomral distribution)
# (1)for each bgen sample from UKB:
# each sample gets the following ethnicity scores:
# Asian | Black | Indian | White
# Black score is based on Africans from Africa and Barbados, 
# with US excluded due to bias towards white ethnicity and to max the 
# UKB african sample
# (2)for each cluster Asian | Indian | White 
# calculate ethicity scores for the respective 
# population
# as tresholds to be the max of the respective ethnicity scores on 1KGP samples 
# (3)the treshold for the black population from UKB, we use
# the minimal ethnicity score for all from all 
# (the maximal Black ethnicity score, for which, when used as a treshold, 
# the intersection of 
# the Black population with the other three populations is empty)
# (4)Add to that all self-declared blacks from UKB whose Black 
# ethicity score is the smallest
# (5)Remove from black UKB population determined from (3) all 
# smaples, who are not self-declared blacks and whose black ethn score is
# not the smallest 

n_PCs=3
ethnic <- c("Asian", 'Black', "Indian_1K", 'White')

UKB_AFR_CAR <- calc_ethn_score(demogr_UKB, demogr_1KGP, ethnic, n_PCs)
trs <- calc_tresh(demogr_1KGP, ethnic, n_PCs)

# calculation of the real treshold for Black score
tb <- UKB_AFR_CAR[Ethn_Score_Asian <= trs['Asian'] |  
                 Ethn_Score_Indian_1K <= trs['Indian_1K'] | 
                 Ethn_Score_White <= trs['White'], ][order(Ethn_Score_Black),][1, Ethn_Score_Black]

UKB_blacks <- c('Black_any', 'Carribbean', 'African', 'Black_other')


UKB_AFR_CAR <- 
  data.table(UKB_AFR_CAR, 
             V=ifelse(as.matrix(UKB_AFR_CAR[, .(Ethn_Score_Black)]) < tb, 1, 0))

colnames(UKB_AFR_CAR)[ncol(UKB_AFR_CAR)] <- 'Black_vs_all'

UKB_AFR_CAR <- determine_ethn(UKB_AFR_CAR, 'Black')

colnames(UKB_AFR_CAR)[ncol(UKB_AFR_CAR)] <- 'Black_min'

UKB_AFR_CAR[,Black_min_score:=Black_min]

UKB_AFR_CAR[!(Ethnicity %in% UKB_blacks), 'Black_min'] <- 0 

UKB_AFR_CAR[, Black:=(Black_vs_all + Black_min - Black_vs_all*Black_min)*Black_min_score]

write.table(UKB_AFR_CAR, 
            "E:\\Genomics\\Resources\\AFR_sample_data\\AFR_PC_ethn_scores.txt", 
            append = FALSE, sep = "\t", quote = FALSE, col.names=TRUE, row.names=FALSE)
