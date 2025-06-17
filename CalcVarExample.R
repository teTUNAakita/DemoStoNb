rm(list=ls())

ExpSigmaDemo2 = function(Nb, NbN){
  tmp_exact = 1 + (1 - NbN) / (Nb/2 - 1)
  tmp_approx = 1 + 1 / (Nb/2 - 1)
  return(list(tmp_exact,tmp_approx))
} # Eq.14&15

CV = function(sigmaE, sigmaD){
  (exp(sigmaE^2 + sigmaD^2) - 1)^0.5
} # Eq.12

Ratio = function(sigmaE, sigmaD){
  (exp(sigmaE^2) - 1) / (exp(sigmaE^2 + sigmaD^2) - 1)
} # Eq.17

res = NULL
res_approx = NULL
sigmaE = 0.4

# Example 1
# Michigan Invasive Red Swamp Crayfish
# 930 or 2262 SNPs
# from https://onlinelibrary.wiley.com/doi/full/10.1111/eva.70007

# Hotel 1 (cohort 1, before chemical eradication)
# sample size: 90
Nb_crayfish_1 = (27.0 + 33) / 2 # arithmetic mean of Nb_LD and Nb_wang
N_crayfish_1 = 74.6 # NS^hat
NbN_crayfish_1 = Nb_crayfish_1 / N_crayfish_1

res$crayfish1[1] = ExpSigmaDemo2(Nb_crayfish_1, NbN_crayfish_1)[[1]]
sigmaD = (log(res$crayfish1[1]))^0.5
res$crayfish1[2] = CV(sigmaE, sigmaD)
res$crayfish1[3] = Ratio(sigmaE, sigmaD)

res_approx$crayfish1[1] = ExpSigmaDemo2(Nb_crayfish_1, NbN_crayfish_1)[[2]]
sigmaD = (log(res_approx$crayfish1[1]))^0.5
res_approx$crayfish1[2] = CV(sigmaE, sigmaD)
res_approx$crayfish1[3] = Ratio(sigmaE, sigmaD)

# Hotel 1 (cohort 2, after chemical eradication)
# sample size: 147
Nb_crayfish_2 = (17.6 + 24) / 2
N_crayfish_2 = 45.1
NbN_crayfish_2 = Nb_crayfish_2 / N_crayfish_2

res$crayfish2[1] = ExpSigmaDemo2(Nb_crayfish_2, NbN_crayfish_2)[[1]]
sigmaD = (log(res$crayfish2[1]))^0.5
res$crayfish2[2] = CV(sigmaE, sigmaD)
res$crayfish2[3] = Ratio(sigmaE, sigmaD)

res_approx$crayfish2[1] = ExpSigmaDemo2(Nb_crayfish_2, NbN_crayfish_2)[[2]]
sigmaD = (log(res_approx$crayfish2[1]))^0.5
res_approx$crayfish2[2] = CV(sigmaE, sigmaD)
res_approx$crayfish2[3] = Ratio(sigmaE, sigmaD)
##########################################################################################

# Example 2
# Reticulated Flatwood Salamander
# 9 microsatellite markers
# from https://zslpublications.onlinelibrary.wiley.com/doi/full/10.1111/acv.12871

# cite 1 2015-16
# sample size: 31 + 92 = 123 (larvae + metamorphs)
Nb_salamander_1 = 37 # weighted mean of Nb_LD and Nb_SF
N_salamander_1 = 180 # 
NbN_salamander_1 = Nb_salamander_1 / N_salamander_1

res$salamander1[1] = ExpSigmaDemo2(Nb_salamander_1, NbN_salamander_1)[[1]]
sigmaD = (log(res$salamander1[1]))^0.5
res$salamander1[2] = CV(sigmaE, sigmaD)
res$salamander1[3] = Ratio(sigmaE, sigmaD)

res_approx$salamander1[1] = ExpSigmaDemo2(Nb_salamander_1, NbN_salamander_1)[[2]]
sigmaD = (log(res_approx$salamander1[1]))^0.5
res_approx$salamander1[2] = CV(sigmaE, sigmaD)
res_approx$salamander1[3] = Ratio(sigmaE, sigmaD)

# cite 2 2015-16
# sample size: 51 + 60 = 111
Nb_salamander_2 = 25
N_salamander_2 = 100 
NbN_salamander_2 = Nb_salamander_2 / N_salamander_2

res$salamander2[1] = ExpSigmaDemo2(Nb_salamander_2, NbN_salamander_2)[[1]]
sigmaD = (log(res$salamander2[1]))^0.5
res$salamander2[2] = CV(sigmaE, sigmaD)
res$salamander2[3] = Ratio(sigmaE, sigmaD)

res_approx$salamander2[1] = ExpSigmaDemo2(Nb_salamander_2, NbN_salamander_2)[[2]]
sigmaD = (log(res_approx$salamander2[1]))^0.5
res_approx$salamander2[2] = CV(sigmaE, sigmaD)
res_approx$salamander2[3] = Ratio(sigmaE, sigmaD)

##########################################################################################

# Example 3
# Brook Trout
# 11 microsatellite markers
# from https://royalsocietypublishing.org/doi/suppl/10.1098/rspb.2015.2601

# stream HE (Healey)
# sample size: 71
Nb_trout_1 = 30.9 # Nb_adj2
N_trout_1 = 1878 # Nc
NbN_trout_1 = Nb_trout_1 / N_trout_1

res$trout1[1] = ExpSigmaDemo2(Nb_trout_1, NbN_trout_1)[[1]]
sigmaD = (log(res$trout1[1]))^0.5
res$trout1[2] = CV(sigmaE, sigmaD)
res$trout1[3] = Ratio(sigmaE, sigmaD)

res_approx$trout1[1] = ExpSigmaDemo2(Nb_trout_1, NbN_trout_1)[[2]]
sigmaD = (log(res_approx$trout1[1]))^0.5
res_approx$trout1[2] = CV(sigmaE, sigmaD)
res_approx$trout1[3] = Ratio(sigmaE, sigmaD)

# stream GK (Gaskill)
# sample size: 78
Nb_trout_2 = 43.5
N_trout_2 = 1299 
NbN_trout_2 = Nb_trout_2 / N_trout_2

res$trout2[1] = ExpSigmaDemo2(Nb_trout_2, NbN_trout_2)[[1]]
sigmaD = (log(res$trout2[1]))^0.5
res$trout2[2] = CV(sigmaE, sigmaD)
res$trout2[3] = Ratio(sigmaE, sigmaD)

res_approx$trout2[1] = ExpSigmaDemo2(Nb_trout_2, NbN_trout_2)[[2]]
sigmaD = (log(res_approx$trout2[1]))^0.5
res_approx$trout2[2] = CV(sigmaE, sigmaD)
res_approx$trout2[3] = Ratio(sigmaE, sigmaD)

res
res_approx