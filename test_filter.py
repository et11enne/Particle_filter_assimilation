#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 10:54:38 2022

@author: cape
"""
import matplotlib.pyplot as plt
import numpy as np
from Filter import ParticleFilter 



Massif="GdesRousses"
Saison="2018-2019"
S1_type="S1_ASC" #choice with [S1_ASC S1_DSC]


'''# S1 Obs prepa
prepaObsS1 = PrepaObs("/home/cape/SXCEN_LOCAL/save_npy/slope/slope_inf20/"+Massif+"/"+S1_type+"/"+Saison+"/",
                      Massif,S1_type,Saison)
obs_frac = prepaObsS1.get_obs_array()


# Pro prepa DISTRIB
proprepa1 = PrepaPro("/home/cape/SXCEN_LOCAL/save_npy/slope/slope_inf20/"+Massif+"/"+S1_type+"/"+Saison+"/",
                     Massif,"ASC", Saison,
                     "/home/cape/OUTPUT_ET_PRO/DISTRIBUTED_"+Massif+"/PRO_2018080106_2019080106.nc","DIST")
list_date_obs_simu = proprepa1.date_fusion()
proprepa1.build_massif_class()
pro_tel = proprepa1.get_prep_array_class()
pro_tel_smaller = proprepa1.get_reduce2obs_date() 
# 

# Pro prepa SEMI
proprepa2 = PrepaPro("/home/cape/SXCEN_LOCAL/save_npy/slope/slope_inf20/"+Massif+"/"+S1_type+"/"+Saison+"/",
                     Massif,"ASC", Saison,
                     "/home/cape/OUTPUT_ET_PRO/PRO_GdesRousses_SEMI/PRO_2018080106_2019080106.nc","SEMI")

proprepa2.load("/home/cape/SXCEN_LOCAL/save_npy/CROCUS/SEMI_DISTRIBUTED/GdesRousses/2018-2019/ASC/sum_output_season.npy")
list_date_obs_simu1 = proprepa2.date_fusion()

pro_tel2 = proprepa2.get_prep_array_class()
pro_tel_smaller2 = proprepa2.get_reduce2obs_date() 
#

plt.imshow((proprepa1.class_tel_massif[:,0,:]).T)
plt.show()

# Correlation test DIST/OBS
corr_test = Correlation(obs_frac,pro_tel_smaller, list_date_obs_simu, Massif)
corr_test.spearman_calc()

#Correlation test SEMI/Obs
corr_test1 = Correlation(obs_frac,pro_tel_smaller2, list_date_obs_simu, Massif)
corr_test1.spearman_calc()

'''

'''
ensemble_test = ProEns("/home/cape/OUTPUT_ET_PRO/DISTRIBUTED_"+Massif+"/PRO_2018080106_2019080106.nc",
                       "/home/cape/save_npy/slope/slope_inf20/"+Massif+"/"+S1_type+"/"+Saison+"/",
                       "/cnrm/cen/users/NO_SAVE/lafaysse/cache/cache/vortex/s2m/grandesrousses/E2@lafaysse/",
                       Saison,Massif,S1_type,"ASC","2018080106_2019080106")

list_date_of_obs = ensemble_test.date_fusion()

dict_ensemble_date_obs = ensemble_test.load_array()
essai = ensemble_test.calc_corr_spearman_day(day=42)
'''

PF = ParticleFilter(Massif,S1_type,Saison)

#ok = PF.weight_normalize(42,'N')




#effweight_test = PF.efficientWeights()
test = PF.find_best_member_ori_day()


#test = PF.find_best_member_ori_day('N')


#test = PF.find_best_member_ori_day('S')
