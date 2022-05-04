#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 10:28:43 2022

@author: cape
"""
from Ensemble import ProEns
import collections
import random
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
###############################################################################################
#                                       TODO
#
#   verifier que S1 et crocus sont dans le même sens  
#
#
###############################################################################################

class ParticleFilter(ProEns,object):

    def __init__(self,Massif,S1_type,Saison):
        '''
        Constructor
        '''
        #self.options = options
        #self.dd = dd  # specific date
        self.load()
        
        self.score_members = dict()
        self.weight_dict = dict()
        self.ensemble_test = ProEns("/home/cape/OUTPUT_ET_PRO/DISTRIBUTED_"+Massif+"/PRO_2018080106_2019080106.nc",
                       "/home/cape/save_npy/slope/slope_inf20/"+Massif+"/"+S1_type+"/"+Saison+"/",
                       "/cnrm/cen/users/NO_SAVE/lafaysse/cache/cache/vortex/s2m/grandesrousses/E2@lafaysse/",
                       Saison,Massif,S1_type,"ASC","2018080106_2019080106") 
        self.list_date_of_obs = self.ensemble_test.date_fusion()

        dict_ensemble_date_obs = self.ensemble_test.load_array()
        
        self.nmembers = self.ensemble_test.nmembers

        
    def weight_normalize(self,ensemble):
        #normalize the result so that the probabilities correctly sum to 1.0. 
        #Normalization is done by dividing each element by the sum of all elements in the list.
        self.weight_dict = ensemble
        
        total = sum(self.weight_dict.values(), 0.0)
        self.weight_dict = {k: v / total for k, v in self.weight_dict.items()}
        
        
        return(self.weight_dict)
        
    def load(self):

    	#self.ens = ProEns(self.options,self.dd)
        #self.obs = PrepaObs()
        print("dvlp")
        
    def efficientWeights(self):
    # reinvented formula 20/03
        
    
        
        effweights = np.around(1. / (self.nmembers * np.mean([w**2 for w in self.weight_dict.values()])))

        return(effweights)
    
    def find_best_member_ori_day(self,ori=None):
    
        list_bmember = []
        list_bmember_N = []
        list_bmember_S = []
        
        for dd in range(0,len(self.list_date_of_obs),1):
            bmember_N = self.weight_normalize(self.ensemble_test.calc_corr_spearman_day(dd,'N',0,10))
            bmember_S = self.weight_normalize(self.ensemble_test.calc_corr_spearman_day(dd,'S',0,10))
            bmember   = self.weight_normalize(self.ensemble_test.calc_corr_spearman_day(dd,ori))
            
            new_data_N = {k: v for k, v in bmember_N.items() if v >= (max(bmember_N.values())  -0.01)}
            new_data_S = {k: v for k, v in bmember_S.items() if v >= (max(bmember_S.values())  -0.01)}
            new_data   = {k: v for k, v in bmember.items()   if v >= (max(bmember.values())    -0.01)}
            
            
            #print(new_data_N)
            #print(new_data_S)
            for key in new_data_S.keys() & new_data_N.keys():
                
                list_bmember.append(key)
            for key in new_data_S.keys() :
                list_bmember_S.append(key)
            for key in new_data_N.keys():   
                list_bmember_N.append(key)
           
            

        occurence_member = Counter(list_bmember)
        occurence_member = collections.OrderedDict(sorted(occurence_member.items()))
        self.barplot_dict(occurence_member,"occurence best member North/South")
        occurence_member = Counter(list_bmember_S)
        occurence_member = collections.OrderedDict(sorted(occurence_member.items()))
        self.barplot_dict(occurence_member,"occurence best member South")
        occurence_member = Counter(list_bmember_N)
        occurence_member = collections.OrderedDict(sorted(occurence_member.items()))
        self.barplot_dict(occurence_member,"occurence member North")
        
        plt.show()
        
        
        
        ### weigth ponderate by altitude on specific day 
        dd=42
        bmember_N_low_alt = self.weight_normalize(self.ensemble_test.calc_corr_spearman_day(dd,'N',0,5))
        bmember_N_high_alt = self.weight_normalize(self.ensemble_test.calc_corr_spearman_day(dd,'N',5,10))
        new_data_N_low_alt = {k: v for k, v in bmember_N_low_alt.items() if v >= (max(bmember_N_low_alt.values())  - 0.01)}
        new_data_N_high_alt = {k: v for k, v in bmember_N_high_alt.items() if v >= (max(bmember_N_high_alt.values())  - 0.01)}
        
        print("##############################################################################")
        print("LOW ALTITUDE : ", new_data_N_low_alt)
        print("##############################################################################")
        print("HIGH ALTITUDE : ", new_data_N_high_alt)
        print("##############################################################################")
              
        a_counter = Counter(new_data_N_low_alt)
        b_counter = Counter(new_data_N_high_alt)
        for k in b_counter.keys():
            b_counter[k] = b_counter[k] * 2
        add_dict = a_counter + b_counter
        bmember_N_fusion = dict(add_dict)
        
        
        bmember_N_fusion_normalize = self.weight_normalize(bmember_N_fusion)
        

        
        self.barplot_dict(bmember_N_fusion_normalize,"score member North ponderate")
        plt.savefig("/home/cape/code_py/Filter/New_Filter_based_Cluzet/fig/alt_ponderate/test_ponderate"+str(dd)+"_N.png")
        plt.show()
        if(ori==None):    
            return(new_data)
        if(ori=='N'):
            return(new_data_N)
        if(ori=='S'):
            return(new_data_S)
            

        #new_data_N = {k: v for k, v in bmember_N.items() if v > max(bmember_N.values)-0.02 }    
       

    def barplot_dict(self,dictio,title=''):
        fig=plt.figure(figsize=(5,6)) 
        names = list(dictio.keys())
        values = list(dictio.values())
        plt.bar(range(len(dictio)), values, tick_label=names)
        plt.title(title)
        plt.xticks(rotation=90)
        
    def resampling(self):
        N=self.nmembers
        weights=[]
        new_particles=[]
        index = int(random.random() * N)
        maxW=0

        for i in range(N):
            weights.append(self.particles[i].weight)

        beta = 0.0
        mw = max(weights)
        for i in range(N):
            beta += random.random() * 2.0 * mw
        		#print "beta =", beta
            while beta > weights[index]:
            		beta -= weights[index]
            		index = (index + 1) % N
            		#print "\tbeta= %f, index = %d, weight = %f" % (beta, index, weights[index])
            new_particles.append(self.particles[index])
        
        # Taken from git code 
        '''
        # resampling (systematic resampling) step
        if self.n_eff < self.n_eff_threshold:
            indices = self.resample_fn(self.weights)
            self.particles = self.particles[indices, :]
            self.weights = np.ones(self.n_particles) / self.n_particles

        # randomly resample some particles from the prior
        if self.resample_proportion > 0:
            random_mask = (
                np.random.random(size=(self.n_particles,)) < self.resample_proportion
            )
            self.resampled_particles = random_mask
            self.init_filter(mask=random_mask)
        '''
    
    '''
    def localisation(self,):
        #Tackle degeneracy by reducing the number of observations simultaneously assimilated
        #Solution
        print("dvlp")
        
    def inflation(self,):
        # Tackle degeneracy by increasing the PF tolerance
        print("dvlp")
        
    def Pf(self,):
        likelihoods = np.array([np.exp(-0.5 * np.matmul(np.matmul((matEns[:, i] - vectObs).T, R_1), (matEns[:, i] - vectObs)))
                                    for i in range(self.pb.nmembers)])
        weights = likelihoods / np.sum(likelihoods)
    
    def Q_Q_score(self):
        
      
    def resampling(self):
        
        #dvlp
        print("dvlp")
    
    def get_score_member(self, dd):
        
        
    def error_add(self):'''