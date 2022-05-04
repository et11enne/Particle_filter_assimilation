#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 11:13:18 2022

@author: cape
"""
from Correlation import Correlation
from snowtools.utils.prosimu import prosimu
from ObsS1 import PrepaObs
from ProPrepa import PrepaPro
import numpy as np
import os
import os.path
import matplotlib.pyplot as plt

from collections import defaultdict

def find_all(name, path):
    result = []
    for root, dirs, files in os.walk(path):
        if name in files:
            result.append(os.path.join(root, name))
    return result



class ProEns(PrepaPro, PrepaObs):
    
    labels = ["N","NE","E","SE","S","SW","W","NW"]
    alt_niv = np.arange(600,3900,300)
    nniveaux = ['600','', '1200',  '', '1800', '', '2400', '','3000', '', '3600']
    
    def __init__(self, path_pro, path_npy , path_ensembliste , Saison , Massif, S1_type,type_mod, Saison_pro):
        
        self.path_file_ensembliste = path_ensembliste
        # Creation dates simu
        self.Saison = Saison
        self.list_date_str_S1 = np.load(path_npy+"dates_S1_str.npy")
        self.nc = prosimu(path_pro)
        self.Type = type_mod
        self.date_list_xticks = []
        self.nmembers = 0
        self.ens = defaultdict(list)
        self.member = []
        self.stack = dict()
        self.isloaded = False
        self.isstacked = False
        self.nbdates = 0
        self.index_obs = []
        self.pro_filename = "PRO_"+Saison_pro+".nc"
        self.type_mod = type_mod
        self.result_pro_loc = []
        self.coefcorr_spearman = dict()
        # S1 Obs prepa
        prepaObsS1 = PrepaObs("/home/cape/save_npy/slope/slope_inf20/"+Massif+"/"+S1_type+"/"+Saison+"/",
                              Massif,S1_type,self.Saison)
        self.obs_array_S1 = prepaObsS1.get_obs_array()
        
        

    def load_array(self):
        

        # Finding folder member name
        self.result_pro_loc = find_all(self.pro_filename,self.path_file_ensembliste)
      
        # looping through every member 
        for i in range(0,len(self.result_pro_loc)):           
           
            self.member.append(self.result_pro_loc[i][83:len(self.result_pro_loc[i])-33:1])
            
            member_uniq=self.member[i]
            
                #variable de chaque membre
            sum_layer_season=np.load("/home/cape/ensembliste/"+self.type_mod+"/"+self.Saison+"/"+member_uniq+"/sum_output_season.npy")
            smaller_sum_layer_season=np.zeros((len(self.index_obs),len(self.alt_niv),len(self.orient)), float)
            smaller_sum_layer_season=sum_layer_season[[self.index_obs],:,:]
            smaller_sum_layer_season=smaller_sum_layer_season[0,:,:,:]
            
            self.ens[member_uniq].append(smaller_sum_layer_season)
            
            
            #Selection of only interesting dates
            #smaller_sum_layer_season = smaller_sum_layer_season[:,alt_min:alt_max,:]
        self.nmembers=len(self.member)
        print("loading arrays successfull ,", len(self.member) , " members found")
        self.coefcorr_spearman = dict.fromkeys(self.member)
        
        return(self.ens)


    def date_fusion(self):
        # Creation dates simu
        
        list_date = self.nc.readtime()
        self.list_date_asc = list_date[0::2]
        self.list_date_dsc = list_date[1::2]
        self.liste_date1 = []
        if(self.Type == "ASC"):
            self.liste_date1 = self.list_date_asc
        elif(self.Type == "DSC"):
            self.liste_date1 = self.list_date_dsc
 
        list_date_str_simu = [None]*(len(self.liste_date1))

        for idx, i in enumerate(self.liste_date1):
            self.list_date_asc[idx] = str(self.list_date_asc[idx])
            list_date_str_simu[idx] = self.list_date_asc[idx][:10]

        list_date_str_simu = [s.replace("-", "") for s in list_date_str_simu]
        
        idx = 0
        self.index_obs = []
        for i in self.list_date_str_S1:

            try:
                while i != list_date_str_simu[idx]:
                    idx = idx + 1
                else:
                    self.index_obs.append(idx)
                    self.date_list_xticks.append(list_date_str_simu[idx])
            except IndexError:
                print("")

        self.nbdates = np.shape(self.date_list_xticks)[0]
    
        print("date to obs finish (ensemble)")
        return(self.date_list_xticks)
        
    def calc_corr_spearman_day(self,day=None,ori=None,alti_min=None,alti_max=None):
        
        corr_per_day = Correlation()
        
        
        for i in range(0,len(self.result_pro_loc)):
            
            member_uniq = self.member[i]
            tempo_member = self.ens[member_uniq] # extracting array [1,dd,ori,alt]
            tempo_member = tempo_member [0] # list is coming from the dict 
            if(day==None):
                for dd in range(0,self.nbdates):
    
                    self.coefcorr_spearman[member_uniq] = corr_per_day.spearman_calc_day(tempo_member[dd].T,self.obs_array_S1[dd,:,:],ori)
                    plt.plot(self.coefcorr_spearman[member_uniq])
                    
            else:
                if((alti_min==None) & (alti_max==None)):
                    self.coefcorr_spearman[member_uniq] = corr_per_day.spearman_calc_day(tempo_member[day].T,self.obs_array_S1[day,:,:],ori)
                else:
                    self.coefcorr_spearman[member_uniq] = corr_per_day.spearman_calc_day(tempo_member[day].T,self.obs_array_S1[day,:,:],ori,alti_min,alti_max)
                #print(self.coefcorr_spearman)
        plt.show()
        
        return(self.coefcorr_spearman)
            
        '''
        for mb in range(1, self.nmembers + 1):
            self.ens[mb].load()
        self.listvar = []
        for var in list(self.ens[1].data.keys()):
            self.listvar.append(var)
        self.pgd = self.ens[1].pgd
        self.isloaded = True
        '''
        
        
        
        
'''
    def stackit(self):
        # for median calc and manipulation, to stack the members.
        if not self.isstacked:
            # self.data = dict()
            if not self.isloaded:
                self.load()

            for var in self.listvar:
                # print ('shape before transposing stack', np.shape(self.ens[1].data[var]))
                self.stack[var] = self.ens[1].data[var]
                for mb in range(2, self.nmembers + 1):
                    self.stack[var] = np.vstack((self.stack[var], self.ens[mb].data[var]))
                self.stack[var] = self.stack[var].T  # change in the dimension of the stack (pts as lines, mbs as cols)
                # self.stack[var] = np.ma.masked_where(self.stack[var] == self.stack[var][:].fill_value, self.stack[var])
        self.isstacked = True

'''