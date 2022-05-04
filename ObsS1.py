#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 11:13:48 2022

@author: cape
"""
import numpy as np 

class PrepaObs(object):

    def __init__(self,path_npy,Massif,S1_type,Saison):
		# Creation dates simu
		
        self.list_date_str_S1 = np.load(path_npy+"dates_S1_str.npy")

   
		#Importing SENTINEL DATA  
        self.AltS1_downgrad = np.load(path_npy+Massif+"downgrad"+"_AltS1a_diag.npy")
        self.AltS1_downgrad_North = self.AltS1_downgrad[:,:1,:]
        self.AltS1_downgrad_South = np.delete(self.AltS1_downgrad, [0,1,2,3,5,6,7], axis=1)
        self.vector1D_S1 = []
    
   
    def get_obs_array(self):
        
        return(self.AltS1_downgrad)
        
    def dim1_vect_day(self):
        #PREPARING SENTINEL1 VECTOR , SAME VECTOR IS GOING TO BE USED FOR EACH MEMBER
        vector1D_S1 = []
        for dd in range(0,nbdat):
                vec_S1_N=np.reshape(self.AltS1_downgrad[dd,:,:],na)
                vector1D_S1.append(vec_S1_N)
	

