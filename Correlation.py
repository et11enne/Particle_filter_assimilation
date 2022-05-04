#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 11:33:33 2022

@author: cape
"""

import netCDF4
import os
import os.path
import numpy as np
from scipy.stats.stats import spearmanr 
import matplotlib.pyplot as plt
from datetime import date, timedelta, datetime
from snowtools.utils.prosimu import prosimu
import numpy as np
from osgeo import gdal



class Correlation(object):
    
    labels = ["N","NE","E","SE","S","SW","W","NW"]
    alt_niv = np.arange(600,3900,300)
    nniveaux = ['600','', '1200',  '', '1800', '', '2400', '','3000', '', '3600']
    orient = np.arange(0,360,45)
    alt_min=0
    alt_max=10 # 10 for max altitude = 3300 +-150
    
    def __init__(self, array_S1 = None, array_PRO=None, list_dates_obs=None, massif=None):
    
        # 3D of same size 
        self.array_S1 = array_S1
        self.array_PRO = array_PRO
        self.dates = list_dates_obs
        #self.nbdates = np.shape(self.dates)[0]
        self.massif = massif
   
    def spearman_calc_season(self):
        
        # signatificativy test
        p_value=[]
        p_value_N=[]
        p_value_S=[]
        
        # Init of list containing Spearman coef
        list_coef_corr_spearman=[]
        list_coef_corr_spearman_N=[]
        list_coef_corr_spearman_S=[]
        
        # I removed alt_min alt_max can be replace with : self.array_S1 = self.array_S1[:,:,alt_min:alt_max]
        # AltS1_downgrad_D = np.delete(AltS1_downgrad_D,delete_extra_date_S1, axis=0) #ATTENTION PROBLABLEMENT DELETE EXTRA EXTRA DATE 
        south_ori_selec=[0,1,2,3,5,6,7]
        
        # S1 array focus ori
        self.array_S1 = self.array_S1[:,:,self.alt_min:self.alt_max] # Useful only for GdesRousses
        array_S1_N = self.array_S1[:,:1,:]
        array_S1_S = np.delete(self.array_S1, south_ori_selec, axis=1)
        
        # Crocus array focus ori 
        self.array_PRO = self.array_PRO[:,:,self.alt_min:self.alt_max] # Useful only for GdesRousses
        array_PRO_N = self.array_PRO[:,:1,:]
        array_PRO_S = np.delete(self.array_PRO, south_ori_selec, axis=1)
        
        naa = np.shape(array_S1_N)[1]
        nab = np.shape(array_S1_N)[2]
        na_N = naa*nab
        naa = np.shape(array_S1_S)[1]
        nab = np.shape(array_S1_S)[2]
        na_S = naa*nab
        naa = np.shape(self.array_S1)[1]
        nab = np.shape(self.array_S1)[2]
        na = naa*nab
        
        for dd in range(0,self.nbdates):
        
            # 1D all ori
            vec_S1 = np.reshape(self.array_S1[dd,:,:],na)
            vec_CROCUS = np.reshape(self.array_PRO[dd,:,:],na)
            
            # 1D vector for North
            vec_S1_N = np.reshape(array_S1_N[dd,:,:],na_N)
            vec_CROCUS_N = np.reshape(array_PRO_N[dd,:,:],na_N)
    
            # 1D vector for south
            vec_S1_S = np.reshape(array_S1_S[dd,:,:],na_S)
            vec_CROCUS_S = np.reshape(array_PRO_S[dd,:,:],na_S)
            
            
            rho, pval = spearmanr(vec_CROCUS,vec_S1)
            list_coef_corr_spearman.append(rho)
            p_value.append(pval)
            rho, pval = spearmanr(vec_CROCUS_N,vec_S1_N)
            list_coef_corr_spearman_N.append(rho)
            p_value_N.append(pval)
            rho, pval = spearmanr(vec_CROCUS_S,vec_S1_S)
            list_coef_corr_spearman_S.append(rho)
            p_value_S.append(pval)
        fig, host = plt.subplots(figsize=(12,5))   

        plt.plot(list_coef_corr_spearman,label="all alt Spearman")
        plt.plot(list_coef_corr_spearman_N,label="all alt Spearman N")
        plt.plot(list_coef_corr_spearman_S,label="all alt Spearman S")
        
        ax1=plt.gca()
        ax1.set_title('Spearman correlations comparisons  , '+self.massif)
        plt.xticks(range(0,len(self.dates)),self.dates,rotation=90)
        fig.tight_layout()
        plt.legend()

    def spearman_calc_day(self,arrayS1_1D,arrayCROCUS_1D,ori=None,alti_min=None,alti_max=None):
        #
        # 0 - 600 ; 1 - 900 ; 2 -1200 ; 3 - 1500 ; 4 - 1800 ; 5 - 2100 ...
        #
        self.array_S1 = arrayS1_1D 
        self.array_PRO = arrayCROCUS_1D
        
    
        
        # I removed alt_min alt_max can be replace with : self.array_S1 = self.array_S1[:,:,alt_min:alt_max]
        # AltS1_downgrad_D = np.delete(AltS1_downgrad_D,delete_extra_date_S1, axis=0) #ATTENTION PROBLABLEMENT DELETE EXTRA EXTRA DATE 
        south_ori_selec=[0,1,2,3,5,6,7]
        
        # S1 array focus ori
        if((alti_min==None) & (alti_max==None)):
            self.array_S1 = self.array_S1[:,self.alt_min:self.alt_max] # Useful only for GdesRousses
            self.array_PRO = self.array_PRO[:,self.alt_min:self.alt_max] # Useful only for GdesRousses
        else:
            self.array_S1 = self.array_S1[:,alti_min:alti_max] # Useful only for GdesRousses
            self.array_PRO = self.array_PRO[:,alti_min:alti_max] # Useful only for GdesRousses
        array_S1_N = self.array_S1[:1,:]
        array_S1_S = np.delete(self.array_S1, south_ori_selec, axis=1)
        
        # Crocus array focus ori 
        
        array_PRO_N = self.array_PRO[:1,:]
        array_PRO_S = np.delete(self.array_PRO, south_ori_selec, axis=1)
        
        naa = np.shape(array_S1_N)[0]
        nab = np.shape(array_S1_N)[1]
        na_N = naa*nab
        naa = np.shape(array_S1_S)[0]
        nab = np.shape(array_S1_S)[1]
        na_S = naa*nab
        naa = np.shape(self.array_S1)[0]
        nab = np.shape(self.array_S1)[1]
        na = naa*nab
        
        
        
        # 1D all ori
        vec_S1 = np.reshape(self.array_S1[:,:],na)
        vec_CROCUS = np.reshape(self.array_PRO[:,:],na)
        
        # 1D vector for North
        vec_S1_N = np.reshape(array_S1_N[:,:],na_N)
        vec_CROCUS_N = np.reshape(array_PRO_N[:,:],na_N)

        # 1D vector for south
        vec_S1_S = np.reshape(array_S1_S[:,:],na_S)
        vec_CROCUS_S = np.reshape(array_PRO_S[:,:],na_S)
        
        
        rho, pval = spearmanr(vec_CROCUS,vec_S1)
        rho_N, pval = spearmanr(vec_CROCUS_N,vec_S1_N)
        rho_S, pval = spearmanr(vec_CROCUS_S,vec_S1_S)
        #rho, pval = spearmanr(vec_CROCUS_N,vec_S1_N)
        #list_coef_corr_spearman_N.append(rho)
        
        #rho, pval = spearmanr(vec_CROCUS_S,vec_S1_S)
        #list_coef_corr_spearman_S.append(rho)
        if(ori==None):
            return(rho)
        elif (ori=='N'):
            return(rho_N)
        elif (ori=="S"):
            return(rho_S)
        
''' 
    def hamming_calc(self):
        #Hamming calc over binary vectors BY CLASS
    def covariance_calc(self):
        #cov
    def global_calc(self):
        # Only 1 score --> 1 obs
    def R_loc_TEL(self):
        # taking pixel b pixel 
        # for each pixel looking for a group of pixel around
    def Q_Q(self):'''