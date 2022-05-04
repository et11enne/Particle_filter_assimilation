#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 12:10:31 2022

@author: cape
"""
import netCDF4
import os
import os.path
import numpy as np
import matplotlib.pyplot as plt
from datetime import date, timedelta, datetime
from snowtools.utils.prosimu import prosimu
import numpy as np
from osgeo import gdal

def read_tif(filename, bd):
    
    driver = gdal.GetDriverByName('GTiff')
    dataset = gdal.Open(filename)
    band = dataset.GetRasterBand(bd)
    cols = dataset.RasterXSize
    rows = dataset.RasterYSize
    transform = dataset.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = -transform[5]
    data = band.ReadAsArray(0, 0, cols, rows)
    return data

class PrepaPro(object):
    
    labels = ["N","NE","E","SE","S","SW","W","NW"]
    alt_niv = np.arange(600,3900,300)
    nniveaux = ['600','', '1200',  '', '1800', '', '2400', '','3000', '', '3600']
    orient = np.arange(0,360,45)
        
    def __init__(self, path_npy, Massif, Type, Saison, path_pro, type_mod):
    
        # Creation dates simu
        self.list_date_str_S1 = np.load(path_npy+"dates_S1_str.npy")
        self.nc = prosimu(path_pro)
        self.Type = Type
        self.date_list_xticks = []

        self.type_mod = type_mod
        self.list_date_asc = []
        self.list_date_dsc = []
        self.liste_date1 = []
        self.class_tel_massif = np.zeros((365,101,143), float)
        self.index_obs = []
        self.vec_binary_crocus = []
         
    def load(self,path_to_load):
        
        self.class_tel_massif = np.load(path_to_load)
        
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

        nbdat = np.shape(self.date_list_xticks)[0]
    
        
        return(self.date_list_xticks)
            
    def TEL_calc(self):
        date_simu=datetime(2018,8,1,18)
        if(self.type_mod == "SEMI"):
            
            output_season=np.zeros((len(self.list_date_asc),len(self.alt_niv),len(self.orient),50), float)
        
            # We take ascendant
            for idx_date,date_value in enumerate(self.liste_date1):
                print("day:",idx_date)
                print(date_value)
                for idx_alt,alt in enumerate(self.alt_niv):
                    for idx_ori,ori in enumerate(self.orient):
                        snow_dz = self.nc.read_var('SNOWDZ',time=self.nc.get_time(date_value),
                                              Number_of_points = self.nc.get_points(massif_num=12,slope=20,ZS=alt,aspect=ori))
                        snow_liq = self.nc.read_var('SNOWLIQ',time=self.nc.get_time(date_value),
                                               Number_of_points = self.nc.get_points(massif_num=12,slope=20,ZS=alt,aspect=ori))
                        for j in range(0,50):
                            output_season[idx_date,idx_alt,idx_ori,j]=snow_dz[j]*snow_liq[j]
                        output_season[np.isnan(output_season)] = 0  
                        
        elif(self.type_mod == "DIST"):
            
            list_date = self.nc.readtime()
            self.list_date_asc = list_date[0::2]
            
            xx = self.nc.read_var('xx',time=self.nc.get_time(date_simu))
            yy = self.nc.read_var('yy',time=self.nc.get_time(date_simu))
            
            
            
            output_season=np.zeros((365,len(yy),len(xx),50), float)
            for idx_date,date_value in enumerate(self.list_date_asc):
                
                snow_dz = self.nc.read_var('SNOWDZ',time=self.nc.get_time(date_value))
                snow_liq = self.nc.read_var('SNOWLIQ',time=self.nc.get_time(date_value))
                print("date",idx_date , " : ", date_value)
                
                for j in range(0,50):
                            output_season[idx_date,:,:,j]=snow_dz[j,:,:]*snow_liq[j,:,:]
            output_season[output_season == 10**20] = 0
            sum_layer_season = np.einsum('ijkl->ijk', output_season)
            return(sum_layer_season)
            
    def build_massif_class(self):
    
        path_MNT = "/home/cape/OUTPUT_ET_PRO/Ange/MNTLouisGRoussecorrected.nc"
        aniveau1=600
        aniveau2=3900
        pas=300
        
        nniv1 = np.arange(aniveau1,aniveau2,pas)
        nbniv = np.shape(nniv1)[0]
        nc_MNT = prosimu(path_MNT)
        zs_mnt = nc_MNT.read_var('ZS')
        
        slope = read_tif("/home/cape/SXCEN_LOCAL/MNT_MASSIF/test_mnt/T2.tif",1)
        aspect = read_tif("/home/cape/SXCEN_LOCAL/MNT_MASSIF/test_mnt/T3.tif",1)

        slope=slope[::-1,:]
        aspect=aspect[::-1,:]
                
    
        # dimension of our vector
        naa = np.shape(zs_mnt)[0]
        nab = np.shape(zs_mnt)[1]
        na = naa*nab
        # vector creation
        vec_zs = np.reshape(zs_mnt.copy(),na)
        vec_slope = np.reshape(slope.copy(),na)
        vec_aspect = np.reshape(aspect.copy(),na)
        
        test_list_orient=[[337.5,22.5],[22.5,67.5],[67.5,112.5],[112.5,157.5],
                          [157.5,202.5],[202.5,247.5],[247.5,292.5],[292.5,337.5]]
        #test_list_orient=[[270,90],[90,270]]
        
        norien=len(test_list_orient)
        # We selct the range of slope
        slo1 = 15
        slo2 = 25
        nniveaux = ['600','', '1200',  '', '1800', '', '2400', '','3000', '', '3600']
        
        crocus_TEL = np.zeros((365,norien,nbniv), float)
        mean_crocus_TEL = np.zeros((365,norien,nbniv), float)
        test_TEL = np.zeros((365,norien,nbniv), float)
        
        #sum_output_season[sum_output_season < 9] = 0
        
        
        sum_layer_season = self.TEL_calc()
        
        for i in range (0,364):
            vec_CROCUS_TEL=np.reshape(sum_layer_season[i,:,:], na)
            for jj in range(0,norien):
        
                test_asp1=test_list_orient[jj][0]
                test_asp2=test_list_orient[jj][1]
                #print("Orientation :[",test_asp1,"][",test_asp2,"]")
                for ii in range(0, nbniv):
                            niv1=nniv1[ii]
                            
                                ######.
                            if(test_asp1==337.5 ):
                               
                                ######################### TEST MEAN
                                
                                test = np.nonzero(((vec_aspect >= test_asp1) | (vec_aspect < test_asp2)) &
                                                (vec_slope < slo2) & (vec_slope >= slo1) & 
                                                (vec_zs >= (niv1-150)) & (vec_zs < (niv1+150)))
                                vec_tempo = (vec_CROCUS_TEL[test])
                                vec_tempo = vec_tempo[vec_tempo!=0]
        
                                
                                vec_CROCUS_TEL_1 = np.mean(vec_tempo)
                                test_TEL[i,jj,ii]=vec_CROCUS_TEL_1
        
                                
                               
                            else:   
                            
                                ######################### TEST MEAN

                              
                                test=np.nonzero(((vec_aspect >= test_asp1) & (vec_aspect < test_asp2)) &
                                (vec_slope < slo2) & (vec_slope >= slo1) & 
                                (vec_zs >= (niv1-150)) & (vec_zs < (niv1+150)))
                                vec_tempo=(vec_CROCUS_TEL[test])
                                
                                vec_tempo=vec_tempo[vec_tempo!=0]
                                
                                
                                vec_CROCUS_TEL_1=np.mean(vec_tempo)
                                test_TEL[i,jj,ii]=vec_CROCUS_TEL_1
        test_TEL[np.isnan(test_TEL)] = 0
        self.class_tel_massif=test_TEL
        
    def binary_with_threshold(self,min_TEL):
        
        image_crocus = self.class_tel_massif[:,:,:] # precise day?
        minimum_TEL = min_TEL
        crocus_binary = np.copy(image_crocus)
        
        crocus_binary[np.copy(crocus_binary) <= minimum_TEL] = 0
        crocus_binary[np.copy(crocus_binary) > minimum_TEL ] = 1
        
        self.vec_binary_crocus = crocus_binary
        
        
    
    
    def get_reduce2obs_date(self):
    
        smaller_class_tel_massif = np.zeros((len(self.index_obs),len(self.alt_niv),len(self.orient)), float)
        smaller_class_tel_massif = self.class_tel_massif[[self.index_obs],:,:]
        smaller_class_tel_massif= smaller_class_tel_massif[0,:,:,:]
        return(smaller_class_tel_massif)
            
    
    def get_prep_array_class(self):
        if(self.type_mod == 'SEMI'):
            print("OKKKK")
            self.class_tel_massif = self.class_tel_massif.transpose(0, 2, 1)

        return(self.class_tel_massif)
   
 