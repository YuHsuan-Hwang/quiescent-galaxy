#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 13:34:58 2019

@author: yuhsuan
"""
# =============================================================================
# packages
# =============================================================================
from __future__ import division
import time
from astropy.table import Table
import matplotlib.pyplot as plt
# =============================================================================
# functions
# =============================================================================
def ReadCatalog(index,catalog):
    data[index] = Table.read(catalog, hdu=1)
    return

def Mask_M(index):
    down = -30
    up = 0
    mask_color1 = (data[index][color1]>down) & (data[index][color1]<up)
    mask_color2 = (data[index][color2]>down) & (data[index][color2]<up)
    mask_color3 = (data[index][color3]>down) & (data[index][color3]<up)
    mask_M = mask_color1 & mask_color2 & mask_color3
    return mask_M

def Mask_error(index):
    mask_V_errorV = (data[index][V_error]<0.1) & (data[index][V_error]>0)
    #Suprime-Cam:i+band
    mask_ip_error = (data[index][ip_error]<0.1) & (data[index][ip_error]>0)
    mask_J_error = (data[index][J_error]<0.1) & (data[index][J_error]>0)
    mask_Ks_error = (data[index][Ks_error]<0.1) & (data[index][Ks_error]>0)
    return mask_V_errorV | mask_ip_error | mask_J_error | mask_Ks_error

def Mask_photoz(index):
    return (data[index][photoz]>0) & (data[index][photoz]<8)

def Mask_class_star(index):
    return (data[index][class_star]==0)
# =============================================================================
# main code
# =============================================================================
time1 = time.time()

###set catalog names
number = 2
catalog = [None]*number
catalog[0] = "COSMOS2015_Laigle+_v1.1_simple.fits"
#"01_COSMOS2015catalog/COSMOS2015/COSMOS2015_Laigle+_v1.1.fits" #COSMOS2015
catalog[1] = "COSMOS+mips24_allmatches.fits" #mips24 in COSMOS2015

###set ID column names
galaxyid = [None]*number
galaxyid[0] = "NUMBER"
galaxyid[1] = "NUMBER"

###set colors
colorname1 = "NUV"
colorname2 = "r"
colorname3 = "J"

###set columns in the main catalog
color1 = "MNUV"
color2 = "MR"
color3 = "MJ"
color = [ color1, color2, color3 ]
photoz = "PHOTOZ"
V_error = "V_MAGERR_APER3"
ip_error = "ip_MAGERR_APER3"
J_error = "J_MAGERR_APER3"
Ks_error = "Ks_MAGERR_APER3"
mass = "MASS_MED"
class_star = "TYPE"
class_SFG = "CLASS"
ra = "ALPHA_J2000"
dec = "DELTA_J2000"

###read catalog
data = [None]*number
for i in range(number):
    ReadCatalog(i,catalog[i])

###identify 24 sources
'''
SOURCE = [False]*len( data[0].filled() )
for i in range(91589,91597,1):#range(len( data[0].filled() )):
    print data[0][galaxyid[0]][i]
    for j in range(len( data[1].filled() )):
        if ( data[0][galaxyid[0]][i] == data[1][galaxyid[1]][j] ):
            SOURCE[i] = True
    print SOURCE[i]
 '''   

mask_data0 = Mask_M(0) & Mask_photoz(0) & Mask_error(0) & Mask_class_star(0)
data0_tmp = data[0][mask_data0]
 
mask_data1 = (data[1][galaxyid[1]]>200000)&(data[1][galaxyid[1]]<996000) \
& Mask_M(1) & Mask_photoz(1) & Mask_error(1) & Mask_class_star(1)
data1_tmp = data[1][mask_data1]

SOURCE = [False]*len( data[0].filled() )
for i in range(200):#range(len( data0_tmp.filled() )):
    #print data0_tmp[galaxyid[0]][i]
    mask = (data0_tmp[galaxyid[0]][i] == data1_tmp[galaxyid[1]])
    data_tmp = data1_tmp[mask]
    if len( data_tmp.filled() )!=0:
        SOURCE[data0_tmp[galaxyid[0]][i]] = True
        #print 'True'
    #print SOURCE[data[0][galaxyid[0]][i]]

'''
SOURCE = [False]*len( data[0].filled() )
for i in range(500):#range(len( data0_tmp.filled() )):
    #print data0_tmp[galaxyid[0]][i]
    mask = (data0_tmp[galaxyid[0]][i] == data[1][galaxyid[1]])
    data_tmp = data[1][mask]
    if len( data_tmp.filled() )!=0:
        SOURCE[data[0][galaxyid[0]][i]] = True
    #print SOURCE[data[0][galaxyid[0]][i]]
'''
data[0]['24SOURCE'] = SOURCE
print data[0]['24SOURCE']
#data[0].write('COSMOS2015_Laigle+_v1.1_simple_24sources.fits')

time2 = time.time()
print 'done! time =', time2-time1 , 'sec'