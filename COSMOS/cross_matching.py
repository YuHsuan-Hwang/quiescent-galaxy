#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  3 16:09:52 2019

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

def MaskData(data1, ra1, dec1, data2, ra2, dec2, i):
    delta_ra = data2[ra2][i]-data1[ra1]
    mask_ra = (delta_ra>=-searching_radius)&(delta_ra<=searching_radius)
    data_tmp1 = data1[mask_ra]
    
    delta_ra = data2[ra2][i]-data_tmp1[ra1]
    delta_dec = data2[dec2][i]-data_tmp1[dec1]
    distance = ( delta_ra**2 + delta_dec**2 )**0.5
    
    mask_dist = distance<=searching_radius
    data_tmp = data_tmp1[mask_dist]
    #print data_tmp
    return data_tmp
    
def Temp():
    totalnumber = [0]*len(data[1].filled())
    for i in range(3):#( len(data[1].filled()) ):
        data_tmp = MaskData(data[0],ra[0],dec[0],data[1],ra[1],dec[1],i)
        data_tmp2 = MaskData(data[2],ra[2],dec[2],data[1],ra[1],dec[1],i)
        print data[1][ra[1]][i]
        print data[1][dec[1]][i]
        print 'len(data_tmp.filled())',len(data_tmp.filled())
        print 'len(data_tmp2.filled())',len(data_tmp2.filled())
        print data_tmp[galaxyid[0]]
        print data_tmp2[galaxyid[0]]
        print data_tmp2[ra[0]]
        print data_tmp2[dec[0]]
        for j in range( len(data_tmp.filled()) ):
            for k in range( len(data_tmp2.filled()) ):
                if data_tmp[galaxyid[0]][j] == data_tmp2[galaxyid[0]][k]:
                    totalnumber[i] += 1
        print 'totalnumber[i]',totalnumber[i]
    return

def Matching_temp(total):
    for i in range( len(data[1].filled()) ):
        data_tmp = MaskData(data[2],ra[2],dec[2],data[1],ra[1],dec[1],i)
        total += len(data_tmp.filled())
    return total

def Matching(total):
    for i in range( len(data[1].filled()) ):
        mask_ra = (data[0][ra[0]]<=data[1][ra[1]][i]+searching_radius) &\
        (data[0][ra[0]]>=data[1][ra[1]][i]-searching_radius)
        data_tmp1 = data[0][mask_ra]
        
        delta_ra = data[1][ra[1]][i]-data_tmp1[ra[0]]
        delta_dec = data[1][dec[1]][i]-data_tmp1[dec[0]]
        distance = ( delta_ra**2 + delta_dec**2 )**0.5
        mask_dist = distance<=searching_radius
        data_tmp = data_tmp1[mask_dist]
        for j in range( len(data_tmp.filled()) ):
            total += 1
    return total

def Matchingprime01(total):
    for i in range( len(data[1].filled()) ):
        mask_ra = (data[0][ra[0]]<data[1][ra[1]][i]+searching_radius) &\
        (data[0][ra[0]]>data[1][ra[1]][i]-searching_radius)
        data_tmp1 = data[0][mask_ra]
        mask_dec = (data_tmp1[dec[0]]<data[1][dec[1]][i]+searching_radius) &\
        (data_tmp1[dec[0]]>data[1][dec[1]][i]-searching_radius)
        data_tmp = data_tmp1[mask_dec]
        for j in range( len(data_tmp.filled()) ):
            delta_ra = data[1][ra[1]][i]-data_tmp[ra[0]][j]
            delta_dec = data[1][dec[1]][i]-data_tmp[dec[0]][j]
            distance = ( delta_ra**2 + delta_dec**2 )**0.5
            if (distance<=searching_radius):
                total+=1
    return total

def Matchingprime02(total):
    sorted_ra = data[0].argsort(keys=ra[0], kind=None)
    for i in range( 15 ):#len(data[1].filled()) ):
        flag=0
        for j in sorted_ra:
            if flag==1:
                if (data[0][ra[0]][j]>data[1][ra[1]][i]+searching_radius):
                    break
                else:
                    delta_ra = data[1][ra[1]][i]-data[0][ra[0]][j]
                    delta_dec = data[1][dec[1]][i]-data[0][dec[0]][j]
                    distance = ( delta_ra**2 + delta_dec**2 )**0.5
                    if (distance<=searching_radius):
                        total+=1
            else:
                if (data[0][ra[0]][j]>data[1][ra[1]][i]-searching_radius):
                    flag=1
    return total

def RegionFile(index,filename,color,size):
    f = open('/Users/yuhsuan/Documents/research/05WH/data/'+filename+'.reg','w')
    f.write('global color='+color+' font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    N = len( data[index].filled() )
    for n in range(N):
        f.write('fk5;circle('+str(data[index][ra[index]][n])+','+str(data[index][dec[index]][n])+','+size+'") # text={'+str(n)+'}\n')
    f.close()
# =============================================================================
# main code
# =============================================================================

time1 = time.time()

###set catalog names
number = 4
catalog = [None]*number
catalog[0] = "COSMOS2015_Laigle+_v1.1_simple.fits" #COSMOS2015
#"COSMOS+wide850_allmatches.fits"
#"COSMOS2015_Laigle+_v1.1_simple.fits" #COSMOS2015
catalog[1] = "04_COSMOS450_850/S2COSMOS/catalog/S2COSMOS_sourcecat850_Simpson18.fits" #wide850
#"COSMOS+wide850_allmatches_simple.fits"
#"04_COSMOS450_850/S2COSMOS/catalog/S2COSMOS_sourcecat850_Simpson18.fits" #wide850
catalog[2] = "COSMOS+mips24_allmatches.fits" #mips24 in COSMOS2015
#"COSMOS+mips24＋wide850_allmatches.fits" #mips24 in COSMOS2015
catalog[3] = "VLA_3GHz_counterpart_array_20170210_paper_smolcic_et_al.fits.txt" #3GHz in COSMOS2015

###set ID column names
galaxyid = [None]*number
galaxyid[0] = "NUMBER"
galaxyid[1] = "Short_ID"
galaxyid[2] = "NUMBER"
galaxyid[3] = "ID_CPT"

##set RA DEC column names
ra = [None]*number
dec = [None]*number
ra[0] = "ALPHA_J2000"
dec[0] = "DELTA_J2000"
ra[1] = "RA_deg"
dec[1] = "DEC_deg"
ra[2] = "ALPHA_J2000"
dec[2] = "DELTA_J2000"

###read catalog
data = [None]*number
for i in range(number):
    ReadCatalog(i,catalog[i])

###matching
searching_radius = 7.0/60/60
total = 0
#print Matching(total)
#print Matching_temp(total)
Temp()

#RegionFile(1,'850sources','green','4.0')

time2 = time.time()
print 'done! time =', time2-time1 , 'sec'