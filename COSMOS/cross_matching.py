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
#import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
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
    mask_V_errorV = (data[index][V_error]<0.2) & (data[index][V_error]>0)
    #Suprime-Cam:i+band
    mask_ip_error = (data[index][ip_error]<0.2) & (data[index][ip_error]>0)
    mask_J_error = (data[index][J_error]<0.2) & (data[index][J_error]>0)
    mask_Ks_error = (data[index][Ks_error]<0.2) & (data[index][Ks_error]>0)
    return mask_V_errorV | mask_ip_error | mask_J_error | mask_Ks_error

def Mask_photoz(index):
    return (data[index][photoz]>0) & (data[index][photoz]<8)

def Mask_class_star(index):
    return (data[index][class_star]==0)

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
number = 5
catalog = [None]*number
catalog[0] = "COSMOS2015_Laigle+_v1.1_simple.fits"
catalog[1] = "04_COSMOS450_850/S2COSMOS/catalog/S2COSMOS_sourcecat850_Simpson18.fits" #wide850
catalog[2] = "COSMOS+mips24_allmatches.fits" #mips24 in COSMOS2015
catalog[3] = "VLA_3GHz_counterpart_array_20170210_paper_smolcic_et_al.fits.txt" #3GHz in COSMOS2015
catalog[4] = "01_COSMOS2015catalog/COSMOS2015/COSMOS2015_Laigle+_v1.1.fits" #COSMOS2015

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
ra[3] = "RA_CPT_J2000"
dec[3] = "DEC_CPT_J2000"

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

###read catalog
data = [None]*number
for i in range(4):
    ReadCatalog(i,catalog[i])

###matching
searching_radius = 7.0
#total = 0
#print Matching(total)

mask_data0 = Mask_M(0) & Mask_photoz(0) & Mask_error(0) & Mask_class_star(0)
data0_tmp = data[0][mask_data0]
print 'len(data0_tmp.filled()) ',len(data0_tmp.filled())

mask_data2 = Mask_M(2) & Mask_photoz(2) & Mask_error(2) & Mask_class_star(2)

#(data[2][galaxyid[2]]>200000)&(data[2][galaxyid[2]]<996000) \
#& Mask_M(2) & Mask_photoz(2) & Mask_error(2) & Mask_class_star(2)
data2_tmp = data[2][mask_data2]

print 'len(data2_tmp.filled()) ',len(data2_tmp.filled())

#mask_data3 = 
#(data[3][galaxyid[3]]>200000)&(data[3][galaxyid[3]]<996000) 
data3_tmp = data[3]#[mask_data3]

print 'len(data3_tmp.filled()) ',len(data3_tmp.filled())

SOURCE = [0]*len( data[0].filled() ) #twomatch:1, onematch:2/3, nomatch:4

print 'start loop'
nomatch_count = 0
nomatch_samplecount = 0
onematch_count = 0
onematch_samplecount = 0
twomatch_count = 0
twomatch_samplecount = 0

c0 = SkyCoord(ra=data0_tmp[ra[0]], dec=data0_tmp[dec[0]])
c1 = SkyCoord(ra=data[1][ra[1]]*u.degree, dec=data[1][dec[1]]*u.degree)
c2 = SkyCoord(ra=data2_tmp[ra[2]], dec=data2_tmp[dec[2]])
c3 = SkyCoord(ra=data3_tmp[ra[3]]*u.degree, dec=data3_tmp[dec[3]]*u.degree)

for i in range( len(data[1].filled()) ):
    sep = c0.separation(c1[i])
    data_tmp = data0_tmp[ sep<=searching_radius*u.arcsec ]
    
    sep = c2.separation(c1[i])
    data_tmp2 = data2_tmp[ sep<=searching_radius*u.arcsec ]
    
    sep = c3.separation(c1[i])
    data_tmp3 = data3_tmp[ sep<=searching_radius*u.arcsec ]   
    
    flag = 0
    predictsource24 = []
    predictsource3 = []
    for j in range( len(data_tmp.filled()) ):
        if data_tmp[galaxyid[0]][j] in data_tmp2[galaxyid[2]]:
            flag_24 = True
            predictsource24.append(data_tmp[galaxyid[0]][j])
        else:
            flag_24 = False
            
        if data_tmp[galaxyid[0]][j] in data_tmp3[galaxyid[3]]:
            flag_3 = True
            predictsource3.append(data_tmp[galaxyid[0]][j])
        else:
            flag_3 = False
        
        if (flag_24)&(flag_3):
            flag = 1
            twomatch_samplecount += 1
            if (SOURCE[ data_tmp[galaxyid[0]][j] ] != 0):
                print "!!!already tagged"
            else:
                SOURCE[ data_tmp[galaxyid[0]][j] ] = 1
            
    if flag==0:
        if ( len( predictsource24 )!=0 )or( len( predictsource3 )!=0 ):
            onematch_count +=1
            onematch_samplecount += len( predictsource24 )
            for k in range( len( predictsource24 ) ):
                if (SOURCE[ predictsource24[k] ] != 0):
                    print "!!!already tagged"
                else:
                        SOURCE[ predictsource24[k] ] = 2
            
            onematch_samplecount += len( predictsource3 )
            for k in range( len( predictsource3 ) ):
                if (SOURCE[ predictsource3[k] ] != 0):
                    print "!!!already tagged"
                else:
                    SOURCE[ predictsource3[k] ] = 3
        else:
            nomatch_count +=1
            nomatch_samplecount += len( data_tmp.filled() )
            for k in range( len( data_tmp.filled() ) ):
                if (SOURCE[ data_tmp[galaxyid[0]][k] ] != 0):
                    print "!!!already tagged"
                else:
                    SOURCE[ data_tmp[galaxyid[0]][k] ] = 4
            
    else:
        twomatch_count +=1
    
# TODO: not only one matched

print 'nomatch_count ',nomatch_count
print 'nomatch_samplecount ',nomatch_samplecount
print 'onematch_count ',onematch_count
print 'onematch_samplecount ',onematch_samplecount
print 'twomatch_count ',twomatch_count
print 'twomatch_samplecount ',twomatch_samplecount

'''
print 'start writing table file'
ReadCatalog(4,catalog[4])              
data[4]['850SOURCE'] = SOURCE
data[4].write('COSMOS2015_Laigle+_v1.1_850sources_newerr.fits')
'''

#RegionFile(1,'850sources','green','4.0')

time2 = time.time()
print 'done! time =', (time2-time1)/60.0 , 'min'