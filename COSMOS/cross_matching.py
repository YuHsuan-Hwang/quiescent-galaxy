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

def MaskPosition(data1, ra1, dec1, data2, ra2, dec2, i):
    mask_ra = (data1[ra1]<=data2[ra2][i]+searching_radius) &\
    (data1[ra1] >=data2[ra2][i]-searching_radius)
    data_tmp1 = data1[mask_ra]
        
    delta_ra = data2[ra2][i]-data_tmp1[ra1]
    delta_dec = data2[dec2][i]-data_tmp1[dec1]
    distance = ( delta_ra**2 + delta_dec**2 )**0.5
    mask_dist = distance<=searching_radius
    data_tmp = data_tmp1[mask_dist]
    return data_tmp

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
catalog[0] = "COSMOS2015_Laigle+_v1.1_simple.fits"
#"01_COSMOS2015catalog/COSMOS2015/COSMOS2015_Laigle+_v1.1.fits" #COSMOS2015
catalog[1] = "04_COSMOS450_850/S2COSMOS/catalog/S2COSMOS_sourcecat850_Simpson18.fits" #wide850
catalog[2] = "COSMOS+mips24_allmatches.fits" #mips24 in COSMOS2015
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
for i in range(number):
    ReadCatalog(i,catalog[i])

###matching
searching_radius = 7.0/60/60
#total = 0
#print Matching(total)

mask_data0 = Mask_M(0) & Mask_photoz(0) & Mask_error(0) & Mask_class_star(0)
data0_tmp = data[0][mask_data0]

mask_data2 = (data[2][galaxyid[2]]>200000)&(data[2][galaxyid[2]]<996000) \
& Mask_M(2) & Mask_photoz(2) & Mask_error(2) & Mask_class_star(2)
data2_tmp = data[2][mask_data2]

mask_data3 = (data[3][galaxyid[3]]>200000)&(data[3][galaxyid[3]]<996000) 
data3_tmp = data[3][mask_data3]

SOURCE = [0]*len( data[0].filled() )
totalnumber = [0]*len(data[1].filled())

print 'start loop'
nomatch_count = 0
onematch_count = 0
onematch_samplecount = 0
twomatch_count = 0
twomatch_samplecount = 0
for i in range(1100,len(data[1].filled()),1):#range( len(data[1].filled()) ):
    #print ''
    #print 'i ', i
    data_tmp = MaskPosition(data0_tmp,ra[0],dec[0],data[1],ra[1],dec[1],i) 
    data_tmp2 = MaskPosition(data2_tmp,ra[2],dec[2],data[1],ra[1],dec[1],i)
    data_tmp3 = MaskPosition(data3_tmp,ra[3],dec[3],data[1],ra[1],dec[1],i)
    #print data_tmp[galaxyid[0]]
    #print data_tmp2[galaxyid[2]]
    #print data_tmp3[galaxyid[3]]
    #print '==='
    #print len(data_tmp.filled())
    #print len(data_tmp2.filled())
    #totalnumber[i] += 1000*len( data_tmp.filled() )
    flag = 0
    predictsource = []
    for j in range( len(data_tmp.filled()) ):
        mask = (data_tmp[galaxyid[0]][j] == data_tmp2[galaxyid[2]])
        data_tmptmp = data_tmp2[mask]
        for k in range( len( data_tmptmp.filled() ) ):
            predictsource.append(data_tmptmp[galaxyid[2]][k])
        mask2 = (data_tmp[galaxyid[0]][j] == data_tmp3[galaxyid[3]])
        data_tmptmp2 = data_tmp3[mask2]
        for k in range( len( data_tmptmp2.filled() ) ):
            predictsource.append(data_tmptmp2[galaxyid[3]][k])
        #totalnumber[i] += 100*len( data_tmptmp.filled() )
        #totalnumber[i] += 10*len( data_tmptmp2.filled() )
        for k in range( len( data_tmptmp.filled() ) ):
            masksource = data_tmptmp[galaxyid[2]][k] == data_tmptmp2[galaxyid[3]]
            datasource = data_tmptmp2[masksource]
            #totalnumber[i] += len( datasource.filled() )
            if len( datasource.filled() )!=0:
                flag = 1
                twomatch_samplecount += len( datasource.filled() )
                #for l in range( len( datasource.filled() ) ):
                    #SOURCE[ datasource[galaxyid[3]][l] ] = 1
                    
                    #print 'True: ',datasource[galaxyid[3]][l]
    if flag==0:
        if ( len( predictsource )!=0 ):
            onematch_count +=1
            onematch_samplecount += len( predictsource )
            #for k in range( len( predictsource ) ):
                #SOURCE[ predictsource[k] ] = 1
                #print 'True: ',predictsource[k]
        else:
            nomatch_count +=1
            #for k in range( len( data_tmp.filled() ) ):
                #SOURCE[ data_tmp[galaxyid[0]][k] ] = 1
            
    else:
        twomatch_count +=1
        print i

print nomatch_count
print onematch_count
print onematch_samplecount
print twomatch_count
print twomatch_samplecount
#print 'start writing table file'              
#data[0]['850SOURCE'] = SOURCE
#data[0].write('COSMOS2015_Laigle+_v1.1_850sources.fits')

#RegionFile(1,'850sources','green','4.0')

time2 = time.time()
print 'done! time =', (time2-time1)/60.0 , 'min'