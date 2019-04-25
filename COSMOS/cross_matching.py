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
    #mask_V_errorV = (data[index][V_error]<0.1) & (data[index][V_error]>0)
    #Suprime-Cam:i+band
    #mask_ip_error = (data[index][ip_error]<0.1) & (data[index][ip_error]>0)
    #mask_J_error = (data[index][J_error]<0.1) & (data[index][J_error]>0)
    #mask_Ks_error = (data[index][Ks_error]<0.1) & (data[index][Ks_error]>0)
    #return mask_V_errorV | mask_ip_error | mask_J_error | mask_Ks_error
    return (data[index][J_error]<0.2) & (data[index][J_error]>0)

def Mask_photoz(index):
    return (data[index][photoz]>0) & (data[index][photoz]<8)

def Mask_class_star(index):
    return (data[index][class_star]==0)

def CrossMatching_farIR():
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
    
    mask_data3 = (data[3]['CAT_CPT']=='COSMOS2015     ') 
    data3_tmp = data[3][mask_data3]
    
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
    
    print 'len(data[1].filled()) ',len(data[1].filled())
    
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
                if (SOURCE[ data_tmp[galaxyid[0]][j] -1] != 0):
                    print "!!!already tagged"
                else:
                    SOURCE[ data_tmp[galaxyid[0]][j] -1] = 1
                
        if flag==0:
            if ( len( predictsource24 )!=0 )or( len( predictsource3 )!=0 ):
                onematch_count +=1
                onematch_samplecount += len( predictsource24 )
                for k in range( len( predictsource24 ) ):
                    if (SOURCE[ predictsource24[k] -1] != 0):
                        print "!!!already tagged"
                    else:
                            SOURCE[ predictsource24[k] -1] = 2
                
                onematch_samplecount += len( predictsource3 )
                for k in range( len( predictsource3 ) ):
                    if (SOURCE[ predictsource3[k] -1] != 0):
                        print "!!!already tagged"
                    else:
                        SOURCE[ predictsource3[k] -1] = 3
            else:
                nomatch_count +=1
                nomatch_samplecount += len( data_tmp.filled() )
                for k in range( len( data_tmp.filled() ) ):
                    if (SOURCE[ data_tmp[galaxyid[0]][k] -1] != 0):
                        print "!!!already tagged"
                    else:
                        SOURCE[ data_tmp[galaxyid[0]][k] - 1] = 4
                
        else:
            twomatch_count +=1
        
    # TODO: not only one matched
    
    print 'nomatch_count ',nomatch_count
    print 'nomatch_samplecount ',nomatch_samplecount
    print 'onematch_count ',onematch_count
    print 'onematch_samplecount ',onematch_samplecount
    print 'twomatch_count ',twomatch_count
    print 'twomatch_samplecount ',twomatch_samplecount
    return SOURCE

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

###set crossmatching catalog


############
# 850 wide #
############
'''
CATALOG = "04_COSMOS450_850/S2COSMOS/catalog/S2COSMOS_sourcecat850_Simpson18.fits"
ID = "Short_ID"
RA = "RA_deg"
DEC = "DEC_deg"
searching_radius = 7.0
SOURCE_COLUMN = "850SOURCE"
INPUT_MAINCATALOG = "01_COSMOS2015catalog/COSMOS2015/COSMOS2015_Laigle+_v1.1.fits" #COSMOS2015
OUTPUT_MAINCATALOG = 'COSMOS2015_Laigle+_v1.1_850wide.fits'
print '850wide'
'''


##############
# 850 narrow #
##############
'''
CATALOG = "04_COSMOS450_850/STUDIES/sources_850.fits"
ID = ""
RA = "RA_850"
DEC = "DEC_850"
searching_radius = 7.0
SOURCE_COLUMN = "850NARROW"
INPUT_MAINCATALOG = "COSMOS2015_Laigle+_v1.1_850wide.fits" #COSMOS2015
OUTPUT_MAINCATALOG = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow.fits'
print '850narrow'
'''


##############
# 450 narrow #
##############
'''
CATALOG = "04_COSMOS450_850/STUDIES/sources_450.fits"
ID = ""
RA = "RA_450"
DEC = "DEC_450"
searching_radius = 4.0
SOURCE_COLUMN = "450NARROW"
INPUT_MAINCATALOG = "COSMOS2015_Laigle+_v1.1_850wide+850narrow.fits" #COSMOS2015
OUTPUT_MAINCATALOG = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow.fits'
print '450narrow'
'''

'''
###set catalog names
number = 5
catalog = [None]*number
catalog[0] = "COSMOS2015_Laigle+_v1.1_simple.fits"
catalog[1] = CATALOG
catalog[2] = "COSMOS+mips24_allmatches.fits" #mips24 in COSMOS2015
catalog[3] = "VLA_3GHz_counterpart_array_20170210_paper_smolcic_et_al.fits.txt" #3GHz in COSMOS2015
catalog[4] = INPUT_MAINCATALOG

###set ID column names
galaxyid = [None]*number
galaxyid[0] = "NUMBER"
galaxyid[1] = ID
galaxyid[2] = "NUMBER"
galaxyid[3] = "ID_CPT"

##set RA DEC column names
ra = [None]*number
dec = [None]*number
ra[0] = "ALPHA_J2000"
dec[0] = "DELTA_J2000"
ra[1] = RA
dec[1] = DEC
ra[2] = "ALPHA_J2000"
dec[2] = "DELTA_J2000"
ra[3] = "RA_CPT_J2000"
dec[3] = "DEC_CPT_J2000"

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
time2 = time.time()
print 'current time :', (time2-time1)/60.0 , 'min'
SOURCE = CrossMatching_farIR()

###output
time2 = time.time()
print 'current time :', (time2-time1)/60.0 , 'min'
print 'start writing table file'
ReadCatalog(4,catalog[4])              
data[4][SOURCE_COLUMN] = SOURCE
data[4].write(OUTPUT_MAINCATALOG)
'''

#############
# 24 micorn #
#############
'''
CATALOG = "COSMOS+mips24_allmatches.fits" #mips24 in COSMOS2015
ID = "NUMBER"
SOURCE_COLUMN = "24MICRON"
INPUT_MAINCATALOG = "COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow.fits" #COSMOS2015
OUTPUT_MAINCATALOG = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron.fits'
print '24micron'

###set catalog names
number = 2
catalog = [None]*number
catalog[0] = CATALOG
catalog[1] = INPUT_MAINCATALOG

###set ID column names
galaxyid = ID

###read catalog
data = [None]*number
for i in range(number):
    ReadCatalog(i,catalog[i])

print 'len(data[0]) ',len(data[0])

time2 = time.time()
print 'current time :', (time2-time1)/60.0 , 'min'
print 'loop start'
SOURCE = [0]*len( data[1].filled() )
for i in range(len(data[0])):
    if SOURCE[ data[0][galaxyid][i] -1 ]==1:
        print "!!!",data[0][galaxyid][i]-1,"!!!already tagged"
    else:
        SOURCE[ data[0][galaxyid][i] -1 ] = 1

print 'start writing table file'            
data[1][SOURCE_COLUMN] = SOURCE
data[1].write(OUTPUT_MAINCATALOG)
'''

#########
# 3 GHz #
#########
'''
CATALOG = "VLA_3GHz_counterpart_array_20170210_paper_smolcic_et_al.fits.txt" #3GHz in COSMOS2015
ID = "ID_CPT"
SOURCE_COLUMN = "3GHZ"
INPUT_MAINCATALOG = "COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron.fits" #COSMOS2015
OUTPUT_MAINCATALOG = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz.fits'
print '3GHz'

###set catalog names
number = 2
catalog = [None]*number
catalog[0] = CATALOG
catalog[1] = INPUT_MAINCATALOG

###set ID column names
galaxyid = ID

###read catalog
data = [None]*number
for i in range(number):
    ReadCatalog(i,catalog[i])

print 'len(data[0]) ',len(data[0])

mask = (data[0]['CAT_CPT']=='COSMOS2015     ')
data0_tmp = data[0][mask]
print 'len(mask_data0) ',len(data0_tmp)

time2 = time.time()
print 'current time :', (time2-time1)/60.0 , 'min'
print 'loop start'
SOURCE = [0]*len( data[1].filled() )
for i in range(len(data0_tmp)):
    if SOURCE[ data0_tmp[galaxyid][i] -1 ]==1:
        print "!!!",data0_tmp[galaxyid][i]-1,"!!!already tagged"
    else:
        SOURCE[ data0_tmp[galaxyid][i] -1 ] = 1

print 'start writing table file'            
data[1][SOURCE_COLUMN] = SOURCE
data[1].write(OUTPUT_MAINCATALOG)
'''

##############
# IR_AGN_ALL #
##############
'''
CATALOG = "ca_all_iragns.fits" #mips24 in COSMOS2015
ID = "IRAGN"
SOURCE_COLUMN = "IR_AGN_ALL"
INPUT_MAINCATALOG = "COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz.fits" #COSMOS2015
OUTPUT_MAINCATALOG = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz+iragnall.fits'
print 'ir_agn_all'

###set catalog names
number = 2
catalog = [None]*number
catalog[0] = CATALOG
catalog[1] = INPUT_MAINCATALOG

###set ID column names
galaxyid = ID

###read catalog
data = [None]*number
for i in range(number):
    ReadCatalog(i,catalog[i])

print 'start writing table file'            
data[1][SOURCE_COLUMN] = data[0][galaxyid]
data[1].write(OUTPUT_MAINCATALOG)
'''

CATALOG = "ca_mir_iragns.fits"
ID = "NUMBER"
SOURCE_COLUMN = "IR_AGN_MIR"
INPUT_MAINCATALOG = "COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz+iragnall.fits" #COSMOS2015
OUTPUT_MAINCATALOG = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz+iragnallmir.fits'
print 'ir_agn_mir'

###set catalog names
number = 2
catalog = [None]*number
catalog[0] = CATALOG
catalog[1] = INPUT_MAINCATALOG

###set ID column names
galaxyid = ID

###read catalog
data = [None]*number
for i in range(number):
    ReadCatalog(i,catalog[i])
print 'len(data[0]) ',len(data[0])

time2 = time.time()
print 'current time :', (time2-time1)/60.0 , 'min'
print 'loop start'
SOURCE = [0]*len( data[1].filled() )
for i in range(len(data[0])):
    if data[0]['IRAGN'][i]==1:
        if SOURCE[ data[0][galaxyid][i] -1 ]==1:
            print "!!!",data[0][galaxyid][i]-1,"!!!already tagged"
        else:
            SOURCE[ data[0][galaxyid][i] -1 ] = 1
        
print 'start writing table file'            
data[1][SOURCE_COLUMN] = SOURCE
data[1].write(OUTPUT_MAINCATALOG)

#RegionFile(1,'850sources','green','4.0')

time2 = time.time()
print 'done! time :', (time2-time1)/60.0 , 'min'