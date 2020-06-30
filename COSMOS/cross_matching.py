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
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from astropy.cosmology import FlatLambdaCDM

# =============================================================================
# functions
# =============================================================================
def ReadCat( index, data,cat ):
    data[index] = Table.read(cat, hdu=1)
    return

def SetupData( inputmask ):
    
    x.append( data[0][color2] - data[0][color3] )
    y.append( data[0][color1] - data[0][color2] )
    mask.append( inputmask )
    x_masked.append( x[0][inputmask] )
    y_masked.append( y[0][inputmask] )
    
    return

def OutputCat( column_name, column_list ):
    t              = Table()
    t[column_name] = column_list
    return t

def Mask_M(index):
    
    down, up    = -30, 0
    mask_color1 = ( data[index][color1]>down ) & ( data[index][color1]<up )
    mask_color2 = ( data[index][color2]>down ) & ( data[index][color2]<up )
    mask_color3 = ( data[index][color3]>down ) & ( data[index][color3]<up )
   
    return mask_color1 & mask_color2 & mask_color3

def Mask_error( mode, up, index ):
    
    if mode==1:
        mask_V_errorV = ( data[index][V_error] <up ) & ( data[index][V_error] >0 )
        mask_ip_error = ( data[index][ip_error]<up ) & ( data[index][ip_error]>0 ) #Suprime-Cam:i+band
        mask_J_error  = ( data[index][J_error] <up ) & ( data[index][J_error] >0 )
        mask_Ks_error = ( data[index][Ks_error]<up ) & ( data[index][Ks_error]>0 )
        return mask_V_errorV | mask_ip_error | mask_J_error | mask_Ks_error
    
    if mode==2:
        return (data[index][V_error]<up) & (data[index][V_error]>0)
    
    if mode==3:
        return (data[index][Ks_mag]<24)
    

def Mask_photoz( index ):
    return (data[index][photoz]>0)

def Mask_class_star( index ):
    return ( data[index][class_star]==0 ) | ( data[index][class_star]==2 )

def MaskAll( index ):
    return Mask_M(index) & Mask_error( 1, 0.1, index ) & Mask_photoz( index ) & Mask_class_star( index )

def Mask_myclassQG ( index ): return ( y[index]> 3.1 ) & ( y[index]> 3.0*x[index]+1.0 )
def Mask_myclassSFG( index ): return ( y[index]<=3.1 ) | ( y[index]<=3.0*x[index]+1.0 )



def Matching_24Micron():
    
    print "Matching 24um sources"
    
    ###set martirials
    path = "/Users/yuhsuan/Documents/research/05WH/data/COSMOS/"
    CAT  = path + "02_mips24/COSMOS+mips24_allmatches.fits" #mips24 in COSMOS2015
    
    ID            = "NUMBER"
    SOURCE_COLUMN = "24MICRON"
    MAIN_CAT      = "COSMOS2015_Laigle+_v1.1_simple.fits" #COSMOS2015
    OUTPUT_CAT    = 'COSMOS2015_24micron.fits'
    
    ###set catalog names
    number      = 2
    cat_name    = [None]*number
    cat_name[0] = MAIN_CAT
    cat_name[1] = CAT
    
    ###set ID column names
    galaxyid = ID
    
    ###read catalog
    print "reading catalog ..."
    data = [None]*number
    for i in range( number ): ReadCat( i, data, cat_name[i] )
    
    time2 = time.time()
    print 'current time :', ( time2-time1 )/60.0 , 'min'
    
    ###matching
    num = 0
    SOURCE = [0]*len( data[0].filled() ) #create an array with length of COSMOS2015
    
    print "matching loop start ..."
    print 'CAT length ',len( data[1] )
    for i in range( len(data[1]) ): #go through every id in CAT
        
        source_id = data[1][galaxyid][i] - 1 #source id start with 0
        
        if SOURCE[ source_id ]==1:
            print "!!!",source_id,"!!!already tagged"
        else:
            SOURCE[ source_id ] = 1
            num +=1
    
    print "source number ",num
    
    ###output catalog
    print "writing table file ..."
    output_cat = OutputCat( SOURCE_COLUMN, SOURCE )  
    output_cat.write( OUTPUT_CAT )
    
    return
   
def Matching_3GHz():
    
    print "Matching 3GHz sources"
    
    ###set martirials
    path = "/Users/yuhsuan/Documents/research/05WH/data/COSMOS/05_3GHz/"
    CAT  = path + "VLA_3GHz_counterpart_array_20170210_paper_smolcic_et_al.fits.txt" #3GHz in COSMOS2015
    
    ID            = "ID_CPT"
    SOURCE_COLUMN = "3GHZ"
    MAIN_CAT      = "COSMOS2015_Laigle+_v1.1_simple.fits" #COSMOS2015
    OUTPUT_CAT    = 'COSMOS2015_3GHz.fits'
    LABEL_COLUMN  = ["Xray_AGN","MIR_AGN","SED_AGN","Quiescent_MLAGN","SFG",
                     "Clean_SFG","HLAGN","MLAGN","Radio_excess"]
    
    ###set catalog names
    number      = 2
    cat_name    = [None]*number
    cat_name[0] = MAIN_CAT
    cat_name[1] = CAT
    
    ###set ID column names
    galaxyid = ID
    
    ###read catalog
    print "reading catalog ..."
    data = [None]*number
    for i in range( number ): ReadCat( i,data,cat_name[i] )
    
    print 'CAT length ',len(data[1])
    
    ###mask cat label
    mask      = ( data[1]['CAT_CPT']=='COSMOS2015     ' )
    data1_tmp = data[1][mask]
    print 'label masked CAT length ',len( data1_tmp )
    
    time2 = time.time()
    print 'current time :', ( time2-time1 )/60.0 , 'min'
    
    ###matching
    num = 0
    SOURCE = [0]*len( data[0].filled() ) #create an array with length of COSMOS2015
    
    print "matching loop start ..."
    for i in range( len(data1_tmp) ): #go through every id in CAT
        
        source_id = data1_tmp[galaxyid][i] - 1
        
        if SOURCE[ source_id ]==1:
            print "!!!",source_id,"!!!already tagged"
        else:
            SOURCE[ source_id ] = 1
            num +=1
    
    print "source number ",num
    
    output_cat                = Table()
    output_cat[SOURCE_COLUMN] = SOURCE
    
    ###matching label
    print "matching AGN labels ..."
    for i in range( len(LABEL_COLUMN) ):
        
        LABEL_SOURCE = [-99]*len( data[0].filled() )
        
        for j in range( len(data1_tmp) ): #go through every id in CAT
            
            source_id = data1_tmp[galaxyid][j] - 1
            
            if data1_tmp[LABEL_COLUMN[i]][j].strip()=="true": #VLA cat shows true
                
                #tag the sample as an AGN
                if LABEL_SOURCE[ source_id ]==1:
                    print "!!!",source_id,"!!!already tagged"
                else:
                    LABEL_SOURCE[ source_id ] = 1
                    
            else: #VLA cat shows false
                
                #tag the sample with zero
                if LABEL_SOURCE[ source_id ]==0:
                    print "!!!",source_id,"!!!already tagged"
                else:
                    LABEL_SOURCE[ source_id ] = 0
                    
        output_cat[LABEL_COLUMN[i]] = LABEL_SOURCE
    
    ###output catalog
    print "writing table file ..."
    #print output_cat
    output_cat.write(OUTPUT_CAT)
    return

def Matching_AS2COSMOS():
    
    print "Matching AS2COSMOS catalog"
    
    ###set martirials
    path = "/Users/yuhsuan/Documents/research/05WH/data/COSMOS/"
    CAT  = path + "07_ALMA/AS2COSMOScatalog_2019-06-01_merged_wJMSphotom_wZspec.fits"
    
    ID            = "ID"
    RA            = "RA_final"
    DEC           = "Dec_final"
    SOURCE_COLUMN = "AS2COSMOS"
    MAIN_CAT      = "COSMOS2015_Laigle+_v1.1_simple.fits" #COSMOS2015
    OUTPUT_CAT    = 'COSMOS2015_AS2COSMOS.fits'
    
    searching_radius = 1.0
    
    ###set catalog names
    number      = 2
    cat_name    = [None]*number
    cat_name[0] = MAIN_CAT
    cat_name[1] = CAT
    
    ###set ID column names
    galaxyid    = [None]*number
    galaxyid[0] = "NUMBER"
    galaxyid[1] = ID
    
    ##set RA DEC column names
    ra     = [None]*number
    dec    = [None]*number
    ra[0]  = "ALPHA_J2000"
    dec[0] = "DELTA_J2000"
    ra[1]  = RA
    dec[1] = DEC
    
    ###read catalog
    print "reading catalog ..."
    data = [None]*number
    for i in range(2): ReadCat( i, data,cat_name[i] )
    
    time2 = time.time()
    print 'current time :', ( time2-time1 )/60.0 , 'min'
    
    #matching
    SOURCE = [0]*len( data[0].filled() ) #create an array with length of COSMOS2015
    c0 = SkyCoord( ra=data[0][ra[0]],          dec=data[0][dec[0]]          )
    c1 = SkyCoord( ra=data[1][ra[1]]*u.degree, dec=data[1][dec[1]]*u.degree )
    
    num = 0
    print "matching loop start ..."
    print 'CAT length ',len(data[1].filled())
    
    for i in range( len(data[1].filled()) ):
        
        sep = c0.separation(c1[i])
        data0_incircle = data[0][ sep<=searching_radius*u.arcsec ]
        
        for j in range( len(data0_incircle.filled()) ):
            source_id = data0_incircle[galaxyid[0]][j] -1
            if ( SOURCE[ source_id ] != 0 ):
                print "!!!",source_id,"!!!already tagged"
            else:
                SOURCE[ source_id ] = 1
                num += 1

    print "source number ",num
    
    time2 = time.time()
    print 'current time :', ( time2-time1 )/60.0 , 'min'
    
    ###output
    print "writing table file ..."
    output_cat = OutputCat( SOURCE_COLUMN, SOURCE )  
    output_cat.write( OUTPUT_CAT )
    return

def Matching_A3COSMOS():
    
    print "Matching A3COSMOS catalog"
    
    ###set martirials
    path = "/Users/yuhsuan/Documents/research/05WH/data/COSMOS/"
    CAT  = path+"07_ALMA/apjsab42da_table4/A-COSMOS_blind.fits"
    
    RA            = "RA"
    DEC           = "DEC"
    SOURCE_COLUMN = "A3COSMOS"
    MAIN_CAT      = "COSMOS2015_Laigle+_v1.1_simple.fits" #COSMOS2015
    OUTPUT_CAT    = 'COSMOS2015_A3COSMOS.fits'
    
    searching_radius = 1.0
    
    ###set catalog names
    number      = 2
    cat_name    = [None]*number
    cat_name[0] = MAIN_CAT
    cat_name[1] = CAT
    
    ###set ID column names
    galaxyid    = [None]*number
    galaxyid[0] = "NUMBER"
    
    ##set RA DEC column names
    ra     = [None]*number
    dec    = [None]*number
    ra[0]  = "ALPHA_J2000"
    dec[0] = "DELTA_J2000"
    ra[1]  = RA
    dec[1] = DEC
    
    ###read catalog
    print "reading catalog ..."
    data = [None]*number
    for i in range(2): ReadCat(i,data,cat_name[i])
    
    time2 = time.time()
    print 'current time :', ( time2-time1 )/60.0 , 'min'
    
    ###matching
    SOURCE = [0]*len( data[0].filled() ) #create an array with length of COSMOS2015
    c0 = SkyCoord(ra=data[0][ra[0]], dec=data[0][dec[0]])
    c1 = SkyCoord(ra=data[1][ra[1]], dec=data[1][dec[1]])
    
    num = 0
    print "matching loop start ..."
    print 'CAT length ',len(data[1].filled())

    multimatched_num = 0
    for i in range( len(data[1].filled()) ):
        
        sep = c0.separation( c1[i] )
        data0_incircle = data[0][ sep<=searching_radius*u.arcsec ]
        
        if ( len(data0_incircle.filled())>1 ): print i, len(data0_incircle.filled())
        for j in range( len(data0_incircle.filled()) ):
            data0_id = data0_incircle[galaxyid[0]][j] -1
            if (SOURCE[ data0_id ] != 0):
                #print "!!!",data0_id,"!!!already tagged"
                multimatched_num +=1
            else:
                SOURCE[ data0_id ] = 1
                num += 1
                
    print "multimatched_num ",multimatched_num
    print "source number ",num
    
    time2 = time.time()
    print 'current time :', ( time2-time1 )/60.0 , 'min'
    
    ###output
    print "writing table file ..."
    output_cat = OutputCat( SOURCE_COLUMN, SOURCE )  
    output_cat.write( OUTPUT_CAT )
    return


def RegionFileUnmask(data,filename,mode,color,size,ra,dec):
    print 'Output region file...'
    f = open(filename+'.reg','w')
    f.write('global color='+color+' font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    
    N = len( data )
    
    print 'source number: ',N
    
    if (mode==0): #no label
        for n in range(N):
            f.write('fk5;circle('+str(data[ra][n])+','+str(data[dec][n])+','+size+'") # text={''}\n')
            #f.write('fk5;circle('+str(data[index][ra][mask[index]][n])+','+str(data[index][dec][mask[index]][n])+','+size+'") # text={''}\n')
    elif (mode==1): # redshift label
        for n in range(N):
            f.write('fk5;circle('+str(data[0][ra][mask[index]][n])+','+str(data[0][dec][mask[index]][n])+','+size+'") # text={'+str(data[0][photoz][mask[index]][n])+'}\n')
    
    #f.write('fk5;circle('+str(data[index][ra][mask[index]][n])+','+str(data[index][dec][mask[index]][n])+','+size+'") # text={''}\n')
    #f.write('fk5;circle('+str(data[index][ra][mask[index]][n])+','+str(data[index][dec][mask[index]][n])+','+size+'") # text={'+str(data[index]["NUMBER"][mask[index]][n])+'}\n')
    f.close()
    return

def A3COSMOS_single():
    
    print "A3COSMOS catalog multisource managing ..."
    
     ###set catalog names
    path = "/Users/yuhsuan/Documents/research/05WH/data/COSMOS/"
    number      = 2
    cat_name    = [None]*number
    cat_name[0] = path+"07_ALMA/apjsab42da_table4/A-COSMOS_blind.fits"
    cat_name[1] = path+"07_ALMA/AS2COSMOScatalog_2019-06-01_merged_wJMSphotom_wZspec.fits"    
    
    ##set RA DEC column names
    ra     = [None]*number
    dec    = [None]*number
    ra[0]  = "RA"
    dec[0] = "DEC"
    ra[1]  = "RA_final"
    dec[1] = "Dec_final" 

    OUTPUT_CAT    = 'A-COSMOS_blind_single.fits'
    
    searching_radius = 1.0
    
    ###read catalog
    print "reading catalog ..."
    data = [None]*2
    for i in range(2): ReadCat( i,data,cat_name[i] )
    
    time2 = time.time()
    print 'current time :', ( time2-time1 )/60.0 , 'min'
    
    c0 = SkyCoord( ra=data[0][ra[0]],          dec=data[0][dec[0]]          )
    c1 = SkyCoord( ra=data[1][ra[1]]*u.degree, dec=data[1][dec[1]]*u.degree )
    
    ###finding multiple source
    data[0]["index"] = np.arange( 0, len(data[0]), 1 )
    data[0]["flag"] = [0]*len( data[0] )
    c0 = SkyCoord(ra=data[0][ra[0]], dec=data[0][dec[0]])
    for i in range( len(data[0]) ):
        if data[0]["flag"][i]==0:
            sep = c0.separation( c0[i] )
            data0_incircle = data[0][ (sep<=searching_radius*u.arcsec)&(sep>0.0) ]
            for j in range( len(data0_incircle) ):
                data[0]["flag"][ data0_incircle["index"][j] ] = 1
            
            sep = c1.separation( c0[i] )
            data1_incircle = data[1][ (sep<=searching_radius*u.arcsec)&(sep>0.0) ]
            if ( len(data1_incircle)>0 ): data[0]["flag"][i] = 1
    
        
    mask = data[0]["flag"]==0 
    data_masked = data[0][mask]
    print len( data_masked )
    #RegionFileUnmask(data_masked,'TestA3Single',0,'red',"10.0",ra,dec)
    ###output
    print "writing table file ..."
    data_masked.write(OUTPUT_CAT)
    return

    
    
    

def Labeling_lensing():
    #id = 659416
    print 'lensing sample'
    ###set martirials
    ID = "NUMBER"
    SOURCE_COLUMN = "LENSING"
    MAIN_CAT = "COSMOS2015_Laigle+_v1.1_simple.fits" #COSMOS2015
    OUTPUT_CAT = 'COSMOS2015_lensing.fits'
    
    ###set catalog names
    number = 1
    cat_name = [None]*number
    cat_name[0] = MAIN_CAT
    
    ###set ID column names
    galaxyid = ID
    
    ###read catalog
    data = [None]*number
    for i in range(number):
        ReadCat(i,data,cat_name[i])
    
    ###label the sample
    SOURCE = [0]*len( data[0].filled() ) #create an array with length of COSMOS2015
    SOURCE[ 659416 - 1 ] = 1
    #print data[0]["PHOTOZ"][659416 - 1]
    
    ###output catalog
    print 'start writing table file'
    output_cat = OutputCat(SOURCE_COLUMN, SOURCE)  
    #print output_cat
    output_cat.write(OUTPUT_CAT)
    return


def TagSource( data_in_circle, j, tag ):
    
    source_id = data_in_circle[galaxyid][j] - 1
    if ( SOURCE[source_id]!=0 )&( SOURCE[source_id]!=7 ):
        print "SOURCE!!!",source_id,"already tagged",SOURCE[source_id],", would like to tag ",tag
    else:
        SOURCE[source_id] = tag
        label_num[tag-1] += 1
        
    return

def TagSource_ALMA( data_in_circle, j ):
    
    source_id = data_in_circle[galaxyid][j] - 1
    if ( SOURCE_DIRECT_ALMA[source_id]!=0 ):
        print "SOURCE_DIRECT_ALMA!!!",source_id,"already tagged",SOURCE_DIRECT_ALMA[source_id],", would like to tag 1"
    else:
        SOURCE_DIRECT_ALMA[source_id] = 1
        label_num2[0] += 1
        
    return

def TagSource_nonALMA( data_in_circle, j ):
    
    source_id = data_in_circle[galaxyid][j] - 1
    if ( SOURCE_DIRECT_NONALMA[source_id]!=0 ):
        print "SOURCE_DIRECT_NONALMA!!!",source_id,"already tagged",SOURCE_DIRECT_NONALMA[source_id],", would like to tag 1"
    else:
        SOURCE_DIRECT_NONALMA[source_id] = 1
        label_num2[1] += 1
        
    return


def MatchingSubmm():
    
    MAIN_CAT         = "COSMOS2015_merged_tmp.fits"
    CAT_24micron_all = path + "02_mips24/mips24_whwang.fits"
    CAT_3GHz_all     = path + "05_3GHz/vla3_cosmos_sources_160321_public5sig.fits.txt"
    CAT_ALMA_all     = path + "07_ALMA/AS2COSMOScatalog_2019-06-01_merged_wJMSphotom_wZspec.fits"
    CAT_A3COSMOS_all = 'A-COSMOS_blind_single.fits'#path + "07_ALMA/apjsab42da_table4/A-COSMOS_blind.fits"#
    
    ###set catalog names
    number      = 6
    cat_name    = [None]*number
    cat_name[0] = MAIN_CAT
    cat_name[1] = CAT
    
    cat_name[2] = CAT_24micron_all
    cat_name[3] = CAT_3GHz_all
    cat_name[4] = CAT_ALMA_all
    cat_name[5] = CAT_A3COSMOS_all
    
    ###set ID column names
    global galaxyid
    galaxyid = "NUMBER"
    
    ##set RA DEC column names
    ra     = [None]*number
    dec    = [None]*number
    ra[0]  = "ALPHA_J2000"
    dec[0] = "DELTA_J2000"
    ra[1]  = RA
    dec[1] = DEC
    ra[2]  = "RA_mips24"
    dec[2] = "DEC_mips24"
    ra[3]  = "ra"
    dec[3] = "dec"
    ra[4]  = "RA_final"
    dec[4] = "Dec_final"
    ra[5]  = "RA"
    dec[5] = "DEC"
    
    ###set columns in the main catalog
    global color1, color2, color3, photoz
    global V_error, ip_error, J_error, Ks_error, class_star, Ks_mag
    color1     = "MNUV"
    color2     = "MR"
    color3     = "MJ"
    photoz     = "REDSHIFT"
    V_error    = "V_MAGERR_APER3"
    ip_error   = "ip_MAGERR_APER3"
    J_error    = "J_MAGERR_APER3"
    Ks_error   = "Ks_MAGERR_APER3"
    class_star = "TYPE"
    Ks_mag     = "Ks_MAG_APER3"
    
    ###read catalog
    global data
    data = [None]*number
    print 'reading catalogs ...'
    for i in range( number ): ReadCat( i,data,cat_name[i] )
    
    time2 = time.time()
    print 'current time :', ( time2-time1 )/60.0 , 'min'
    
    print "len( data[5] )", len( data[5] )
    
    
    ###mask cat           
    data0_masked     = data[0][MaskAll(0)]
    data0_masked_24  = data[0][MaskAll(0)&(data[0]["24MICRON"]==1 )]
    data0_masked_3   = data[0][MaskAll(0)&(data[0]["3GHZ"]==1     )]
    data0_masked_AS2 = data[0][MaskAll(0)&(data[0]["AS2COSMOS"]==1)]
    data0_masked_A3  = data[0][MaskAll(0)&(data[0]["A3COSMOS"]==1 )]
    print "data0_masked     ",len(data0_masked)
    print "data0_masked_24  ",len(data0_masked_24)
    print "data0_masked_3   ",len(data0_masked_3)
    print "data0_masked_AS2 ",len(data0_masked_AS2)
    print "data0_masked_A3  ",len(data0_masked_A3)

    print 'masked MAIN_CAT length ',len( data0_masked.filled() )
    
    ###calculate coordinate
    c0 = SkyCoord( ra=data0_masked[ra[0]],     dec=data0_masked[dec[0]]     ) #MAIN_CAT
    c1 = SkyCoord( ra=data[1][ra[1]]*u.degree, dec=data[1][dec[1]]*u.degree ) #CAT
    c2 = SkyCoord( ra=data[2][ra[2]]*u.degree, dec=data[2][dec[2]]*u.degree ) #24
    c3 = SkyCoord( ra=data[3][ra[3]]         , dec=data[3][dec[3]]          ) #3
    c4 = SkyCoord( ra=data[4][ra[4]]*u.degree, dec=data[4][dec[4]]*u.degree ) #AS2COSMOS
    c5 = SkyCoord( ra=data[5][ra[5]]         , dec=data[5][dec[5]]          ) #A3COSMOS
    
    c6 = SkyCoord( ra=data0_masked_24[ra[0]] , dec=data0_masked_24[dec[0]]  ) #24
    c7 = SkyCoord( ra=data0_masked_3[ra[0]]  , dec=data0_masked_3[dec[0]]   ) #3
    c8 = SkyCoord( ra=data0_masked_AS2[ra[0]], dec=data0_masked_AS2[dec[0]] ) #AS2COSMOS_opt
    c9 = SkyCoord( ra=data0_masked_A3[ra[0]] , dec=data0_masked_A3[dec[0]]  ) #A3COSMOS_opt
    
    ###matching
    
    global SOURCE, SOURCE_DIRECT_ALMA, SOURCE_DIRECT_NONALMA, label_num, label_num2
    
    SOURCE                = [0]*len( data[0].filled() ) #create an array with length of COSMOS2015
    SOURCE_DIRECT_ALMA    = [0]*len( data[0].filled() )
    SOURCE_DIRECT_NONALMA = [0]*len( data[0].filled() )
    
    #all AS2COSMO S850 sources = ALMA_detect_Num + MIPSandVLA_detect_Num + MIPS_detect_Num
    #                            + VLA_detect_Num + MIPSorVLA_detect_no_opt_Num + non_detect_Num
    label_num = [ 0,0,0,0,0,0,0 ]

    ALMA_detect_Num              = 0 #with ALMA observation
    
    A3andAS2_detect_Num          = 0 #with A3 and AS2 observation  label1
    A3COSMOS_detect_Num          = 0 #only with A3 observation     label2
    AS2COSMOS_detect_Num         = 0 #only with AS2 observation    label3
    A3orAS2_detect_no_opt_Num    = 0 #but without opt samples
    
    MIPSorVLA_detect_Num         = 0
    
    MIPSandVLA_detect_Num        = 0 #with 24 and 3 observation    label4
    MIPS_detect_Num              = 0 #only with 24 observation     label5
    VLA_detect_Num               = 0 #only with 3 observation      label6
    MIPSorVLA_detect_no_opt_Num  = 0 #but without opt samples
    
    non_detect_Num               = 0 #no any other observation     label7: every opt sample within the radius
    non_detect_no_opt_Num        = 0 #even without any opt sample (outside COSMOS region)
    
    label_num2 = [ 0,0 ]
    
    ALMA_num = 0
    
    print "matching loop start ..."
    print 'CAT length ',len(data[1].filled())
    
    for i in range( len(data[1].filled()) ): #go through every id in CAT
        
        sep4            = c4.separation( c1[i] )
        data4_in_circle = data[4][ sep4<=searching_radius*u.arcsec ]
        sep5            = c5.separation( c1[i] )
        data5_in_circle = data[5][ sep5<=searching_radius*u.arcsec ]
        

        
        if ( len(data5_in_circle)>0 )|( len(data4_in_circle)>0 ): #ALMA detected
            
            #print len(data5_in_circle), len(data4_in_circle)
            ALMA_num += len(data5_in_circle)+len(data4_in_circle)
            
            # tag every opt sample within the radius
            sep0            = c0.separation( c1[i] )
            data0_in_circle = data0_masked[ sep0<=searching_radius*u.arcsec ]
                
            if ( len(data0_in_circle)>0 ):
                for j in range ( len(data0_in_circle) ): # go through every opt sources
                    TagSource_ALMA( data0_in_circle, j )
            
            
            
            ALMA_detect_Num +=1
            
            sep8            = c8.separation( c1[i] )
            data8_in_circle = data0_masked_AS2[ sep8<=searching_radius*u.arcsec ]
            sep9            = c9.separation( c1[i] )
            data9_in_circle = data0_masked_A3[ sep9<=searching_radius*u.arcsec ]
            
            if ( len(data9_in_circle)>0 )&( len(data8_in_circle)>0 ):
                
                label_flag = 0
                for j in range ( len(data9_in_circle) ): # go through every 24 opt sources
                    for k in range( len(data8_in_circle) ):
                        if ( data9_in_circle[galaxyid][j]==data8_in_circle[galaxyid][k] ):
                                
                            TagSource( data9_in_circle, j, 1 )
                            label_flag = 1

                if label_flag==1 :
                    A3andAS2_detect_Num  += 1
                    continue               
            
            if ( len(data9_in_circle)>0 )&( len(data8_in_circle)==0 ): #A3COSMOS detected
            
                A3COSMOS_detect_Num += 1         
                for j in range( len(data9_in_circle) ):
                    TagSource( data9_in_circle, j, 2 )
                            
            elif ( len(data9_in_circle)==0 )&( len(data8_in_circle)>0 ): #AS2COSMOS detected
                
                AS2COSMOS_detect_Num += 1           
                for j in range( len(data8_in_circle) ):
                    TagSource( data8_in_circle, j, 3 )
                            
            else: #no opt counterpart
                A3orAS2_detect_no_opt_Num += 1
            
      
            
        else: #no ALMA observation
            
            
            # tag every opt sample within the radius
            sep0            = c0.separation( c1[i] )
            data0_in_circle = data0_masked[ sep0<=searching_radius*u.arcsec ]
            
            if ( len(data0_in_circle)>0 ):
                for j in range ( len(data0_in_circle) ): # go through every opt sources
                    TagSource_nonALMA( data0_in_circle, j )
        
            MIPSorVLA_detect_Num += 1
            
            sep2            = c2.separation( c1[i] )
            data2_in_circle = data[2][ sep2<=searching_radius*u.arcsec ]
            sep3            = c3.separation( c1[i] )
            data3_in_circle = data[3][ sep3<=searching_radius*u.arcsec ]  
            
            if ( len(data2_in_circle)>0 )|( len(data3_in_circle)>0 ): # with 24 or 3 detected
                
                sep6            = c6.separation( c1[i] )
                data6_in_circle = data0_masked_24[ sep6<=searching_radius*u.arcsec ]
                sep7            = c7.separation( c1[i] )
                data7_in_circle = data0_masked_3[ sep7<=searching_radius*u.arcsec ]  
            
                if ( len(data6_in_circle)>0 )&( len(data7_in_circle)>0 ): # with 24 or 3 opt sources
                    
                    label_flag = 0
                    
                    for j in range ( len(data6_in_circle) ): # go through every 24 opt sources
                        for k in range( len(data7_in_circle) ):
                            if ( data6_in_circle[galaxyid][j]==data7_in_circle[galaxyid][k] ):
                                
                                TagSource( data6_in_circle, j, 4 )
                                label_flag = 1

                    if label_flag==1 :
                        MIPSandVLA_detect_Num += 1
                        continue
                                    
                if ( len(data6_in_circle)>0 )&( len(data7_in_circle)==0 ): # with 24 opt sources
                    
                    MIPS_detect_Num += 1                    
                    for j in range ( len(data6_in_circle) ): # go through every 24 opt sources                      
                        TagSource( data6_in_circle, j, 5 )
                            
                elif ( len(data6_in_circle)==0 )&( len(data7_in_circle)>0 ): # with 3 opt sources
                    
                    VLA_detect_Num += 1                    
                    for j in range ( len(data7_in_circle) ): # go through every 3 opt sources                       
                        TagSource( data7_in_circle, j, 6 )
                            
                else: # with 24 or 3 detected but no opt sample
                    
                    MIPSorVLA_detect_no_opt_Num += 1
                
              
            
            
            else: # no ALMA, 24, or 3 detection
            
                non_detect_Num += 1
                MIPSorVLA_detect_Num -= 1
                
                if ( len(data0_in_circle)>0 ):
                    for j in range ( len(data0_in_circle) ): # go through every opt sources
                        TagSource(  data0_in_circle, j, 7 ) 

                else:                   
                    non_detect_no_opt_Num += 1
                            
    print
    print "ALMA_detect_Num             ", ALMA_detect_Num
    print
    print "A3andAS2_detect_Num         ", A3andAS2_detect_Num
    print "opt sample num (1)          ", label_num[0]
    print "A3COSMOS_detect_Num         ", A3COSMOS_detect_Num
    print "opt sample num (2)          ", label_num[1]
    print "AS2COSMOS_detect_Num        ", AS2COSMOS_detect_Num
    print "opt sample num (3)          ", label_num[2]
    print "A3orAS2_detect_no_opt_Num   ", A3orAS2_detect_no_opt_Num
    print
    print "MIPSorVLA_detect_Num        ", MIPSorVLA_detect_Num
    print
    print "MIPSandVLA_detect_Num       ", MIPSandVLA_detect_Num
    print "opt sample num (4)          ", label_num[3]
    print "MIPS_detect_Num             ", MIPS_detect_Num
    print "opt sample num (5)          ", label_num[4]
    print "VLA_detect_Num              ", VLA_detect_Num
    print "opt sample num (6)          ", label_num[5]
    print "MIPSorVLA_detect_no_opt_Num ", MIPSorVLA_detect_no_opt_Num
    print
    print "non_detect_Num              ", non_detect_Num
    print "non_detect_no_opt_Num       ", non_detect_no_opt_Num
    print "opt sample num (7)          ", label_num[6]
    print
    
    print
    print "direct matching"
    print "opt sample num (ALMA)       ", label_num2[0]
    print "opt sample num (nonALMA)    ", label_num2[1]
    print
    print "ALMA num                    ", ALMA_num
    print "average ALMA num            ", ( "%.2f" % (ALMA_num/ALMA_detect_Num) )
    print 
    ###output
    print 'writing table file ...'
    output_cat = OutputCat( SOURCE_COLUMN, SOURCE ) 
    
    output_cat[SOURCE_COLUMN+"_DIRECT_ALMA"   ]  = SOURCE_DIRECT_ALMA
    output_cat[SOURCE_COLUMN+"_DIRECT_NONALMA"]  = SOURCE_DIRECT_NONALMA
    
    
    output_cat.write( OUTPUT_CAT )
    print
    return

def Matching_850wide():
    
    print "Matching 850 wide sources"
    
    global path, CAT, RA, DEC, SOURCE_COLUMN, MAIN_CAT, OUTPUT_CAT, searching_radius
    
    ###set martirials
    path = "/Users/yuhsuan/Documents/research/05WH/data/COSMOS/"
    
    CAT           = path + "04_COSMOS450_850/S2COSMOS/catalog/S2COSMOS_sourcecat850_Simpson18.fits"
    RA            = "RA_deg"
    DEC           = "DEC_deg"
    SOURCE_COLUMN = "850WIDE"
    
    OUTPUT_CAT       = "COSMOS2015_850wide.fits"    
    searching_radius = 7.0

    MatchingSubmm()
    
    return


def Matching_450narrow():
    
    print "Matching 450 narrow sources"
    
    global path, CAT, RA, DEC, SOURCE_COLUMN, MAIN_CAT, OUTPUT_CAT, searching_radius
    
    ###set martirials
    path = "/Users/yuhsuan/Documents/research/05WH/data/COSMOS/"
    
    CAT           = path + "04_COSMOS450_850/STUDIES/sources_450.fits"
    RA            = "RA_450"
    DEC           = "DEC_450"
    SOURCE_COLUMN = "450NARROW"
    
    OUTPUT_CAT       = "COSMOS2015_450narrow.fits"
    searching_radius = 4.0
    
    MatchingSubmm()
    
    return

def Mask_450Region( data, ra, dec ): 
    
    c0     = SkyCoord( ra=data[ra], dec=data[dec] )#*u.degree
    center = SkyCoord( '10h00m25.0s', '2d24m22.0s', frame='fk5' )
    radius = 0.2*u.degree
    sep    = center.separation( c0 )
    
    return sep<=radius

def cat_outer():
    
    ###set catalog names
    number      = 2
    cat_name    = [None]*number
    cat_name[0] = MAIN_CAT
    cat_name[1] = CAT
    
    ###set ID column names
    galaxyid    = [None]*number
    galaxyid[0] = "NUMBER"
    
    ##set RA DEC column names
    ra     = [None]*number
    dec    = [None]*number
    ra[0]  = "ALPHA_J2000"
    dec[0] = "DELTA_J2000"
    ra[1]  = RA
    dec[1] = DEC
    
    ###read catalog
    print "reading catalog ..."
    data = [None]*number
    for i in range(2): ReadCat(i,data,cat_name[i])
    
    time2 = time.time()
    print 'current time :', ( time2-time1 )/60.0 , 'min'
    
    mask0 = ((data[0]["FLAG_HJMCC"]==0) | (data[0]["FLAG_HJMCC"]==2)) &(data[0]["FLAG_COSMOS"]==1)
    data0_masked = data[0][mask0]
    mask1 = Mask_450Region( data[1], ra[1], dec[1] )
    data1_masked = data[1][mask1]
    
    ###matching
    c0 = SkyCoord( ra=data0_masked[ra[0]],          dec=data0_masked[dec[0]]  )
    c1 = SkyCoord( ra=data1_masked[ra[1]], dec=data1_masked[dec[1]] ) #*u.degree
    
    num = 0
    print "matching loop start ..."
    print 'CAT length ',len( data1_masked )

    for i in range( len( data1_masked ) ):
        
        sep = c0.separation( c1[i] )
        data0_incircle = data0_masked[ sep<=searching_radius*u.arcsec ]
        
        if ( len(data0_incircle.filled())==0 ): num+=1
                
    print "num ",num
    
    time2 = time.time()
    print 'current time :', ( time2-time1 )/60.0 , 'min'
    
    return

def S2COSMOS_outer():
    
    print "Number of S2COSMOS in the outer region"
    
    ###set martirials
    global CAT, RA, DEC, MAIN_CAT, searching_radius
    path = "/Users/yuhsuan/Documents/research/05WH/data/COSMOS/"
    CAT  = path + "04_COSMOS450_850/S2COSMOS/catalog/S2COSMOS_sourcecat850_Simpson18.fits"
    RA            = "RA_deg"
    DEC           = "DEC_deg"
    MAIN_CAT      = "COSMOS2015_Laigle+_v1.1_simple.fits" #COSMOS2015
    searching_radius = 7.0
    
    cat_outer()
    
    return

def STUDIES_outer():
    
    print "Number of STUDIES in the outer region"
    
    ###set martirials
    global CAT, RA, DEC, MAIN_CAT, searching_radius
    path = "/Users/yuhsuan/Documents/research/05WH/data/COSMOS/"
    CAT  = path + "04_COSMOS450_850/STUDIES/sources_450.fits"
    RA            = "RA_450"
    DEC           = "DEC_450"
    MAIN_CAT      = "COSMOS2015_Laigle+_v1.1_simple.fits" #COSMOS2015
    searching_radius = 7.0
    
    cat_outer()
    
    return

def A3COSMOS_outer():

    print "Number of A3COSMOS in the outer region"
    
    ###set martirials
    global CAT, RA, DEC, MAIN_CAT, searching_radius
    CAT  = "A-COSMOS_blind_single.fits"
    RA            = "RA"
    DEC           = "DEC"
    MAIN_CAT      = "COSMOS2015_Laigle+_v1.1_simple.fits" #COSMOS2015
    searching_radius = 7.0
    
    cat_outer()
    
    return

def AS2COSMOS_outer():
    
    print "Number of AS2COSMOS in the outer region"
    
    ###set martirials
    global CAT, RA, DEC, MAIN_CAT, searching_radius
    path = "/Users/yuhsuan/Documents/research/05WH/data/COSMOS/"
    CAT  = path + "07_ALMA/AS2COSMOScatalog_2019-06-01_merged_wJMSphotom_wZspec.fits"
    RA            = "RA_final"
    DEC           = "Dec_final"
    MAIN_CAT      = "COSMOS2015_Laigle+_v1.1_simple.fits" #COSMOS2015
    searching_radius = 7.0
    
    cat_outer()
    
    return



def Random_matching():
    
    print "Random matching "
    
    ALMA_MODE = 0  # 0: 4 and 7 radius, 1: 1 radius
    MODE = 0       # 0: 450um,          1: 850um
    
    ###set martirials
    CAT  = "COSMOS2015_Laigle+_v1.1_simple_z.fits" #COSMOS2015
    
    if MODE==0:
        print "450 um"
        searching_radius = 4.0/60.0/60.0
    else:
        print "850 um"
        searching_radius = 7.0/60.0/60.0
        
    area_searching_radius = 7.0/60.0/60.0
    
    if ALMA_MODE==1: 
        searching_radius = 1.0/60.0/60.0
        print "ALMA"
        
    print "searching radius ", searching_radius*60.0*60.0
    
    ra = "ALPHA_J2000"
    dec = "DELTA_J2000"
    
    ###set columns in the main catalog
    global color1, color2, color3, photoz
    global V_error, ip_error, J_error, Ks_error, class_star, Ks_mag
    color1     = "MNUV"
    color2     = "MR"
    color3     = "MJ"
    photoz     = "REDSHIFT"
    V_error    = "V_MAGERR_APER3"
    ip_error   = "ip_MAGERR_APER3"
    J_error    = "J_MAGERR_APER3"
    Ks_error   = "Ks_MAGERR_APER3"
    class_star = "TYPE"
    Ks_mag     = "Ks_MAG_APER3"
    
    ###read catalog
    print "reading catalog ..."
    global data,x,y,mask,x_masked,y_masked
    data = [None]*1
    x,y,mask,x_masked,y_masked = [],[],[],[],[]
    
    for i in range(1): ReadCat(i,data,CAT)
    
    time2 = time.time()
    print 'current time :', ( time2-time1 )/60.0 , 'min'
    
    SetupData( MaskAll(0) )
    SetupData( ((data[0]["FLAG_HJMCC"]==0)|(data[0]["FLAG_HJMCC"]==2)) & (data[0]["FLAG_COSMOS"]==1) )
    SetupData( MaskAll(0) & Mask_myclassQG(0)  )
    SetupData( MaskAll(0) & Mask_myclassSFG(0) )


    c_area      = np.array( [data[0][mask[1]][ra],data[0][mask[1]][dec]] )
    c_area      = c_area.T
    c_area_tree = cKDTree(c_area)

    c_QG        = np.array( [data[0][mask[2]][ra],data[0][mask[2]][dec]] )
    c_QG        = c_QG.T
    c_QG_tree   = cKDTree(c_QG)
    
    c_SFG       = np.array( [data[0][mask[3]][ra],data[0][mask[3]][dec]] )
    c_SFG       = c_SFG.T
    c_SFG_tree  = cKDTree(c_SFG)
    
    if MODE==1 :
        ra_min  = data[0][ra ][mask[0]].min()
        ra_max  = data[0][ra ][mask[0]].max()
        dec_min = data[0][dec][mask[0]].min()
        dec_max = data[0][dec][mask[0]].max()
    else:
        center450 = SkyCoord('10h00m25.0s', '2d24m22.0s', frame='fk5')
    
    if MODE==0 : pos_num_list = [353,78,275]
    else: pos_num_list = [981,391,590]   
    
    if ALMA_MODE==1:
        if MODE==0 : pos_num_list = [85]
        else: pos_num_list = [452]   
    
    for n in range( len(pos_num_list) ):
        
        test_num = 1000
        matchedQG = np.zeros(test_num)
        matchedSFG = np.zeros(test_num)
        pos_num = pos_num_list[n]
        print
        print "pos_num ", pos_num
        
        zero_alert = 0
        
        for i in range(test_num):
            
            if MODE==1 :
                random1 = np.random.uniform(0,1,pos_num)
                random2 = np.random.uniform(0,1,pos_num)
                random_ra  = random1*(ra_max-ra_min) + ra_min
                random_dec = random2*(dec_max-dec_min) + dec_min
            else:
                random_ang = np.random.uniform(0,2*np.pi,pos_num)
                random_dis = 0.2*np.sqrt(np.random.uniform(0,1,pos_num))
                random_ra  = center450.ra.value  + random_dis*np.cos(random_ang)
                random_dec = center450.dec.value + random_dis*np.sin(random_ang)
            
            random_pos = np.array( [random_ra,random_dec] )
            random_pos = random_pos.T
            
            #print random_pos
            
            random_pos_tree = cKDTree(random_pos)
            
            area_in_circle = random_pos_tree.query_ball_tree(c_area_tree, r=area_searching_radius)
            
            for j in range( len(area_in_circle) ):
                #print len(area_in_circle[i])
                onearea_in_circle = len(area_in_circle[j])
                
                while(1):
                    if ( onearea_in_circle==0 ):
                        
                        #print random_pos[j][0], random_pos[j][1]
                        
                        zero_alert +=1
                        
                        if MODE==1 :
                            onerandom1 = np.random.uniform(0,1,1)
                            onerandom2 = np.random.uniform(0,1,1)
                            onerandom_ra  = onerandom1*(ra_max-ra_min) + ra_min
                            onerandom_dec = onerandom2*(dec_max-dec_min) + dec_min
                        else:
                            onerandom_ang = np.random.uniform(0,2*np.pi,1)
                            onerandom_dis = 0.2*np.sqrt(np.random.uniform(0,1,1))
                            onerandom_ra  = center450.ra.value  + onerandom_dis*np.cos(onerandom_ang)
                            onerandom_dec = center450.dec.value + onerandom_dis*np.sin(onerandom_ang)

                        random_pos[j][0] = onerandom_ra
                        random_pos[j][1] = onerandom_dec
                        
                        #nerandom_pos = np.array(random_pos[i])
                        onerandom_pos_tree = cKDTree( np.array( [random_pos[j]] ) )
                        onearea_in_circle = onerandom_pos_tree.query_ball_tree(c_area_tree, r=area_searching_radius)
                        
                    else: break
                        
            #print
            #random_pos_tree = cKDTree(random_pos)
            #area_in_circle = random_pos_tree.query_ball_tree(c_area_tree, r=area_searching_radius)
            #for i in range( len(area_in_circle) ):
            #    print len(area_in_circle[i])
            
            random_pos_tree = cKDTree(random_pos)
            
            QG_in_circle  = random_pos_tree.query_ball_tree(c_QG_tree,  r=searching_radius)
            SFG_in_circle = random_pos_tree.query_ball_tree(c_SFG_tree, r=searching_radius)
            
            for j in range( len(QG_in_circle) ): matchedQG[i] += len(QG_in_circle[j])
            for j in range( len(SFG_in_circle) ): matchedSFG[i] += len(SFG_in_circle[j])
            

        
        if MODE==0:
            actual_num_QG  = [29,5,24]
            actual_num_SFG = [420,83,337]
            
        else:
            actual_num_QG  = [191,77,114]
            actual_num_SFG = [1993,801,1192]
        
        if ALMA_MODE==1:
            if MODE==0:
                actual_num_QG  = [2]
                actual_num_SFG = [47]
            else:
                actual_num_QG  = [9]
                actual_num_SFG = [223]
        
        
        print "zero alert num ",zero_alert
        print
        print  "QG  actual ", actual_num_QG[n]
        print  "QG  mean ",np.mean(matchedQG) 
        print  "QG  std  ",( "%.2f" % np.std(matchedQG) )
        print  "QG  16%  ",np.percentile(matchedQG, 16)
        print  "QG  84%  ",np.percentile(matchedQG, 84)
        print  "QG  number      ", np.count_nonzero(matchedQG == actual_num_QG[n])
        print  "QG  possibility ",float( np.count_nonzero(matchedQG == actual_num_QG[n]) )/1000.0
        
        print
        print  "SFG actual ", actual_num_SFG[n]
        print  "SFG mean ",np.mean(matchedSFG) 
        print  "SFG std  ",( "%.2f" % np.std(matchedSFG) )
        print  "SFG 16%  ",np.percentile(matchedSFG, 16)
        print  "SFG 84%  ",np.percentile(matchedSFG, 84)
        print  "SFG number      ", np.count_nonzero(matchedSFG == actual_num_SFG[n])
        print  "SFG possibility ",float( np.count_nonzero(matchedSFG == actual_num_SFG[n]) )/1000.0
        
        print
        #print np.mean(matchedQG),"$\pm$",np.std(matchedQG)," & ",np.mean(matchedSFG),"$\pm$",np.std(matchedSFG)
        print
        
        
        #from scipy.optimize import curve_fit
        #from scipy.stats import poisson
 
        #def fit_function(k, lamb, scale):
        #    return poisson.pmf(k, lamb) * scale    
        
        
        
        plt.figure()
        
        if ALMA_MODE==1: input_bin = [0,1,2,3,4,5,6]
        else: input_bin = range(int(min(matchedQG)),int(max(matchedQG))+2,2)
        
        data_entries, bins, _ = plt.hist( matchedQG, bins=input_bin, label="expect", density=True )
        #binscenters = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
        #popt, pcov = curve_fit(fit_function, xdata=binscenters, ydata=data_entries)
        #print(popt)
        #print data_entries
        #print bins
        
        #plt.plot( bins, fit_function(bins, *popt), label="poisson, lambda = "+ ( "%.1f"%popt[0] ) )
        
        if actual_num_QG[n]<max(matchedQG): plt.axvline( actual_num_QG[n], color='k', linestyle='--' )
        plt.axvline( np.mean(matchedQG), color='k' )
        
        plt.title("QG, actual num: "+str(actual_num_QG[n]))
        plt.xlabel("matched QGs")
        plt.ylabel("number")
        plt.legend()

        
        plt.figure()
        
        input_bin = range(int(min(matchedSFG)),int(max(matchedSFG))+10,10)
        
        data_entries, bins, _ = plt.hist( matchedSFG, bins=input_bin, label="expect", density=True )
        #binscenters = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
        #popt, pcov = curve_fit(fit_function, xdata=binscenters, ydata=data_entries)
        #print(popt)

        
        #plt.plot( bins, fit_function(bins, *popt), label="poisson, lambda = "+ ( "%.1f"%popt[0] ) )


        if actual_num_SFG[n]<max(matchedSFG): plt.axvline( actual_num_SFG[n], color='k', linestyle='--' )
        plt.axvline( np.mean(matchedSFG), color='k' )
        
        plt.title("SFG, actual num: "+str(actual_num_SFG[n]))
        plt.xlabel("matched SFGs")
        plt.ylabel("number")
        plt.legend()
         
    return

def CalDist():
    
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    dist_arcmin = cosmo.kpc_proper_per_arcmin(1)
    print dist_arcmin*7.0/60.0
    print dist_arcmin*4.0/60.0
    print dist_arcmin*1.0/60.0
    
    return

def Plot_Random_matching():
    
    label_list = [ "all", "with ALMA", "without ALMA", "ALMA source"]
    
    QG_actual_450 = np.array( [29,5,24,2] )
    QG_mean_450 = np.array( [20.02,4.30,15.53,0.28] )
    QG_16p_450 = np.array( [15.0,2.0,12.0,0.0] )
    QG_84p_450 = np.array( [25.0,6.0,19.0,1.0] )
    QG_std_450 = np.array( [4.76,2.08,3.97,0.51] )
    
    print QG_mean_450
    print QG_84p_450 - QG_mean_450
    print QG_mean_450 - QG_16p_450
    print np.array( [(QG_mean_450-QG_16p_450)/QG_mean_450 *100,(QG_84p_450-QG_mean_450)/QG_mean_450 *100]  )
    print
    
    SFG_actual_450 = np.array( [420,83,337,47] )
    SFG_mean_450 = np.array( [163.05,36.19,127.88 ,2.50] )
    SFG_16p_450 = np.array( [149.84,30.0,116.0 ,1.0] )
    SFG_84p_450 = np.array( [177.0,42.0,139.0,4.0] )
    SFG_std_450 = np.array( [13.22,6.15,11.71,1.52] )
    
    print SFG_mean_450
    print SFG_84p_450 - SFG_mean_450
    print SFG_mean_450 - SFG_16p_450
    print np.array( [(SFG_mean_450-SFG_16p_450)/SFG_mean_450 *100,(SFG_84p_450-SFG_mean_450)/SFG_mean_450 *100]  )
    print
    
    QG_actual_850 = np.array( [191,77,114,9] )
    QG_mean_850 = np.array( [117.36,47.21,71.09,1.096] )
    QG_16p_850 = np.array( [105.0,40.0,62.0,0.0] )
    QG_84p_850 = np.array( [129.0,54.0,80.0,2.0] )
    QG_std_850 = np.array( [12.16,7.17,9.3,1.06] )

    print QG_mean_850
    print QG_84p_850 - QG_mean_850
    print QG_mean_850 - QG_16p_850
    print np.array( [(QG_mean_850-QG_16p_850)/QG_mean_850 *100,(QG_84p_850-QG_mean_850)/QG_mean_850 *100]  )
    print

    SFG_actual_850 = np.array( [1993,801,1192,223] )
    SFG_mean_850 = np.array( [1050.89,420.42,630.64,9.78] )
    SFG_16p_850 = np.array( [1014.68,398.0,600.0,7.0] )
    SFG_84p_850 = np.array( [1086.16,443.16,660.0,13.0] )
    SFG_std_850 = np.array( [36.43,23.00,29.72,3.11] )
    
    print SFG_mean_850
    print SFG_84p_850 - SFG_mean_850
    print SFG_mean_850 - SFG_16p_850
    print np.array( [(SFG_mean_850-SFG_16p_850)/SFG_mean_850 *100,(SFG_84p_850-SFG_mean_850)/SFG_mean_850 *100]  )
    print
    
    '''
    fig, axes = plt.subplots()
    
    plt.errorbar( [1,2,3,4], (QG_actual_450-QG_mean_450)/QG_mean_450 *100,
                 np.array( [(QG_mean_450-QG_16p_450)/QG_mean_450 *100,(QG_84p_450-QG_mean_450)/QG_mean_450 *100]  ),
                 color='r', fmt='o', capsize=5 ,alpha=0.8 )
    plt.errorbar( [1,2,3,4], (SFG_actual_450-SFG_mean_450)/SFG_mean_450 *100,
                 np.array( [(SFG_mean_450-SFG_16p_450)/SFG_mean_450 *100,(SFG_84p_450-SFG_mean_450)/SFG_mean_450 *100]  ),
                 color='b', fmt='D', capsize=5 ,alpha=0.8 )
    
    patch1 = mpatches.Patch( color='r', label='QG' )
    patch2 = mpatches.Patch( color='b', label='SFG' )
    plt.legend( handles=[patch1,patch2],loc=4 )
    
    plt.setp( axes, xticks=[y_axes+1 for y_axes in range(4)], xticklabels=label_list )
    plt.ylabel('difference (%)', fontdict = {'fontsize' : 14})
    
    plt.axis([0.2,3.8,-50,200])
    fig.savefig('spatial_4.png', bbox_inches = 'tight', format='png', dpi=400)
    
    fig, axes = plt.subplots()
    
    plt.errorbar( [1,2,3,4], (QG_actual_850-QG_mean_850)/QG_mean_850 *100,
                 np.array( [(QG_mean_850-QG_16p_850)/QG_mean_850 *100,(QG_84p_850-QG_mean_850)/QG_mean_850 *100]  ),
                 color='r', fmt='o', capsize=5 ,alpha=0.8 )
    plt.errorbar( [1,2,3,4], (SFG_actual_850-SFG_mean_850)/SFG_mean_850 *100,
                 np.array( [(SFG_mean_850-SFG_16p_850)/SFG_mean_850 *100,(SFG_84p_850-SFG_mean_850)/SFG_mean_850 *100]  ),
                 color='b', fmt='D', capsize=5 ,alpha=0.8 )
    
    patch1 = mpatches.Patch( color='r', label='QG' )
    patch2 = mpatches.Patch( color='b', label='SFG' )
    plt.legend( handles=[patch1,patch2],loc=4 )
    
    plt.setp( axes, xticks=[y_axes+1 for y_axes in range(4)], xticklabels=label_list )
    plt.ylabel('difference (%)', fontdict = {'fontsize' : 14})
    
    plt.axis([0.2,3.8,0,100])
    fig.savefig('spatial_7.png', bbox_inches = 'tight', format='png', dpi=400)
    
    #print (QG_actual_450[3]-QG_mean_450[3])/QG_mean_450[3] *100
    
    fig, axes = plt.subplots()
    
    plt.errorbar( 1, (QG_actual_450[3]-QG_mean_450[3])/QG_mean_450[3] *100,
                 np.array([[(QG_mean_450[3]-QG_16p_450[3])/QG_mean_450[3] *100,(QG_84p_450[3]-QG_mean_450[3])/QG_mean_450[3] *100]]  ).T,
                 color='r', fmt='o', capsize=5 ,alpha=0.8 )
    plt.errorbar( 1, (SFG_actual_450[3]-SFG_mean_450[3])/SFG_mean_450[3] *100,
                 np.array([[(SFG_mean_450[3]-SFG_16p_450[3])/SFG_mean_450[3] *100,(SFG_84p_450[3]-SFG_mean_450[3])/SFG_mean_450[3] *100]]  ).T,
                 color='b', fmt='D', capsize=5 ,alpha=0.8 )
    plt.errorbar( 2, (QG_actual_850[3]-QG_mean_850[3])/QG_mean_850[3] *100,
                 np.array([[(QG_mean_850[3]-QG_16p_850[3])/QG_mean_850[3] *100,(QG_84p_850[3]-QG_mean_850[3])/QG_mean_850[3] *100]]  ).T,
                 color='r', fmt='o', capsize=5 ,alpha=0.8 )
    plt.errorbar( 2, (SFG_actual_850[3]-SFG_mean_850[3])/SFG_mean_850[3] *100,
                 np.array([[(SFG_mean_850[3]-SFG_16p_850[3])/SFG_mean_850[3] *100,(SFG_84p_850[3]-SFG_mean_850[3])/SFG_mean_850[3] *100] ] ).T,
                 color='b', fmt='D', capsize=5 ,alpha=0.8 )

    patch1 = mpatches.Patch( color='r', label='QG' )
    patch2 = mpatches.Patch( color='b', label='SFG' )
    plt.legend( handles=[patch1,patch2],loc=4 )
    
    plt.setp( axes, xticks=[y_axes+1 for y_axes in range(2)], xticklabels=[450,850] )
    plt.ylabel('difference (%)', fontdict = {'fontsize' : 14})
    
    plt.axis([0.2,2.8,250,2500])
    
    fig.savefig('spatial_1.png', bbox_inches = 'tight', format='png', dpi=400)
    
    '''
    
    '''
    fig, axes = plt.subplots()
    
    plt.errorbar( [1,2,3,4], (QG_actual_450-QG_mean_450)/QG_mean_450 *100,
                 QG_std_450/QG_mean_450 *100,
                 color='r', fmt='o', capsize=5 ,alpha=0.8 )
    plt.errorbar( [1,2,3,4], (SFG_actual_450-SFG_mean_450)/SFG_mean_450 *100,
                 SFG_std_450/SFG_mean_450 *100,
                 color='b', fmt='D', capsize=5 ,alpha=0.8 )
    
    
    plt.setp( axes, xticks=[y_axes+1 for y_axes in range(4)], xticklabels=label_list )
    plt.ylabel('difference (%)', fontdict = {'fontsize' : 14})
    
    plt.axis([0.2,3.8,-50,200])
    
    
    fig, axes = plt.subplots()
    
    plt.errorbar( [1,2,3,4], (QG_actual_850-QG_mean_850)/QG_mean_850 *100,
                 QG_std_850/QG_mean_850 *100,
                 color='r', fmt='o', capsize=5 ,alpha=0.8 )
    plt.errorbar( [1,2,3,4], (SFG_actual_850-SFG_mean_850)/SFG_mean_850 *100,
                 SFG_std_850/SFG_mean_850 *100,
                 color='b', fmt='D', capsize=5 ,alpha=0.8 )
    
    
    plt.setp( axes, xticks=[y_axes+1 for y_axes in range(4)], xticklabels=label_list )
    plt.ylabel('difference (%)', fontdict = {'fontsize' : 14})
    
    plt.axis([0.2,3.8,0,100])
    
    print (QG_actual_450[3]-QG_mean_450[3])/QG_mean_450[3] *100
    
    fig, axes = plt.subplots()
    
    plt.errorbar( 1, (QG_actual_450[3]-QG_mean_450[3])/QG_mean_450[3] *100,
                 QG_std_450[3]/QG_mean_450[3] *100,
                 color='r', fmt='o', capsize=5 ,alpha=0.8 )
    plt.errorbar( 1, (SFG_actual_450[3]-SFG_mean_450[3])/SFG_mean_450[3] *100,
                 SFG_std_450[3]/SFG_mean_450[3] *100,
                 color='b', fmt='D', capsize=5 ,alpha=0.8 )
    plt.errorbar( 2, (QG_actual_850[3]-QG_mean_850[3])/QG_mean_850[3] *100,
                 QG_std_850[3]/QG_mean_850[3] *100,
                 color='r', fmt='o', capsize=5 ,alpha=0.8 )
    plt.errorbar( 2, (SFG_actual_850[3]-SFG_mean_850[3])/SFG_mean_850[3] *100,
                 SFG_std_850[3]/SFG_mean_850[3] *100,
                 color='b', fmt='D', capsize=5 ,alpha=0.8 )
    
    plt.setp( axes, xticks=[y_axes+1 for y_axes in range(2)], xticklabels=[450,850] )
    plt.ylabel('difference (%)', fontdict = {'fontsize' : 14})
    
    plt.axis([0.2,2.8,250,2500])





    fig, axes = plt.subplots()
    
    plt.errorbar( [1,2,3,4], (QG_actual_450-QG_mean_450)/QG_mean_450 *100,
                 np.sqrt(QG_mean_450)/QG_mean_450 *100,
                 color='r', fmt='o', capsize=5 ,alpha=0.8 )
    plt.errorbar( [1,2,3,4], (SFG_actual_450-SFG_mean_450)/SFG_mean_450 *100,
                 np.sqrt(SFG_mean_450)/SFG_mean_450 *100,
                 color='b', fmt='D', capsize=5 ,alpha=0.8 )
    
    
    plt.setp( axes, xticks=[y_axes+1 for y_axes in range(4)], xticklabels=label_list )
    plt.ylabel('difference (%)', fontdict = {'fontsize' : 14})
    
    plt.axis([0.2,3.8,-50,200])
    
    
    fig, axes = plt.subplots()
    
    plt.errorbar( [1,2,3,4], (QG_actual_850-QG_mean_850)/QG_mean_850 *100,
                 np.sqrt(QG_mean_850)/QG_mean_850 *100,
                 color='r', fmt='o', capsize=5 ,alpha=0.8 )
    plt.errorbar( [1,2,3,4], (SFG_actual_850-SFG_mean_850)/SFG_mean_850 *100,
                 np.sqrt(SFG_mean_850)/SFG_mean_850 *100,
                 color='b', fmt='D', capsize=5 ,alpha=0.8 )
    
    
    plt.setp( axes, xticks=[y_axes+1 for y_axes in range(4)], xticklabels=label_list )
    plt.ylabel('difference (%)', fontdict = {'fontsize' : 14})
    
    plt.axis([0.2,3.8,0,100])
    
    print (QG_actual_450[3]-QG_mean_450[3])/QG_mean_450[3] *100
    
    fig, axes = plt.subplots()
    
    plt.errorbar( 1, (QG_actual_450[3]-QG_mean_450[3])/QG_mean_450[3] *100,
                 np.sqrt(QG_mean_450[3])/QG_mean_450[3] *100,
                 color='r', fmt='o', capsize=5 ,alpha=0.8 )
    plt.errorbar( 1, (SFG_actual_450[3]-SFG_mean_450[3])/SFG_mean_450[3] *100,
                np.sqrt(SFG_mean_450[3])/SFG_mean_450[3] *100,
                 color='b', fmt='D', capsize=5 ,alpha=0.8 )
    plt.errorbar( 2, (QG_actual_850[3]-QG_mean_850[3])/QG_mean_850[3] *100,
                 np.sqrt(QG_mean_850[3])/QG_mean_850[3] *100,
                 color='r', fmt='o', capsize=5 ,alpha=0.8 )
    plt.errorbar( 2, (SFG_actual_850[3]-SFG_mean_850[3])/SFG_mean_850[3] *100,
                 np.sqrt(SFG_mean_850[3])/SFG_mean_850[3] *100,
                 color='b', fmt='D', capsize=5 ,alpha=0.8 )
    
    plt.setp( axes, xticks=[y_axes+1 for y_axes in range(2)], xticklabels=[450,850] )
    plt.ylabel('difference (%)', fontdict = {'fontsize' : 14})
    
    plt.axis([0.2,2.8,250,2500])
    
    '''
    
    return
    

# =============================================================================
# main code
# =============================================================================

time1 = time.time()


# ===== 24 micorn =====
#Matching_24Micron()
#print ''

# ===== 3 GHz =====
#Matching_3GHz()
#print ''

# ===== 450 narrow - direct =====
# ===== 850 narrow - direct =====

# ===== ALMA =====
#Matching_AS2COSMOS()

# ===== ALMA =====
#Matching_A3COSMOS()

# ===== lensing =====
#id = 659416
#Labeling_lensing()

# ===== A3COSMOS multisource issue =====
#A3COSMOS_single()

# ===== 850 wide =====
#COSMOS2015_24micron,3GHz,ALMA needed
#Matching_850wide()

# ===== 450 narrow =====
#COSMOS2015_24micron,3GHz,ALMA needed
#Matching_450narrow()

#S2COSMOS_outer()
#STUDIES_outer()
#A3COSMOS_outer()
#AS2COSMOS_outer()
#A3COSMOS_outer()

Random_matching()
#Plot_Random_matching()

#CalDist()

time2 = time.time()
print 'done! time :', (time2-time1)/60.0 , 'min'