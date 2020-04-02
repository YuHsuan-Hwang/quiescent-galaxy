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
def ReadCat(index,data,cat):
    data[index] = Table.read(cat, hdu=1)
    return

def OutputCat(column_name, column_list):
    t = Table()
    t[column_name] = column_list
    return t

def Mask_M(index,data,color1,color2,color3):
    down = -30
    up = 0
    mask_color1 = (data[index][color1]>down) & (data[index][color1]<up)
    mask_color2 = (data[index][color2]>down) & (data[index][color2]<up)
    mask_color3 = (data[index][color3]>down) & (data[index][color3]<up)
    mask_M = mask_color1 & mask_color2 & mask_color3
    return mask_M

def Mask_error(index,data,V_error,ip_error,J_error,Ks_error):
    mask_V_error = (data[index][V_error]<0.1) & (data[index][V_error]>0)
    #Suprime-Cam:i+band
    mask_ip_error = (data[index][ip_error]<0.1) & (data[index][ip_error]>0)
    mask_J_error = (data[index][J_error]<0.1) & (data[index][J_error]>0)
    mask_Ks_error = (data[index][Ks_error]<0.1) & (data[index][Ks_error]>0)
    return mask_V_error | mask_ip_error | mask_J_error | mask_Ks_error

def Mask_photoz(index,data,photoz):
    return (data[index][photoz]>0) & (data[index][photoz]<8)

def Mask_class_star(index,data,class_star):
    return (data[index][class_star]==0)  


def Matching_24Micron():
    print '24micron'
    ###set martirials
    CAT = "COSMOS+mips24_allmatches.fits" #mips24 in COSMOS2015
    ID = "NUMBER"
    SOURCE_COLUMN = "24MICRON"
    MAIN_CAT = "COSMOS2015_Laigle+_v1.1_simple.fits" #COSMOS2015
    OUTPUT_CAT = 'COSMOS2015_24micron.fits'
    
    ###set catalog names
    number = 2
    cat_name = [None]*number
    cat_name[0] = MAIN_CAT
    cat_name[1] = CAT
    
    ###set ID column names
    galaxyid = ID
    
    ###read catalog
    data = [None]*number
    for i in range(number):
        ReadCat(i,data,cat_name[i])
    
    print 'CAT length ',len(data[1])
    
    time2 = time.time()
    print 'current time :', (time2-time1)/60.0 , 'min'
    
    ###matching
    num = 0
    SOURCE = [0]*len( data[0].filled() ) #create an array with length of COSMOS2015
    print 'loop start'
    for i in range(len(data[1])): #go through every id in CAT
        data1_id = data[1][galaxyid][i] - 1
        if SOURCE[ data1_id ]==1:
            print "!!!",data1_id,"!!!already tagged"
        else:
            SOURCE[ data1_id ] = 1
            num +=1
    print "source number ",num
    
    ###output catalog
    print 'start writing table file'
    output_cat = OutputCat(SOURCE_COLUMN, SOURCE)  
    #print output_cat
    output_cat.write(OUTPUT_CAT)
    return
   
def Matching_3GHz():
    print '3GHz'
    ###set martirials
    path = "/Users/yuhsuan/Documents/research/05WH/data/COSMOS/05_3GHz/"
    #CAT = path+"vla3_cosmos_sources_160321_public5sig.fits.txt"
    CAT = path+"VLA_3GHz_counterpart_array_20170210_paper_smolcic_et_al.fits.txt" #3GHz in COSMOS2015
    ID = "ID_CPT"
    SOURCE_COLUMN = "3GHZ"
    MAIN_CAT = "COSMOS2015_Laigle+_v1.1_simple.fits" #COSMOS2015
    OUTPUT_CAT = 'COSMOS2015_3GHz.fits'
    LABEL_COLUMN = ["Xray_AGN","MIR_AGN","SED_AGN","Quiescent_MLAGN","SFG","Clean_SFG",\
                 "HLAGN","MLAGN","Radio_excess"]
    
    ###set catalog names
    number = 2
    cat_name = [None]*number
    cat_name[0] = MAIN_CAT
    cat_name[1] = CAT
    
    ###set ID column names
    galaxyid = ID
    
    ###read catalog
    data = [None]*number
    for i in range(number):
        ReadCat(i,data,cat_name[i])
    
    print 'CAT length ',len(data[1])
    
    ###mask cat label
    mask = (data[1]['CAT_CPT']=='COSMOS2015     ')
    data1_tmp = data[1][mask]
    print 'label masked CAT length ',len(data1_tmp)
    
    time2 = time.time()
    print 'current time :', (time2-time1)/60.0 , 'min'
    
    ###matching
    num = 0
    SOURCE = [0]*len( data[0].filled() ) #create an array with length of COSMOS2015
    print 'loop start'
    for i in range(len(data1_tmp)): #go through every id in CAT
        data1_id = data1_tmp[galaxyid][i] - 1
        if SOURCE[ data1_id ]==1:
            print "!!!",data1_id,"!!!already tagged"
        else:
            SOURCE[ data1_id ] = 1
            num +=1
    print "source number ",num
    
    
    output_cat = Table()
    output_cat[SOURCE_COLUMN] = SOURCE
    
    ###matching label
    for i in range(len(LABEL_COLUMN)):
        LABEL_SOURCE = [-99]*len( data[0].filled() )
        for j in range(len(data1_tmp)): #go through every id in CAT
            data1_id = data1_tmp[galaxyid][j] - 1
            if data1_tmp[LABEL_COLUMN[i]][j].strip()=="true":
                if LABEL_SOURCE[ data1_id ]==1:
                    print "!!!",data1_id,"!!!already tagged"
                else:
                    LABEL_SOURCE[ data1_id ] = 1
                    
            else:
                if LABEL_SOURCE[ data1_id ]==0:
                    print "!!!",data1_id,"!!!already tagged"
                else:
                    LABEL_SOURCE[ data1_id ] = 0
        output_cat[LABEL_COLUMN[i]] = LABEL_SOURCE
    
    ###output catalog
    print 'start writing table file'
    #print output_cat
    output_cat.write(OUTPUT_CAT)
    return

def Matching_ALMA():
    print 'ALMA catalog'
    ###set martirials
    CAT = "06_ALMA/AS2COSMOScatalog_2019-06-01_merged_wJMSphotom_wZspec.fits"
    ID = "ID"
    RA = "RA_final"
    DEC = "Dec_final"
    SOURCE_COLUMN = "ALMA_10"
    MAIN_CAT = "COSMOS2015_Laigle+_v1.1_simple.fits" #COSMOS2015
    OUTPUT_CAT = 'COSMOS2015_ALMA10.fits'
    searching_radius = 1.0
    
    ###set catalog names
    number = 2
    cat_name = [None]*number
    cat_name[0] = MAIN_CAT
    cat_name[1] = CAT
    
    ###set ID column names
    galaxyid = [None]*number
    galaxyid[0] = "NUMBER"
    galaxyid[1] = ID
    
    ##set RA DEC column names
    ra = [None]*number
    dec = [None]*number
    ra[0] = "ALPHA_J2000"
    dec[0] = "DELTA_J2000"
    ra[1] = RA
    dec[1] = DEC
    
    ###read catalog
    data = [None]*number
    for i in range(2):
        ReadCat(i,data,cat_name[i])
    
    time2 = time.time()
    print 'current time :', (time2-time1)/60.0 , 'min'
    
    #matching
    SOURCE = [0]*len( data[0].filled() ) #create an array with length of COSMOS2015
    c0 = SkyCoord(ra=data[0][ra[0]], dec=data[0][dec[0]])
    c1 = SkyCoord(ra=data[1][ra[1]]*u.degree, dec=data[1][dec[1]]*u.degree)
    
    num = 0
    print 'start loop'
    print 'CAT length ',len(data[1].filled())
    #unmatched_num = 0
    for i in range( len(data[1].filled()) ):
        sep = c0.separation(c1[i])
        data0_incircle = data[0][ sep<=searching_radius*u.arcsec ]
        
        for j in range( len(data0_incircle.filled()) ):
            data0_id = data0_incircle[galaxyid[0]][j] -1
            if (SOURCE[ data0_id ] != 0):
                print "!!!",data0_id,"!!!already tagged"
            else:
                SOURCE[ data0_id ] = 1
                num += 1
                #unmatched_num += 1
                #break
    #print unmatched_num
    print "source number ",num
    
    time2 = time.time()
    print 'current time :', (time2-time1)/60.0 , 'min'
    
    ###output
    print 'start writing table file'
    output_cat = OutputCat(SOURCE_COLUMN, SOURCE)  
    #print output_cat
    output_cat.write(OUTPUT_CAT)
    return

def Matching_A3COSMOS():
    print 'A3COSMOS catalog'
    ###set martirials
    path = "/Users/yuhsuan/Documents/research/05WH/data/COSMOS/"
    CAT = path+"07_ALMA/apjsab42da_table4/A-COSMOS_blind.fits"
    
    RA = "RA"
    DEC = "DEC"
    SOURCE_COLUMN = "A3COSMOS"
    MAIN_CAT = "COSMOS2015_Laigle+_v1.1_simple.fits" #COSMOS2015
    OUTPUT_CAT = 'COSMOS2015_A3COSMOS10.fits'
    searching_radius = 1.0
    
    ###set catalog names
    number = 2
    cat_name = [None]*number
    cat_name[0] = MAIN_CAT
    cat_name[1] = CAT
    
    ###set ID column names
    galaxyid = [None]*number
    galaxyid[0] = "NUMBER"
    
    ##set RA DEC column names
    ra = [None]*number
    dec = [None]*number
    ra[0] = "ALPHA_J2000"
    dec[0] = "DELTA_J2000"
    ra[1] = RA
    dec[1] = DEC
    
    ###read catalog
    data = [None]*number
    for i in range(2):
        ReadCat(i,data,cat_name[i])
    
    time2 = time.time()
    print 'current time :', (time2-time1)/60.0 , 'min'
    
    ###matching
    SOURCE = [0]*len( data[0].filled() ) #create an array with length of COSMOS2015
    c0 = SkyCoord(ra=data[0][ra[0]], dec=data[0][dec[0]])
    c1 = SkyCoord(ra=data[1][ra[1]], dec=data[1][dec[1]])
    
    num = 0
    print 'start loop'
    print 'CAT length ',len(data[1].filled())
    #unmatched_num = 0
    multimatched_num = 0
    for i in range( len(data[1].filled()) ):
        sep = c0.separation(c1[i])
        data0_incircle = data[0][ sep<=searching_radius*u.arcsec ]
        
        if (len(data0_incircle.filled())>1): print i, len(data0_incircle.filled())
        for j in range( len(data0_incircle.filled()) ):
            data0_id = data0_incircle[galaxyid[0]][j] -1
            if (SOURCE[ data0_id ] != 0):
                #print "!!!",data0_id,"!!!already tagged"
                multimatched_num +=1
            else:
                SOURCE[ data0_id ] = 1
                num += 1
                #unmatched_num += 1
                #break
    #print unmatched_num
    print "multimatched_num ",multimatched_num
    print "source number ",num
    
    time2 = time.time()
    print 'current time :', (time2-time1)/60.0 , 'min'
    
    ###output
    print 'start writing table file'
    output_cat = OutputCat(SOURCE_COLUMN, SOURCE)  
    #print output_cat
    output_cat.write(OUTPUT_CAT)
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


def Matching_850wide():
    print '850 wide'
    ###set martirials
    CAT = "04_COSMOS450_850/S2COSMOS/catalog/S2COSMOS_sourcecat850_Simpson18.fits"
    RA = "RA_deg"
    DEC = "DEC_deg"
    SOURCE_COLUMN = "850WIDE"
    MAIN_CAT = "COSMOS2015_Laigle+_v1.1_simple.fits"
    
    CAT_24micron_all = "02_mips24/mips24_whwang.fits"
    CAT_24micron_cpt = "COSMOS2015_24micron.fits"
    CAT_3GHz_all = "05_3GHz/vla3_cosmos_sources_160321_public5sig.fits.txt"
    CAT_3GHz_cpt = "COSMOS2015_3GHz.fits"
    CAT_ALMA_all = "06_ALMA/AS2COSMOScatalog_2019-06-01_merged_wJMSphotom_wZspec.fits"
    CAT_ALMA_cpt = "COSMOS2015_ALMA10.fits"
    
    OUTPUT_CAT = 'COSMOS2015_850wide.fits'
    searching_radius_850 = 7.0
    
    
    ###set catalog names
    number = 8
    cat_name = [None]*number
    cat_name[0] = MAIN_CAT
    cat_name[1] = CAT
    cat_name[2] = CAT_24micron_all
    cat_name[3] = CAT_3GHz_all
    cat_name[4] = CAT_ALMA_all
    cat_name[5] = CAT_24micron_cpt
    cat_name[6] = CAT_3GHz_cpt
    cat_name[7] = CAT_ALMA_cpt
    
    ###set ID column names
    galaxyid = [None]*number
    galaxyid[0] = "NUMBER"
    
    ##set RA DEC column names
    ra = [None]*number
    dec = [None]*number
    ra[0] = "ALPHA_J2000"
    dec[0] = "DELTA_J2000"
    ra[1] = RA
    dec[1] = DEC
    ra[2] = "RA_mips24"
    dec[2] = "DEC_mips24"
    ra[3] = "ra"
    dec[3] = "dec"
    ra[4] = "RA_final"
    dec[4] = "Dec_final"
    
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
    
    ###set columns for mutiple detection
    column = [None]*number
    column[5] = "24MICRON"
    column[6] = "3GHZ"
    column[7] = "ALMA_10"
    
    ###read catalog
    data = [None]*number
    for i in range(number):
        ReadCat(i,data,cat_name[i])
    
    time2 = time.time()
    print 'current time :', (time2-time1)/60.0 , 'min'
    
    ###mask cat
    mask_data0 = Mask_M(0,data,color1,color2,color3) \
    & Mask_error(0,data,V_error,ip_error,J_error,Ks_error) \
    & Mask_photoz(0,data,photoz) & Mask_class_star(0,data,class_star)
    data0_tmp = data[0][mask_data0]
    print 'masked MAIN_CAT length ',len(data0_tmp.filled())
    
    ###calculate coordinate
    c0 = SkyCoord(ra=data0_tmp[ra[0]], dec=data0_tmp[dec[0]]) #MAIN_CAT
    c1 = SkyCoord(ra=data[1][ra[1]]*u.degree, dec=data[1][dec[1]]*u.degree) #CAT
    c2 = SkyCoord(ra=data[2][ra[2]]*u.degree, dec=data[2][dec[2]]*u.degree) #24
    c3 = SkyCoord(ra=data[3][ra[3]], dec=data[3][dec[3]]) #3
    c4 = SkyCoord(ra=data[4][ra[4]]*u.degree, dec=data[4][dec[4]]*u.degree) #ALMA
    
    ###matching
    
    SOURCE = [0]*len( data[0].filled() ) #create an array with length of COSMOS2015
    
    typeA_num = 0 #24 3 both detected
    label1_num = 0
    
    typeB_num = 0 #24 or 3 detected
    label2_num = 0 #24
    label3_num = 0 #3
    
    typeC_num = 0 #no opt counterpart
    
    typeD_num = 0 #all selected
    label4_num = 0
    
    typeE_num = 0 #ALMA detected, no opt couterpart
    ALMA_no_opt_num = 0
    
    typeF_num = 0 #ALMA detected
    label5_num = 0
    
    print 'loop start'
    print 'CAT length ',len(data[1].filled())
    for i in range( len(data[1].filled()) ): #go through every id in CAT
        sep0 = c0.separation(c1[i])
        data0_in_circle = data0_tmp[ sep0<=searching_radius_850*u.arcsec ]
        sep4 = c4.separation(c1[i])
        data4_in_circle = data[4][ sep4<=searching_radius_850*u.arcsec ]
     
        if len(data4_in_circle) >0:
            #ALMA follow-ups
            flag = 0
            opt_num = 0
            for j in range(len(data0_in_circle)):
                
                
                data0_id = data0_in_circle[galaxyid[0]][j] - 1
                if data[7][column[7]][data0_id]==1: #ALMA detected
                    if SOURCE[data0_id]!=0:
                        print "!!!",data0_id,"!!!already tagged",SOURCE[data0_id],", would like to tag 5"
                    else:
                        SOURCE[data0_id] = 5
                        opt_num += 1
                        flag = 1
            if flag == 1:
                typeF_num += 1
                label5_num += opt_num
                ALMA_no_opt_num += len(data4_in_circle) - opt_num
            else:
                typeE_num += 1
                ALMA_no_opt_num += len(data4_in_circle)
            
            
        else:
            #not ALMA follow-ups
            
            sep2 = c2.separation(c1[i])
            data2_in_circle = data[2][ sep2<=searching_radius_850*u.arcsec ]
            sep3 = c3.separation(c1[i])
            data3_in_circle = data[3][ sep3<=searching_radius_850*u.arcsec ]
            
            if len(data2_in_circle)+len(data3_in_circle)>0:
                #sources with clues, no need to select all
                
                flag = 0
                for j in range(len(data0_in_circle)):
                    
                    
                    data0_id = data0_in_circle[galaxyid[0]][j] - 1
                    if data[5][column[5]][data0_id]==1: #24 micron detected
                        if data[6][column[6]][data0_id]==1: #3 GHz detected
                            if SOURCE[data0_id]!=0:
                                print "!!!",data0_id,"!!!already tagged",SOURCE[data0_id],", would like to tag 1"
                            else:
                                SOURCE[data0_id] = 1
                                label1_num += 1
                                flag = 1
                                
                                
                if flag == 1:
                    typeA_num += 1
    
                else:
                    for j in range(len(data0_in_circle)):
                        
                        
                        data0_id = data0_in_circle[galaxyid[0]][j] - 1
                        if data[5][column[5]][data0_id]==1:  #24
                            if SOURCE[data0_id]!=0:
                                print "!!!",data0_id,"!!!already tagged",SOURCE[data0_id],", would like to tag 2"
                            else:
                                SOURCE[data0_id] = 2
                                label2_num += 1
                                flag = 1
                        if data[6][column[6]][data0_id]==1:  #3
                            if SOURCE[data0_id]!=0:
                                print "!!!",data0_id,"!!!already tagged",SOURCE[data0_id],", would like to tag 3"
                            else:
                                SOURCE[data0_id] = 2
                                label3_num += 1
                                flag = 1
                            
                            
                    if flag == 1:
                        typeB_num += 1
                    else:
                        #clues given, but no opt sample found
                        typeC_num += 1
                    
            else:
                #sources with no clues, select all
                for j in range(len(data0_in_circle)):
                    
                    data0_id = data0_in_circle[galaxyid[0]][j] - 1
                    if SOURCE[data0_id]!=0:
                        print "!!!",data0_id,"!!!already tagged",SOURCE[data0_id],", would like to tag 4"
                    else:
                        SOURCE[data0_id] = 4
                        label4_num += 1
                typeD_num += 1
    
    print 'typeA_num ', typeA_num #24 3 both detected
    print 'label1_num ', label1_num
    
    print 'typeB_num ', typeB_num #24 or 3 detected
    print 'label2_num ', label2_num #24
    print 'label3_num ', label3_num #3
    
    print 'typeC_num ', typeC_num #no opt counterpart
    
    print 'typeD_num ', typeD_num #all selected
    print 'label4_num ', label4_num
    
    print 'typeE_num ', typeE_num #ALMA detected, no opt couterpart
    print 'ALMA_no_opt_num ', ALMA_no_opt_num
    
    print 'typeF_num ', typeF_num #ALMA detected
    print 'label5_num ', label5_num
    
    ###output
    print 'start writing table file'
    output_cat = OutputCat(SOURCE_COLUMN, SOURCE)  
    print output_cat
    #output_cat.write(OUTPUT_CAT)
    return



def Matching_450narrow():
    print '450 narrow'
    ###set martirials
    CAT = "04_COSMOS450_850/STUDIES/sources_450.fits"
    RA = "RA_450"
    DEC = "DEC_450"
    SOURCE_COLUMN = "450NARROW"
    MAIN_CAT = "COSMOS2015_Laigle+_v1.1_simple.fits"
    
    CAT_24micron_all = "02_mips24/mips24_whwang.fits"
    CAT_24micron_cpt = "COSMOS2015_24micron.fits"
    CAT_3GHz_all = "05_3GHz/vla3_cosmos_sources_160321_public5sig.fits.txt"
    CAT_3GHz_cpt = "COSMOS2015_3GHz.fits"
    CAT_ALMA_all = "06_ALMA/AS2COSMOScatalog_2019-06-01_merged_wJMSphotom_wZspec.fits"
    CAT_ALMA_cpt = "COSMOS2015_ALMA10.fits"
    
    OUTPUT_CAT = 'COSMOS2015_450narrow.fits'
    searching_radius_450 = 4.0
    
    
    ###set catalog names
    number = 8
    cat_name = [None]*number
    cat_name[0] = MAIN_CAT
    cat_name[1] = CAT
    cat_name[2] = CAT_24micron_all
    cat_name[3] = CAT_3GHz_all
    cat_name[4] = CAT_ALMA_all
    cat_name[5] = CAT_24micron_cpt
    cat_name[6] = CAT_3GHz_cpt
    cat_name[7] = CAT_ALMA_cpt
    
    ###set ID column names
    galaxyid = [None]*number
    galaxyid[0] = "NUMBER"
    
    ##set RA DEC column names
    ra = [None]*number
    dec = [None]*number
    ra[0] = "ALPHA_J2000"
    dec[0] = "DELTA_J2000"
    ra[1] = RA
    dec[1] = DEC
    ra[2] = "RA_mips24"
    dec[2] = "DEC_mips24"
    ra[3] = "ra"
    dec[3] = "dec"
    ra[4] = "RA_final"
    dec[4] = "Dec_final"
    
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
    
    ###set columns for mutiple detection
    column = [None]*number
    column[5] = "24MICRON"
    column[6] = "3GHZ"
    column[7] = "ALMA_10"
    
    ###read catalog
    data = [None]*number
    for i in range(number):
        ReadCat(i,data,cat_name[i])
    
    time2 = time.time()
    print 'current time :', (time2-time1)/60.0 , 'min'
    
    ###mask cat
    mask_data0 = Mask_M(0,data,color1,color2,color3) \
    & Mask_error(0,data,V_error,ip_error,J_error,Ks_error) \
    & Mask_photoz(0,data,photoz) & Mask_class_star(0,data,class_star)
    data0_tmp = data[0][mask_data0]
    print 'masked MAIN_CAT length ',len(data0_tmp.filled())
    
    ###calculate coordinate
    c0 = SkyCoord(ra=data0_tmp[ra[0]], dec=data0_tmp[dec[0]]) #MAIN_CAT
    c1 = SkyCoord(ra=data[1][ra[1]]*u.degree, dec=data[1][dec[1]]*u.degree) #CAT
    c2 = SkyCoord(ra=data[2][ra[2]]*u.degree, dec=data[2][dec[2]]*u.degree) #24
    c3 = SkyCoord(ra=data[3][ra[3]], dec=data[3][dec[3]]) #3
    c4 = SkyCoord(ra=data[4][ra[4]]*u.degree, dec=data[4][dec[4]]*u.degree) #ALMA
    
    ###matching
    
    SOURCE = [0]*len( data[0].filled() ) #create an array with length of COSMOS2015
    
    typeA_num = 0 #24 3 both detected
    label1_num = 0
    
    typeB_num = 0 #24 or 3 detected
    label2_num = 0 #24
    label3_num = 0 #3
    
    typeC_num = 0 #no opt counterpart
    
    typeD_num = 0 #all selected
    label4_num = 0
    
    typeE_num = 0 #ALMA detected, no opt couterpart
    ALMA_no_opt_num = 0
    
    typeF_num = 0 #ALMA detected
    label5_num = 0
    
    print 'loop start'
    print 'CAT length ',len(data[1].filled())
    for i in range( len(data[1].filled()) ): #go through every id in CAT
        sep0 = c0.separation(c1[i])
        data0_in_circle = data0_tmp[ sep0<=searching_radius_450*u.arcsec ]
        sep4 = c4.separation(c1[i])
        data4_in_circle = data[4][ sep4<=searching_radius_450*u.arcsec ]
     
        if len(data4_in_circle) >0:
            #ALMA follow-ups
            flag = 0
            opt_num = 0
            for j in range(len(data0_in_circle)):
                
                
                data0_id = data0_in_circle[galaxyid[0]][j] - 1
                if data[7][column[7]][data0_id]==1: #ALMA detected
                    if SOURCE[data0_id]!=0:
                        print "!!!",data0_id,"!!!already tagged",SOURCE[data0_id],", would like to tag 5"
                    else:
                        SOURCE[data0_id] = 5
                        opt_num += 1
                        flag = 1
            if flag == 1:
                typeF_num += 1
                label5_num += opt_num
                ALMA_no_opt_num += len(data4_in_circle) - opt_num
            else:
                typeE_num += 1
                ALMA_no_opt_num += len(data4_in_circle)
            
            
        else:
            #not ALMA follow-ups
            
            sep2 = c2.separation(c1[i])
            data2_in_circle = data[2][ sep2<=searching_radius_450*u.arcsec ]
            sep3 = c3.separation(c1[i])
            data3_in_circle = data[3][ sep3<=searching_radius_450*u.arcsec ]
            
            if len(data2_in_circle)+len(data3_in_circle)>0:
                #sources with clues, no need to select all
                
                flag = 0
                for j in range(len(data0_in_circle)):
                    
                    
                    data0_id = data0_in_circle[galaxyid[0]][j] - 1
                    if data[5][column[5]][data0_id]==1: #24 micron detected
                        if data[6][column[6]][data0_id]==1: #3 GHz detected
                            if SOURCE[data0_id]!=0:
                                print "!!!",data0_id,"!!!already tagged",SOURCE[data0_id],", would like to tag 1"
                            else:
                                SOURCE[data0_id] = 1
                                label1_num += 1
                                flag = 1
                                
                                
                if flag == 1:
                    typeA_num += 1
    
                else:
                    for j in range(len(data0_in_circle)):
                        
                        
                        data0_id = data0_in_circle[galaxyid[0]][j] - 1
                        if data[5][column[5]][data0_id]==1:  #24
                            if SOURCE[data0_id]!=0:
                                print "!!!",data0_id,"!!!already tagged",SOURCE[data0_id],", would like to tag 2"
                            else:
                                SOURCE[data0_id] = 2
                                label2_num += 1
                                flag = 1
                        if data[6][column[6]][data0_id]==1:  #3
                            if SOURCE[data0_id]!=0:
                                print "!!!",data0_id,"!!!already tagged",SOURCE[data0_id],", would like to tag 3"
                            else:
                                SOURCE[data0_id] = 2
                                label3_num += 1
                                flag = 1
                            
                            
                    if flag == 1:
                        typeB_num += 1
                    else:
                        #clues given, but no opt sample found
                        typeC_num += 1
                    
            else:
                #sources with no clues, select all
                for j in range(len(data0_in_circle)):
                    
                    data0_id = data0_in_circle[galaxyid[0]][j] - 1
                    if SOURCE[data0_id]!=0:
                        print "!!!",data0_id,"!!!already tagged",SOURCE[data0_id],", would like to tag 4"
                    else:
                        SOURCE[data0_id] = 4
                        label4_num += 1
                typeD_num += 1
    
    print 'typeA_num ', typeA_num #24 3 both detected
    print 'label1_num ', label1_num
    
    print 'typeB_num ', typeB_num #24 or 3 detected
    print 'label2_num ', label2_num #24
    print 'label3_num ', label3_num #3
    
    print 'typeC_num ', typeC_num #no opt counterpart
    
    print 'typeD_num ', typeD_num #all selected
    print 'label4_num ', label4_num
    
    print 'typeE_num ', typeE_num #ALMA detected, no opt couterpart
    print 'ALMA_no_opt_num ', ALMA_no_opt_num
    
    print 'typeF_num ', typeF_num #ALMA detected
    print 'label5_num ', label5_num
    
    ###output
    print 'start writing table file'
    output_cat = OutputCat(SOURCE_COLUMN, SOURCE)  
    #print output_cat
    output_cat.write(OUTPUT_CAT)
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
#Matching_ALMA()

# ===== ALMA =====
Matching_A3COSMOS()

# ===== lensing =====
#id = 659416
#Labeling_lensing()


# ===== 850 wide =====
#COSMOS2015_24micron,3GHz,ALMA needed
#Matching_850wide()

# ===== 450 narrow =====
#COSMOS2015_24micron,3GHz,ALMA needed
#Matching_450narrow()


time2 = time.time()
print 'done! time :', (time2-time1)/60.0 , 'min'