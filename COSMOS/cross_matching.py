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
def ReadCat( index, data,cat ):
    data[index] = Table.read(cat, hdu=1)
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
    if SOURCE[source_id]!=0:
        print "!!!",source_id,"already tagged",SOURCE[source_id],", would like to tag ",tag
    else:
        SOURCE[source_id] = tag
        label_num[tag-1] += 1
        
    return


def MatchingSubmm():
    
    MAIN_CAT         = "COSMOS2015_merged_tmp.fits"
    CAT_24micron_all = path + "02_mips24/mips24_whwang.fits"
    CAT_3GHz_all     = path + "05_3GHz/vla3_cosmos_sources_160321_public5sig.fits.txt"
    CAT_ALMA_all     = path + "07_ALMA/AS2COSMOScatalog_2019-06-01_merged_wJMSphotom_wZspec.fits"
    CAT_A3COSMOS_all = path + "07_ALMA/apjsab42da_table4/A-COSMOS_blind.fits"
    
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
    
    global SOURCE, label_num
    
    SOURCE = [0]*len( data[0].filled() ) #create an array with length of COSMOS2015
    
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
    
    print "matching loop start ..."
    print 'CAT length ',len(data[1].filled())
    
    for i in range( len(data[1].filled()) ): #go through every id in CAT
        
        sep4            = c4.separation( c1[i] )
        data4_in_circle = data[4][ sep4<=searching_radius*u.arcsec ]
        sep5            = c5.separation( c1[i] )
        data5_in_circle = data[5][ sep5<=searching_radius*u.arcsec ]
    
        
        if ( len(data5_in_circle)>0 )|( len(data4_in_circle)>0 ): #ALMA detected
            
            ALMA_detect_Num +=1
            
            sep8            = c8.separation( c1[i] )
            data8_in_circle = data0_masked_AS2[ sep8<=searching_radius*u.arcsec ]
            sep9            = c9.separation( c1[i] )
            data9_in_circle = data0_masked_A3[ sep9<=searching_radius*u.arcsec ]
            
            if ( len(data9_in_circle)>0 )|( len(data8_in_circle)>0 ):
                
                label_flag = 0
                for j in range ( len(data9_in_circle) ): # go through every 24 opt sources
                    for k in range( len(data8_in_circle) ):
                        if ( data9_in_circle[galaxyid][j]==data8_in_circle[galaxyid][k] ):
                                
                            TagSource( data9_in_circle, j, 1 )
                            label_flag = 1

                if label_flag==1 :
                    A3andAS2_detect_Num  += 1
                    continue               
            
            if ( len(data9_in_circle)>0 )|( len(data8_in_circle)==0 ): #A3COSMOS detected
            
                A3COSMOS_detect_Num += 1         
                for j in range( len(data9_in_circle) ):
                    TagSource( data9_in_circle, j, 2 )
                            
            elif ( len(data9_in_circle)==0 )|( len(data8_in_circle)>0 ): #AS2COSMOS detected
                
                AS2COSMOS_detect_Num += 1           
                for j in range( len(data8_in_circle) ):
                    TagSource( data8_in_circle, j, 3 )
                            
            else: #no opt counterpart
                   A3orAS2_detect_no_opt_Num += 1
            
        else: #no ALMA observation
            
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
            
                if ( len(data6_in_circle)>0 )|( len(data7_in_circle)>0 ): # with 24 or 3 opt sources
                    
                    label_flag = 0
                    
                    for j in range ( len(data6_in_circle) ): # go through every 24 opt sources
                        for k in range( len(data7_in_circle) ):
                            if ( data6_in_circle[galaxyid][j]==data7_in_circle[galaxyid][k] ):
                                
                                TagSource( data6_in_circle, j, 4 )
                                label_flag = 1

                    if label_flag==1 :
                        MIPSandVLA_detect_Num += 1
                        continue
                                    
                if ( len(data6_in_circle)>0 )|( len(data7_in_circle)==0 ): # with 24 opt sources
                    
                    MIPS_detect_Num += 1                    
                    for j in range ( len(data6_in_circle) ): # go through every 24 opt sources                      
                        TagSource( data6_in_circle, j, 5 )
                            
                elif ( len(data6_in_circle)==0 )|( len(data7_in_circle)>0 ): # with 3 opt sources
                    
                    VLA_detect_Num += 1                    
                    for j in range ( len(data7_in_circle) ): # go through every 3 opt sources                       
                        TagSource( data7_in_circle, j, 6 )
                            
                else: # with 24 or 3 detected but no opt sample
                    
                    MIPSorVLA_detect_no_opt_Num += 1
            
            else: # no ALMA, 24, or 3 detection
            
                non_detect_Num += 1
                
                sep0            = c0.separation( c1[i] )
                data0_in_circle = data0_masked[ sep0<=searching_radius*u.arcsec ]
                
                if ( len(data0_in_circle)>0 ):
                    for j in range ( len(data0_in_circle) ): # go through every opt sources
                        TagSource( data0_in_circle, j, 7 )    

                else:                   
                    non_detect_no_opt_Num += 1
                            
    
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
    
    ###output
    print 'writing table file ...'
    output_cat = OutputCat(SOURCE_COLUMN, SOURCE)  
    #output_cat.write(OUTPUT_CAT)
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


# ===== 850 wide =====
#COSMOS2015_24micron,3GHz,ALMA needed
Matching_850wide()

# ===== 450 narrow =====
#COSMOS2015_24micron,3GHz,ALMA needed
Matching_450narrow()


time2 = time.time()
print 'done! time :', (time2-time1)/60.0 , 'min'