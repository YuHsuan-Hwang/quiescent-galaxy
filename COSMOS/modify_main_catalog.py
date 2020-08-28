#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 13:34:58 2019

@author: yuhsuan
"""
# =============================================================================
# packages
# =============================================================================
import time
from astropy.table import Table, hstack

# =============================================================================
# functions
# =============================================================================
def ReadCatalog(catalog):
    return Table.read(catalog, hdu=1)

def AddColumns(t,data,columns):
    for i in range(len(columns)):
        t[columns[i]] = data[columns[i]]

def SimplifyMainCat():
    ###set catalog name
    path = "/Users/yuhsuan/Documents/research/05WH/data/COSMOS/01_COSMOS2015catalog/COSMOS2015/"
    inputcatalog = "COSMOS2015_Laigle+_v1.1"
    catalog = path+inputcatalog+".fits"
    
    ###set columns
    columns = [ "ALPHA_J2000", "DELTA_J2000", "NUMBER", "MNUV", "MR", "MJ", "MU", "MV", "PHOTOZ","ZQ",
                "V_MAGERR_APER3", "ip_MAGERR_APER3", "J_MAGERR_APER3", "Ks_MAGERR_APER3",
                "V_MAG_APER3", "ip_MAG_APER3", "J_MAG_APER3", "Ks_MAG_APER3",
                "MASS_MED","MASS_MED_MIN68","MASS_MED_MAX68", "TYPE", "CLASS",
                "MASS_BEST","FLAG_HJMCC","FLAG_PETER","FLAG_COSMOS","FLAG_DEEP","FLAG_SHALLOW"  ]
    
    ###read catalog
    data = ReadCatalog(catalog)
    
    ###add columns
    t = Table()
    AddColumns(t,data,columns)
    
    ###output fits file
    t.write(inputcatalog+'_simple.fits')
    return

def SimplifyAGNCat():
    ###set catalog name
    inputcatalog = "cos_agn"
    catalog = inputcatalog+".fits"
    
    ###set columns
    columns = [ "agn_c17b","agn_xxx" ]
    
    ###read catalog
    data = ReadCatalog(catalog)
    
    ###add columns
    t = Table()
    AddColumns(t,data,columns)
    
    ###output fits file
    t.write(inputcatalog+'_simple.fits')
    return

def SimplifyAGNcaCat():
    ###set catalog name
    inputcatalog = "ca_all_iragns"
    catalog = inputcatalog+".fits"
    
    ###set columns
    columns = [ "IRAGN" ]
    
    ###read catalog
    data = ReadCatalog(catalog)
    
    ###add columns
    t = Table()
    AddColumns(t,data,columns)
    
    ###output fits file
    t.write("COSMOS2015_"+inputcatalog+'_simple.fits')
    return

def MergingMainCatTmp():
    
    ###set catalog name
    inputcatalog = "COSMOS2015_"
    cats = ["Laigle+_v1.1_simple_z",
            "24micron","3GHz","AS2COSMOS","A3COSMOS","lensing"]
    
    ###read tables
    tables = []
    for i in range(len(cats)):
        tables.append(Table.read(inputcatalog+cats[i]+".fits",hdu=1))
    
    ###merage tables
    output_cat = hstack([tables[0],tables[1]], join_type='exact')
    for i in range(len(cats)-2):
        output_cat = hstack([output_cat,tables[i+2]], join_type='exact')
    
    ###output
    output_cat.write("COSMOS2015_merged_tmp.fits")
    return


def MergingMainCat():
    
    ###set catalog name
    inputcatalog = "COSMOS2015_"
    cats = ["Laigle+_v1.1_simple_z",
            "24micron","3GHz","AS2COSMOS","A3COSMOS","lensing","850wide","450narrow",
            "850wide_oneband","450narrow_oneband",
            "cos_agn_simple","ca_all_iragns_simple"]
    
    ###read tables
    tables = []
    for i in range(len(cats)):
        tables.append(Table.read(inputcatalog+cats[i]+".fits",hdu=1))
    
    ###merage tables
    output_cat = hstack([tables[0],tables[1]], join_type='exact')
    for i in range(len(cats)-2):
        output_cat = hstack([output_cat,tables[i+2]], join_type='exact')
    
    ###output
    output_cat.write("COSMOS2015_merged.fits")
    return

def MergingColumnRedshift():
    
    ###set catalog name
    inputcatalog = "COSMOS2015_Laigle+_v1.1_simple"
    catalog = inputcatalog+".fits"
    
    SOURCE_COLUMN = "REDSHIFT"
    
    ###read catalog
    data = ReadCatalog(catalog)
    
    SOURCE = [-999.99]*len( data )
    
    print len(data)
    
    print "loop start"
    for i in range( len(data) ):
        if data["PHOTOZ"][i]>9.9:
            SOURCE[i] = data["ZQ"][i]
        else:
            SOURCE[i] = data["PHOTOZ"][i]
        if i==10000:
            time2 = time.time()
            print 'current time :', (time2-time1)/60.0 , 'min'
        
    data[SOURCE_COLUMN] = SOURCE
    
    ###output
    data.write("COSMOS2015_Laigle+_v1.1_simple_z.fits")
    #print SOURCE
    return



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
        return (data[index][Ks_error]<up) & (data[index][Ks_error]>0)
    
    if mode==3:
        return (data[index][Ks_mag]<24)
    

def Mask_photoz( index ):
    return (data[index][photoz]>0)

def Mask_class_star( index ):
    return ( data[index][class_star]==0 ) | ( data[index][class_star]==2 )

def MaskAll( index ):
    return Mask_M( index ) & Mask_error( 1, 0.1, index ) & Mask_photoz( index ) & Mask_class_star( index )

def MaskMainCat():
    
    mask = MaskAll(0)
    data_masked = data[0][mask]
    
    ###output
    data_masked.write("COSMOS2015_Laigle+_v1.1_simple_z_masked.fits")
    #print SOURCE
    return    
    
    
# =============================================================================
# main code
# =============================================================================

time1 = time.time()

# ===== Simplify =====
#SimplifyMainCat()

# ===== Simplify AGN catalog from yu-yen =====
#cos_agn.fits
#SimplifyAGNCat()

#ca_all_iragns.fits
#SimplifyAGNcaCat()

# ===== Merging two redshift in COSMOS2015 =====
#MergingColumnRedshift()

# ===== Merging =====
#MergingMainCatTmp()
#MergingMainCat()


# ===== Masking =====
###set colors
colorname1 = "NUV"
colorname2 = "r"
colorname3 = "J"

###set columns in the main catalog

color1 = "MNUV"
color2 = "MR"
color3 = "MJ"
color  = [ color1, color2, color3 ]

photoz     = "REDSHIFT"
V_error    = "V_MAGERR_APER3"
ip_error   = "ip_MAGERR_APER3"
J_error    = "J_MAGERR_APER3"
Ks_error   = "Ks_MAGERR_APER3"
V_mag      = "V_MAG_APER3"
ip_mag     = "ip_MAG_APER3"
J_mag      = "J_MAG_APER3"
Ks_mag     = "Ks_MAG_APER3"
#mass       = "MASS_MED"
class_star = "TYPE"
class_SFG  = "CLASS"
ra         = "ALPHA_J2000"
dec        = "DELTA_J2000"
mass       = "MASS_BEST"

set_xlable = '$M_{'+colorname2+'}-M_{'+colorname3+'}$'
set_ylable = '$M_{'+colorname1+'}-M_{'+colorname2+'}$' 

###read catalog
catalog    = [None]*1
catalog[0] = "COSMOS2015_Laigle+_v1.1_simple_z.fits"
data       = [None]*1
x,y,mask,x_masked,y_masked = [],[],[],[],[]

data[0] = ReadCatalog( catalog[0] )
MaskMainCat()

time2 = time.time()
print 'done! time =', (time2-time1)/60.0 , 'min'