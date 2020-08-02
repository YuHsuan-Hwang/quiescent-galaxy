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
MergingMainCat()

time2 = time.time()
print 'done! time =', (time2-time1)/60.0 , 'min'