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
from astropy.table import Table

# =============================================================================
# functions
# =============================================================================
def ReadCatalog(catalog):
    return Table.read(catalog, hdu=1)

def AddColumns():
    for i in range(len(columns)):
        t[columns[i]] = data[columns[i]]
        
# =============================================================================
# main code
# =============================================================================
time1 = time.time()

###set catalog name
inputcatalog = "COSMOS2015_Laigle+_v1.1_5+2band_2agn_9cat_ALMA4"
#"COSMOS2015_Laigle+_v1.1_850sources"
catalog = inputcatalog+".fits"

###set columns
columns = [ "ALPHA_J2000", "DELTA_J2000", "NUMBER", "MNUV", "MR", "MJ", "PHOTOZ",\
           "V_MAGERR_APER3", "ip_MAGERR_APER3", "J_MAGERR_APER3", "Ks_MAGERR_APER3",\
           "MASS_MED", "TYPE", "CLASS", "Ks_MAG_APER3",\
           "850SOURCE","850NARROW","450NARROW","24MICRON","3GHZ", \
           "IR_AGN_ALL","IR_AGN_MIR",\
           "Xray_AGN","MIR_AGN","SED_AGN","Quiescent_MLAGN","SFG","Clean_SFG",\
           "HLAGN","MLAGN","Radio_excess","450NARROW_SIMPLE","850WIDE_SIMPLE",\
           "MASS_BEST","ALMA_05","ALMA_10","850WIDE_ALMA","850WIDE_ALMA_10",\
           "FLAG_HJMCC","FLAG_PETER","FLAG_COSMOS","FLAG_DEEP","FLAG_SHALLOW"]

###read catalog
data = ReadCatalog(catalog)

###add columns
t = Table()
AddColumns()

###output fits file
t.write(inputcatalog+'_simple.fits')

time2 = time.time()
print 'done! time =', time2-time1 , 'sec'