#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 23:51:12 2019

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
    x[index] = data[index][color2]-data[index][color3]
    y[index] = data[index][color1]-data[index][color2]
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

def Mask_classQG(index):
    return (data[index][class_SFG]==0)
def Mask_classSFG(index):
    return (data[index][class_SFG]==1)
def Mask_myclassQG(index):
    return (y[index]>3.1) & (y[index]>3.0*x[index]+1.0)
def Mask_myclassSFG(index):
    return (y[index]<=3.1) | (y[index]<=3.0*x[index]+1.0)

#def Mask_mass(index,low,high):
#    return (data[index][mass]<high) & (data[index][mass]>low)
    
def PlotHist_photoz(index):
    NBINS = 50
    plt.hist(data[index][photoz][mask[index]], NBINS, color='b', alpha=0.7)
    plt.title('Histogram photoz')
    plt.ylabel('Number')
    plt.xlabel('z')
    plt.show()
    return

def PlotHist_M(index):
    NBINS = 70
    histcolor = [ 'b', 'g', 'r' ]
    for i in range(3):
        plt.subplot(3, 1, i+1)
        histdata = data[index][color[i]][mask[index]]
        plt.hist(histdata, NBINS, color=histcolor[i], alpha=0.7,label=color[i])
        if (i==0):
            plt.title('Histogram')
        plt.ylabel('Number')
        if (i==2):
            plt.xlabel('magnitude')
        plt.legend(frameon=False,loc=1)
    plt.show()
    return

def DrawLine():
    plt.plot([-1.5,(3.1-1.0)/3.0],[3.1,3.1],'k-',lw=0.7)
    plt.plot([(3.1-1.0)/3.0,2.0],[3.1,3.0*2.0+1.0],'k-',lw=0.7)
    return

def Plot(index, scale, struc, limit, line):
    if (struc==1):
        set_s = 0.5
        set_alpha = 0.1
    else:
        set_s = 1
        set_alpha = 1
    plt.scatter( x_masked[index], y_masked[index], s=set_s, alpha=set_alpha)
    plt.title(colorname1+colorname2+colorname3)
    plt.xlabel(set_xlable)
    plt.ylabel(set_ylable)
    if (limit==1):
        plt.axis([-1.5,2,-1,7])
    if (scale==1):
        plt.axis('scaled')
    if (line==1):
        DrawLine()
    plt.show()
    return


def Construct_mask_zbin(index):
    mask_zbin = []
    mask_zbin.append((data[index][photoz]>0) & (data[index][photoz]<0.5))
    mask_zbin.append((data[index][photoz]>0.5) & (data[index][photoz]<1))
    mask_zbin.append((data[index][photoz]>1) & (data[index][photoz]<1.5))
    mask_zbin.append((data[index][photoz]>1.5) & (data[index][photoz]<2))
    mask_zbin.append((data[index][photoz]>2) & (data[index][photoz]<2.5))
    mask_zbin.append((data[index][photoz]>2.5) & (data[index][photoz]<3))
    mask_zbin.append((data[index][photoz]>3) & (data[index][photoz]<3.5))
    mask_zbin.append((data[index][photoz]>3.5) & (data[index][photoz]<4))
    #mask_zbin.append((data[index][photoz]>3) & (data[index][photoz]<4))
    return mask_zbin

def Construct_zbin_title():
    zbin_title = []
    zbin_title.append('$0<z<0.5$')
    zbin_title.append('$0.5<z<1$')
    zbin_title.append('$1<z<1.5$')
    zbin_title.append('$1.5<z<2$')
    zbin_title.append('$2<z<2.5$')
    zbin_title.append('$2.5<z<3$')
    zbin_title.append('$3<z<3.5$')
    zbin_title.append('$3.5<z<4$')
    return zbin_title

def Plot_zbin(index, scale, struc, limit, line):
    if (struc==1):
        set_s = 0.5
        set_alpha = 0.1
    else:
        set_s = 1
        set_alpha = 1
    mask_zbin = Construct_mask_zbin(index)
    zbin_title = Construct_zbin_title()

    for i in range(8):
        plt.subplot(2, 4, i+1)
        datamask = mask[index] & mask_zbin[i]
        datax = x[index][datamask]
        datay = y[index][datamask]
        plt.scatter( datax, datay, s=set_s, alpha=set_alpha)
        plt.title(zbin_title[i])
        if ((i==4) | (i==5) | (i==6) | (i==7)):
            plt.xlabel(set_xlable)
        if ((i==0) | (i==4)):
            plt.ylabel(set_ylable)
        if (limit==1):
            plt.axis([-1.5,2,-1,7])
        if (scale==1):
            plt.axis('scaled')
        if (line==1):
            DrawLine()
    plt.show()
    return

def PrintNum_zbin(index):
    mask_zbin = Construct_mask_zbin(index)
    #print len( data[index].filled() )
    print len( data[index].filled()[mask[index]] )
    for i in range(8):
        print len( data[index].filled()[mask[index] & mask_zbin[i]] )
    return

def PrintNumdiff_zbin(index1, index2):
    mask_zbin1 = Construct_mask_zbin(index1)
    mask_zbin2 = Construct_mask_zbin(index2)
    print len( data[index1].filled()[mask[index1]] )-len( data[index2].filled()[mask[index2]] )
    for i in range(8):
        len1 = len( data[index1].filled()[mask[index1] & mask_zbin1[i]] )
        len2 = len( data[index2].filled()[mask[index2] & mask_zbin2[i]] )
        print len1-len2
    return

def PrintPercentage_zbin(index1, index2):
    mask_zbin1 = Construct_mask_zbin(index1)
    mask_zbin2 = Construct_mask_zbin(index2)
    for i in range(8):
        len1 = len( data[index1].filled()[mask[index1] & mask_zbin1[i]] )
        if (len1==0.0):
            print 'none'
        else:
            len2 = len( data[index2].filled()[mask[index2] & mask_zbin2[i]] )
            print "%.4f" % ((100*len2)/len1)
    return

def RegionFile(index,filename,color,size,galaxyid):
    f = open('/Users/yuhsuan/Documents/research/05WH/data/'+filename+'.reg','w')
    f.write('global color='+color+' font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    N = len( data[index].filled()[mask[index]] )
    for n in range(N):
        f.write('fk5;circle('+str(data[index][ra][mask[index]][n])+','+str(data[index][dec][mask[index]][n])+','+size+'") # text={'+str(data[index][galaxyid][mask[index]][n])+'}\n')
    f.close()
    
# =============================================================================
# main code
# =============================================================================

time1 = time.time()

###set catalogs
number = 5
catalog = [None]*number
catalog[0] = "COSMOS2015_Laigle+_v1.1_simple.fits"
#"01_COSMOS2015catalog/COSMOS2015/COSMOS2015_Laigle+_v1.1.fits"
catalog[1] = "COSMOS+mips24_allmatches.fits"
catalog[2] = "COSMOS+wide850_allmatches.fits"
catalog[3] = "COSMOS+mips24＋wide850_allmatches.fits"
catalog[4] = "COSMOS+wide850_bestmatchfor850.fits"


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
ra = "ALPHA_J2000"
dec = "DELTA_J2000"

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[i])

###mask data
mask = []
x_masked = []
y_masked = []

for i in range(number):
    mask.append( Mask_M(i) & Mask_photoz(i) & Mask_error(i) & Mask_class_star(i) )
    #& Mask_myclassQG(i)
    #& Mask_mass(i,7.0,9.8)
    #& Mask_myclassQG(i) & Mask_classSFG(i)
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])
    
#(data24['flux_24_2']>80)

###plot
        
#PlotHist_photoz(0)    
#PlotHist_M(0)

set_xlable = '$M_{'+colorname2+'}-M_{'+colorname3+'}$'
set_ylable = '$M_{'+colorname1+'}-M_{'+colorname2+'}$'

###NUVrJ_unscaled_struc_minus1point5to2_minus1to7
#Plot(0,0,1,1,1)#Plot_zbin(index, scale, struc, limit, line)

###NUVrJ_zbin_unscaled_struc_minus1point5to2_minus1to7
#Plot_zbin(0,0,1,1,1) #Plot_zbin(index, scale, struc, limit, line)
#Plot_zbin(1,0,1,1,1)
###NUVrJ_zbin_unscaled_minus1point5to2_minus1to7
#Plot_zbin(0,0,1,1,1)
#Plot_zbin(1,0,1,1,1)

'''
PrintNum_zbin(0)
print "---"
PrintNum_zbin(3)
print "---"
PrintNumdiff_zbin(0,2)
print "---"
PrintPercentage_zbin(0,2)
'''

#RegionFile(2, 'COSMOS850sources', 'red','4.0','NUMBER')
#RegionFile(1, 'COSMOS24sources', 'blue','3.0','NUMBER')
#RegionFile(3, 'COSMOS24850sources', 'yellow','2.0','NUMBER')
#RegionFile(0, 'COSMOSsources', 'pink','1.0','NUMBER')

data = data[0][mask[0]]

data.write('COSMOS2015_Laigle+_v1.1_simple_masked.fits')


#PlotHist_massmed(data24,mask24)

time2 = time.time()
print 'done! time =', time2-time1 , 'sec'