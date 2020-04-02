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
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
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

def Mask_error(mode,inputup,index):
    up = inputup
    if mode==1:
        mask_V_errorV = (data[index][V_error]<up) & (data[index][V_error]>0)
        #Suprime-Cam:i+band
        mask_ip_error = (data[index][ip_error]<up) & (data[index][ip_error]>0)
        mask_J_error = (data[index][J_error]<up) & (data[index][J_error]>0)
        mask_Ks_error = (data[index][Ks_error]<up) & (data[index][Ks_error]>0)
        return mask_V_errorV | mask_ip_error | mask_J_error | mask_Ks_error
    if mode==2:
        return (data[index][V_error]<up) & (data[index][V_error]>0)
    if mode==3:
        return (data[index][Ks_mag]<24)
    

def Mask_photoz(index):
    return (data[index][photoz]>0) #& (data[index][photoz]<8)

def Mask_class_star(index):
    return (data[index][class_star]==0) | (data[index][class_star]==2)###

def Mask_classQG(index):
    return (data[index][class_SFG]==0)
def Mask_classSFG(index):
    return (data[index][class_SFG]==1)
def Mask_myclassQG(index):
    return (y[index]>3.1) & (y[index]>3.0*x[index]+1.0)
def Mask_myclassSFG(index):
    return (y[index]<=3.1) | (y[index]<=3.0*x[index]+1.0)

def Mask_mass(index):
    return (data[index][mass]>11)

#def Mask_mass(index,low,high):
#    return (data[index][mass]<high) & (data[index][mass]>low)
    
def PlotHist_photoz(index,inputcolor,inputbins,labelname):
    NBINS = inputbins
    plt.hist(data[index][photoz][mask[index]], NBINS, color=inputcolor, alpha=0.8, label = labelname)
    plt.title('Histogram photoz')
    plt.ylabel('Number')
    plt.xlabel('z')
    plt.legend()
    plt.show()
    return

def PlotHist_photoz_para(indexlist,inputcolorlist,inputlabel,inputbins):
    NBINS = inputbins
    plotdata = []
    for i in range(len(indexlist)):
        plotdata.append(data[indexlist[i]][photoz][mask[indexlist[i]]])
    plt.hist(plotdata, NBINS, color=inputcolorlist, alpha=0.7,label=inputlabel)
    plt.title('COSMOS2015')
    plt.ylabel('Number', fontdict = {'fontsize' : 14})
    plt.xlabel('z', fontdict = {'fontsize' : 14})
    plt.xlim(0,3.5)
    plt.legend()
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

def PlotHist_mass(index,inputcolor,inputbins,labelname):
    NBINS = inputbins
    plt.hist(data[index][mass][mask[index]], NBINS, color=inputcolor, alpha=0.8, label = labelname)
    plt.title('Histogram mass')
    plt.ylabel('Number')
    plt.xlabel('log(mass)')
    plt.legend()
    plt.show()
    return

def PlotHist_mass_para(indexlist,inputcolorlist,inputlabel,inputbins):
    NBINS = inputbins
    plotdata = []
    for i in range(len(indexlist)):
        plotdata.append(data[indexlist[i]][mass][mask[indexlist[i]]])
    plt.hist(plotdata, NBINS, color=inputcolorlist, alpha=0.7,label=inputlabel)
    plt.title('COSMOS2015')
    plt.ylabel('Number', fontdict = {'fontsize' : 14})
    plt.xlabel('log(mass)', fontdict = {'fontsize' : 14})
    plt.legend()
    plt.show()
    return

def DrawLine():
    plt.plot([-1.5,(3.1-1.0)/3.0],[3.1,3.1],'k-',lw=0.7)
    plt.plot([(3.1-1.0)/3.0,2.0],[3.1,3.0*2.0+1.0],'k-',lw=0.7)
    return

def Plot(index, scale, struc, limit, line, inputcolor, label, labelname):
    if (struc==1):
        set_s = 0.5
        set_alpha = 0.1
    else:
        set_s = 1.5
        set_alpha = 1
    plt.scatter( x_masked[index], y_masked[index], s=set_s, alpha=set_alpha,color=inputcolor, label=labelname)
    
    maskQG = mask[index] & Mask_myclassQG(index)
    maskSFG = mask[index] & Mask_myclassSFG(index)
    plt.title(colorname1+colorname2+colorname3, fontdict = {'fontsize' : 16})
    
    plt.title(colorname1+colorname2+colorname3+'\n'\
              +str(len(data[index].filled()[maskQG]))+' QGs ('+\
              str( "%.2f" % (len(data[index].filled()[maskQG])*100/Num_zbin_total(index)) )+'%), '\
              +str(len(data[index].filled()[maskSFG]))+' SFGs('+
              str( "%.2f" % (len(data[index].filled()[maskSFG])*100/Num_zbin_total(index)) )+'%)')
    
    plt.xlabel(set_xlable, fontdict = {'fontsize' : 14})
    plt.ylabel(set_ylable, fontdict = {'fontsize' : 14})
    if (limit==1):
        plt.axis([-1.5,2,-1,7])
        #plt.axis([-0.5,2,-1,7])
    if (scale==1):
        plt.axis('scaled')
    if (line==1):
        DrawLine()
    if (label==1):
        plt.legend()
    plt.show()
    return

def PrintQGfraction(index):
    maskQG = mask[index] & Mask_myclassQG(index)
    all_num = Num_zbin_total(index)
    QG_num = len(data[index].filled()[maskQG])
    print all_num, QG_num
    fraction = QG_num*100.0 / all_num
    error =  np.sqrt(QG_num)*100.0 / all_num
    #error = np.sqrt( (1/QG_num)*100.0/all_num )
    #error = 100.0 / all_num * np.sqrt( (1/QG_num)+(QG_num**2/all_num**3) )
    return fraction, error


def PrintAGNfraction(index):
    maskQG = mask[1]
    maskAGN = mask[index] & Mask_myclassQG(index)
    AGN_num = len(data[index].filled()[maskAGN])
    QG_all_num = len(data[index].filled()[maskQG])
    print AGN_num, QG_all_num
    fraction = AGN_num*100.0 / QG_all_num
    error =  np.sqrt(AGN_num)*100.0 / QG_all_num
    #error = np.sqrt( (1/QG_num)*100.0/all_num )
    #error = 100.0 / all_num * np.sqrt( (1/QG_num)+(QG_num**2/all_num**3) )
    return fraction, error


def Construct_mask_zbin(index):
    mask_zbin = []
    mask_zbin.append((data[index][photoz]>0) & (data[index][photoz]<0.5))
    mask_zbin.append((data[index][photoz]>0.5) & (data[index][photoz]<1))
    mask_zbin.append((data[index][photoz]>1) & (data[index][photoz]<1.5))
    mask_zbin.append((data[index][photoz]>1.5) & (data[index][photoz]<2))
    mask_zbin.append((data[index][photoz]>2) & (data[index][photoz]<2.5))
    mask_zbin.append((data[index][photoz]>2.5) & (data[index][photoz]<3))
    mask_zbin.append((data[index][photoz]>3) & (data[index][photoz]<3.5))
    mask_zbin.append((data[index][photoz]>3.5))
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
    zbin_title.append('$z>3.5$')
    return zbin_title

def Plot_zbin(index, scale, struc, limit, line, inputcolor,labelname):
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
        datamaskQG = mask[index] & mask_zbin[i] & Mask_myclassQG(index)
        datamask1 = mask[1] & mask_zbin[i]
        datamaskQG1 = mask[1] & mask_zbin[i] & Mask_myclassQG(1)
        datax = x[index][datamask]
        datay = y[index][datamask]
        plt.scatter( datax, datay, s=set_s, alpha=set_alpha,color=inputcolor, label=labelname)
        #print len(x[index][datamaskQG].filled())
        if len(x[index][datamask].filled())==0:
            plt.title(zbin_title[i]+', 0.00%')
        else:
            plt.title(zbin_title[i]+', '+str(len(x[index][datamask].filled()))+', '+\
                      str( "%.2f" % ((len(x[index][datamaskQG].filled()))*100.0/len(x[index][datamask].filled())) )+'%')
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
        plt.legend()
    plt.show()
    return

def Print_zbin_QGfraction(index):
    fraction = []
    error = []
    mask_zbin = Construct_mask_zbin(index)
    for i in range(8):
        datamask = mask[index] & mask_zbin[i]
        datamaskQG = mask[index] & mask_zbin[i] & Mask_myclassQG(index)
        all_num = len(x[index][datamask].filled())
        QG_num = len(x[index][datamaskQG].filled())
        
        if (all_num==0)|(QG_num==0):
            fraction.append(0.0)
            error.append(0.0)
        else:
            fraction.append(QG_num*100.0/all_num)
            error.append( np.sqrt(QG_num)*100.0/all_num )
            #error.append(100.0 / all_num * np.sqrt( (1/QG_num)+(QG_num**2/all_num**3) ))
            print  all_num, QG_num, error[i]
    return fraction, error

def PrintNum_zbin(index):
    mask_zbin = Construct_mask_zbin(index)
    #print len( data[index].filled() )
    print len( data[index].filled()[mask[index]] )
    for i in range(8):
        print len( data[index].filled()[mask[index] & mask_zbin[i]] )
    return

def PrintNum_zbin_total(index):
    print len( data[index].filled()[mask[index]] )
    return

def Num_zbin_total(index): 
    return len( data[index].filled()[mask[index]] )

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

def RegionFile(index,filename,color,size):
    f = open(filename+'.reg','w')
    f.write('global color='+color+' font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    N = len( data[index].filled()[mask[index]] )
    print N
    for n in range(N):
        #f.write('fk5;circle('+str(data[index][ra][mask[index]][n])+','+str(data[index][dec][mask[index]][n])+','+size+'") # text={'+str(n)+'}\n')
        #f.write('fk5;circle('+str(data[index][ra][mask[index]][n])+','+str(data[index][dec][mask[index]][n])+','+size+'") # text={''}\n')
        f.write('fk5;circle('+str(data[index][ra][mask[index]][n])+','+str(data[index][dec][mask[index]][n])+','+size+'") # text={'+str(data[index][photoz][mask[index]][n])+'}\n')
        #f.write('fk5;circle('+str(data[index][ra][mask[index]][n])+','+str(data[index][dec][mask[index]][n])+','+size+'") # text={'+str(data[index]["NUMBER"][mask[index]][n])+'}\n')
    f.close()

def PlotMatchedResult(matchedband):
    
    #twomatch:1, onematch24/3:2/3, nomatch:4 rgcb
    if matchedband=='850ALL':
        mask.append(((data[0]['850NARROW']==1) | (data[0]['850SOURCE']==1)) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#& Mask_myclassQG(0))
        mask.append(((data[0]['850NARROW']==2) | (data[0]['850SOURCE']==2)) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#& Mask_myclassQG(0))
        mask.append(((data[0]['850NARROW']==3) | (data[0]['850SOURCE']==3)) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#& Mask_myclassQG(0))
        mask.append(((data[0]['850NARROW']==4) | (data[0]['850SOURCE']==4)) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#& Mask_myclassQG(0))
    else:
        mask.append((data[0][matchedband]==1) & Mask_M(0) & Mask_photoz(0) & Mask_error(2,0.2,0) & Mask_class_star(0))#& Mask_myclassQG(0))
        mask.append((data[0][matchedband]==2) & Mask_M(0) & Mask_photoz(0) & Mask_error(2,0.2,0) & Mask_class_star(0))#& Mask_myclassQG(0))
        mask.append((data[0][matchedband]==3) & Mask_M(0) & Mask_photoz(0) & Mask_error(2,0.2,0) & Mask_class_star(0))#& Mask_myclassQG(0))
        mask.append((data[0][matchedband]==4) & Mask_M(0) & Mask_photoz(0) & Mask_error(2,0.2,0) & Mask_class_star(0))#& Mask_myclassQG(0))
    for i in range(4):
        x_masked.append(x[i][mask[i]])
        y_masked.append(y[i][mask[i]])
    
    plt.figure(1)
    labellist = ['both','24 micron','3 GHz','none']
    PlotHist_photoz_para([0,1,2,3],['r','g','c','b'],labellist,10)
    
    plt.figure(2)
    Plot(0,0,0,1,1,'r',1,'both matched')
    Plot(1,0,0,1,1,'g',1,'24 micron matched')
    Plot(2,0,0,1,1,'c',1,'3GHz matched')
    Plot(3,0,0,1,1,'b',1,'none matched')
    plt.title(colorname1+colorname2+colorname3)
    
    plt.figure(3)
    plt.figure(figsize=(10,20))
    plt.subplot(2, 2, 2)
    PlotHist_photoz(0,'r',10,'both matched')
    plt.subplot(2, 2, 4)
    PlotHist_photoz(1,'g',10,'24 micron matched')
    plt.subplot(2, 2, 1)
    PlotHist_photoz(2,'c',10,'3GHz matched') 
    plt.subplot(2, 2, 3)
    PlotHist_photoz(3,'b',10,'none matched') 
    plt.tight_layout()
    
    plt.figure(4)
    plt.figure(figsize=(10,20))
    plt.subplot(2, 2, 2)
    Plot(0,0,0,1,1,'r',1,'both matched, '+str(Num_zbin_total(0))+' samples')
    plt.subplot(2, 2, 4)
    Plot(1,0,0,1,1,'g',1,'24 micron matched, '+str(Num_zbin_total(1))+' samples')
    plt.subplot(2, 2, 1)
    Plot(2,0,0,1,1,'c',1,'3GHz matched, '+str(Num_zbin_total(2))+' samples')
    plt.subplot(2, 2, 3)
    Plot(3,0,0,1,1,'b',1,'none matched, '+str(Num_zbin_total(3))+' samples')
    plt.tight_layout()
    return

# =============================================================================
# main code
# =============================================================================

time1 = time.time()

###set colors
colorname1 = "NUV"
colorname2 = "r"
colorname3 = "J"

###set columns in the main catalog
color1 = "MNUV"
color2 = "MR"
color3 = "MJ"
color = [ color1, color2, color3 ]
photoz = "REDSHIFT"
V_error = "V_MAGERR_APER3"
ip_error = "ip_MAGERR_APER3"
J_error = "J_MAGERR_APER3"
Ks_error = "Ks_MAGERR_APER3"
Ks_mag = "Ks_MAG_APER3"
mass = "MASS_MED"
class_star = "TYPE"
class_SFG = "CLASS"
ra = "ALPHA_J2000"
dec = "DELTA_J2000"
mass = "MASS_BEST"

set_xlable = '$M_{'+colorname2+'}-M_{'+colorname3+'}$'
set_ylable = '$M_{'+colorname1+'}-M_{'+colorname2+'}$' 

########
# main #
########

'''
###set catalogs
number = 5
catalog = [None]*number
catalog[0] = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow_simple.fits'
#"COSMOS2015_Laigle+_v1.1_850sources_simple.fits"
#"COSMOS2015_Laigle+_v1.1_850sources_1.fits"
#"01_COSMOS2015catalog/COSMOS2015/COSMOS2015_Laigle+_v1.1.fits"
catalog[1] = "COSMOS+mips24_allmatches_simple.fits"
catalog[2] = "COSMOS+wide850_allmatches_simple.fits"
catalog[3] = "COSMOS+mips24＋wide850_allmatches.fits"
catalog[4] = "COSMOS+wide850_bestmatchfor850.fits"

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
    #mask.append(data[0]['850SOURCE']==4)
    mask.append( Mask_M(i) & Mask_photoz(i) & Mask_error(2,0.1,i) & Mask_class_star(i) )
    #& (data[0]['850SOURCE']!=0)
    #twomatch:1, onematch:2/3, nomatch:4 rgcb
    #& Mask_myclassQG(i)
    #& Mask_mass(i,7.0,9.8)
    #& Mask_myclassQG(i) & Mask_classSFG(i)
    #x_masked.append(x[i][mask[i]])
    #y_masked.append(y[i][mask[i]])

#(data24['flux_24_2']>80)

###plot
#PlotHist_photoz(0,'C0',50)
#PlotHist_M(0)

###NUVrJ_unscaled_struc_minus1point5to2_minus1to7
#Plot(0,0,1,1,1,'C0') #Plot_zbin(index, scale, struc, limit, line)

###NUVrJ_zbin_unscaled_struc_minus1point5to2_minus1to7
#Plot_zbin(0,0,0,1,1) #Plot_zbin(index, scale, struc, limit, line)
#Plot_zbin(1,0,1,1,1)
###NUVrJ_zbin_unscaled_minus1point5to2_minus1to7
#Plot_zbin(0,0,1,1,1)
#Plot_zbin(1,0,1,1,1)

#PrintNum_zbin_total(0)
#PrintNum_zbin_total(1)
#PrintNum_zbin_total(2)
#PrintNum_zbin_total(3)
#PrintNum_zbin(0)
    
print "---"
PrintNum_zbin(3)
print "---"
PrintNumdiff_zbin(0,2)
print "---"
PrintPercentage_zbin(0,2)


#RegionFile(2, 'COSMOS850sources', 'red','5.0')
#RegionFile(1, 'COSMOS24sources', 'blue','4.0')
#RegionFile(0, 'COSMOS24sources', 'blue','1.0')

#PlotHist_massmed(data24,mask24)
'''

############
# 850 wide #
############
'''
number = 5
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []

PlotMatchedResult('850SOURCE')
'''

##############
# 850 narrow #
##############
'''
number = 5
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []

PlotMatchedResult('850NARROW')
'''

######################
# 850 wide or narrow #
######################
'''
number = 5
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []

PlotMatchedResult('850ALL')
'''

##############
# 450 narrow #
##############
'''
number = 5
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []

PlotMatchedResult('450NARROW')
'''

#############
# 24 micorn #
#############
'''
number = 2
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []

mask.append( Mask_M(0) & Mask_photoz(0) & Mask_error(2,0.2,0) & Mask_class_star(0))#\
        #& Mask_myclassQG(0))
mask.append( (data[0]['24MICRON']==1)\
        & Mask_M(0) & Mask_photoz(0) & Mask_error(2,0.2,0) & Mask_class_star(0))#\
        #& Mask_myclassQG(0))

for i in range(2):
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])

plt.figure(1)
labellist = ['all','24 micron source']
PlotHist_photoz_para([0,1],['C0','C1'],labellist,20)

plt.figure(2)
Plot(0,0,1,1,1,'C0',1,'all, '+str(Num_zbin_total(0))+' samples')
Plot(1,0,1,1,1,'C1',1,'24 micron source, '+str(Num_zbin_total(1))+' samples')
'''

#########
# 3 GHz #
#########
'''
number = 2
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []

mask.append( Mask_M(0) & Mask_photoz(0) & Mask_error(2,0.2,0) & Mask_class_star(0))#\
        #& Mask_myclassQG(0))
mask.append( (data[0]['3GHZ']==1)\
        & Mask_M(0) & Mask_photoz(0) & Mask_error(2,0.2,0) & Mask_class_star(0))#\
        #& Mask_myclassQG(0))

for i in range(2):
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])

plt.figure(1)
labellist = ['all','3 GHz source']
PlotHist_photoz_para([0,1],['C0','C1'],labellist,20)

plt.figure(2)#Plot_zbin(index, scale, struc, limit, line)
Plot(0,0,1,1,1,'C0',1,'all, '+str(Num_zbin_total(0))+' samples')
Plot(1,0,0,1,1,'C1',1,'3 GHz source, '+str(Num_zbin_total(1))+' samples')
'''

######################
# 3 GHz or 24 micron #
######################
'''
number = 2
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []

mask.append( Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#\
        #& Mask_myclassQG(0))
mask.append( ( (data[0]['3GHZ']==1) | (data[0]['24MICRON']==1) )\
        & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#\
        #& Mask_myclassQG(0))

for i in range(2):
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])

plt.figure(1)
labellist = ['all','3 GHz or 24 micron source']
PlotHist_photoz_para([0,1],['C0','C1'],labellist,20)

plt.figure(2)#Plot_zbin(index, scale, struc, limit, line)
Plot(0,0,1,1,1,'C0',1,'all, '+str(Num_zbin_total(0))+' samples')
Plot(1,0,0,1,1,'C1',1,'3 GHz or 24 micron source, '+str(Num_zbin_total(1))+' samples')
'''

######################
# all dusty galaxies #
######################
'''
number = 2
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []

mask.append( Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#\
        #& Mask_myclassQG(0))
mask.append(  (  (data[0]['24MICRON']==1) | (data[0]['3GHZ']==1) | (data[0]['850SOURCE']!=0)\
               |(data[0]['850NARROW']!=0)|(data[0]['450NARROW']!=0) )
        & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0)    )#\
        #& Mask_myclassQG(0))

for i in range(2):
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])

plt.figure(1)
labellist = ['all','IR and radio source']
PlotHist_photoz_para([0,1],['C0','C1'],labellist,20)

plt.figure(2)#Plot_zbin(index, scale, struc, limit, line)
Plot(0,0,1,1,1,'C0',1,'all, '+str(Num_zbin_total(0))+' samples')
Plot(1,0,0,1,1,'C1',1,'IR and radio source, '+str(Num_zbin_total(1))+' samples')
'''


###test different selection
'''
number = 5
catalog = [None]*number
catalog[0] = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []

mask.append((data[0]['850SOURCE']==4) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))
mask.append( Mask_M(0) & Mask_photoz(0) & Mask_error(2,0.2,0) & Mask_class_star(0))#& Mask_myclassQG(0))
mask.append( Mask_M(0) & Mask_photoz(0) & Mask_error(2,0.1,0) & Mask_class_star(0))#& Mask_myclassQG(0))
mask.append( Mask_M(0) & Mask_photoz(0) & Mask_error(3,0,0) & Mask_class_star(0))#& Mask_myclassQG(0))

for i in range(4):
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])

plt.figure(1)
PlotHist_photoz(0,'C0',50)

plt.figure(1)
PlotHist_photoz(1,'palevioletred',50)
plt.figure(2)
PlotHist_photoz(2,'coral',50) 
plt.figure(1)
PlotHist_photoz(3,'peru',50) 

plt.figure(1)
PlotHist_photoz(1,'seagreen',50)
plt.figure(2)
PlotHist_photoz(2,'yellowgreen',50) 

plt.figure(2)
labellist = ['4 bands error < 0.1','K band error < 0.2','K band error < 0.1','K band < 24']
PlotHist_photoz_para([0,1,2,3],['C0','palevioletred','coral','peru'],labellist)

plt.figure(3)
labellist = ['4 bands error < 0.1','V band error < 0.2','V band error < 0.1']
PlotHist_photoz_para([0,1,2],['C0','seagreen','yellowgreen'],labellist)

plt.figure(4)
Plot(0,0,1,1,1,'C0')

plt.figure(5)
Plot(1,0,1,1,1,'palevioletred')
plt.figure(6)
Plot(2,0,1,1,1,'coral')
plt.figure(3)
Plot(3,0,1,1,1,'peru')

plt.figure(4)
Plot(1,0,1,1,1,'seagreen')
plt.figure(5)
Plot(2,0,1,1,1,'yellowgreen')
'''

# radio galaxies without 24 detection #
'''
number = 2
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
mask.append((data[0]['3GHZ']==1) & (data[0]['24MICRON']==0) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0) & Mask_myclassQG(0))
mask.append((data[0]['3GHZ']==1) & (data[0]['24MICRON']==0) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0) )

RegionFile(0, 'COSMOS3without24QG', 'orange','1.5')
RegionFile(1, 'COSMOS3without24', 'green','3.0')
'''

# 450 micorn stacking #
'''
number = 2
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
mask.append((data[0]['450NARROW']!=0) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0) & Mask_myclassQG(0))
mask.append( Mask_M(1) & Mask_photoz(1) & Mask_error(1,0.1,1) & Mask_class_star(1) & Mask_myclassQG(1))

RegionFile(0, 'COSMOS450QG', 'green','1.5')
RegionFile(1, 'COSMOSQG', 'blue','2.0')
'''
'''
ra = "RA_450"
dec = "DEC_450"
number = 1
catalog = [None]*1
catalog[0] = "04_COSMOS450_850/STUDIES/sources_450.fits" 

###read catalog
data = [None]*number
data[0] = Table.read(catalog[0], hdu=1)

###mask data
mask = []
mask.append(data[0][ra]>0)

RegionFile(0, '450source', 'yellow','4.0')
'''


# redshift clump #
'''
number = 1
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
mask1 = data[0]['450NARROW']==1
mask2 = Mask_M(0) & Mask_error(1,0.1,0) & Mask_class_star(0)
mask3 =(data[0][photoz]>0.5) & (data[0][photoz]<1.0)
mask.append(mask1 & mask2 & mask3)

PrintNum_zbin_total(0)

RegionFile(0, 'COSMOS450_both_zclump', 'red','5.0')
'''

# 3GHz without 24 micron #
'''
number = 2
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz+iragnall_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []

mask.append( Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#\
        #& Mask_myclassQG(0))
mask.append( (data[0]['3GHZ']==1) & (data[0]['24MICRON']==0) & (data[0]['IR_AGN_ALL']==1)\
        & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#\
        #& Mask_myclassQG(0))

for i in range(2):
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])

plt.figure(1)
labellist = ['all','3 GHz without 24 micron detection']
PlotHist_photoz_para([0,1],['C0','C1'],labellist,20)

plt.figure(2)#Plot_zbin(index, scale, struc, limit, line)
Plot(0,0,1,1,1,'C0',1,'all, '+str(Num_zbin_total(0))+' samples')
Plot(1,0,0,1,1,'C1',1,'3 GHz without 24 micron detection, '+str(Num_zbin_total(1))+' samples')
'''

# ir agn #
'''
number = 2
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz+iragnallmir_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []

mask.append( Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#\
        #& Mask_myclassQG(0))
mask.append(  (data[0]['IR_AGN_MIR']==1)&(data[0]['IR_AGN_ALL']==1)\
        & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#\
        #& Mask_myclassQG(0))

for i in range(2):
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])

plt.figure(1)
labellist = ['all','IR AGN (mir)']
PlotHist_photoz_para([0,1],['C0','C1'],labellist,20)

plt.figure(2)#Plot_zbin(index, scale, struc, limit, line)
Plot(0,0,1,1,1,'C0',1,'all, '+str(Num_zbin_total(0))+' samples')
Plot(1,0,0,1,1,'C1',1,'IR AGN (mir), '+str(Num_zbin_total(1))+' samples')

plt.figure(3)
Plot_zbin(1, 0, 0, 1, 1)
'''
'''
number = 2
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz+iragnallmir_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []

mask.append( (data[0]['IR_AGN_MIR']==1)\
        & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#\
        #& Mask_myclassQG(0))
mask.append(  (data[0]['IR_AGN_MIR']==1)&(data[0]['3GHZ']==1)&(data[0]['24MICRON']==0)\
        & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#\
        #& Mask_myclassQG(0))

for i in range(2):
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])

plt.figure(1)
labellist = ['IR AGN (mir)','with 3 GHz but without 24 micron']
PlotHist_photoz_para([0,1],['C1','b'],labellist,20)

plt.figure(2)#Plot_zbin(index, scale, struc, limit, line)
Plot(0,0,0,1,1,'C1',1,'IR AGN (mir), '+str(Num_zbin_total(0))+' samples')
Plot(1,0,0,1,1,'b',1,'with 3 GHz but without 24 micron, '+str(Num_zbin_total(1))+' samples')
'''
'''
number = 2
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz+iragnallmir_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []

mask.append( (data[0]['IR_AGN_MIR']==1)\
        & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#\
        #& Mask_myclassQG(0))
mask.append(  (data[0]['IR_AGN_MIR']==1)&(data[0]['24MICRON']==1)\
        & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#\
        #& Mask_myclassQG(0))

for i in range(2):
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])

plt.figure(1)
labellist = ['IR AGN (mir)','IR AGN with 24 micron detection']
PlotHist_photoz_para([0,1],['C1','g'],labellist,20)

plt.figure(2)#Plot_zbin(index, scale, struc, limit, line)
Plot(0,0,0,1,1,'C1',1,'IR AGN (mir), '+str(Num_zbin_total(0))+' samples')
Plot(1,0,0,1,1,'g',1,'IR AGN with 24 micron detection, '+str(Num_zbin_total(1))+' samples')
'''

# 3 GHz AGN #
'''
number = 2
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_5band_2agn_9cat_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []

mask.append( (data[0]['3GHZ']==1)\
        & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#\
        #& Mask_myclassQG(0))
mask.append(  (data[0]['HLAGN']==1)\
        & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#\
        #& Mask_myclassQG(0))

for i in range(number):
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])

plt.figure(1)
labellist = ['3 GHz','HLAGN',]
PlotHist_photoz_para([0,1],['C1','r'],labellist,20)

plt.figure(2)#Plot_zbin(index, scale, struc, limit, line)
Plot(0,0,1,1,1,'C1',1,'3GHz, '+str(Num_zbin_total(0))+' samples')
Plot(1,0,0,1,1,'r',1,'HLAGN, '+str(Num_zbin_total(1))+' samples')
'''

# poster #
'''
number = 1
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_5band_2agn_9cat_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []

mask.append( Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))

for i in range(number):
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])

fig = plt.figure(1)
labellist = ['COSMOS2015']
PlotHist_photoz_para([0],['C0'],labellist,20)
#fig.savefig('COSMOS2015_hist.png', format='png', dpi=1200)

fig = plt.figure(2)#Plot_zbin(index, scale, struc, limit, line)
Plot(0,0,1,1,1,'C0',0,'COSMOS2015, '+str(Num_zbin_total(0))+' samples')
#fig.savefig('COSMOS2015.png', format='png', dpi=1200)

'''
'''
number = 3
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_5band_2agn_9cat_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []
mask.append( Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))
mask.append( (data[0]['24MICRON']==1)\
        & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))
mask.append( (data[0]['3GHZ']==1)\
        #((data[0]['450NARROW']==1)|(data[0]['24MICRON']==1)|(data[0]['SFG']==1))&(data[0]['HLAGN']!=1)
        & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))

for i in range(number):
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])

fig = plt.figure(1)
labellist = ['COSMOS2015','24 micron',"3 GHz"]
PlotHist_photoz_para([0,1,2],['C0','C1','r'],labellist,20)
#fig.savefig('COSMOS2015_3_hist.png', format='png', dpi=1200)

fig = plt.figure(2)#Plot_zbin(index, scale, struc, limit, line)
Plot(0,0,1,1,1,'C0',0,'COSMOS2015, '+str(Num_zbin_total(0))+' samples')
Plot(1,0,1,1,1,'C1',0,'24 micron, '+str(Num_zbin_total(1))+' samples')
Plot(2,0,1,1,1,'r',0,"3GHz, "+str(Num_zbin_total(2))+' samples')
#fig.savefig('COSMOS2015_3.png', format='png', dpi=1200)
'''

# 450 detected QGs case study #
'''
number = 2
catalog = [None]*1
catalog[0] = 'COSMOS2015_merged.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

CENTER_RA,CENTER_DEC = '10h00m25.0s', '2d24m22.0s'
RADIUS = 0.2*u.degree
def MaskRegion(index,inputra,inputdec,inputradius):
    c0 = SkyCoord(ra=data[index][ra], dec=data[index][dec])
    center = SkyCoord(inputra, inputdec, frame='fk5')
    radius = inputradius
    sep = center.separation(c0)
    return sep<=radius

###mask data
mask = []
x_masked = []
y_masked = []

#mask_tmp = (data[0]['FLAG_HJMCC']==0) &(data[0]['FLAG_PETER']==0)& (data[0]['FLAG_COSMOS']==1)

#Mask_M(0) & Mask_photoz(0) & Mask_class_star(0) & Mask_error(1,0.1,0)
mask.append( Mask_M(0) & Mask_photoz(0) & Mask_class_star(0) & Mask_error(1,0.1,0))
mask.append(Mask_M(0) & Mask_photoz(0) & Mask_class_star(0) & Mask_error(1,0.1,0) & (data[0]['agn_c17b']==True))

#(data[0]['3GHZ']==1) & (data[0]['Clean_SFG']!=1)

#'10h01m42.2187s','1d40m35.795s',0.05*u.degree
#1 '10h01m32.1382s','2d04m28.354s'
#2 '9h58m51.6261s','2d03m51.510s'
#3 '10h02m24.0136s','2d38m23.892s'

for i in range(number):
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])

#plt.figure(1)
#Plot(0, 0, 1, 1, 1, 'C0', 1,'IRAGN, '+str(Num_zbin_total(0))+' samples')
#Plot(1, 0, 0, 1, 1, 'C1', 1,'ALMA, '+str(Num_zbin_total(1))+' samples')
PlotHist_mass(1,"C1",20,"")

#Plot(2, 0, 0, 1, 1, 'darkseagreen', 1,'no counterpart, '+str(Num_zbin_total(2))+' samples')

#Plot(1, 0, 0, 1, 1, 'black', 1,'850 detected, '+str(Num_zbin_total(1))+' samples')
#plt.figure(1)
#labellist = ['3 with 24','3 without 24']
#PlotHist_photoz_para([0,1],['C0','g'],labellist,10)

#RegionFile(0, 'COSMOS_HLAGN_QGonlyone', 'yellow','10.0')
#RegionFile(1, 'COSMOS_SFG_Kselected_z', 'orange','1.0')
#RegionFile(2, 'COSMOS_all_SFG_z', 'pink','0.8')
#RegionFile(3, 'COSMOS_450QG_simple', 'yellow','0.8')
'''



'''
number = 1
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_5band_2agn_9cat_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []
mask.append( (data[0]['450NARROW']!=0)\
            &Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0)\
            &Mask_myclassQG(0))

for i in range(number):
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])

fig = plt.figure(1)
labellist = ['450QG']
PlotHist_photoz_para([0],['g'],labellist,5)

fig = plt.figure(2)#Plot_zbin(index, scale, struc, limit, line)
Plot(0,0,0,1,1,'g',0,'450QG, '+str(Num_zbin_total(0))+' samples')
num_array = range( Num_zbin_total(0) )
for i in num_array:
    plt.text(x_masked[0][i]-0.05, y_masked[0][i], i, fontsize=6)
'''
'''
number = 2
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_5band_2agn_9cat_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []
mask.append( (data[0]['HLAGN']==1)\
            & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0)\
            & Mask_myclassQG(0) )
mask.append( (data[0]['HLAGN']==1)\
            & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))

for i in range(number):
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])

fig = plt.figure(1)
labellist = ['3GHZ QG, '+str(Num_zbin_total(0))+' samples']
PlotHist_photoz_para([0],['r'],labellist,20)

fig = plt.figure(2)
Plot_zbin(1, 0, 0, 1, 1, 'r')
'''
'''
number = 5
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_5+2band_2agn_9cat_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []
mask.append(  Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0)\
            & Mask_mass(0) )#& Mask_myclassQG(0))
mask.append( Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0)\
            & (data[0]['3GHZ']==1) & (data[0]['Clean_SFG']!=1) & Mask_mass(0))
mask.append( Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0)\
            & (data[0]['3GHZ']==1) & (data[0]['HLAGN']==1) & Mask_mass(0))#& Mask_myclassQG(0))
mask.append( Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0)\
            & (data[0]['3GHZ']==1) & (data[0]['MLAGN']==1) & Mask_mass(0))#& Mask_myclassQG(0))
mask.append( Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0)\
            & (data[0]['3GHZ']==1) & Mask_mass(0))

for i in range(number):
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])

fig = plt.figure(1)
#PlotHist_mass(1,'C0',30,'COSMOS2015')
#PlotHist_mass_para([2,3],['r','g'],['HLAGN','MLAGN'],30)

#fig = plt.figure(2)
#Plot(0, 0, 1, 1, 1, 'C0',1,'all, '+str(Num_zbin_total(0))+' samples')
#Plot(1, 0, 0, 1, 1, 'k',1,'AGN, '+str(Num_zbin_total(1))+' samples')
#Plot(2, 0, 0, 1, 1, 'r',1,'HLAGN, '+str(Num_zbin_total(2))+' samples')
#Plot(3, 0, 0, 1, 1, 'darkgreen',1,'MLAGN, '+str(Num_zbin_total(3))+' samples')
#Plot(4, 0, 0, 1, 1, 'C1',1,'3 GHz, '+str(Num_zbin_total(4))+' samples')

#fig = plt.figure(3)
#Plot_zbin(0, 0, 0, 1, 1, 'C0','all')
#Plot_zbin(2, 0, 0, 1, 1, 'r','HLAGN')
#fig = plt.figure(4)
#Plot_zbin(3, 0, 0, 1, 1, 'g','MLAGN')
'''

# 850 region file #
'''
number = 1
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_5+1band_2agn_9cat_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
mask.append( (data[0]['850SOURCE']!=0)&Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0) & Mask_myclassQG(0))

RegionFile(0, 'COSMOS_850wideQG', 'yellow','15.0')
'''
'''
number = 1
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_5+2band_2agn_9cat_simple.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

    
for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []
mask.append( (data[0]['850WIDE_SIMPLE']==1)&Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0) & Mask_myclassQG(0))
for i in range(number):
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])
    
Plot(0,0,0,1,1,'magenta',0,'850QG, '+str(Num_zbin_total(0))+' samples')
print Num_zbin_total(0)
'''

number = 6
catalog = [None]*1
catalog[0] = 'COSMOS2015_merged.fits' 

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []
mask_all = Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0) & Mask_mass(0)
mask.append( mask_all )

mask.append( mask_all & Mask_myclassQG(0))
mask.append( mask_all &\
            ( (data[0]['Radio_excess']==1)|(data[0]['agn_c17b']==True)|(data[0]['agn_xxx']==True) ))
mask.append( mask_all & (data[0]['Radio_excess']==1) )
mask.append( mask_all & (data[0]['agn_c17b']==True) )
mask.append( mask_all & (data[0]['agn_xxx']==True) )

#& Mask_mass(0)

for i in range(number):
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])


'''
plt.figure(1)
Plot(0, 0, 1, 1, 1, 'C0', 1,'all, '+str(Num_zbin_total(0))+' samples')
plt.figure(2)
Plot(2, 0, 0, 1, 1, 'C1', 1,'radio AGN, '+str(Num_zbin_total(2))+' samples')
plt.figure(3)
Plot(3, 0, 0, 1, 1, 'g', 1,'IR AGN, '+str(Num_zbin_total(3))+' samples')
plt.figure(4)
Plot(4, 0, 0, 1, 1, 'r', 1,'X-ray AGN, '+str(Num_zbin_total(4))+' samples')
'''


#Plot(5, 0, 0, 1, 1, 'r', 1,'HLAGN, '+str(Num_zbin_total(5))+' samples')
'''
plt.figure(1)
PlotHist_mass(2,'C1',30,'radio AGN')
plt.figure(2)
PlotHist_mass(3,'g',30,'IR AGN')
plt.figure(3)
PlotHist_mass(4,'r',30,'X-ray')
'''
'''
plt.figure(1)
Plot_zbin(0, 0, 0, 1, 1, 'C0','all')
print "all"
#plt.figure(2)
#Plot_zbin(1, 0, 0, 1, 1, 'k','QG')
#print "QG"
#PrintNum_zbin(1)
plt.figure(3)
Plot_zbin(2, 0, 0, 1, 1, 'C1','radio AGN')
print "radio AGN"
#PrintNum_zbin(2)
plt.figure(4)
Plot_zbin(3, 0, 0, 1, 1, 'g','IR AGN')
#print "IR AGN"
#PrintNum_zbin(5)
plt.figure(5)
Plot_zbin(4, 0, 0, 1, 1, 'r','X-ray AGN')
#print "X-ray AGN"
#PrintNum_zbin(4)
'''

'''
plt.figure(1,figsize=(10,6))
set_s = 0.1
set_alpha = 1
#plt.scatter( data[0][mask[0]][photoz], data[0][mask[0]][mass], s=set_s, alpha=set_alpha,color='C0', label='all')
plt.scatter( data[1][mask[1]][photoz], data[1][mask[1]][mass], s=0.3, alpha=set_alpha,color='k', label='QG')
plt.scatter( data[2][mask[2]][photoz], data[2][mask[2]][mass], s=0.3, alpha=set_alpha,color='C1', label='radio AGN')
#plt.scatter( data[3][mask[3]][photoz], data[3][mask[3]][mass], s=0.3, alpha=set_alpha,color='C1', label='IR AGN')
#plt.scatter( data[4][mask[4]][photoz], data[4][mask[4]][mass], s=0.3, alpha=set_alpha,color='C1', label='X-ray AGN')

zbin1 = [0.175,0.5,0.8,1.125,1.475,2.0,2.5,3.0,3.75]
zbin2 = [0.175,0.5,0.8,1.125,1.475,2.0,2.5,3.0,3.75,4.4]
m_lim_QG = [8.4,9.0,9.4,9.6,9.9,10.1,10.3,10.4,10.5]
m_lim_all = [8.1,8.7,9.1,9.3,9.7,9.9,10.0,10.1,10.1,10.1]
#plt.plot( zbin2, m_lim_all, 'o-C1', label='all mass limit Laigle(2016)', mfc='none')
plt.plot( zbin1, m_lim_QG, 'o-C3', label='QG mass limit Laigle(2016)', mfc='none')


plt.title('mass v.s. redshift', fontdict = {'fontsize' : 16})
plt.xlabel('z', fontdict = {'fontsize' : 14})
plt.ylabel('log(mass)', fontdict = {'fontsize' : 14})
plt.legend()
plt.show()
'''


fig, axes = plt.subplots()

plt.setp(axes, xticks=[y_axes+1 for y_axes in range(number)], xticklabels=['COSMOS2015','all AGN', 'radio AGN', 'IR AGN', 'X-ray AGN'])

fraction, error = PrintQGfraction(0)
print fraction, error
plt.errorbar(1,fraction,error,color='k', fmt='o', capsize=5)
plt.text(1.0-0.2, fraction-error-5.5, str("%.2f" % fraction), fontsize=10)
plt.hlines(fraction,0,6, colors='k', linestyles='dotted')

fraction, error = PrintQGfraction(2)
print fraction, error
plt.errorbar(2,fraction,error,color='C0', fmt='o', capsize=5)
plt.text(2.0-0.2, fraction-error-5.5, str("%.2f" % fraction), fontsize=10)

fraction, error = PrintQGfraction(3)
print fraction, error
plt.errorbar(3,fraction,error,color='r', fmt='o', capsize=5)
plt.text(3.0-0.2, fraction-error-5.5, str("%.2f" % fraction), fontsize=10)

fraction, error = PrintQGfraction(4)
print fraction, error
plt.errorbar(4,fraction,error,color='C1', fmt='o', capsize=5)
plt.text(4.0-0.2, fraction-error-5.5, str("%.2f" % fraction), fontsize=10)

fraction, error = PrintQGfraction(5)
print fraction, error
plt.errorbar(5,fraction,error,color='g', fmt='o', capsize=5)
plt.text(5.0-0.2, fraction-error-5.5, str("%.2f" % fraction), fontsize=10)


#axes.plot([1,2,3,4],[51.51,57.65,20.61,31.33],'o',['C0','C1','g','r'])
plt.axis([0,6,0,70])
#plt.title('QG percentage, mass>10^11', fontdict = {'fontsize' : 16})
plt.ylabel('QG fraction(%)', fontdict = {'fontsize' : 18})#14})
plt.show()


'''
fig, axes = plt.subplots(figsize=(8,6))

zbin = [0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75]

fraction, error = Print_zbin_QGfraction(0)
plt.errorbar(zbin,fraction,error,fmt='o-',color='k',label='COSMOS2015', capsize=5)

fraction, error = Print_zbin_QGfraction(2)
a=0.01
zbin = [0.25+a,0.75+a,1.25+a,1.75+a,2.25+a,2.75+a,3.25+a,3.75+a]
plt.errorbar(zbin,fraction,error,fmt='o-',color='C0',label='all AGN', capsize=5)

fraction, error = Print_zbin_QGfraction(3)
a=0.02
zbin = [0.25+a,0.75+a,1.25+a,1.75+a,2.25+a,2.75+a,3.25+a,3.75+a]
plt.errorbar(zbin,fraction,error,fmt='o-',color='r',label='radio AGN', capsize=5)

fraction, error = Print_zbin_QGfraction(4)
a=0.03
zbin = [0.25+a,0.75+a,1.25+a,1.75+a,2.25+a,2.75+a,3.25+a,3.75+a]
plt.errorbar(zbin,fraction,error,fmt='o-',color='C1',label='IR AGN', capsize=5)

fraction, error = Print_zbin_QGfraction(5)
a=0.04
zbin = [0.25+a,0.75+a,1.25+a,1.75+a,2.25+a,2.75+a,3.25+a,3.75+a]
plt.errorbar(zbin,fraction,error,fmt='o-',color='g',label='X-Ray AGN', capsize=5)

#plt.title('QG percentage, mass>10^11', fontdict = {'fontsize' : 20})
plt.xlabel('z', fontdict = {'fontsize' : 16})
plt.ylabel('QG fraction (%)', fontdict = {'fontsize' : 16})
plt.legend()
plt.show()

'''

'''
fig, axes = plt.subplots()

plt.setp(axes, xticks=[y_axes+1 for y_axes in range(5)], xticklabels=['all', 'radio AGN', 'IR AGN', 'X-ray AGN'])

fraction, error = PrintQGfraction(0)
print fraction, error
plt.errorbar(1,fraction,error,color='C0', fmt='o', capsize=5)
plt.text(1.0-0.2, fraction-error-5.5, str("%.2f" % fraction), fontsize=10)

fraction, error = PrintQGfraction(2)
print fraction, error
plt.errorbar(2,fraction,error,color='C1', fmt='o', capsize=5)
plt.text(2.0-0.2, fraction-error-5.5, str("%.2f" % fraction), fontsize=10)

fraction, error = PrintQGfraction(3)
print fraction, error
plt.errorbar(3,fraction,error,color='g', fmt='o', capsize=5)
plt.text(3.0-0.2, fraction-error-5.5, str("%.2f" % fraction), fontsize=10)

fraction, error = PrintQGfraction(4)
print fraction, error
plt.errorbar(4,fraction,error,color='r', fmt='o', capsize=5)
plt.text(4.0-0.2, fraction-error-5.5, str("%.2f" % fraction), fontsize=10)


#axes.plot([1,2,3,4],[51.51,57.65,20.61,31.33],'o',['C0','C1','g','r'])
plt.axis([0,5,0,70])
plt.title('QG percentage, mass>10^11', fontdict = {'fontsize' : 16})
plt.ylabel('%', fontdict = {'fontsize' : 14})
plt.show()
'''
'''
fig, axes = plt.subplots()

plt.setp(axes, xticks=[y_axes+1 for y_axes in range(5)], xticklabels=['all AGN', 'radio AGN', 'IR AGN', 'X-ray AGN'])

fraction, error = PrintAGNfraction(0)
print fraction, error
plt.errorbar(1,fraction,error,color='C0', fmt='o', capsize=5)
plt.text(1.0-0.2, fraction+error+0.5, str("%.2f" % fraction), fontsize=10)

fraction, error = PrintAGNfraction(2)
print fraction, error
plt.errorbar(2,fraction,error,color='C1', fmt='o', capsize=5)
plt.text(2.0-0.2, fraction+error+0.5, str("%.2f" % fraction), fontsize=10)

fraction, error = PrintAGNfraction(3)
print fraction, error
plt.errorbar(3,fraction,error,color='g', fmt='o', capsize=5)
plt.text(3.0-0.2, fraction+error+0.5, str("%.2f" % fraction), fontsize=10)

fraction, error = PrintAGNfraction(4)
print fraction, error
plt.errorbar(4,fraction,error,color='r', fmt='o', capsize=5)
plt.text(4.0-0.2, fraction+error+0.5, str("%.2f" % fraction), fontsize=10)


#axes.plot([1,2,3,4],[51.51,57.65,20.61,31.33],'o',['C0','C1','g','r'])
plt.axis([0,5,0,30])
plt.title('AGN fraction for QGs, mass>10^11', fontdict = {'fontsize' : 16})
plt.ylabel('%', fontdict = {'fontsize' : 14})
plt.show()
'''


'''
fig, axes = plt.subplots()

plt.setp(axes, xticks=[y_axes+1 for y_axes in range(5)], xticklabels=['all AGN', 'radio AGN', 'IR AGN', 'X-ray AGN'])

fraction, error = PrintAGNfraction(0)
print fraction, error
plt.errorbar(1,fraction,error,color='C0', fmt='o', capsize=5)
plt.text(1.0-0.2, fraction+error+0.5, str("%.2f" % fraction), fontsize=10)

fraction, error = PrintAGNfraction(2)
print fraction, error
plt.errorbar(2,fraction,error,color='C1', fmt='o', capsize=5)
plt.text(2.0-0.2, fraction+error+0.5, str("%.2f" % fraction), fontsize=10)

fraction, error = PrintAGNfraction(3)
print fraction, error
plt.errorbar(3,fraction,error,color='g', fmt='o', capsize=5)
plt.text(3.0-0.2, fraction+error+0.5, str("%.2f" % fraction), fontsize=10)

fraction, error = PrintAGNfraction(4)
print fraction, error
plt.errorbar(4,fraction,error,color='r', fmt='o', capsize=5)
plt.text(4.0-0.2, fraction+error+0.5, str("%.2f" % fraction), fontsize=10)


#axes.plot([1,2,3,4],[51.51,57.65,20.61,31.33],'o',['C0','C1','g','r'])
plt.axis([0,5,0,30])
plt.title('AGN fraction for SFGs, mass>10^11', fontdict = {'fontsize' : 16})
plt.ylabel('%', fontdict = {'fontsize' : 14})
plt.show()
'''
time2 = time.time()
print 'done! time =', time2-time1 , 'sec'