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
        return (data[index][Ks_error]<up) & (data[index][Ks_error]>0)
    if mode==3:
        return (data[index][Ks_mag]<24)
    

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
    plt.title('Histogram photoz')
    plt.ylabel('Number')
    plt.xlabel('z')
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

def DrawLine():
    plt.plot([-1.5,(3.1-1.0)/3.0],[3.1,3.1],'k-',lw=0.7)
    plt.plot([(3.1-1.0)/3.0,2.0],[3.1,3.0*2.0+1.0],'k-',lw=0.7)
    return

def MultiplePlot(index, scale, struc, limit, line, inputcolor, labelname):
    if (struc==1):
        set_s = 0.5
        set_alpha = 0.1
    else:
        set_s = 1
        set_alpha = 1
    plt.scatter( x_masked[index], y_masked[index], s=set_s, alpha=set_alpha,color=inputcolor, label=labelname)
    plt.title(colorname1+colorname2+colorname3)
    plt.xlabel(set_xlable)
    plt.ylabel(set_ylable)
    if (limit==1):
        plt.axis([-1.5,2,-1,7])
    if (scale==1):
        plt.axis('scaled')
    if (line==1):
        DrawLine()
    return

def Plot(index, scale, struc, limit, line, inputcolor, label, labelname):
    if (struc==1):
        set_s = 0.5
        set_alpha = 0.1
    else:
        set_s = 1
        set_alpha = 1
    plt.scatter( x_masked[index], y_masked[index], s=set_s, alpha=set_alpha,color=inputcolor, label=labelname)
    plt.title(colorname1+colorname2+colorname3)
    plt.xlabel(set_xlable)
    plt.ylabel(set_ylable)
    if (limit==1):
        plt.axis([-1.5,2,-1,7])
    if (scale==1):
        plt.axis('scaled')
    if (line==1):
        DrawLine()
    if (line==1):
        plt.legend()
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
    f = open('/Users/yuhsuan/Documents/research/05WH/data/'+filename+'.reg','w')
    f.write('global color='+color+' font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    N = len( data[index].filled()[mask[index]] )
    for n in range(N):
        f.write('fk5;circle('+str(data[index][ra][mask[index]][n])+','+str(data[index][dec][mask[index]][n])+','+size+'") # text={'+str(n)+'}\n')
    f.close()

def PlotMatchedResult(matchedband):
    #twomatch:1, onematch24/3:2/3, nomatch:4 rgcb
    mask.append((data[0][matchedband]==1) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#& Mask_myclassQG(0))
    mask.append((data[0][matchedband]==2) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#& Mask_myclassQG(0))
    mask.append((data[0][matchedband]==3) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#& Mask_myclassQG(0))
    mask.append((data[0][matchedband]==4) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#& Mask_myclassQG(0))
    
    for i in range(4):
        x_masked.append(x[i][mask[i]])
        y_masked.append(y[i][mask[i]])
    
    
    
    plt.figure(1)
    labellist = ['both','24 micron','3 GHz','none']
    PlotHist_photoz_para([0,1,2,3],['r','g','c','b'],labellist,10)
    
    plt.figure(2)
    MultiplePlot(0,0,0,1,1,'r','both matched')
    MultiplePlot(1,0,0,1,1,'g','24 micron matched')
    MultiplePlot(2,0,0,1,1,'c','3GHz matched')
    MultiplePlot(3,0,0,1,1,'b','none matched')
    #plt.legend()
    plt.show()
    
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
# =============================================================================
# main code
# =============================================================================

time1 = time.time()

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
Ks_mag = "Ks_MAG_APER3"
mass = "MASS_MED"
class_star = "TYPE"
class_SFG = "CLASS"
ra = "ALPHA_J2000"
dec = "DELTA_J2000"

set_xlable = '$M_{'+colorname2+'}-M_{'+colorname3+'}$'
set_ylable = '$M_{'+colorname1+'}-M_{'+colorname2+'}$'  

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

#for i in range(number):
#    ReadCatalog(i,catalog[i])

for i in range(number):
    ReadCatalog(i,catalog[0])

###mask data
mask = []
x_masked = []
y_masked = []

#for i in range(number):
    #mask.append(data[0]['850SOURCE']==4)
    #mask.append( Mask_M(i) & Mask_photoz(i) & Mask_error(2,0.1,i) & Mask_class_star(i) )
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


###850 sources
'''
#twomatch:1, onematch24/3:2/3, nomatch:4 rgcb
mask.append((data[0]['850SOURCE']==1) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0)& Mask_myclassQG(0))
mask.append((data[0]['850SOURCE']==2) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0)& Mask_myclassQG(0))
mask.append((data[0]['850SOURCE']==3) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0)& Mask_myclassQG(0))
mask.append((data[0]['850SOURCE']==4) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0)& Mask_myclassQG(0))

for i in range(4):
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])

plt.figure(1)
PlotHist_photoz(0,'r',20)
plt.figure(2)
PlotHist_photoz(1,'g',20)
plt.figure(3)
PlotHist_photoz(2,'c',20) 
plt.figure(4)
PlotHist_photoz(3,'b',20) 

#plt.figure(5)
labellist = ['both','24 micron','3 GHz','none']
PlotHist_photoz_para([0,1,2,3],['r','g','c','b'],labellist,10)

plt.figure(6)
Plot(0,0,0,1,1,'r')
plt.figure(7)
Plot(1,0,0,1,1,'g')
plt.figure(8)
Plot(2,0,0,1,1,'c')
plt.figure(9)
Plot(3,0,0,1,1,'b')
'''

###850 narrow
'''
#twomatch:1, onematch24/3:2/3, nomatch:4 rgcb
mask.append((data[0]["850NARROW"]==1) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#& Mask_myclassQG(0))
mask.append((data[0]['850NARROW']==2) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#& Mask_myclassQG(0))
mask.append((data[0]['850NARROW']==3) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#& Mask_myclassQG(0))
mask.append((data[0]['850NARROW']==4) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0))#& Mask_myclassQG(0))

for i in range(4):
    x_masked.append(x[i][mask[i]])
    y_masked.append(y[i][mask[i]])



plt.figure(1)
labellist = ['both','24 micron','3 GHz','none']
PlotHist_photoz_para([0,1,2,3],['r','g','c','b'],labellist,10)

plt.figure(2)
Plot(0,0,0,1,1,'r')
Plot(1,0,0,1,1,'g')
Plot(2,0,0,1,1,'c')
Plot(3,0,0,1,1,'b') 

plt.figure(3)
plt.figure(figsize=(10,20))
plt.subplot(2, 2, 2)
PlotHist_photoz(0,'r',10)
plt.subplot(2, 2, 4)
PlotHist_photoz(1,'g',10)
plt.subplot(2, 2, 1)
PlotHist_photoz(2,'c',10) 
plt.subplot(2, 2, 3)
PlotHist_photoz(3,'b',10) 
plt.tight_layout()

plt.figure(4)
plt.figure(figsize=(10,20))
plt.subplot(2, 2, 2)
Plot(0,0,0,1,1,'r')
plt.subplot(2, 2, 4)
Plot(1,0,0,1,1,'g')
plt.subplot(2, 2, 1)
Plot(2,0,0,1,1,'c')
plt.subplot(2, 2, 3)
Plot(3,0,0,1,1,'b')
plt.tight_layout()
'''

###450 narrow

PlotMatchedResult('850NARROW')




###test different selection
'''
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






#PrintNum_zbin_total(0)
#PrintNum_zbin_total(1)
#PrintNum_zbin_total(2)
#PrintNum_zbin_total(3)
#PrintNum_zbin(0)



'''
print "---"
PrintNum_zbin(3)
print "---"
PrintNumdiff_zbin(0,2)
print "---"
PrintPercentage_zbin(0,2)
'''

#RegionFile(2, 'COSMOS850sources', 'red','5.0')
#RegionFile(1, 'COSMOS24sources', 'blue','4.0')
#RegionFile(0, 'COSMOS24sources', 'blue','1.0')



#PlotHist_massmed(data24,mask24)

time2 = time.time()
print 'done! time =', time2-time1 , 'sec'