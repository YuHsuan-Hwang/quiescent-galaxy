#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 15:48:33 2019

@author: yuhsuan
"""
# =============================================================================
# packages
# =============================================================================
from __future__ import division
import time
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from astropy.coordinates import SkyCoord
from astropy import units as u
# =============================================================================
# functions
# =============================================================================
def ReadCatalog(index,catalog):
    data[index] = Table.read(catalog, hdu=1)
    return
def ReadXY(index,data_index):
    x[index] = data[data_index][color2]-data[data_index][color3]
    y[index] = data[data_index][color1]-data[data_index][color2]
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
    return (data[index][photoz]>0) ###& (data[index][photoz]<8)

def Mask_class_star(index):
    return (data[index][class_star]==0) | (data[index][class_star]==2)###

#def Mask_classQG(index):
#    return (data[index][class_SFG]==0)
#def Mask_classSFG(index):
#    return (data[index][class_SFG]==1)
def Mask_myclassQG(index):
    mask_nondustySFG = ((data[0]['850WIDE']==0)|(data[0]['850WIDE']==4)|(data[0]['LENSING']==1))\
            &(data[0]['24MICRON']==0)&(data[0]['Radio_excess']!=0)
    #mask_nonSMG = ((data[0]['850WIDE']==0)|(data[0]['850WIDE']==4)|(data[0]['LENSING']==1))
    return (y[index]>3.1) & (y[index]>3.0*x[index]+1.0) #& mask_nondustySFG
def Mask_myclassSFG(index):
    return (y[index]<=3.1) | (y[index]<=3.0*x[index]+1.0)

def Mask_masscut(index):
    return (data[index][mass]>11)

#def Mask_mass(index,low,high):
#    return (data[index][mass]<high) & (data[index][mass]>low)

def PlotHist_photoz(index,inputcolor,inputbins,labelname):
    NBINS = inputbins
    plt.hist(data[0][photoz][mask[index]], NBINS, color=inputcolor, alpha=0.8, label = labelname)
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
        plotdata.append(data[0][photoz][mask[indexlist[i]]])
    plt.hist(plotdata, NBINS, color=inputcolorlist, alpha=0.7,label=inputlabel)
    #plt.title('COSMOS2015')
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

def Num_total(index): 
    return len( data[0].filled()[mask[index]] )

def Plot(index, scale, struc, limit, line, inputcolor, label, labelname,fill):
    if (struc==1): #small
        set_s = 0.5
        set_alpha = 0.1
    elif (struc==2): #very small
        set_s = 0.01
        set_alpha = 0.05
    elif (struc==0.5):
        set_s = 1.5
        set_alpha = 0.2
    elif (struc==0): #big
        set_s = 1.5
        set_alpha = 1
    else: #struc==-1 very big
        set_s = 15
        set_alpha = 1
    
    if fill==0:
        plt.scatter( x_masked[index], y_masked[index], s=set_s, alpha=set_alpha,color=inputcolor, label=labelname, facecolor='none')
    else:
        plt.scatter( x_masked[index], y_masked[index], s=set_s, alpha=set_alpha,color=inputcolor, label=labelname)
    maskQG = mask[index] & Mask_myclassQG(index)
    maskSFG = mask[index] & Mask_myclassSFG(index)
    #plt.title(colorname1+colorname2+colorname3, fontdict = {'fontsize' : 18})#16})
    '''
    plt.title(colorname1+colorname2+colorname3+'\n'\
              +str(len(data[0].filled()[maskQG]))+' QGs ('+\
              str( "%.2f" % (len(data[0].filled()[maskQG])*100/Num_total(index)) )+'%), '\
              +str(len(data[0].filled()[maskSFG]))+' SFGs('+
              str( "%.2f" % (len(data[0].filled()[maskSFG])*100/Num_total(index)) )+'%)')
    '''
    plt.xlabel(set_xlable, fontdict = {'fontsize' : 14})#14 20})
    plt.ylabel(set_ylable, fontdict = {'fontsize' : 14})#14 20})
    if (limit==1):
        plt.axis([-1.5,2,-1,7])
        #plt.axis([-0.5,2,-1,7])
    if (scale==1):
        plt.axis('scaled')
    if (line==1):
        DrawLine()
    
    plt.legend()

    plt.text(-1.4, 3.3, 'Quiescent Galaxies', fontsize=10)
    plt.text(-1.4, 2.7, 'Star-Forming Galaxies', fontsize=10)
    plt.show()
    return

def PrintQGfraction(index):
    maskQG = mask[index] & Mask_myclassQG(index)
    all_num = Num_total(index)
    QG_num = len(data[0].filled()[maskQG])
    print all_num, QG_num
    fraction = QG_num*100.0 / all_num
    error =  np.sqrt(QG_num)*100.0 / all_num
    if (fraction+error>100):
        up_err = 100.0-fraction
    else:
        up_err = error
    if (fraction-error<0.0):
        low_err = fraction
    else:
        low_err = error
    return fraction, np.array([[low_err,up_err]]).T


def PrintAGNfraction(index):
    maskQG = mask[1]
    maskAGN = mask[index] & Mask_myclassQG(index)
    AGN_num = len(data[0].filled()[maskAGN])
    QG_all_num = len(data[0].filled()[maskQG])
    print AGN_num, QG_all_num
    fraction = AGN_num*100.0 / QG_all_num
    error =  np.sqrt(AGN_num)*100.0 / QG_all_num
    return fraction, error

def Construct_mask_zbin(index):
    mask_zbin = []
    mask_zbin.append((data[0][photoz]>0) & (data[0][photoz]<0.5))
    mask_zbin.append((data[0][photoz]>0.5) & (data[0][photoz]<1))
    mask_zbin.append((data[0][photoz]>1) & (data[0][photoz]<1.5))
    mask_zbin.append((data[0][photoz]>1.5) & (data[0][photoz]<2))
    mask_zbin.append((data[0][photoz]>2) & (data[0][photoz]<2.5))
    mask_zbin.append((data[0][photoz]>2.5) & (data[0][photoz]<3))
    mask_zbin.append((data[0][photoz]>3) & (data[0][photoz]<3.5))
    mask_zbin.append((data[0][photoz]>3.5))
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
                      str( "%.2f" % ((len(x[index][datamaskQG].filled()))*100.0\
                                     /len(x[index][datamask].filled())) )+'%')
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
            error.append(np.array([0.0,0.0]))
        else:
            fraction_tmp = QG_num*100.0/all_num
            fraction.append(fraction_tmp)
            error_tmp = np.sqrt(QG_num)*100.0/all_num
            if (fraction_tmp+error_tmp>100.0):
                up_error = 100.0-fraction_tmp
            else:
                up_error = error_tmp
            if (fraction_tmp-error_tmp<0.0):
                low_error = fraction_tmp
            else:
                low_error = error_tmp
            error.append(np.array([low_error,up_error]))
            print  all_num, QG_num, up_error, low_error
    return fraction, np.array(error).T
    
def PlotMatchedResult(matchedband):
    
    #twomatch:1, onematch24/3:2/3, nomatch:4 rgcb
    mask_all =  Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0) #& Mask_myclassQG(0)
    if matchedband=='850ALL':
        mask.append( mask_all& ((data[0]['850NARROW']==1) | (data[0]['850SOURCE']==1)) )
        mask.append( mask_all& ((data[0]['850NARROW']==2) | (data[0]['850SOURCE']==2)) )
        mask.append( mask_all& ((data[0]['850NARROW']==3) | (data[0]['850SOURCE']==3)) )
        mask.append( mask_all& ((data[0]['850NARROW']==4) | (data[0]['850SOURCE']==4)) )
    else:
        mask.append( mask_all& (data[0][matchedband]==1) )
        mask.append( mask_all& (data[0][matchedband]==2) )
        mask.append( mask_all& (data[0][matchedband]==3) )
        mask.append( mask_all& (data[0][matchedband]==4) )
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
    Plot(0,0,0,1,1,'r',1,'both matched, '+str(Num_total(0))+' samples')
    plt.subplot(2, 2, 4)
    Plot(1,0,0,1,1,'g',1,'24 micron matched, '+str(Num_total(1))+' samples')
    plt.subplot(2, 2, 1)
    Plot(2,0,0,1,1,'c',1,'3GHz matched, '+str(Num_total(2))+' samples')
    plt.subplot(2, 2, 3)
    Plot(3,0,0,1,1,'b',1,'none matched, '+str(Num_total(3))+' samples')
    plt.tight_layout()
    return

def PlotColorColor():
    for i in range(number):
        ReadXY(i,0)
    mask_all =  Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0)# & Mask_myclassQG(0)
    mask.append(mask_all)
    mask.append(mask_all&Mask_myclassQG(0))
    mask.append( (data[0]['24MICRON']==1) & mask_all )
    mask.append( (data[0]['3GHZ']==1) & mask_all)
    
    for i in range(number):
        x_masked.append(x[i][mask[i]])
        y_masked.append(y[i][mask[i]])
    
    fig = plt.figure(1)
    '''
    labellist = ['COSMOS2015','24 micron',"3 GHz"]
    PlotHist_photoz_para([0,1,2],['C0','C1','r'],labellist,20)
    fig.savefig('hist.png', format='png', dpi=1200)
    
    plt.figure(2)#Plot_zbin(index, scale, struc, limit, line)
    '''
    Plot(0,0,1,1,1,'C0',0,'COSMOS2015',1)

    patch1 = mpatches.Patch(color='C0', label='COSMOS2015')
    #patch2 = mpatches.Patch(color='C1', label='24 micron')
    #patch3 = mpatches.Patch(color='r', label='3GHz')
    #plt.legend(handles=[patch1,patch2,patch3])
    plt.legend(handles=[patch1])
    
    
    #fill_x = [-2,(3.1-1.0)/3.0,2.0,-2]
    #fill_y = [3.1,3.1,7,7]
    #plt.fill(fill_x,fill_y,'k',alpha=0.3)
    
    fig.savefig('NUVrJ_1.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    fig = plt.figure(2)
    PlotHist_photoz(0,'C0',20,'COSMOS2015,QG')
    fig.savefig('hist.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    fig = plt.figure(3)
    PlotHist_photoz(1,'C0',20,'COSMOS2015,QG')
    fig.savefig('hist_QG.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    fig = plt.figure(4)
    Plot(2,0,0.5,1,1,'C1',0,'24 micron',1)
    patch = mpatches.Patch(color='C1', label='24 micron')
    plt.legend(handles=[patch])
    fig.savefig('NUVrJ_24.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    fig = plt.figure(5)
    Plot(3,0,0.5,1,1,'r',0,'3 GHz',1)
    patch = mpatches.Patch(color='r', label='3 GHz')
    plt.legend(handles=[patch])
    fig.savefig('NUVrJ_3.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    return

def PlotColorColor_Submm850():
    for i in range(number):
        ReadXY(i,0)
    mask_all =  Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0) 
    mask.append(mask_all)
    mask.append( ((data[0]['850WIDE']==2)|(data[0]['850WIDE']==3)) & mask_all )
    mask.append( (data[0]['850WIDE']==5) & mask_all )
    #mask.append( (data[0]['850WIDE']==5)&(data[0]['LENSING']==0) & mask_all )
    #mask.append( (data[0]['LENSING']==0) & mask_all )
    #mask.append( ((data[0]['450NARROW']==1)|(data[0]['450NARROW']==2)|(data[0]['450NARROW']==3)) & mask_all)
    
    for i in range(number):
        x_masked.append(x[i][mask[i]])
        y_masked.append(y[i][mask[i]])
    
    fig = plt.figure(1)
    
    #Plot_zbin(index, scale, struc, limit, line)
    
    Plot(0,0,2,1,1,'C0',0,'COSMOS2015',1)
    Plot(1,0,-1,1,1,'C4',1,'24 micron or 3 GHz counterpart',0)
    Plot(2,0,-1,1,1,'C4',1,'ALMA 870 micron counterpart',1)
    #Plot(2,0,0,1,1,'C4',1,'ALMA 870 micron counterpart',1)
    '''
    patch1 = mpatches.Patch(color='C0', label='COSMOS2015')
    patch2 = mpatches.Patch(color='C4', label='850 micron with 24 micron or 3 GHz')
    patch3 = mpatches.Patch(color='C4', label='850 micron with ALMA')
    
    plt.legend(handles=[patch1,patch2,patch3])
    '''
    
    #fill_x = [-2,(3.1-1.0)/3.0,2.0,-2]
    #fill_y = [3.1,3.1,7,7]
    #plt.fill(fill_x,fill_y,'k',alpha=0.3)
    
    fig.savefig('NUVrJ_850.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    return

def PlotColorColor_Submm450():
    for i in range(number):
        ReadXY(i,0)
    mask_all =  Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0) 
    mask.append(mask_all)
    mask.append( ((data[0]['450NARROW']==2)|(data[0]['450NARROW']==3)) & mask_all )
    mask.append( ((data[0]['450NARROW']==5)) & mask_all )
    #mask.append( ((data[0]['450NARROW']==1)|(data[0]['450NARROW']==2)|(data[0]['450NARROW']==3)) & mask_all)
    
    for i in range(number):
        x_masked.append(x[i][mask[i]])
        y_masked.append(y[i][mask[i]])
    
    fig = plt.figure(1)
    
    #Plot_zbin(index, scale, struc, limit, line)
    
    Plot(0,0,2,1,1,'C0',0,'COSMOS2015',1)
    Plot(1,0,-1,1,1,'C6',1,'24 micron or 3 GHz counterpart',0)
    Plot(2,0,-1,1,1,'C6',1,'ALMA 870 micron counterpart',1)
    '''
    patch1 = mpatches.Patch(color='C0', label='COSMOS2015')
    patch2 = mpatches.Patch(color='C4', label='850 micron with 24 micron or 3 GHz')
    patch3 = mpatches.Patch(color='C4', label='850 micron with ALMA')
    
    plt.legend(handles=[patch1,patch2,patch3])
    '''
    
    #fill_x = [-2,(3.1-1.0)/3.0,2.0,-2]
    #fill_y = [3.1,3.1,7,7]
    #plt.fill(fill_x,fill_y,'k',alpha=0.3)
    
    fig.savefig('NUVrJ_450.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    return

def PlotColorColor_radioAGN():
    for i in range(number):
        ReadXY(i,0)
    mask_all =  Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0) 
    mask.append( (data[0]['Radio_excess']==1)&(Mask_masscut(0)) & mask_all )
    
    for i in range(number):
        x_masked.append(x[i][mask[i]])
        y_masked.append(y[i][mask[i]])
    
    fig = plt.figure(1)
    
    #Plot_zbin(index, scale, struc, limit, line)
    
    Plot(0,0,0,1,1,'r',0,'radio AGN',1)

    '''
    patch1 = mpatches.Patch(color='C0', label='COSMOS2015')
    patch2 = mpatches.Patch(color='C4', label='850 micron with 24 micron or 3 GHz')
    patch3 = mpatches.Patch(color='C4', label='850 micron with ALMA')
    
    plt.legend(handles=[patch1,patch2,patch3])
    '''
    
    #fill_x = [-2,(3.1-1.0)/3.0,2.0,-2]
    #fill_y = [3.1,3.1,7,7]
    #plt.fill(fill_x,fill_y,'k',alpha=0.3)
    
    #fig.savefig('NUVrJradioAGN.png', format='png', dpi=1200)
    
    return

def PlotColorColor_IRAGN():
    for i in range(number):
        ReadXY(i,0)
    mask_all =  Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0) 
    mask.append( (data[0]['agn_c17b']==True)&(Mask_masscut(0)) & mask_all )
    
    for i in range(number):
        x_masked.append(x[i][mask[i]])
        y_masked.append(y[i][mask[i]])
    
    fig = plt.figure(1)
    
    #Plot_zbin(index, scale, struc, limit, line)
    
    Plot(0,0,0,1,1,'C1',0,'IR AGN',1)

    '''
    patch1 = mpatches.Patch(color='C0', label='COSMOS2015')
    patch2 = mpatches.Patch(color='C4', label='850 micron with 24 micron or 3 GHz')
    patch3 = mpatches.Patch(color='C4', label='850 micron with ALMA')
    
    plt.legend(handles=[patch1,patch2,patch3])
    '''
    
    #fill_x = [-2,(3.1-1.0)/3.0,2.0,-2]
    #fill_y = [3.1,3.1,7,7]
    #plt.fill(fill_x,fill_y,'k',alpha=0.3)
    
    #fig.savefig('NUVrJradioAGN.png', format='png', dpi=1200)
    
    return


def PlotColorColor_XrayAGN():
    for i in range(number):
        ReadXY(i,0)
    mask_all =  Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0) 
    mask.append( (data[0]['agn_xxx']==True)&(Mask_masscut(0)) & mask_all )
    
    for i in range(number):
        x_masked.append(x[i][mask[i]])
        y_masked.append(y[i][mask[i]])
    
    fig = plt.figure(1)
    
    #Plot_zbin(index, scale, struc, limit, line)
    
    Plot(0,0,0,1,1,'g',0,'X-ray AGN',1)

    '''
    patch1 = mpatches.Patch(color='C0', label='COSMOS2015')
    patch2 = mpatches.Patch(color='C4', label='850 micron with 24 micron or 3 GHz')
    patch3 = mpatches.Patch(color='C4', label='850 micron with ALMA')
    
    plt.legend(handles=[patch1,patch2,patch3])
    '''
    
    #fill_x = [-2,(3.1-1.0)/3.0,2.0,-2]
    #fill_y = [3.1,3.1,7,7]
    #plt.fill(fill_x,fill_y,'k',alpha=0.3)
    
    #fig.savefig('NUVrJradioAGN.png', format='png', dpi=1200)
    
    return
    

def Plot_nondustySFG():
    for i in range(number):
        ReadXY(i,0)
    mask_all =  Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0)# & Mask_myclassQG(0)
    mask_nondustySFG = ((data[0]['850WIDE']==0)|(data[0]['850WIDE']==4)|(data[0]['LENSING']==1))\
            &(data[0]['24MICRON']==0)&(data[0]['Radio_excess']!=0)
    mask.append(mask_all & mask_nondustySFG & Mask_myclassQG(0))
    
    for i in range(number):
        x_masked.append(x[i][mask[i]])
        y_masked.append(y[i][mask[i]])
    
    fig = plt.figure(1)
    Plot(0,0,1,1,1,'C0',0,'COSMOS2015',1)

    patch1 = mpatches.Patch(color='C0', label='COSMOS2015')
    plt.legend(handles=[patch1])
    fig.savefig('NUVrJ_nondustySFG.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    fig = plt.figure(2)
    PlotHist_photoz(0,'C0',20,'COSMOS2015')
    fig.savefig('hist_nondustySFG.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    fig = plt.figure(3)
    PlotHist_mass(0,'C0',20,'COSMOS2015')
    fig.savefig('masshist_nondustySFG.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    return

def PrintNum():
    for i in range(number):
        ReadXY(i,0)
    mask_all =  Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0) & Mask_myclassQG(0)
    mask_850 = (data[0]['850WIDE']==4)
    mask.append(mask_all &mask_850)
    print len(data[0][mask[0]])
    return

def Plot_QGfraction():
    for i in range(number):
        ReadXY(i,0)

    mask_all = Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0) & Mask_masscut(0)
    mask.append( mask_all )
    mask.append( mask_all & Mask_myclassQG(0))
    mask.append( mask_all &\
                ( (data[0]['Radio_excess']==1)|(data[0]['agn_c17b']==True)|(data[0]['agn_xxx']==True) ))
    mask.append( mask_all & (data[0]['Radio_excess']==1) )
    mask.append( mask_all & (data[0]['agn_c17b']==True) )
    mask.append( mask_all & (data[0]['agn_xxx']==True) )
    
    
    fig, axes = plt.subplots()

    plt.setp(axes, xticks=[y_axes+1 for y_axes in range(number)], xticklabels=['COSMOS2015','all AGN', 'radio AGN', 'IR AGN', 'X-ray AGN'])
    
    fraction, error = PrintQGfraction(0)
    print fraction, error
    plt.errorbar(1,fraction,error,color='k', fmt='o', capsize=5)
    plt.text(1.0-0.2, fraction-error[0][0]-5.5, str("%.2f" % fraction), fontsize=10)
    plt.hlines(fraction,0,6, colors='k', linestyles='dotted')
    
    fraction, error = PrintQGfraction(2)
    print fraction, error
    plt.errorbar(2,fraction,error,color='C0', fmt='o', capsize=5)
    plt.text(2.0-0.2, fraction-error[0][0]-5.5, str("%.2f" % fraction), fontsize=10)
    
    fraction, error = PrintQGfraction(3)
    print fraction, error
    plt.errorbar(3,fraction,error,color='r', fmt='o', capsize=5)
    plt.text(3.0-0.2, fraction-error[0][0]-5.5, str("%.2f" % fraction), fontsize=10)
    
    fraction, error = PrintQGfraction(4)
    print fraction, error
    plt.errorbar(4,fraction,error,color='C1', fmt='o', capsize=5)
    plt.text(4.0-0.2, fraction+error[0][0]+5.5, str("%.2f" % fraction), fontsize=10)
    
    fraction, error = PrintQGfraction(5)
    print fraction, error
    plt.errorbar(5,fraction,error,color='g', fmt='o', capsize=5)
    plt.text(5.0-0.2, fraction-error[0][0]-5.5, str("%.2f" % fraction), fontsize=10)
    
    
    #axes.plot([1,2,3,4],[51.51,57.65,20.61,31.33],'o',['C0','C1','g','r'])
    plt.axis([0,6,0,70])
    #plt.title('QG percentage, mass>10^11', fontdict = {'fontsize' : 16})
    plt.ylabel('QG fraction(%)', fontdict = {'fontsize' : 18})#14})
    plt.show()
    
    #fig.savefig('QGfraction.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    
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
    
    #fig.savefig('QGfraction_zbin.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    
def RegionFile(index,filename,mode,color,size):
    print 'Output region file...'
    f = open(filename+'.reg','w')
    f.write('global color='+color+' font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    
    N = len( data[0].filled()[mask[0]][mask[index]] )
    #N = len( data[index] )
    
    print 'source number: ',N
    
    if (mode==0): #no label
        for n in range(N):
            f.write('fk5;circle('+str(data[0][ra][mask[0]][mask[index]][n])+','+str(data[0][dec][mask[0]][mask[index]][n])+','+size+'") # text={''}\n')
            #f.write('fk5;circle('+str(data[index][ra][mask[index]][n])+','+str(data[index][dec][mask[index]][n])+','+size+'") # text={''}\n')
    elif (mode==1): # redshift label
        for n in range(N):
            f.write('fk5;circle('+str(data[0][ra][mask[0]][mask[index]][n])+','+str(data[0][dec][mask[0]][mask[index]][n])+','+size+'") # text={'+str(data[0][photoz][mask[0]][mask[index]][n])+'}\n')
    
    #f.write('fk5;circle('+str(data[index][ra][mask[index]][n])+','+str(data[index][dec][mask[index]][n])+','+size+'") # text={''}\n')
    #f.write('fk5;circle('+str(data[index][ra][mask[index]][n])+','+str(data[index][dec][mask[index]][n])+','+size+'") # text={'+str(data[index]["NUMBER"][mask[index]][n])+'}\n')
    f.close()
    return

def RegionFileUnmask(index,filename,mode,color,size,ra,dec):
    print 'Output region file...'
    f = open(filename+'.reg','w')
    f.write('global color='+color+' font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    
    N = len( data[index] )
    
    print 'source number: ',N
    
    if (mode==0): #no label
        for n in range(N):
            f.write('fk5;circle('+str(data[0][ra][n])+','+str(data[0][dec][n])+','+size+'") # text={''}\n')
            #f.write('fk5;circle('+str(data[index][ra][mask[index]][n])+','+str(data[index][dec][mask[index]][n])+','+size+'") # text={''}\n')
    elif (mode==1): # redshift label
        for n in range(N):
            f.write('fk5;circle('+str(data[0][ra][mask[index]][n])+','+str(data[0][dec][mask[index]][n])+','+size+'") # text={'+str(data[0][photoz][mask[index]][n])+'}\n')
    
    #f.write('fk5;circle('+str(data[index][ra][mask[index]][n])+','+str(data[index][dec][mask[index]][n])+','+size+'") # text={''}\n')
    #f.write('fk5;circle('+str(data[index][ra][mask[index]][n])+','+str(data[index][dec][mask[index]][n])+','+size+'") # text={'+str(data[index]["NUMBER"][mask[index]][n])+'}\n')
    f.close()
    return

def taskRegionFile():
    path = '/Users/yuhsuan/Documents/research/05WH/data/COSMOS/'
    #catalog = path+'05_3GHz/vla3_cosmos_sources_160321_public5sig.fits.txt'
    catalog = path+"07_ALMA/apjsab42da_table4/A-COSMOS_blind.fits"
    ReadCatalog(0,catalog)
    #RegionFileUnmask(0,'COSMOS_3GHz_5sig',0,'red','2.0')
    RegionFileUnmask(0,'COSMOS_A3COSMOS',0,'yellow','0.5',"RA","DEC")
    return


def MaskRegion(inputra,inputdec,inputradius):
    c0 = SkyCoord(ra=data[0][mask[0]][ra], dec=data[0][mask[0]][dec])
    center = SkyCoord(inputra, inputdec, frame='fk5')
    radius = inputradius
    sep = center.separation(c0)
    return sep<=radius


def test():
    for i in range(number):
        ReadXY(i,0)
    
    
    mask_all =  Mask_myclassSFG(0)#Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0) & Mask_myclassSFG(0)
    mask_24_3 = ((data[0]['850WIDE']==2)|(data[0]['850WIDE']==3)) & mask_all
    #mask_24_3 = ((data[0]['450NARROW']==2)|(data[0]['450NARROW']==3))& mask_all
    mask_870 = (data[0]['850WIDE']==5) & mask_all
    #mask_870 = (data[0]['450NARROW']==5) & mask_all
    
    
    mask.append(mask_all)
    mask.append( mask_24_3  )
    mask.append( mask_870 )
    
    mask.append((data[0]['Radio_excess']==1)&mask_24_3)
    mask.append((data[0]['Radio_excess']==1)&mask_870 )
    mask.append((data[0]['agn_c17b']==True)&mask_24_3)
    mask.append((data[0]['agn_c17b']==True)&mask_870 )
    mask.append((data[0]['agn_xxx']==True)&mask_24_3)
    mask.append((data[0]['agn_xxx']==True)&mask_870 )
    #mask.append( MaskRegion('10h01m47.0s', '2d3m54.0s',0.3*u.arcmin) )
    mask.append( MaskRegion('10h00m14.0s', '1d56m41.0s',1.0*u.arcmin) )
    for i in range(9):
        x_masked.append(x[i][mask[i]])
        y_masked.append(y[i][mask[i]])
    
    plt.figure(1)
    
    #Plot_zbin(index, scale, struc, limit, line)
    
    Plot(0,0,2,1,1,'C0',0,'COSMOS2015',1)
    Plot(1,0,-1,1,1,'C6',1,'24 micron or 3 GHz counterpart',0)
    Plot(2,0,-1,1,1,'C6',1,'ALMA 870 micron counterpart',1)
    Plot(3,0,-1,1,1,'r',1,'',0)
    Plot(4,0,-1,1,1,'r',1,'',1)
    Plot(5,0,-1,1,1,'C1',1,'',0)
    Plot(6,0,-1,1,1,'C1',1,'',1)
    Plot(7,0,-1,1,1,'g',1,'',0)
    Plot(8,0,-1,1,1,'g',1,'',1)
    
    print len(data[0][mask[0]])
    print len(data[0][mask[1]])
    print len(data[0][mask[2]])
    print len(data[0][mask[3]])
    print len(data[0][mask[4]])
    print len(data[0][mask[5]])
    print len(data[0][mask[6]])
    print len(data[0][mask[7]])
    print len(data[0][mask[8]])
    #print len(data[0][mask[9]])
    
    #print data[0][mask[2]][photoz]
    
    #RegionFile(3,'COSMOS_450_23_rAGN',0,'white','20.0')
    #RegionFile(4,'COSMOS_850_A_rAGN',0,'pink','10.0')
    #RegionFile(9,'COSMOS_850_23_toptop',0,'white','20.0')
    
    #RegionFile(3,'COSMOS_rAGN',0,'red','1.2')
    #RegionFile(5,'COSMOS_IRrAGN',0,'yellow','1.2')
    #RegionFile(7,'COSMOS_XAGN',0,'green','1.2')
    
    #RegionFile(1,'COSMOS_850_23_highz_top',0,'black','20.0')
    
    RegionFile(9,'COSMOS_A3multi_allSFG2',0,'orange','1.0')
    
    return

def test2():
    for i in range(number):
        ReadXY(i,0)
    mask_all =  Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0) & Mask_myclassQG(0)
    mask_QG_850 = ((data[0]['850WIDE']==0)|(data[0]['850WIDE']==4))\
            &(data[0]['24MICRON']==0)&(data[0]['Radio_excess']!=0)
    mask_QG_450 = (data[0]['450NARROW']==0)&(data[0]['24MICRON']==0)&(data[0]['Radio_excess']!=0)
    
    
    mask_IRBQG_850 = ((data[0]['850WIDE']==0)|(data[0]['850WIDE']==4))\
            &( (data[0]['24MICRON']==1)|(data[0]['Radio_excess']==0) )
    mask_IRBQG_450 = (data[0]['450NARROW']==0)&( (data[0]['24MICRON']==1)|(data[0]['Radio_excess']==0) )

    mask.append(mask_all & mask_QG_850)
    mask.append(mask_all & mask_IRBQG_850)
    mask.append(mask_all & mask_QG_450)
    mask.append(mask_all & mask_IRBQG_450)
    plt.figure(1)
    PlotHist_photoz_para([0,2],['C4','C6'],['850','450'],20)
    plt.figure(2)
    PlotHist_photoz_para([1,3],['C4','C6'],['850','450'],20)
    
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

###read catalog
catalog = [None]*1
catalog[0] = 'COSMOS2015_merged.fits' 
data = [None]*1
ReadCatalog(0,catalog[0])



# ===== color-color diagram =====
'''
number = 4
x = [None]*number
y = [None]*number
mask = []
x_masked = []
y_masked = []

PlotColorColor()
'''

# ===== color-color diagram =====
'''
number = 3
x = [None]*number
y = [None]*number
mask = []
x_masked = []
y_masked = []

#PlotColorColor_Submm850()
#PlotColorColor_Submm450()
'''

# ===== color-color diagram =====
'''
number = 1
x = [None]*number
y = [None]*number
mask = []
x_masked = []
y_masked = []

#PlotColorColor_radioAGN()
#PlotColorColor_IRAGN()
PlotColorColor_XrayAGN()
'''

# ===== nondustySFG =====
'''
number = 1
x = [None]*number
y = [None]*number
mask = []
x_masked = []
y_masked = []

Plot_nondustySFG()
'''

# ===== AGN QG correlation =====
'''
number = 6
x = [None]*number
y = [None]*number
mask = []
x_masked = []
y_masked = []

Plot_QGfraction()
'''

# ===== RegionFile =====
'''
number = 1
x = [None]*number
y = [None]*number
mask = []
x_masked = []
y_masked = []

taskRegionFile()
'''

# ===== PrintNum =====
'''
number = 1
x = [None]*number
y = [None]*number
mask = []
x_masked = []
y_masked = []

PrintNum()
'''
# ===== test =====

number = 10
x = [None]*number
y = [None]*number
mask = []
x_masked = []
y_masked = []

test()



time2 = time.time()
print 'done! time =', time2-time1 , 'sec'