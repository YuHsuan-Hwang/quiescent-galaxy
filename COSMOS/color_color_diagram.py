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

def ReadCatalog( index, catalog ):
    
    data[index] = Table.read( catalog, hdu=1 )
    
    return

def SetupData( inputmask ):
    
    x.append( data[0][color2] - data[0][color3] )
    y.append( data[0][color1] - data[0][color2] )
    mask.append( inputmask )
    x_masked.append( x[0][inputmask] )
    y_masked.append( y[0][inputmask] )
    
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
        return (data[index][V_error]<up) & (data[index][V_error]>0)
    
    if mode==3:
        return (data[index][Ks_mag]<24)
    

def Mask_photoz( index ):
    return (data[index][photoz]>0)

def Mask_class_star( index ):
    return ( data[index][class_star]==0 ) | ( data[index][class_star]==2 )

def MaskAll( index ):
    return Mask_M(index) & Mask_error( 1, 0.1, index ) & Mask_photoz( index ) & Mask_class_star( index )


def Mask_classQG(    index ): return ( data[index][class_SFG]==0 )
def Mask_classSFG(   index ): return ( data[index][class_SFG]==1 )
def Mask_myclassQG(  index ): return ( y[index]> 3.1 ) & ( y[index]> 3.0*x[index]+1.0 )
def Mask_myclassSFG( index ): return ( y[index]<=3.1 ) | ( y[index]<=3.0*x[index]+1.0 )

def Mask_nonSMG( index ): return ( data[index]['850WIDE'] ==0 ) | ( data[index]['850WIDE']     ==7 )
def Mask_nonIRB( index ): return ( data[index]['24MICRON']==0 ) & ( data[index]['Radio_excess']!=0 ) #radio excess = -99 or 1

def Mask_masscut(index):
    return (data[index][mass]>11)

#def Mask_mass(index,low,high):
#    return (data[index][mass]<high) & (data[index][mass]>low)

def PlotHist_photoz( index, inputcolor, inputbins, labelname ):
    NBINS = inputbins
    plt.hist( data[0][photoz][mask[index]], NBINS,
             color=inputcolor, alpha=0.8, label = labelname) 
    plt.title( 'Histogram photoz' )
    plt.ylabel( 'Number' )
    plt.xlabel( 'z'      )
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
    plt.hist(data[0][mass][mask[index]], NBINS, color=inputcolor, alpha=0.8, label = labelname)
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



def Plot( index, scale, struc, limit, line, inputcolor,
          label, labelname, fill, inputmarker ):
    
    print "plotting ",labelname,"..."
    
    if   ( struc==1    ): set_s, set_alpha = 0.5 , 0.1   #small
    elif ( struc==2    ): set_s, set_alpha = 0.01, 0.05  #very small
    elif ( struc==0.5  ): set_s, set_alpha = 1.5 , 0.2  
    elif ( struc==0    ): set_s, set_alpha = 1.5 , 1     #big
    elif ( struc==-0.5 ): set_s, set_alpha = 10  , 1     # very big
    elif ( struc==-1   ): set_s, set_alpha = 15  , 1     # very big
    else :                set_s, set_alpha = 30  , 1     #struc==-2 very big
    
    if ( label==1 ):
        if fill==0:
            plt.scatter( x_masked[index], y_masked[index], marker=inputmarker,
                         s=set_s, alpha=set_alpha, color=inputcolor, label=labelname, facecolor='none' )
        else:
            plt.scatter( x_masked[index], y_masked[index], marker=inputmarker,
                         s=set_s, alpha=set_alpha, color=inputcolor, label=labelname )
    else:
        if fill==0:
            plt.scatter( x_masked[index], y_masked[index], marker=inputmarker,
                         s=set_s, alpha=set_alpha, color=inputcolor, facecolor='none' )
        else:
            plt.scatter( x_masked[index], y_masked[index], marker=inputmarker,
                         s=set_s, alpha=set_alpha, color=inputcolor )        
    
    #plt.title(colorname1+colorname2+colorname3, fontdict = {'fontsize' : 18})#16})
    plt.xlabel( set_xlable, fontdict = {'fontsize' : 14} )#14 20})
    plt.ylabel( set_ylable, fontdict = {'fontsize' : 14} )#14 20})
    
    if ( limit==1 ): plt.axis( [-1.5,2,-1,7] )  #plt.axis([-0.5,2,-1,7])
    if ( scale==1 ): plt.axis( 'scaled' )
    if ( line ==1 ): DrawLine()
    if ( label==1 ): plt.legend()
    
    plt.text( -1.4, 3.3, 'Quiescent Galaxies',    fontsize=10 )
    plt.text( -1.4, 2.7, 'Star-Forming Galaxies', fontsize=10 )
    plt.show()
    
    maskQG    = mask[index] & Mask_myclassQG( index )
    maskSFG   = mask[index] & Mask_myclassSFG(index )
    total_num = len( data[0].filled()[mask[index]]  )
    QG_num    = len( data[0].filled()[maskQG]       )
    SFG_num   = len( data[0].filled()[maskSFG]      )
    
    print "total sample size ", total_num
    print "QG    sample size ", QG_num
    print "QG    percentage  ", ( "%.2f" % (QG_num/total_num*100)  ),"%"
    print "SFG   sample size ", SFG_num
    print "SFG   percentage  ", ( "%.2f" % (SFG_num/total_num*100) ),"%"
    print
    
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
    


def PlotColorColor():
    
    SetupData( MaskAll(0)                            )
    SetupData( MaskAll(0) & Mask_myclassQG(0)        )
    SetupData( MaskAll(0) & (data[0]['24MICRON']==1) )
    SetupData( MaskAll(0) & (data[0]['3GHZ']    ==1) )
    
    fig = plt.figure(1)
    Plot( 0,0,1,1,1,'C0',0,'COSMOS2015',1,'o' )
    patch = mpatches.Patch( color='C0', label='COSMOS2015' )
    plt.legend( handles=[patch] )
    #fill_x = [-2,(3.1-1.0)/3.0,2.0,-2]
    #fill_y = [3.1,3.1,7,7]
    #plt.fill(fill_x,fill_y,'k',alpha=0.3)
    #fig.savefig('NUVrJ_1.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    plt.figure(2)
    PlotHist_photoz( 0,'C0',20,'COSMOS2015'    )
    #fig.savefig('hist.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    plt.figure(3)
    PlotHist_photoz( 1,'C0',20,'COSMOS2015,QG' )
    #fig.savefig('hist_QG.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    plt.figure(4)
    Plot( 2,0,0.5,1,1,'C1',0,'24 micron',1,'o' )
    patch = mpatches.Patch(color='C1', label='24 micron')
    plt.legend(handles=[patch])
    #fig.savefig('NUVrJ_24.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    plt.figure(5)
    Plot( 3,0,0.5,1,1,'r',0,'3 GHz',1,'o' )
    patch = mpatches.Patch(color='r', label='3 GHz')
    plt.legend(handles=[patch])
    #fig.savefig('NUVrJ_3.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    plt.show()
    return

def PlotColorColor_Submm850():
    
    mask_ALMAcpt = ( (data[0]['850WIDE']==1)|(data[0]['850WIDE']==2)|(data[0]['850WIDE']==3) ) & (data[0]['LENSING']==0)
    mask_IRcpt   =   (data[0]['850WIDE']==4)|(data[0]['850WIDE']==5)|(data[0]['850WIDE']==6)
    
    SetupData( MaskAll(0)                           )
    SetupData( MaskAll(0) &  mask_IRcpt             )
    SetupData( MaskAll(0) &  mask_ALMAcpt           )
    SetupData( MaskAll(0) & (data[0]['LENSING']==1) )
    
    plt.figure(1)  
    
    plt.scatter(-1000,-1000,s=1.5,color='C0',alpha=1,label='COSMOS2015')
    
    #Plot_zbin(index, scale, struc, limit, line)  
    Plot( 1,0,-0.5,1,1, 'C4', 1,'24um or 3GHz counterpart', 0,'o' )
    Plot( 2,0,-0.5,1,1, 'C4', 1,'ALMA counterpart',         1,'o' )
    Plot( 3,0,-2,  1,1, 'C4', 1,'lensing case',             1,'*' )
    Plot( 0,0, 2,  1,1, 'C0', 0,'COSMOS2015',               1,'o' )
    
    #fill_x = [-2,(3.1-1.0)/3.0,2.0,-2]
    #fill_y = [3.1,3.1,7,7]
    #plt.fill(fill_x,fill_y,'k',alpha=0.3)  

    #fig.savefig('NUVrJ_850.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    return

def PlotColorColor_Submm450():
    
    mask_ALMAcpt = ( (data[0]['450NARROW']==1)|(data[0]['450NARROW']==2)|(data[0]['450NARROW']==3) ) & (data[0]['LENSING']==0)
    mask_IRcpt   =   (data[0]['450NARROW']==4)|(data[0]['450NARROW']==5)|(data[0]['450NARROW']==6)
    
    SetupData( MaskAll(0)                           )
    SetupData( MaskAll(0) &  mask_IRcpt             )
    SetupData( MaskAll(0) &  mask_ALMAcpt           )
    SetupData( MaskAll(0) & (data[0]['LENSING']==1) )
    
    plt.figure(1)  
    
    plt.scatter(-1000,-1000,s=1.5,color='C0',alpha=1,label='COSMOS2015')
    
    #Plot_zbin(index, scale, struc, limit, line)  
    Plot( 1,0,-0.5,1,1, 'C6', 1,'24um or 3GHz counterpart', 0,'o' )
    Plot( 2,0,-0.5,1,1, 'C6', 1,'ALMA counterpart',         1,'o' )
    Plot( 3,0,-2,  1,1, 'C6', 1,'lensing case',             1,'*' )
    Plot( 0,0, 2,  1,1, 'C0', 0,'COSMOS2015',               1,'o' )
    
    #fill_x = [-2,(3.1-1.0)/3.0,2.0,-2]
    #fill_y = [3.1,3.1,7,7]
    #plt.fill(fill_x,fill_y,'k',alpha=0.3)  

    #fig.savefig('NUVrJ_450.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    return


def PlotColorColor_radioAGN():
    
    SetupData( MaskAll(0) & (data[0]['Radio_excess']==1) & (Mask_masscut(0)) )   
    plt.figure(1)    
    Plot(0,0,0,1,1,'r',0,'radio AGN',1,'o')   
    #fig.savefig('NUVrJradioAGN.png', format='png', dpi=1200)
    
    return

def PlotColorColor_IRAGN():

    SetupData( MaskAll(0) & (data[0]['agn_c17b']==True) & (Mask_masscut(0)) )   
    plt.figure(1)
    Plot(0,0,0,1,1,'C1',0,'IR AGN',1,'o')
    #fig.savefig('NUVrJradioAGN.png', format='png', dpi=1200)
    
    return


def PlotColorColor_XrayAGN():

    SetupData( MaskAll(0) & (data[0]['agn_xxx']==True) & (Mask_masscut(0)) )    
    plt.figure(1)
    Plot(0,0,0,1,1,'g',0,'X-ray AGN',1,'o')
    #fig.savefig('NUVrJradioAGN.png', format='png', dpi=1200)
    
    return
    

def Plot_nonIRQG():

    SetupData( MaskAll(0) & Mask_nonSMG(0) & Mask_nonIRB(0)                     )
    SetupData( MaskAll(0) & Mask_nonSMG(0) & Mask_nonIRB(0) & Mask_myclassQG(0) )
    
    plt.figure(1)
    Plot( 0,0,1,1,1,'C0',0,'COSMOS2015',1,'o' )

    patch1 = mpatches.Patch( color='C0', label='COSMOS2015, nonIR-QG' )
    plt.legend(handles=[patch1])
    #fig.savefig('NUVrJ_nonIRQG.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    plt.figure(2)
    PlotHist_photoz( 1,'C0',20,'COSMOS2015, nonIR-QG' )
    #fig.savefig('hist_nonIRQG.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    plt.figure(3)
    PlotHist_mass( 1,'C0',20,'COSMOS2015, nonIR-QG' )
    #fig.savefig('masshist_nonIRQG.png', bbox_inches = 'tight', format='png', dpi=1200)
    
    return

def Plot_QGfraction():
    
    SetupData( MaskAll(0) & Mask_masscut(0) )
    SetupData( MaskAll(0) & Mask_masscut(0) & Mask_myclassQG(0) )
    SetupData( MaskAll(0) & Mask_masscut(0) &\
              ((data[0]['Radio_excess']==1)|(data[0]['agn_c17b']==True)|(data[0]['agn_xxx']==True)) )
    SetupData( MaskAll(0) & Mask_masscut(0) & (data[0]['Radio_excess']==1   ) )
    SetupData( MaskAll(0) & Mask_masscut(0) & (data[0]['agn_c17b']    ==True) )
    SetupData( MaskAll(0) & Mask_masscut(0) & (data[0]['agn_xxx']     ==True) )
    
    label_list = ['COSMOS2015','all AGN', 'radio AGN', 'IR AGN', 'X-ray AGN']
    color_list = ['C0','k','r','C1','g']
    index_list = [0,2,3,4,5]
    
    fig, axes = plt.subplots()

    plt.setp( axes, xticks=[y_axes+1 for y_axes in range(5)], xticklabels=label_list )
    
    for i in range(5):
        print "plotting ", label_list[i], "..."
        fraction, error = PrintQGfraction( index_list[i] )
        print "QG fraction ", fraction
        print "error ", error[0],error[1]
        print
        plt.errorbar( i+1, fraction, error, color=color_list[i], fmt='o', capsize=5 )
        plt.text( i+0.8, fraction-error[0][0]-5.5, str("%.2f" % fraction), fontsize=10 )
        if ( i==0 ): plt.hlines( fraction,0,6, colors=color_list[i], linestyles='dotted' )

    plt.axis([0,6,0,70])
    #plt.title('QG percentage, mass>10^11', fontdict = {'fontsize' : 16})
    plt.ylabel('QG fraction(%)', fontdict = {'fontsize' : 18})#14})
    plt.show()
    
    #fig.savefig('QGfraction.png', bbox_inches = 'tight', format='png', dpi=1200)
 


    fig, axes = plt.subplots(figsize=(8,6))
    
    for i in range(5):
        print "plotting ", label_list[i], "..."
        fraction, error = Print_zbin_QGfraction( index_list[i] )
        print
        a    = i*0.01
        zbin = [0.25+a,0.75+a,1.25+a,1.75+a,2.25+a,2.75+a,3.25+a,3.75+a]
        plt.errorbar( zbin, fraction, error, fmt='o-',
                      color=color_list[i], label=label_list[i], capsize=5 )

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
color  = [ color1, color2, color3 ]

photoz     = "REDSHIFT"
V_error    = "V_MAGERR_APER3"
ip_error   = "ip_MAGERR_APER3"
J_error    = "J_MAGERR_APER3"
Ks_error   = "Ks_MAGERR_APER3"
Ks_mag     = "Ks_MAG_APER3"
mass       = "MASS_MED"
class_star = "TYPE"
class_SFG  = "CLASS"
ra         = "ALPHA_J2000"
dec        = "DELTA_J2000"
mass       = "MASS_BEST"

set_xlable = '$M_{'+colorname2+'}-M_{'+colorname3+'}$'
set_ylable = '$M_{'+colorname1+'}-M_{'+colorname2+'}$' 

###read catalog
catalog    = [None]*1
catalog[0] = 'COSMOS2015_merged.fits' 
data       = [None]*1
x,y,mask,x_masked,y_masked = [],[],[],[],[]

ReadCatalog( 0,catalog[0] )


# ===== color-color diagram =====

#PlotColorColor()

#PlotColorColor_Submm850()
#PlotColorColor_Submm450()

#PlotColorColor_radioAGN()
#PlotColorColor_IRAGN()
#PlotColorColor_XrayAGN()


# ===== nonIR-QG =====

#Plot_nonIRQG()


# ===== AGN QG correlation =====

#Plot_QGfraction()


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

# ===== test =====
'''
number = 10
x = [None]*number
y = [None]*number
mask = []
x_masked = []
y_masked = []

test()
'''


time2 = time.time()
print 'done! time =', time2-time1 , 'sec'