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
        return (data[index][Ks_error]<up) & (data[index][Ks_error]>0)
    
    if mode==3:
        return (data[index][Ks_mag]<24)
    

def Mask_photoz( index ):
    return (data[index][photoz]>0)

def Mask_class_star( index ):
    return ( data[index][class_star]==0 ) | ( data[index][class_star]==2 )

def MaskAll( index ):
    return Mask_M( index ) & Mask_error( 1, 0.1, index ) & Mask_photoz( index ) & Mask_class_star( index )
           # &((data[index]["FLAG_HJMCC"]==0) | (data[index]["FLAG_HJMCC"]==2)) &(data[index]["FLAG_COSMOS"]==1)
           # &(data[index]["FLAG_HJMCC"]==0) &(data[index]["FLAG_COSMOS"]==1 )& (data[index]["FLAG_PETER"]==0)
def MaskAll_oneband( index ):
    return Mask_M( index ) & Mask_error( 2, 0.2, index ) & Mask_photoz( index ) & Mask_class_star( index )

def Mask_classQG(    index ): return ( data[index][class_SFG]==0 )
def Mask_classSFG(   index ): return ( data[index][class_SFG]==1 )
def Mask_myclassQG(  index ): return ( y[index]> 3.1 ) & ( y[index]> 3.0*x[index]+1.0 )
def Mask_myclassSFG( index ): return ( y[index]<=3.1 ) | ( y[index]<=3.0*x[index]+1.0 )

def Mask_myclassQG_UVJ(  index,delta ): return ( y[index]> 1.0 ) & ( x[index]< 1.6 ) & ( y[index]> 0.88*x[index]+delta-0.3 )
def Mask_myclassSFG_UVJ(  index,delta ): return ( y[index]<= 1.0 ) | ( x[index]>= 1.6 ) | ( y[index]<= 0.88*x[index]+delta-0.3 )

def Mask_nonSMG( index ): return ( data[index]['850WIDE'] ==0 ) | ( data[index]['850WIDE']     ==7 )
def Mask_nonIRB( index ): return ( data[index]['24MICRON']==0 ) & ( data[index]['Radio_excess']!=0 ) #radio excess = -99 or 1

def Mask_masscut(index):
    return (data[index][mass]>11)

#def Mask_mass(index,low,high):
#    return (data[index][mass]<high) & (data[index][mass]>low)

def PlotHist_photoz( index, inputcolor, inputbins, labelname ):
    NBINS = inputbins
    plt.hist( data[0][photoz][mask[index]], NBINS,
             color=inputcolor, alpha=0.7, label = labelname) 
    #plt.title( 'Histogram photoz' )
    plt.ylabel( 'Number' , fontdict = {'fontsize' : 14} )
    plt.xlabel( 'z'      , fontdict = {'fontsize' : 14} )
    if (index==0): plt.xlim(0,5)
    else: plt.xlim(0,4)
    plt.legend( fontsize=12, frameon=False )
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
    #plt.xlim(0,3.5)
    plt.legend()
    plt.show()
    return

def PlotHist_M(index):
    NBINS = 70
    histcolor = [ 'b', 'g', 'r' ]
    for i in range(3):
        plt.subplot(3, 1, i+1)
        histdata = data[0][color[i]][mask[index]]
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
    
    plt.text( -1.4, 3.3, 'Quiescent Galaxies',    fontsize=12 )
    plt.text( -1.4, 2.7, 'Star-Forming Galaxies', fontsize=12 )
    plt.show()
    
    maskQG    = mask[index] & Mask_myclassQG( index )
    maskSFG   = mask[index] & Mask_myclassSFG(index )
    total_num = len( data[0].filled()[mask[index]]  )
    QG_num    = len( data[0].filled()[maskQG]       )
    SFG_num   = len( data[0].filled()[maskSFG]      )
    
    print "total sample size ", total_num
    print "QG    sample size ", QG_num
    print "QG    percentage  ", ( "%.1f" % (QG_num/total_num*100)  ),"$\pm$", ( "%.1f" % (np.sqrt(QG_num)/total_num*100)  ),"%"
    print "SFG   sample size ", SFG_num
    print "SFG   percentage  ", ( "%.1f" % (SFG_num/total_num*100) ),"%"
    print
    
    return

def PlotMvsZ( index, inputcolor, inputlabel ):
    
    plt.scatter(data[0][mask[index]][photoz],data[0][mask[index]][mass],s=1,alpha=0.1,color=inputcolor,label=inputlabel)
    plt.xlabel( "z", fontdict = {'fontsize' : 14} )
    plt.ylabel( "log( stellar mass($M_{\odot}$) )", fontdict = {'fontsize' : 14} )
    plt.legend()
    
    return

def PrintQGfraction(index):
    maskQG = mask[index] & Mask_myclassQG(index)
    all_num = len(data[0][mask[index]])#Num_total(index)
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

def PrintAGNinAllfraction(index):
    AGN_num = len(data[0][mask[index]])
    All_num = len(data[0][mask[0]])
    print AGN_num, All_num
    fraction = AGN_num*100.0 / All_num
    error =  np.sqrt(AGN_num)*100.0 / All_num
    if (fraction+error>100):
        up_err = 100.0-fraction
    else:
        up_err = error
    if (fraction-error<0.0):
        low_err = fraction
    else:
        low_err = error
    return fraction, np.array([[low_err,up_err]]).T

def PrintAGNinQGfraction(index):
    AGN_num = len(data[0][mask[index]&Mask_myclassQG(0)])
    QG_all_num = len(data[0][mask[1]])
    print AGN_num, QG_all_num
    fraction = AGN_num*100.0 / QG_all_num
    error =  np.sqrt(AGN_num)*100.0 / QG_all_num
    if (fraction+error>100):
        up_err = 100.0-fraction
    else:
        up_err = error
    if (fraction-error<0.0):
        low_err = fraction
    else:
        low_err = error
    return fraction, np.array([[low_err,up_err]]).T

def PrintAGNinSFGfraction(index):
    AGN_num = len(data[0][mask[index]&Mask_myclassSFG(0)])
    SFG_all_num = len(data[0][mask[2]])
    print AGN_num, SFG_all_num
    fraction = AGN_num*100.0 / SFG_all_num
    error =  np.sqrt(AGN_num)*100.0 / SFG_all_num
    if (fraction+error>100):
        up_err = 100.0-fraction
    else:
        up_err = error
    if (fraction-error<0.0):
        low_err = fraction
    else:
        low_err = error
    return fraction, np.array([[low_err,up_err]]).T

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
        #datamask1 = mask[1] & mask_zbin[i]
        #datamaskQG1 = mask[1] & mask_zbin[i] & Mask_myclassQG(1)
        datax = x[index][datamask]
        datay = y[index][datamask]
        plt.scatter( datax, datay, s=set_s, alpha=set_alpha,color=inputcolor, label=labelname)
        #print len(x[index][datamaskQG].filled())
        if len(x[index][datamask].filled())==0:
            plt.title(zbin_title[i]+', 0.00%')
        else:
            plt.title(zbin_title[i]+', '+str(len(x[index][datamaskQG].filled()))+'/'+str(len(x[index][datamask].filled()))+', '+\
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
        #plt.legend()
    plt.show()
    return

def DrawLine_UVJ_orig(delta):
    plt.plot( [1.6,1.6], [0.88*1.6+delta,3], 'k-' )
    plt.plot( [-3,(1.3-delta)/0.88], [1.3,1.3], 'k-' )
    plt.plot( [(1.3-delta)/0.88,1.6], [1.3,0.88*1.6+delta], 'k-' )
    return

def DrawLine_UVJ(delta):
    plt.plot( [1.6,1.6], [0.88*1.6+delta-0.3,3], 'k-' )
    plt.plot( [-3,(1.3-delta)/0.88], [1.3-0.3,1.3-0.3], 'k-' )
    plt.plot( [(1.3-delta)/0.88,1.6], [1.3-0.3,0.88*1.6+delta-0.3], 'k-' )
    return

def Plot_zbin_UVJ(index, scale, struc, limit, line, inputcolor,labelname):
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
        if (i==0):
            datamask = mask[index] & mask_zbin[i]
            datamaskQG = mask[index] & mask_zbin[i] & Mask_myclassQG_UVJ(index,0.69)
            datamask1 = mask[1] & mask_zbin[i]
            datamaskQG1 = mask[1] & mask_zbin[i] & Mask_myclassQG_UVJ(1,0.69)
        elif (i==1):
            datamask = mask[index] & mask_zbin[i]
            datamaskQG = mask[index] & mask_zbin[i] & Mask_myclassQG_UVJ(index,0.59)
            datamask1 = mask[1] & mask_zbin[i]
            datamaskQG1 = mask[1] & mask_zbin[i] & Mask_myclassQG_UVJ(1,0.59)
        else:
            datamask = mask[index] & mask_zbin[i]
            datamaskQG = mask[index] & mask_zbin[i] & Mask_myclassQG_UVJ(index,0.49)
            datamask1 = mask[1] & mask_zbin[i]
            datamaskQG1 = mask[1] & mask_zbin[i] & Mask_myclassQG_UVJ(1,0.49)
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
            plt.axis([-0.7,2.3,-0.5,2.3])
        if (scale==1):
            plt.axis('scaled')
        if (line==1):
            if (i==0):
                DrawLine_UVJ(0.69)
            elif (i==1):
                DrawLine_UVJ(0.59)
            else:
                DrawLine_UVJ(0.49)
        #plt.legend()
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
    





def PlotColorColor( MASKALL ):
    
    SetupData( MASKALL                          )
    SetupData( MASKALL & Mask_myclassQG(0)        )
    SetupData( MASKALL & (data[0]['24MICRON']==1) )
    SetupData( MASKALL & (data[0]['3GHZ']    ==1) )
    
    #print np.ma.median(data[0][mask[0]][photoz])
    #PlotHist_M(1)
    
    
    fig = plt.figure(1)
    Plot( 0,0,1,1,1,'C0',0,'COSMOS2015',1,'o' )
    #Plot_zbin(0, 0, 1, 1, 1, 'C0',"")
    
    patch = mpatches.Patch( color='C0', label='COSMOS2015' )
    plt.legend( handles=[patch], fontsize=12, frameon=False, loc=3 )
    #fill_x = [-2,(3.1-1.0)/3.0,2.0,-2]
    #fill_y = [3.1,3.1,7,7]
    #plt.fill(fill_x,fill_y,'k',alpha=0.3)
    fig.savefig('NUVrJ_1.png', bbox_inches = 'tight', format='png', dpi=400)


    
    fig = plt.figure(2)
    Plot( 2,0,0.5,1,1,'C1',0,'24 micron',1,'o' )
    patch = mpatches.Patch(color='C1', label='24 micron')
    plt.legend(handles=[patch], fontsize=12, frameon=False, loc=3 )
    plt.text( -1.3, 6.5, "(a)", ha='center', va='center', fontsize=14 )
    fig.savefig('NUVrJ_24.png', bbox_inches = 'tight', format='png', dpi=400)
    
    fig = plt.figure(3)
    Plot( 3,0,0.5,1,1,'r',0,'3 GHz',1,'o' )
    patch = mpatches.Patch(color='r', label='3 GHz')
    plt.legend(handles=[patch] , fontsize=12, frameon=False, loc=3 )
    plt.text( -1.3, 6.5, "(b)", ha='center', va='center', fontsize=14 )
    fig.savefig('NUVrJ_3.png', bbox_inches = 'tight', format='png', dpi=400)
    
    
    

    fig = plt.figure(4)    
    PlotHist_photoz( 0,'C0',50,"SFG+QG"    )
    plt.text( 0.3, 15000, "(a)", ha='center', va='center', fontsize=14 ) #17500
    fig.savefig('hist_all.png', bbox_inches = 'tight', format='png', dpi=400)
    
    fig = plt.figure(5)
    PlotHist_photoz( 1,'k',40,"QG"    )
    plt.text( 0.3, 2750, "(b)", ha='center', va='center', fontsize=14 ) #2600
    fig.savefig('hist_QG.png', bbox_inches = 'tight', format='png', dpi=400)
    
    
    fig = plt.figure(6)
    PlotMvsZ( 0, 'C0', "SFG+QG" )
    PlotMvsZ( 1, 'k', "QG" )
    patch1 = mpatches.Patch( color='C0', alpha=0.7, label="SFG+QG" )
    patch2 = mpatches.Patch( color='k', alpha=0.7,  label="QG"     )
    plt.legend( handles=[patch1,patch2], fontsize=12, frameon=False )
    plt.text( 0.1, 12.2, "(c)", ha='center', va='center', fontsize=14 )
    fig.savefig('mass_z.png', bbox_inches = 'tight', format='png', dpi=400)
        
    
    plt.show()
    return

def PlotColorColor_fourband():
    PlotColorColor( MaskAll(0) )
    return

def PlotColorColor_oneband():
    PlotColorColor( MaskAll_oneband(0) )
    return

def PlotColorColor_fouronebandcompare():
    
    SetupData( MaskAll(0) )
    SetupData( MaskAll_oneband(0) )
    SetupData( MaskAll(0) & Mask_classQG(0) )
    SetupData( MaskAll_oneband(0) & Mask_classQG(0) )
    
    plt.figure(1)    
    PlotHist_photoz_para( [0,1],['C0','C1'],['4band','1band'],20 )
    
    plt.figure(2)    
    PlotHist_photoz_para( [2,3],['C0','C1'],['4band','1band'],20 )
    
    return

def PlotColorColor_UVJ():

    
    global colorname1, colorname2, colorname3, color1, color2, color3, color, set_xlable, set_ylable 
    
    ###set colors
    colorname1 = "U"
    colorname2 = "V"
    colorname3 = "J"
    
    ###set columns in the main catalog
    color1 = "MU"
    color2 = "MV"
    color3 = "MJ"
    color  = [ color1, color2, color3 ]

    set_xlable = '$M_{'+colorname2+'}-M_{'+colorname3+'}$'
    set_ylable = '$M_{'+colorname1+'}-M_{'+colorname2+'}$' 

    
    SetupData( MaskAll(0)                    )
    SetupData( MaskAll(0) & Mask_classSFG(0) )
    SetupData( MaskAll(0) & Mask_classQG(0)  )
    
    fig = plt.figure(1)
    #patch = mpatches.Patch( color='C0', label='COSMOS2015' )
    #plt.legend( handles=[patch] )
    #fill_x = [-2,(3.1-1.0)/3.0,2.0,-2]
    #fill_y = [3.1,3.1,7,7]
    #plt.fill(fill_x,fill_y,'k',alpha=0.3)
    #fig.savefig('NUVrJ_1.png', bbox_inches = 'tight', format='png', dpi=400)
    #Plot( 1,0,1,1,1,'k',0,'COSMOS2015',1,'o' )
    
    Plot_zbin_UVJ(1, 0, 1, 1, 1, 'C0',"SFG")
    #Plot_zbin_UVJ(2, 0, 1, 1, 1, 'k',"QG")
    
    return

def PlotColorColor_Submm850():
    
    mask_ALMAcpt = ( (data[0]['850WIDE']==1)|(data[0]['850WIDE']==2)|(data[0]['850WIDE']==3) ) & (data[0]['LENSING']==0)
    mask_IRcpt   =   (data[0]['850WIDE']==4)|(data[0]['850WIDE']==5)|(data[0]['850WIDE']==6)
    
    SetupData( MaskAll(0)                           )
    SetupData( MaskAll(0) &  mask_IRcpt             )
    SetupData( MaskAll(0) &  mask_ALMAcpt           )
    SetupData( MaskAll(0) & (data[0]['LENSING']==1) )
    SetupData( MaskAll(0) & (data[0]['850WIDE']!=0)& (data[0]['850WIDE']!=7) & Mask_myclassQG(0) )
    
    #for i in range( len(data[0][mask[4]]) ):
     #   print  x_masked[4][i], y_masked[4][i], data[0][mask[4]][photoz][i]
    
    fig = plt.figure(1)  
    
    #for i in range( len(data[0][mask[4]]) ):
     #   plt.text( x_masked[4][i], y_masked[4][i], data[0][mask[4]][photoz][i], fontsize=6 )
    
    plt.scatter(-1000,-1000,s=1.5,color='C0',alpha=1,label='COSMOS2015')
    
    #Plot_zbin(index, scale, struc, limit, line)  
    Plot( 1,0,-0.5,1,1, 'C4', 1,'24um or 3GHz counterpart', 0,'o' )
    Plot( 2,0,-0.5,1,1, 'C4', 1,'ALMA counterpart',         1,'o' )
    Plot( 3,0,-2,  1,1, 'C4', 1,'lensed system',            1,'*' )
    Plot( 0,0, 2,  1,1, 'C0', 0,'COSMOS2015',               1,'o' )
    
    plt.text( -1.3, -0.6, "(b)", ha='center', va='center', fontsize=14 )
    
    #Plot( 4,0, 2,  1,1, 'C0', 0,'850',               1,'o' )
    
    #fill_x = [-2,(3.1-1.0)/3.0,2.0,-2]
    #fill_y = [3.1,3.1,7,7]
    #plt.fill(fill_x,fill_y,'k',alpha=0.3)  

    fig.savefig('NUVrJ_850.png', bbox_inches = 'tight', format='png', dpi=400)
    
    return

def PlotColorColor_Submm850_oneband():
    
    mask_ALMAcpt = ( (data[0]['850WIDE_ONEBAND']==1)|(data[0]['850WIDE_ONEBAND']==2)|(data[0]['850WIDE_ONEBAND']==3) ) 
    mask_IRcpt   =   (data[0]['850WIDE_ONEBAND']==4)|(data[0]['850WIDE_ONEBAND']==5)|(data[0]['850WIDE_ONEBAND']==6)
    
    SetupData( MaskAll_oneband(0) )
    SetupData( MaskAll_oneband(0) &  mask_IRcpt )
    SetupData( MaskAll_oneband(0) &  mask_ALMAcpt & (data[0]['LENSING']==0) )
    SetupData( MaskAll_oneband(0) &  mask_ALMAcpt & (data[0]['LENSING']==1) )
    SetupData( MaskAll_oneband(0) & (data[0]['850WIDE_ONEBAND']!=0)& (data[0]['850WIDE_ONEBAND']!=7) & Mask_myclassQG(0) )
    
    #for i in range( len(data[0][mask[4]]) ):
    #    print  x_masked[4][i], y_masked[4][i], data[0][mask[4]][photoz][i]
    
    fig = plt.figure(1)  
    
    #for i in range( len(data[0][mask[4]]) ):
    #    plt.text( x_masked[4][i], y_masked[4][i], data[0][mask[4]][photoz][i], fontsize=6 )
    
    plt.scatter(-1000,-1000,s=1.5,color='C0',alpha=1,label='COSMOS2015')
    
    #Plot_zbin(index, scale, struc, limit, line)  
    Plot( 1,0,-0.5,1,1, 'C4', 1,'24um or 3GHz counterpart', 0,'o' )
    Plot( 2,0,-0.5,1,1, 'C4', 1,'ALMA counterpart',         1,'o' )
    Plot( 3,0,-2,  1,1, 'C4', 1,'lensed system',            1,'*' )
    Plot( 0,0, 2,  1,1, 'C0', 0,'COSMOS2015',               1,'o' )
    
    plt.text( -1.3, -0.6, "(b)", ha='center', va='center', fontsize=14 )
    
    #Plot( 4,0, 2,  1,1, 'C0', 0,'850',               1,'o' )
    
    #fill_x = [-2,(3.1-1.0)/3.0,2.0,-2]
    #fill_y = [3.1,3.1,7,7]
    #plt.fill(fill_x,fill_y,'k',alpha=0.3)  

    fig.savefig('NUVrJ_850.png', bbox_inches = 'tight', format='png', dpi=400)
    
    return

def PlotColorColor_Submm450():
    
    mask_ALMAcpt = (data[0]['450NARROW']==1)|(data[0]['450NARROW']==2)|(data[0]['450NARROW']==3)
    mask_IRcpt   = (data[0]['450NARROW']==4)|(data[0]['450NARROW']==5)|(data[0]['450NARROW']==6)
    
    SetupData( MaskAll(0)                           )
    SetupData( MaskAll(0) &  mask_IRcpt             )
    SetupData( MaskAll(0) &  mask_ALMAcpt & (data[0]['LENSING']==0) )
    SetupData( MaskAll(0) &  mask_ALMAcpt & (data[0]['LENSING']==1) )
    SetupData( MaskAll(0) & (data[0]['450NARROW']!=0)& (data[0]['450NARROW']!=7) & Mask_myclassQG(0) )
    
    fig = plt.figure(1)  
    
    #for i in range( len(data[0][mask[4]]) ):
        #plt.text( x_masked[4][i], y_masked[4][i], data[0][mask[4]][photoz][i], fontsize=6 )
    
    plt.scatter(-1000,-1000,s=1.5,color='C0',alpha=1,label='COSMOS2015')
    
    #Plot_zbin(index, scale, struc, limit, line)  
    Plot( 1,0,-0.5,1,1, 'C6', 1,'24um or 3GHz counterpart', 0,'o' )
    Plot( 2,0,-0.5,1,1, 'C6', 1,'ALMA counterpart',         1,'o' )
    Plot( 3,0,-2,  1,1, 'C6', 1,'lensed system',            1,'*' )
    Plot( 0,0, 2,  1,1, 'C0', 0,'COSMOS2015',               1,'o' )
    
    plt.text( -1.3, -0.6, "(a)", ha='center', va='center', fontsize=14 )
    
    #Plot( 4,0, 2,  1,1, 'C0', 0,'450',               1,'o' )
    
    #fill_x = [-2,(3.1-1.0)/3.0,2.0,-2]
    #fill_y = [3.1,3.1,7,7]
    #plt.fill(fill_x,fill_y,'k',alpha=0.3)  

    fig.savefig('NUVrJ_450.png', bbox_inches = 'tight', format='png', dpi=400)
    
    return

def PlotColorColor_Submm450_oneband():
    
    mask_ALMAcpt = (data[0]['450NARROW_ONEBAND']==1)|(data[0]['450NARROW_ONEBAND']==2)|(data[0]['450NARROW_ONEBAND']==3)
    mask_IRcpt   = (data[0]['450NARROW_ONEBAND']==4)|(data[0]['450NARROW_ONEBAND']==5)|(data[0]['450NARROW_ONEBAND']==6)
    
    SetupData( MaskAll_oneband(0)                           )
    SetupData( MaskAll_oneband(0) &  mask_IRcpt             )
    SetupData( MaskAll_oneband(0) &  mask_ALMAcpt & (data[0]['LENSING']==0) )
    SetupData( MaskAll_oneband(0) &  mask_ALMAcpt & (data[0]['LENSING']==1) )
    SetupData( MaskAll_oneband(0) & (data[0]['450NARROW_ONEBAND']!=0)& (data[0]['450NARROW_ONEBAND']!=7) & Mask_myclassQG(0) )
    
    fig = plt.figure(1)  
    
    #for i in range( len(data[0][mask[4]]) ):
    #    plt.text( x_masked[4][i], y_masked[4][i], data[0][mask[4]][photoz][i], fontsize=6 )
    
    plt.scatter(-1000,-1000,s=1.5,color='C0',alpha=1,label='COSMOS2015')
    
    #Plot_zbin(index, scale, struc, limit, line)  
    Plot( 1,0,-0.5,1,1, 'C6', 1,'24um or 3GHz counterpart', 0,'o' )
    Plot( 2,0,-0.5,1,1, 'C6', 1,'ALMA counterpart',         1,'o' )
    Plot( 3,0,-2,  1,1, 'C6', 1,'lensed system',            1,'*' )
    Plot( 0,0, 2,  1,1, 'C0', 0,'COSMOS2015',               1,'o' )
    
    plt.text( -1.3, -0.6, "(a)", ha='center', va='center', fontsize=14 )
    
    #Plot( 4,0, 2,  1,1, 'C0', 0,'450',               1,'o' )
    
    #fill_x = [-2,(3.1-1.0)/3.0,2.0,-2]
    #fill_y = [3.1,3.1,7,7]
    #plt.fill(fill_x,fill_y,'k',alpha=0.3)  

    fig.savefig('NUVrJ_450.png', bbox_inches = 'tight', format='png', dpi=400)
    
    return

def PrintRow( index, title, index_total ):
    
    num_all    = len( data[0][mask[index]] )
    num_QG     = len( data[0][mask[index] & Mask_myclassQG(0)] )
    num_QG_all = len( data[0][mask[index_total] & Mask_myclassQG(0)] )
    
    #print title," ",num_all,"    ",num_QG,"   ",( "%.1f" % (num_QG/num_all*100) ),"+/-",( "%.1f" % (np.sqrt(num_QG)/num_all*100) )
    print title,"&",num_all,"&",num_QG,\
          "&",( "%.1f" % (num_QG/num_all*100) ),"$\pm$",( "%.1f" % (np.sqrt(num_QG)/num_all*100) ),"\%",\
          "&",( "%.2f" % (num_QG/num_QG_all*100) ),"$\pm$",( "%.2f" % (np.sqrt(num_QG)/num_QG_all*100) ),"\%",\
          "\\\\"

    return

def Mask_Region( inputra,inputdec,inputradius ): #mask is not list
    
    c0     = SkyCoord( ra=data[0][ra], dec=data[0][dec] )
    center = SkyCoord( inputra, inputdec, frame='fk5' )
    radius = inputradius
    sep    = center.separation( c0 )
    
    return sep<=radius

def PrintMatchNumTable( MASKALL, wide850, narrow450 ):
    
    print "printing Matching number table ..."
    print
    
    SetupData( MASKALL                            )
    SetupData( MASKALL & (data[0]['24MICRON']==1) )
    SetupData( MASKALL & (data[0]['3GHZ']    ==1) )
    SetupData( MASKALL & (data[0][wide850]  !=0)& (data[0][wide850]  !=7) )
    
    SetupData( MASKALL & Mask_Region( '10h00m25.0s', '2d24m22.0s', 0.2*u.degree ) )
    SetupData( MASKALL & (data[0][narrow450]!=0)& (data[0][narrow450]!=7)\
               & Mask_Region( '10h00m25.0s', '2d24m22.0s', 0.2*u.degree ) )
    
    print "       SFGs+QGs     QGs      QGs/(SFGs+QGs)      QGs/(all QGs)"
    PrintRow( 0, "total in the COSMOS field",0 )
    PrintRow( 1, "24 um detected",0 )
    PrintRow( 2, "3 GHz detected",0 )
    PrintRow( 3, "850 um detected",0 )
    print "\hline"
    PrintRow( 4, "total in the STUDIES field",4 )
    PrintRow( 5, "450 um detected",4 )
    
    return

def PrintMatchNumTable_fourband():
    PrintMatchNumTable( MaskAll(0),'850WIDE','450NARROW' )
    return

def PrintMatchNumTable_oneband():
    PrintMatchNumTable( MaskAll_oneband(0),'850WIDE_ONEBAND','450NARROW_ONEBAND' )
    return

def PlotColorColor_radioAGN():
    
    SetupData( MaskAll(0) & (data[0]['Radio_excess']==1) & (Mask_masscut(0)) )   
    fig = plt.figure(1)    
    Plot(0,0,0,1,1,'r',0,'radio AGN',1,'o')
    patch = mpatches.Patch(color='r', label='radio AGN')
    plt.legend(handles=[patch] , fontsize=12, frameon=False, loc=3 )
    plt.text( -1.3, 6.5, "(a)", ha='center', va='center', fontsize=14 )
    #PlotHist_photoz( 0,'r',10,"radio"    )
    fig.savefig('NUVrJradioAGN.png', bbox_inches = 'tight', format='png', dpi=400)
    
    return

def PlotColorColor_IRAGN():

    SetupData( MaskAll(0) & (data[0]['agn_c17b']==True) & (Mask_masscut(0)) )   
    fig = plt.figure(1)
    Plot(0,0,0,1,1,'C1',0,'mid-IR AGN',1,'o')
    patch = mpatches.Patch(color='C1', label='mid-IR AGN')
    plt.legend(handles=[patch] , fontsize=12, frameon=False, loc=3 )
    plt.text( -1.3, 6.5, "(b)", ha='center', va='center', fontsize=14 )
    #PlotHist_photoz( 0,'C1',10,"IR"    )
    fig.savefig('NUVrJIRAGN.png', bbox_inches = 'tight', format='png', dpi=400)
    
    return

def PlotColorColor_XrayAGN():

    SetupData( MaskAll(0) & (data[0]['agn_xxx']==True) & (Mask_masscut(0)) )    
    fig = plt.figure(1)
    Plot(0,0,0,1,1,'g',0,'X-ray AGN',1,'o')
    patch = mpatches.Patch(color='g', label='X-ray AGN')
    plt.legend(handles=[patch] , fontsize=12, frameon=False, loc=3 )
    plt.text( -1.3, 6.5, "(c)", ha='center', va='center', fontsize=14 )
    #PlotHist_photoz( 0,'g',10,"X"    )
    fig.savefig('NUVrJXrayAGN.png', bbox_inches = 'tight', format='png', dpi=400)
    
    return



def PlotColorColor_radioAGN_oneband():
    
    SetupData( MaskAll_oneband(0) & (data[0]['Radio_excess']==1) & (Mask_masscut(0)) )   
    fig = plt.figure(1)    
    Plot(0,0,0,1,1,'r',0,'radio AGN',1,'o')
    patch = mpatches.Patch(color='r', label='radio AGN')
    plt.legend(handles=[patch] , fontsize=12, frameon=False, loc=3 )
    plt.text( -1.3, 6.5, "(a)", ha='center', va='center', fontsize=14 )
    #PlotHist_photoz( 0,'r',10,"radio"    )
    fig.savefig('NUVrJradioAGN.png', bbox_inches = 'tight', format='png', dpi=400)
    
    return

def PlotColorColor_IRAGN_oneband():

    SetupData( MaskAll_oneband(0) & (data[0]['agn_c17b']==True) & (Mask_masscut(0)) )   
    fig = plt.figure(1)
    Plot(0,0,0,1,1,'C1',0,'mid-IR AGN',1,'o')
    patch = mpatches.Patch(color='C1', label='mid-IR AGN')
    plt.legend(handles=[patch] , fontsize=12, frameon=False, loc=3 )
    plt.text( -1.3, 6.5, "(b)", ha='center', va='center', fontsize=14 )
    #PlotHist_photoz( 0,'C1',10,"IR"    )
    fig.savefig('NUVrJIRAGN.png', bbox_inches = 'tight', format='png', dpi=400)
    
    return

def PlotColorColor_XrayAGN_oneband():

    SetupData( MaskAll_oneband(0) & (data[0]['agn_xxx']==True) & (Mask_masscut(0)) )    
    fig = plt.figure(1)
    Plot(0,0,0,1,1,'g',0,'X-ray AGN',1,'o')
    patch = mpatches.Patch(color='g', label='X-ray AGN')
    plt.legend(handles=[patch] , fontsize=12, frameon=False, loc=3 )
    plt.text( -1.3, 6.5, "(c)", ha='center', va='center', fontsize=14 )
    #PlotHist_photoz( 0,'g',10,"X"    )
    fig.savefig('NUVrJXrayAGN.png', bbox_inches = 'tight', format='png', dpi=400)
    
    return

def Plot_QGfraction( MASKALL ):
    
    SetupData( MASKALL & Mask_masscut(0) )
    SetupData( MASKALL & Mask_masscut(0) & Mask_myclassQG(0) )
    SetupData( MASKALL & Mask_masscut(0) &\
              ((data[0]['Radio_excess']==1)|(data[0]['agn_c17b']==True)|(data[0]['agn_xxx']==True)) )
    SetupData( MASKALL & Mask_masscut(0) & (data[0]['Radio_excess']==1   ) )
    SetupData( MASKALL & Mask_masscut(0) & (data[0]['agn_c17b']    ==True) )
    SetupData( MASKALL & Mask_masscut(0) & (data[0]['agn_xxx']     ==True) )
    
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
    plt.ylabel('QG fraction(%)', fontdict = {'fontsize' : 16})#14})
    plt.text( 0.5, 65, "(a)", ha='center', va='center', fontsize=14 )
    plt.show()
    
    fig.savefig('QGfraction.png', bbox_inches = 'tight', format='png', dpi=400)
 


    fig, axes = plt.subplots()
    
    for i in range(5):
        print "plotting ", label_list[i], "..."
        fraction, error = Print_zbin_QGfraction( index_list[i] )
        print
        a    = i*0.01
        zbin = [0.25+a,0.75+a,1.25+a,1.75+a,2.25+a,2.75+a,3.25+a,3.75+a]
        plt.errorbar( zbin, fraction, error, fmt='o-',
                      color=color_list[i], label=label_list[i], capsize=5 )

    #plt.title('QG percentage, mass>10^11', fontdict = {'fontsize' : 20})
    plt.xlabel( 'z',               fontdict = {'fontsize' : 16} )
    plt.ylabel( 'QG fraction (%)', fontdict = {'fontsize' : 16} )
    plt.axis([0,4,-5,120])
    plt.legend( fontsize=12, frameon=False )
    plt.text( 0.3, 110, "(b)", ha='center', va='center', fontsize=14 )
    plt.show()
    
    fig.savefig('QGfraction_zbin.png', bbox_inches = 'tight', format='png', dpi=400)
    
    return

def Plot_QGfraction_fourband():
    Plot_QGfraction( MaskAll(0) )
    return

def Plot_QGfraction_oneband():
    Plot_QGfraction( MaskAll_oneband(0) )
    return

def Plot_AGNfraction( MASKALL ):
    
    SetupData( MASKALL & Mask_masscut(0)  )
    SetupData( MASKALL & Mask_masscut(0) & Mask_myclassQG(0)  )
    SetupData( MASKALL & Mask_masscut(0) & Mask_myclassSFG(0)  )
    SetupData( MASKALL & Mask_masscut(0)  &\
              ((data[0]['Radio_excess']==1)|(data[0]['agn_c17b']==True)|(data[0]['agn_xxx']==True)) )
    SetupData( MASKALL & Mask_masscut(0) & (data[0]['Radio_excess']==1   ) )
    SetupData( MASKALL & Mask_masscut(0) & (data[0]['agn_c17b']    ==True) )
    SetupData( MASKALL & Mask_masscut(0) & (data[0]['agn_xxx']     ==True) )
    
    label_list = ['all AGN', 'radio AGN', 'IR AGN', 'X-ray AGN']
    index_list = [3,4,5,6]
    
    fig, axes = plt.subplots()

    plt.setp( axes, xticks=[y_axes+1 for y_axes in range(len(label_list))], xticklabels=label_list )
    
    for i in range(len(label_list)):
        print "plotting ", label_list[i], "/ All ..."
        fraction, error = PrintAGNinAllfraction( index_list[i] )
        print "QG fraction ", fraction
        print "error ", error[0],error[1]
        print
        plt.errorbar( i+1, fraction, error, color='k', fmt='o', capsize=5 ) #color_list[i]
        #plt.text( i+0.8, fraction-error[0][0]-1.5, str("%.2f" % fraction), fontsize=10 )

        print "plotting ", label_list[i], "/ QG ..."
        fraction, error = PrintAGNinQGfraction( index_list[i] )
        print "QG fraction ", fraction
        print "error ", error[0],error[1]
        print
        plt.errorbar( i+1, fraction, error, color='r', fmt='o', capsize=5 ) #color_list[i]
        #plt.text( i+0.8, fraction-error[0][0]-1.5, str("%.2f" % fraction), fontsize=10 )

        print "plotting ", label_list[i], "/ SFG ..."
        fraction, error = PrintAGNinSFGfraction( index_list[i] )
        print "QG fraction ", fraction
        print "error ", error[0],error[1]
        print
        plt.errorbar( i+1, fraction, error, color='b', fmt='o', capsize=5 ) #color_list[i]
        #plt.text( i+0.8, fraction-error[0][0]-1.5, str("%.2f" % fraction), fontsize=10 )

    #plt.axis([0,6,0,70])
    #plt.title('QG percentage, mass>10^11', fontdict = {'fontsize' : 16})
    plt.ylabel('AGN fraction(%)', fontdict = {'fontsize' : 16})#14})
    #plt.text( 0.5, 65, "(a)", ha='center', va='center', fontsize=14 )
    plt.show()
    
    patch1 = mpatches.Patch(color='k', label='QG+SFG')
    patch2 = mpatches.Patch(color='r', label='QG')
    patch3 = mpatches.Patch(color='b', label='SFG')
    plt.legend(handles=[patch1,patch2,patch3] , fontsize=12, frameon=False, loc=1 )
    
    fig.savefig('AGNfraction.png', bbox_inches = 'tight', format='png', dpi=400)
 

    '''
    fig, axes = plt.subplots()
    
    for i in range(5):
        print "plotting ", label_list[i], "..."
        fraction, error = Print_zbin_QGfraction( index_list[i] )
        print
        a    = i*0.01
        zbin = [0.25+a,0.75+a,1.25+a,1.75+a,2.25+a,2.75+a,3.25+a,3.75+a]
        plt.errorbar( zbin, fraction, error, fmt='o-',
                      color=color_list[i], label=label_list[i], capsize=5 )

    #plt.title('QG percentage, mass>10^11', fontdict = {'fontsize' : 20})
    plt.xlabel( 'z',               fontdict = {'fontsize' : 16} )
    plt.ylabel( 'QG fraction (%)', fontdict = {'fontsize' : 16} )
    plt.axis([0,4,-5,120])
    plt.legend( fontsize=12, frameon=False )
    plt.text( 0.3, 110, "(b)", ha='center', va='center', fontsize=14 )
    plt.show()
    
    fig.savefig('QGfraction_zbin.png', bbox_inches = 'tight', format='png', dpi=400)
    '''
    return

def Plot_AGNfraction_fourband():
    Plot_AGNfraction( MaskAll(0) )
    return

def Plot_AGNfraction_oneband():
    Plot_AGNfraction( MaskAll_oneband(0) )
    return

def RegionFile(index,filename,mode,color,size):
    print 'Output region file...'
    f = open(filename+'.reg','w')
    f.write('global color='+color+' font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    
    N = len( data[0].filled()[mask[index]] )
    #N = len( data[index] )
    
    print 'source number: ',N
    
    if (mode==0): #no label
        for n in range(N):
            f.write('fk5;circle('+str(data[0][ra][mask[index]][n])+','+str(data[0][dec][mask[index]][n])+','+size+'") # text={''}\n')
            #f.write('fk5;circle('+str(data[index][ra][mask[index]][n])+','+str(data[index][dec][mask[index]][n])+','+size+'") # text={''}\n')
    elif (mode==1): # redshift label
        for n in range(N):
            f.write('fk5;circle('+str(data[0][ra][mask[index]][n])+','+str(data[0][dec][mask[index]][n])+','+size+'") # text={'+str(data[0][photoz][mask[index]][n])+'}\n')
    
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

def taskRegionFile_Unmask():
    path = '/Users/yuhsuan/Documents/research/05WH/data/COSMOS/'
    #catalog = path+'05_3GHz/vla3_cosmos_sources_160321_public5sig.fits.txt'
    catalog = path+"07_ALMA/apjsab42da_table4/A-COSMOS_blind.fits"
    ReadCatalog(0,catalog)
    #RegionFileUnmask(0,'COSMOS_3GHz_5sig',0,'red','2.0')
    RegionFileUnmask(0,'COSMOS_A3COSMOS',0,'yellow','0.5',"RA","DEC")
    return

'''
def MaskRegion(inputra,inputdec,inputradius):
    c0 = SkyCoord(ra=data[0][mask[0]][ra], dec=data[0][mask[0]][dec])
    center = SkyCoord(inputra, inputdec, frame='fk5')
    radius = inputradius
    sep = center.separation(c0)
    return sep<=radius
'''

def taskRegionFile():
    
    SetupData( MaskAll_oneband(0) )
    SetupData( MaskAll_oneband(0) & Mask_myclassQG(0) )
    
    alma_mask = (data[0]['850WIDE_ONEBAND']==1)|(data[0]['850WIDE_ONEBAND']==2)|(data[0]['850WIDE_ONEBAND']==3)
    SetupData( MaskAll_oneband(0) & alma_mask & Mask_myclassQG(0) )
    
    
    #RegionFile(1,"COSMOS_QG_1band_z",1,'black','1.0')
    RegionFile(2,"COSMOS_QG_ALMA850detected_1band_z",1,'yellow','15.0')

    return

def PrintLensing_z():
    print 0.98*(2.2/14.5e-3)**0.26-1
    print 0.98*(8.84/49.7e-3)**0.26-1
    return


def stacking_sfr_tmp():
    L_IR450 = [ 2.5573646e09, 3.4028356e10, 3.5214292e10, 6.4974896e10,
                1.4658295e10, 6.3934115e08, 3.4028356e10 ]
    L_IR850 = [ 1.0393994e10, 1.1283071e11, 6.7596636e10, 1.3854581e11,
                1.6799213e09, 6.7196851e09, 0.0000000,     6.8162994e09,
                1.5864483e09, 8.0078843e09, -9.4678762e09, 3.5750400e10,
                3.5387951e10, 3.7364623e10, 6.8862904e10, 1.2218141e11,
                1.1661365e11, 6.6354851e10, 2.0047333e10, 1.0498772e11,
                3.1772969e10, 1.6801486e11, 2.4162867e11, 1.0533251e12,
                1.0103987e12, 1.0587750e12, 3.8567891e11]
    for i in range(len(L_IR450)):
        print L_IR450[i]*1.7e-10
    print
    for i in range(len(L_IR850)):
        print L_IR850[i]*1.7e-10
    return

def PlotMagerrMag( band, mag, error, index ):
    
    plt.scatter(data[0][mask[0]][error],data[0][mask[0]][mag],s=2,alpha=0.1)
    plt.xlabel( band+" error", fontdict = {'fontsize' : 14} )
    plt.ylabel( band,          fontdict = {'fontsize' : 14} )
    plt.gca().invert_yaxis()
    plt.xlim(-0.1,0.8)
    plt.ylim(29,17)
    print

    d = 0.001
    SetupData( MaskAll(0) & (data[0][error]<(0.2+d)) & (data[0][error]>(0.2-d)) )
    print len( data[0][mask[index]][mag] )
    mean = np.log10( np.mean( 10**data[0][mask[index]][mag] ) )
    median = np.median( data[0][mask[index]][mag] )
    print 'mean   of err=0.2: ', mean
    print 'median of err=0.2: ', median
    plt.axhline( mean,   linestyle='dotted', color='g' )
    plt.axhline( median, linestyle='dotted', color='r' )  
    plt.axvline( 0.1, color='k' )
    
    return

def PlotLimitMag():
    
    SetupData( MaskAll(0) )
    #SetupData( MaskAll(0) & (data[0][V_error]<0.0) )
    #SetupData( MaskAll(0) & (data[0][V_error]<-20.0) )
    
    print len( data[0][mask[0]] )
    #print len( data[0][mask[1]] )
    #print len( data[0][mask[2]] )
    
    fig = plt.subplot(2,2,1)
    PlotMagerrMag( "V", V_mag, V_error, 1 )
       
    
    fig = plt.subplot(2,2,2)
    PlotMagerrMag( "ip", ip_mag, ip_error, 2 )
    
    fig = plt.subplot(2,2,3)
    PlotMagerrMag( "J", J_mag, J_error, 3 )
    
    fig = plt.subplot(2,2,4)
    PlotMagerrMag( "Ks", Ks_mag, Ks_error, 4 )

    #fig.savefig('1.pdf', bbox_inches = 'tight', format='png', dpi=1200)
    
    return

def PrintBSMGfraction():
    print 14764+1821, 1821.0/(14764.0+1821.0),np.sqrt(1821)/(14764.0+1821.0)
    print 8.0/1821.0, np.sqrt(8)/1821.0
    print 27.0/17811.0, np.sqrt(27)/17811.0
    return
    
def PrintIRBQGfraction():
    print str("%.1f" % (1769.0/18304.0*100.0))+"$\pm$"+str("%.1f" % (np.sqrt(1769.0)/18304.0*100.0))
    print str("%.1f" % (179.0/1846.0*100.0))+"$\pm$"+str("%.1f" % (np.sqrt(179.0)/1846.0*100.0))
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
catalog[0] = 'COSMOS2015_merged.fits' 
data       = [None]*1
x,y,mask,x_masked,y_masked = [],[],[],[],[]

ReadCatalog( 0,catalog[0] )


# ===== color-color diagram =====

#PlotColorColor_fourband()
#PlotColorColor_UVJ()

#PlotColorColor_Submm850()
#PlotColorColor_Submm450()

#PrintMatchNumTable_fourband()

#PlotColorColor_radioAGN()
#PlotColorColor_IRAGN()
#PlotColorColor_XrayAGN()


# ===== AGN QG correlation =====

#Plot_QGfraction_fourband()

#PrintSpatial()

#PrintSpatialSFG()

#PrintLensing_z()


#stacking_sfr_tmp()

#PlotLimitMag()


#PrintBSMGfraction()

#test()


# ===== one band version =====

#PlotColorColor_oneband()
#PlotColorColor_fouronebandcompare()

#PlotColorColor_Submm850_oneband()
#PlotColorColor_Submm450_oneband()

#PrintMatchNumTable_oneband()

#PlotColorColor_radioAGN_oneband()
#PlotColorColor_IRAGN_oneband()
#PlotColorColor_XrayAGN_oneband()

#Plot_QGfraction_oneband()

Plot_AGNfraction_oneband()

# ===== RegionFile =====

#taskRegionFile()

#PrintIRBQGfraction()

time2 = time.time()
print 'done! time =', time2-time1 , 'sec'