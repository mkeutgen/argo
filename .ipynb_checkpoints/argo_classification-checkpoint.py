#!/usr/bin/python

import time
from IPython import display

import numpy as np
import pandas as pd
import xarray as xr

import matplotlib.pyplot as plt
import argo_processing as ap

#####################################################################################################

# Write info in a CSV
def write_in_csv(profile_name, variable, variable_name, csv_file):
    """
    Input:
    profile_name  - str,  filename to modify
    variable      - list, Variable values
    variable_name - list, Variable names
    csv_file      - str,  Modified dataset location
    
    Return:
    """

    df = pd.read_csv(csv_file).iloc[:,1:]               # Open file
    mask = df['file'] == profile_name                   # Create mask for line to modify
    index = df[mask].index                              # Get index of line to modify
    for var,var_name in zip(variable, variable_name):   # For all variables in variable file
        df[var_name][index] = int(var)                  # Write variable
    df.to_csv(csv_file)                                 # Save modified dataset

# Read info in a CSV   
def read_in_csv(prof, csv_file):
    """
    Input:
    profile_name  - str,  filename to read
    csv_file      - str,  Dataset location
    
    Return:
    df - Dictionnary of all the variables
    """
    
    df = pd.read_csv(csv_file).iloc[:,1:]
    mask = np.array(df['file'] == prof)
    df = df[mask].to_dict()
    
    return df

#####################################################################################################

def trainset_analysis(ARGO_path, Analyzed_profiles, Index_file):
    n = len(Analyzed_profiles)

    #%matplotlib inline
    
    # SETUP FIGURE
    fig = plt.figure(figsize=(10,10), dpi=200)

        # Set title
    title = fig.suptitle('', y= 1.02)
    
    # Set panels
    ax1  = plt.subplot(231)
    ax2  = ax1.twiny()
    ax3  = plt.subplot(232)
    ax4  = ax3.twiny()
    ax5  = plt.subplot(233)
    ax7  = plt.subplot(234)
    ax8  = ax7.twiny()
    ax9  = plt.subplot(235)
    ax11 = plt.subplot(236)
    
    # Set lines
    line1,  = ax1.plot([0,1],[0,1000] , '-o', markersize=1, c = 'C0')  # Temperature
    line2,  = ax2.plot([0,1],[0,500]  , '-o', markersize=1, c = 'C1')  # Salinity
    line3,  = ax3.plot([0,1],[0,500]  , '-o', markersize=1, c = 'C4')  # N2
    line4,  = ax4.plot([0,1],[0,1000] , '-o', markersize=1, c = 'k')   # Density 
    line5,  = ax5.plot([0,1],[0,1000] , '-o', markersize=1, c = 'C7')  # Spice    
    line7,  = ax7.plot([0,1],[0,1000] , '-o', markersize=1, c = 'C3')  # O2
    line8,  = ax8.plot([0,1],[0,500]  , '-o', markersize=1, c = 'C8')  # -AOU
    line9,  = ax9.plot([0,1],[0,1000] , '-o', markersize=1, c = 'C5')  # POC
    line10, = ax9.plot([0,1],[0,500]  , '-o', markersize=1, c = 'k')   # POC flt
    line11, = ax11.plot([0,1],[0,1000], '-o', markersize=1, c = 'C2')  # Chla 
    line12, = ax11.plot([0,1],[0,500] , '-o', markersize=1, c = 'k')   # Chla flt
    
    # Set labels
    ax3.ticklabel_format(style = 'sci', axis = 'x', scilimits=(0,0))
    
    ax1.set_xlabel('Temperature [°C]', color = 'C0')
    ax2.set_xlabel('Salinity [psu]', color = 'C1')
    ax3.set_xlabel('N$^{2}$ [s$^{-2}$]', color = 'C4')
    ax4.set_xlabel('Density [kg.m$^{-3}$]', color = 'k')
    ax5.set_xlabel('Spice [kg.m$^{-3}$]', color = 'C7')
    ax7.set_xlabel('Oxygen [umol.kg$^{-1}$]', color = 'C3')
    ax8.set_xlabel('- AOU [umol.kg$^{-1}$]', color = 'C8')
    ax9.set_xlabel('POC [umol.kg$^{-1}$]', color = 'C5')
    ax11.set_xlabel('Chlorophyll [mg.m$^{-3}$]', color = 'C2')
    
    # LOOP for all profiles analyzed
    for i in range(n):

        argo_prof = xr.open_dataset(ARGO_path+np.array(Analyzed_profiles['file'])[i])
        argo_prof = ap.format_profil(argo_prof)
        
        date = str(argo_prof.time.data)[0:10]
        lon = round(float(argo_prof.lon.data),2)
        lat = round(float(argo_prof.lat.data),2)

        title.set_text('Profile '+str(i+1)+'/'+str(n)+' | '+date+' | Lon ='+str(lon)+', Lat='+str(lat))
        
        ## Plot variables
        ii = argo_prof.depth < 1000
        ii_mid = argo_prof.depth_mid < 1000
        
        # Temperature
        line1.set_xdata(argo_prof.ctem),line1.set_ydata(argo_prof.depth)
        xmin,xmax = np.min(argo_prof.ctem[ii]), np.max(argo_prof.ctem[ii])
        ax1.set_xlim(xmin-0.05*(xmax-xmin), xmax+0.05*(xmax-xmin))
        
        # Salinity
        line2.set_xdata(argo_prof.psal),line2.set_ydata(argo_prof.depth)
        xmin,xmax = np.min(argo_prof.psal[ii]), np.max(argo_prof.psal[ii])
        ax2.set_xlim(xmin-0.05*(xmax-xmin), xmax+0.05*(xmax-xmin))
        
        # N2
        line3.set_xdata(argo_prof.N2),line3.set_ydata(argo_prof.depth_mid)
        xmin,xmax = np.min(argo_prof.N2[ii_mid]), np.max(argo_prof.N2[ii_mid])
        ax3.set_xlim(xmin-0.05*(xmax-xmin), xmax+0.05*(xmax-xmin))
        
        # Density
        line4.set_xdata(argo_prof.sig0),line4.set_ydata(argo_prof.depth)
        xmin,xmax = np.min(argo_prof.sig0[ii]), np.max(argo_prof.sig0[ii])
        ax4.set_xlim(xmin-0.05*(xmax-xmin), xmax+0.05*(xmax-xmin))
        
        # Spice    
        line5.set_xdata(argo_prof.spic),line5.set_ydata(argo_prof.depth)
        xmin,xmax = np.min(argo_prof.spic[ii]), np.max(argo_prof.spic[ii])
        ax5.set_xlim(xmin-0.05*(xmax-xmin), xmax+0.05*(xmax-xmin))
        
        # O2 and -AOU
        line7.set_xdata(argo_prof.doxy),line7.set_ydata(argo_prof.depth)
        line8.set_xdata(-argo_prof.aou),line8.set_ydata(argo_prof.depth)
        
        xmin1,xmax1 = np.min(argo_prof.doxy[ii]), np.max(argo_prof.doxy[ii])
        xmin2,xmax2 = np.min(-argo_prof.aou[ii]), np.max(-argo_prof.aou[ii])
        
        dm1,dm2 = xmax1-xmin1,xmax2-xmin2
        
        if dm1>dm2:
            ax7.set_xlim(xmin1-0.05*dm1, xmax1+0.05*dm1)
            d_dm = (dm1-dm2)/2
            ax8.set_xlim(xmin2-d_dm-0.05*dm1, xmax2+d_dm+0.05*dm1)
        else:
            ax8.set_xlim(xmin2-0.05*dm2, xmax2+0.05*dm2)
            d_dm = (dm2-dm1)/2
            ax7.set_xlim(xmin1-d_dm-0.05*dm2, xmax1+d_dm+0.05*dm2)
        
        # POC and POC flt
        line9.set_xdata(argo_prof.poc),line9.set_ydata(argo_prof.depth)
        xmin,xmax = np.min(argo_prof.poc[ii]), np.max(argo_prof.poc[ii])
        ax9.set_xlim(xmin-0.05*(xmax-xmin), xmax+0.05*(xmax-xmin))
        
        line10.set_xdata(argo_prof.poc_flt),line10.set_ydata(argo_prof.depth)
        
        # Chla and Chla flt
        line11.set_xdata(argo_prof.chla),line11.set_ydata(argo_prof.depth)
        xmin,xmax = np.min(argo_prof.chla[ii]), np.max(argo_prof.chla[ii])
        ax11.set_xlim(xmin-0.05*(xmax-xmin), xmax+0.05*(xmax-xmin))
        
        line12.set_xdata(argo_prof.chla_flt),line12.set_ydata(argo_prof.depth)
        
        ## Set limits and grid
        ax1.set_ylim(1000,0)
        ax1.set_yticks(np.array(np.linspace(0,1000,21), dtype = int))
        ax1.yaxis.grid(color='gray', linestyle='dashed')

        ax3.set_ylim(1000,0)
        ax3.set_yticks(np.array(np.linspace(0,1000,21), dtype = int))
        ax3.yaxis.grid(color='gray', linestyle='dashed')

        ax5.set_ylim(1000,0)
        ax5.set_yticks(np.array(np.linspace(0,1000,21), dtype = int))
        ax5.yaxis.grid(color='gray', linestyle='dashed')

        ax7.set_ylim(1000,0)
        ax7.set_yticks(np.array(np.linspace(0,1000,21), dtype = int))
        ax7.yaxis.grid(color='gray', linestyle='dashed')

        ax9.set_ylim(1000,0)
        ax9.set_yticks(np.array(np.linspace(0,1000,21), dtype = int))
        ax9.yaxis.grid(color='gray', linestyle='dashed')

        ax11.set_ylim(1000,0)
        ax11.set_yticks(np.array(np.linspace(0,1000,21), dtype = int))
        ax11.yaxis.grid(color='gray', linestyle='dashed')
        
        
        plt.tight_layout()

        fig.canvas.draw()
        fig.canvas.flush_events()

        display.clear_output(wait=True)
        display.display(plt.gcf())
        time.sleep(0.01)

        # Save input
        ## Quality
        valid_input = False
        while not valid_input:
            print("########### NEW INPUT ###########")

            ## Subduction
            wrong = True
            while wrong:
                try:
                    # No signature, Filament+Organic Mater, Filament, Uncertain, Bad quality
                    sub_bool = int(input("Subduction filament ? [0: No, 1: Fil+MO, 2: Fil, 3: Idk, -1: Bad QC]")) 
                except ValueError:
                    wrong = True
                else:
                    wrong = False
                    break
            
            ## Validation
            wrong = True
            while wrong:
                try:
                    valid_input = float(input("Valid input ? [0: no, 1: yes] "))
                except ValueError:
                    wrong = True
                else:
                    wrong = False
                    break 

            write_in_csv(Analyzed_profiles['file'][i], [sub_bool], ['Subduction'], Index_file)
    
    display.clear_output(wait=True)