#!/usr/bin/python

import time
import datetime as dt
from IPython import display

import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import interp1d

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean as cmo
import matplotlib.pyplot as plt
import gsw

data_path  = "/data/GLOBARGO/float_data"
index_file_path = "argo_synthetic-profile_index_upgrade.csv"

#########################################################################################################
# Select all profiles containing certain variables
def profile_selection(req_param,index_profile):
    """
    Select all argo profiles containing required parameters among a argo filename list
    > INPUTS:
    req_param = List of required parameters (ex: ['DOXY', 'BBP700']) 
    index_profile = Index of argo filename
    
    > OUTPUTS:
    sel_profs: names of selected profiles
    mask: mask of selected profiles
    """
    
    sel_profs = []
    mask = []
    for param, prof in zip(index_profile['parameters'], index_profile['file']):
        sel = True 
        for i in req_param:
            if not i in param.split(' '):
                sel = False
        if sel:
            sel_profs.append(prof)
        mask.append(sel)
    return np.array(sel_profs), np.array(mask)

#########################################################################################################

def hampel(input_series, window_size, n_sigmas=3):    
    n = len(input_series)
    new_series = input_series.copy()
    k = 1.4826 # scale factor for Gaussian distribution
    
    indices = []
    
    # possibly use np.nanmedian
    ### JCG - NaNs are excluded from filter so not relevant
    for i in range((window_size),(n - window_size)):
        x0 = np.median(input_series[(i - window_size):(i + window_size)])
        S0 = k * np.median(np.abs(input_series[(i - window_size):(i + window_size)] - x0))
        if (np.abs(input_series[i] - x0) > n_sigmas * S0):
            new_series[i] = x0
            indices.append(i)
    
    return new_series, indices


# Format an xarray profile to work with function from this lib

def format_profil(Profil_test, with_bbp=False, n_prof = 0):
    time = Profil_test.JULD[n_prof].data
    lon = Profil_test.LONGITUDE[n_prof]
    lat = Profil_test.LATITUDE[n_prof]

    temp_a = Profil_test.TEMP_ADJUSTED[n_prof] # °C
    psal_a = Profil_test.PSAL_ADJUSTED[n_prof] # 0.001
    pres_a = Profil_test.PRES_ADJUSTED[n_prof] # dbar
    depth  = -gsw.z_from_p(pres_a, lat)

    SA_a = gsw.SA_from_SP(SP = psal_a, p = pres_a, lon = lon, lat = lat) # g/kg
    CT_a = gsw.CT_from_t(SA = SA_a, t = temp_a, p = pres_a) # °C
    SIG0_a = gsw.density.sigma0(SA = SA_a, CT = CT_a) # kg/m3

    SPI_a = gsw.spiciness0(SA_a, CT_a) # kg/m3
    O2sol_a = gsw.O2sol(SA = SA_a, CT = CT_a, p = pres_a, lon= lon, lat = lat) # umol/kg
    
    N2_a, p_mid = gsw.stability.Nsquared(SA_a, CT_a, pres_a, lat=lat)
    depth_mid = -gsw.z_from_p(p_mid, lat)
    
    ds = xr.Dataset(
        data_vars=dict(
            time = time,
            lon = lon.data,
            lat = lat.data,
            pres=(["depth"], pres_a.data),
            temp=(["depth"], temp_a.data),
            ctem=(["depth"], CT_a.data),
            psal=(["depth"], psal_a.data),
            asal=(["depth"], SA_a.data),
            sig0=(["depth"], SIG0_a.data),
            spic=(["depth"], SPI_a.data),
            osat=(["depth"], O2sol_a.data),
            N2  =(["depth_mid"], N2_a.data),
            ),
            coords={'depth': (['depth'], depth.data),
                    'depth_mid': (['depth_mid'], depth_mid.data)
                   })
    
    var_list = ['CHLA_ADJUSTED','BBP700_ADJUSTED','BBP532_ADJUSTED',
                'CDOM_ADJUSTED','CP660_ADJUSTED','NITRATE_ADJUSTED']
    var_name = ['chla','bbp7','bbp5','cdom','cp660','nitr']
    
    for var,name in zip(var_list,var_name):
        if var in list(Profil_test.variables):
            tmp = Profil_test[var][n_prof]
            ds[name]=(["depth"],  tmp.data)
    
    if 'CHLA_ADJUSTED' in list(Profil_test.variables):
        chla_fl_a, chla_fl_indices = hampel(ds.chla, 3)
        ds['chla_flt']=(["depth"],  chla_fl_a.data)
        
    if 'BBP700_ADJUSTED' in list(Profil_test.variables):
        bbp7_fl_a, bbp7_fl_indices = hampel(ds.bbp7, 3)
        ds['bbp7_flt']=(["depth"],  bbp7_fl_a.data)
        
        poc = 3.12*10**4*ds.bbp7+3.0
        ds['poc']=(["depth"],  poc.data)
        
        poc_flt = 3.12*10**4*ds.bbp7_flt+3.0
        ds['poc_flt']=(["depth"],  poc_flt.data)
        
    if 'DOXY_ADJUSTED' in list(Profil_test.variables):
        doxy_a = Profil_test.DOXY_ADJUSTED[n_prof] # umol/kg
        AOU_a = O2sol_a - doxy_a # umol/kg
        
        ds['doxy']=(["depth"], doxy_a.data)
        ds['aou'] =(["depth"], AOU_a.data)
        
    return ds

#########################################################################################################
def adjust_doxy_xlim(doxy, maou):
    daou = np.nanmax(maou)-np.nanmin(maou)
    ddox = np.nanmax(doxy)-np.nanmin(doxy)
    delta = abs(daou-ddox)/2
    
    if daou > ddox:
        xmin_aou, xmax_aou = np.nanmin(maou), np.nanmax(maou)
        xmin_dox, xmax_dox = np.nanmin(doxy)-delta, np.nanmax(doxy)+delta
        
    else:
        xmin_dox, xmax_dox = np.nanmin(doxy), np.nanmax(doxy)
        xmin_aou, xmax_aou = np.nanmin(maou)-delta, np.nanmax(maou)+delta
        
    return xmin_aou, xmax_aou, xmin_dox, xmax_dox

# Plot main varialbe profiles of an ARGO_profile

def profile_plot(prof, z_max = 1000):
    print(z_max)
    var_list = list(prof.variables)

    fig = plt.figure(figsize=(10,4), dpi = 200)
    
    date = str(prof.time.data)[0:10]
    lon = round(float(prof.lon.data),2)
    lat = round(float(prof.lat.data),2)
    
    fig.suptitle(date+' | Lon='+str(lon)+', Lat='+str(lat))
             
    ###
    ax1 = fig.add_subplot(141)
    if 'temp' in var_list:
        ax1.plot(prof.temp,prof.depth, '-o', c = 'peru', markersize=1)
        ax1.set_xlabel('Temperature [°C]')
        ax1.xaxis.label.set_color('peru')
        ax1.tick_params(axis='x', colors='peru')

    ax2 = ax1.twiny()
    if 'psal' in var_list:
        ax2.plot(prof.psal,prof.depth, '-o', c = 'teal', markersize=1)
        ax2.set_xlabel('Salinity [psu]')
        ax2.xaxis.label.set_color('teal')
        ax2.tick_params(axis='x', colors='teal')

    ax1.set_yticks(np.array(np.linspace(0,1000,21), dtype = int))
    ax1.set_ylim(z_max,0)
    ax1.yaxis.grid(color='gray', linestyle='dashed')

    ###
    ax3 = fig.add_subplot(142)
    if 'sig0' in var_list:
        ax3.plot(prof.sig0,prof.depth, '-o', c = 'k', markersize=1)
        ax3.set_xlabel('Density [kg/m3]')
        ax3.xaxis.label.set_color('k')
        ax3.tick_params(axis='x', colors='k')

    ax4 = ax3.twiny()
    if 'spic' in var_list:
        ax4.plot(prof.spic,prof.depth, '-o', c = 'grey', markersize=1)
        ax4.set_xlabel('Spiciness [kg/m3]')
        ax4.xaxis.label.set_color('grey')
        ax4.tick_params(axis='x', colors='grey')

    ax3.set_yticks(np.array(np.linspace(0,1000,21), dtype = int))
    ax3.set_ylim(z_max,0)
    ax3.yaxis.grid(color='gray', linestyle='dashed')

    ###
    ax5 = fig.add_subplot(143)
    if 'doxy' in var_list:
        ax5.plot(prof.doxy,prof.depth, '-o', c = 'red', markersize=1)
        #ax5.plot(prof.doxy+prof.aou,prof.depth, '-o', c = 'k', markersize=1)
        ax5.set_xlabel('Oxygen [umol/kg]')
        ax5.xaxis.label.set_color('red')
        ax5.tick_params(axis='x', colors='red')

    ax6 = ax5.twiny()
    if 'aou' in var_list:
        ax6.plot(-prof.aou, prof.depth, '-o', c = 'blue', markersize=1)
        ax6.set_xlabel('-AOU [umol/kg]')
        ax6.xaxis.label.set_color('blue')
        ax6.tick_params(axis='x', colors='blue')
    
    aou_xmin,aou_max,doxy_xmin,doxy_xmax = adjust_doxy_xlim(prof.doxy, -prof.aou)
    ax5.set_xlim(doxy_xmin,doxy_xmax), ax6.set_xlim(aou_xmin,aou_max)
    
    ax5.set_yticks(np.array(np.linspace(0,1000,21), dtype = int))
    ax5.set_ylim(z_max,0)
    ax5.yaxis.grid(color='gray', linestyle='dashed')

    ###
    ax7 = fig.add_subplot(144)
    if 'chla' in var_list:
        ax7.plot(prof.chla,prof.depth, '-o', c = 'green', markersize=1)
        ax7.set_xlabel('Chlorophyll [mg/m3]')
        ax7.xaxis.label.set_color('green')
        ax7.tick_params(axis='x', colors='green')

    ax8 = ax7.twiny()
    if 'bbp7' in var_list:
        ax8.plot(prof.bbp7,prof.depth, '-o', c = 'brown', markersize=1)
        ax8.set_xlabel('Backscatter [/m]')
        ax8.xaxis.label.set_color('brown')
        ax8.tick_params(axis='x', colors='brown')
    
    
    ax7.set_yticks(np.array(np.linspace(0,1000,21), dtype = int))
    ax7.set_ylim(z_max,0)
    ax7.yaxis.grid(color='gray', linestyle='dashed')

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    
    fig.show()

#########################################################################################################
# List all files from same float as input file

def prof_from_float(file):    
    prof = xr.open_dataset(data_path+file)
    float_nbr = int(prof.PLATFORM_NUMBER)
    float_numbers = pd.read_csv(index_file_path)['float_number']
    mask = np.array(float_numbers == float_nbr)
    float_files = np.array(pd.read_csv(index_file_path)['file'][mask])
    
    return float_files

#########################################################################################################
# Interpolate profile

def interp_var_1d(var_raw, depth_raw, depth_int):
    mask = np.isfinite(var_raw)
    var_raw = np.array(var_raw)[mask]
    depth_raw = np.array(depth_raw)[mask]
    if len(var_raw)>10:
        f_interp = interp1d(depth_raw, var_raw, bounds_error = False, fill_value = np.nan)
        prof_interp = f_interp(depth_int)
    else:
        prof_interp = np.full(len(depth_int), np.nan)
    return prof_interp


#########################################################################################################
# Create an array of all float data

def array_from_float(file):
    
    # Select all profiles from same float
    float_files = prof_from_float(file)
    
    float_prof = []
    for f in float_files:
        prof = xr.open_dataset(data_path+f)
        prof = format_profil(prof)
        float_prof.append(prof)
    
    # Interpolate all profiles on a same grid
    lats, lons, times, press = [], [], [], []
    temps, psals, sig0s, spics = [], [], [], []
    doxys, aous, chlas, bbp7s = [], [], [], []
    nitrs = []

    depth_int = np.linspace(1,1000,1000)

    for prof in float_prof:
        lats.append(prof.lat), lons.append(prof.lon)
        times.append(pd.to_datetime(np.array(prof.time)))

        depth_tmp = np.array(prof.depth)

        var_list = list(prof.variables)

        pres_int = interp_var_1d(prof.pres, prof.depth, depth_int)
        temp_int = interp_var_1d(prof.temp, prof.depth, depth_int)
        psal_int = interp_var_1d(prof.psal, prof.depth, depth_int)
        sig0_int = interp_var_1d(prof.sig0, prof.depth, depth_int)
        spic_int = interp_var_1d(prof.spic, prof.depth, depth_int)

        if 'doxy' in var_list:
            doxy_int = interp_var_1d(prof.doxy, prof.depth, depth_int)
        else:
            doxy_int = np.full(len(depth_int), np.nan)
        if 'aou' in var_list:
            aou_int = interp_var_1d(prof.aou, prof.depth, depth_int)
        else:
            aou_int = np.full(len(depth_int), np.nan)
        if 'chla' in var_list:
            chla_int = interp_var_1d(prof.chla, prof.depth, depth_int)
        else:
            chla_int = np.full(len(depth_int), np.nan)
        if 'bbp7' in var_list:
            bbp7_int = interp_var_1d(prof.bbp7, prof.depth, depth_int)    
        else:
            bbp7_int = np.full(len(depth_int), np.nan)
        if 'nitr' in var_list:
            nitr_int = interp_var_1d(prof.nitr, prof.depth, depth_int)
        else:
            nitr_int = np.full(len(depth_int), np.nan)
            
        press.append(pres_int), temps.append(temp_int), psals.append(psal_int), 
        sig0s.append(sig0_int), spics.append(spic_int), doxys.append(doxy_int), 
        aous.append(aou_int), chlas.append(chla_int), bbp7s.append(bbp7_int),
        nitrs.append(nitr_int)
    
    print(np.shape(np.array(nitrs)))
    # Store interpolated profiles in a dataset
    ds = xr.Dataset(
            data_vars=dict( lon = lons, lat = lats,
                pres=(["time","depth"], np.array(press),{'vname':'Pressure', 'units':'bdar'}),
                temp=(["time","depth"], np.array(temps),{'vname':'Temperature', 'units':'°C'}),
                psal=(["time","depth"], np.array(psals),{'vname':'Pra. Salinity', 'units':'psu'}),
                sig0=(["time","depth"], np.array(sig0s),{'vname':'Pot. Density', 'units':'kg/m3'}),
                spic=(["time","depth"], np.array(spics),{'vname':'Spiciness', 'units':'kg/m3'}),
                nitr=(["time","depth"], np.array(nitrs),{'vname':'Nitrate', 'units':'umol/kg'}),
                doxy=(["time","depth"], np.array(doxys),{'vname':'Diss. Oxygen', 'units':'umol/kg'}),
                aou =(["time","depth"], np.array(aous) ,{'vname':'Ap Oxy Util', 'units':'umol/kg'}),
                chla=(["time","depth"], np.array(chlas),{'vname':'Chlorophyll', 'units':'mg/m3'}),
                bbp7=(["time","depth"], np.array(bbp7s),{'vname':'Backscatter', 'units':'/m'}), ),
            coords={'depth': (['depth'], depth_int), 'time': (['time'], times)})
    
    return ds

##############################################################################################
# Plot float trajectory

def traj_from_float(file):
    # Select all profiles from same float
    float_files = prof_from_float(file)
    
    times, lons, lats = [], [], []
    for f in float_files:
        prof = xr.open_dataset(data_path+f)
        lat  = float(prof.LATITUDE[0].data)
        lon  = float(prof.LONGITUDE[0].data)
        time = prof.JULD[0].data
        times.append(time), lons.append(lon), lats.append(lat)
    times = pd.to_datetime(np.array(times))
    
    # Plot trajectory
    fig = plt.figure(figsize=(10,3), dpi = 200)

    #  Global plot
    fig_crs  = ccrs.PlateCarree(central_longitude = 0)
    ax = fig.add_subplot(121, projection=fig_crs)

    ax.add_feature(cfeature.LAND, zorder = 1),ax.add_feature(cfeature.COASTLINE, zorder = 1)
    ax.set_extent([1, 360, -90, 90], crs=fig_crs)

    plt.plot(lons,lats, c = 'k', linewidth = 0.5, zorder = 1, transform = fig_crs)
    c = plt.scatter(lons, lats, c = 'r', s= 1, zorder = 2, transform = fig_crs)

    ax.set_xticks([0, 60, 120, 180, 240, 300, 360], crs=ccrs.PlateCarree())
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())

    # Local plot
    ax = fig.add_subplot(122, projection=fig_crs)
    ax.add_feature(cfeature.LAND, zorder = 1),ax.add_feature(cfeature.COASTLINE, zorder = 1)

    x_min,x_max,y_min,y_max = int(np.nanmin(lons))-1, int(np.nanmax(lons))+1, int(np.nanmin(lats))-1, int(np.nanmax(lats))+1
    ax.set_extent([x_min,x_max,y_min,y_max], crs=fig_crs)

    ax.set_xticks(np.linspace(x_min, x_max, x_max-x_min+1), crs=ccrs.PlateCarree())
    ax.set_xticklabels(np.array(np.linspace(x_min, x_max, x_max-x_min+1), dtype=int), rotation = 90)
    ax.set_yticks(np.linspace(y_min, y_max, y_max-y_min+1), crs=ccrs.PlateCarree())


    plt.plot(lons,lats, c = 'k', linewidth = 1, zorder = 1, transform = fig_crs)
    c = plt.scatter(lons, lats, c = times, s = 3, zorder = 2, transform = fig_crs)

    cb = plt.colorbar(c)
    cb_tick = cb.get_ticks()

    ticklabels = []
    for t in cb_tick:
        t = t*10**-9
        y = str(dt.datetime.fromtimestamp(t).year)
        m = str(dt.datetime.fromtimestamp(t).month)
        m = (2-len(m))*'0'+m
        d = str(dt.datetime.fromtimestamp(t).day)
        d = (2-len(d))*'0'+d
        ticklabels.append(y+'-'+m+'-'+d)
    ticklabels = np.array(ticklabels)

    cb.set_ticks(cb_tick)
    cb.set_ticklabels(ticklabels)
    
    plt.tight_layout()
    
####################################################################################################
# Plot field of float containing

def fields_from_float(file, z_min = 0, z_max = 1000, vmax_chla = 1):
    ds = array_from_float(file)
    
    prof = xr.open_dataset(data_path+file)
    prof = format_profil(prof)
    time = pd.to_datetime(np.array(prof.time))
    
    var_list = list(ds.variables)
    for i in ['lon','lat','time','pres','depth']:
        var_list.remove(i)

    

    col_map = {'temp':plt.cm.Spectral_r, 'psal':plt.cm.nipy_spectral, 'sig0':cmo.cm.haline, 'spic':plt.cm.gist_rainbow,
               'doxy':plt.cm.gnuplot2, 'aou': plt.cm.gnuplot2, 'chla':plt.cm.gist_ncar_r, 'nitr': plt.cm.rainbow, 'bbp7':plt.cm.gist_ncar_r}

    mask_color = np.logical_and(ds.depth>z_min, ds.depth<z_max)
    
    for var in var_list:
        if np.isfinite(ds[var]).sum()>1:

            plt.figure(figsize=(10,2.5), dpi = 200)

            field = ds[var].T
            mask_nan = abs(field).sum('depth')>0
            xx,yy = np.meshgrid(ds.time[mask_nan], ds.depth)
            
            if var in ['chla']:
                #field = np.log10(field)
                vmax  = np.log10(np.nanmean(field) + 2*np.nanstd(field))
                c = plt.pcolormesh(ds.time[mask_nan], ds.depth, field[:,mask_nan], vmin = 0, 
                                   #vmax = (np.max(field)+np.mean(field))/2,
                                   vmax = vmax_chla,
                                   cmap = col_map[var])
            elif var in ['bbp7']:
                #field = np.log10(field)
                
                vmax  = np.log10(np.nanmean(field) + 2*np.nanstd(field))
                c = plt.pcolormesh(ds.time[mask_nan], ds.depth, field[:,mask_nan], vmin = 0, 
                                   vmax = 3*np.std(field)+np.mean(field), 
                                   cmap = col_map[var])
            else:
                c = plt.pcolormesh(ds.time[mask_nan], ds.depth, field[:,mask_nan], 
                                   vmin = np.nanmin(field[mask_color]), 
                                   vmax = np.nanmax(field[mask_color]), 
                                   cmap = col_map[var])

            plt.colorbar(c, label = ds[var].vname+' ['+ds[var].units+']')
            if not var in ['chla','bbp7']:
                CS = plt.contour(xx,yy,(ds[var].T)[:,mask_nan], colors = 'k', linewidths = 0.5)
                plt.clabel(CS, inline=1, fontsize=5, fmt = '%.1f')
            plt.xticks(rotation = 15), plt.ylim(z_max,z_min)
            
            plt.scatter(time, 0, c = 'r')
            plt.show()