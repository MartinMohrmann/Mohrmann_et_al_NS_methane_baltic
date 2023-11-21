import cmocean
import matplotlib
import cartopy
import pandas
import cartopy.crs as ccrs
import numpy as np
from matplotlib import pyplot as plt
import glidertools as gt
import xarray
import dictionaries
import bathy_maps.bathyutils as bathyutils
import dictionaries
import matplotlib.patheffects as pe


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def plot_concentrations_map(fig, ax, bathy=True, cbar=True, labels=True, zoomin=False, gtrajectory=True,
    gridlines=True, biglabels=True):
    import matplotlib.ticker as mticker

    if zoomin:
        S_lim = 54.7
        N_lim = 56.6
        W_lim = 12.5
        E_lim = 17.5
    else:
        S_lim = dictionaries.S_lim
        N_lim = 56.6
        #N_lim = dictionaries.N_lim#57 # dictionaries.N_lim # 57
        W_lim = dictionaries.W_lim
        E_lim = dictionaries.E_lim

    leak1 = dictionaries.leak1
    leak2 = dictionaries.leak2
    leak3 = dictionaries.leak3
    leak4 = dictionaries.leak4

    pc = cartopy.crs.PlateCarree()
    ax.set_facecolor('#00022e')
    ax.set_extent((W_lim, E_lim, S_lim, N_lim), crs=pc)

    feature = cartopy.feature.NaturalEarthFeature(name='land', category='physical',
                                           scale='10m', edgecolor='black', facecolor='burlywood',
                                           zorder=0)


    for leak in [leak1, leak2, leak3, leak4]:
        ax.scatter(leak['lon'], leak['lat'], color='orange', transform=pc, s=150, zorder=20,
                    edgecolor='k', marker="X")


    if bathy:
        emodnet_path = "/home/coffee/VOTO/data/bathymetry/emodnet_baltic/"
        extent = (S_lim, 60, W_lim-1, E_lim+1)
        ds_bathy = bathyutils.emod_subset(extent, emodnet_path)

        stride = 10
        lon = ds_bathy.lon[::stride]
        lat = ds_bathy.lat[::stride]
        bathy = ds_bathy.elevation[::stride, ::stride].values
        ds_bathy.close()
        # set all land to 0 m, we only want bathymetry
        bathy[bathy>0] = 0
        bathy = bathy.clip(-100)

        levels = np.arange(-100, 10, 10)
        bathy_c = ax.contourf(lon, lat, bathy, transform=pc,
                            cmap=truncate_colormap(cmocean.cm.deep_r, 0.0, 0.6),
                            alpha=1, zorder=1,
                            levels=levels,
                            linestyles='-')
        if cbar:
        	fig.colorbar(bathy_c,ax=ax,
                orientation='horizontal',
                label='Bathymetry (m)',
                shrink=0.6,
                extend='min')

    ax.add_feature(feature)
    ax.add_feature(cartopy.feature.LAKES, alpha=1, zorder=11)
    gl = ax.gridlines(draw_labels=True,
                      linewidth=1, alpha=1, color='darkgrey', linestyle='--', zorder=40)
    gl.top_labels = None
    gl.right_labels = None
    if not labels:
    	gl.left_labels = None

    if gtrajectory:
	    ax.plot([15.7,16.3], [55.66, 55.55],
		    transform=ccrs.PlateCarree(),
		    color='white',
		    zorder=5,
		    lw=4,
		    path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()],
		    ls='-')

    if zoomin:
        gl.xlocator = mticker.FixedLocator([13, 17])
    else:
        gl.xlocator = mticker.FixedLocator([10, 20])

    if biglabels:
        gl.ylabel_style = {'size': 18, 'color':'k'}
        gl.xlabels_top = True
        gl.xlabels_bottom = False
        gl.xlabel_style = {'size':18, 'color':'k', 'rotation':0}

    if gridlines:
        gl2 = ax.gridlines(draw_labels=False,
                       linewidth=1, alpha=1, color='darkgrey', linestyle='--', zorder=40)

    return pc, ax, gl


def plot_helcom(ax, df_helcom):
    """"""
    acolors = []
    for area in range(0,len(df_helcom)):
        acolors.append(np.random.rand(3,))

    crs = ccrs.AzimuthalEquidistant()
    crs_proj4 = crs.proj4_init
    df_helcom = df_helcom.to_crs(crs_proj4)
    ax.add_geometries(df_helcom.sort_values('Shape_Area', ascending=False).geometry,
                      crs=crs,
                      zorder=100,
                      edgecolor='k',
                      linewidth=2,
                      #color=acolors,
                      facecolor=(1,1,1,0))
    return ax

def plot_pipelines(ax):
    import geopandas

    nybb = geopandas.read_file('/home/coffee/VOTO/data/Nordstream/_ags_Pipelines_2018/Pipelines_2018.shp')

    alltraj = nybb.to_crs({'proj':'longlat', 'ellps':'WGS84', 'datum':'WGS84'}).geometry
    alltraj = alltraj.to_list()
    Nordstream1 = alltraj[4]
    Nordstream2 = alltraj[5]
    N2lons = np.array(Nordstream2.coords.xy[0])
    N2lats = np.array(Nordstream2.coords.xy[1])
    mask = ((N2lons<14.6) | (N2lons>15.6))
    N1lons = Nordstream1.geoms[1].coords.xy[0]
    N1lats = Nordstream1.geoms[1].coords.xy[1]

    ax.plot(np.where(mask, N2lons, np.nan),
            np.where(mask, N2lats, np.nan),
            transform=ccrs.PlateCarree(),
            color='darkgrey',
            zorder=1,
            lw=5,
            alpha=1)

    ax.plot([15.775656, 15.80 ,15.41, 15.3, 14.5],
            [55.506924, 55.3,54.876667, 54.8, 54.564],
            transform=ccrs.PlateCarree(),
            color='darkgrey',
            zorder=1,
            lw=7,
            alpha=1)

    ax.plot(N1lons, N1lats,
            transform=ccrs.PlateCarree(),
            color='darkgrey',
            zorder=1,
            lw=7,
            alpha=1)


def plot_glider_data_tail(ax, means, enddate, vmin, vmax):
    means4 = means.where((means.time<enddate) & (means.time>enddate-np.timedelta64(6,'D')))

    if not np.isnan(means4['latitude']).all():
        mcolor = ax.scatter(
            means4['longitude'],
            means4['latitude'],
            c='k',
            transform=ccrs.PlateCarree(),
            s=20,
            zorder=16,
            cmap='Reds')

        mcolor = ax.scatter(
            means4['longitude'],
            means4['latitude'],
            c=means4['Methane_tempcorr']*1e-3,
            transform=ccrs.PlateCarree(),
            norm=matplotlib.colors.LogNorm(
               vmin=vmin,
               vmax=vmax),
            s=15,
            zorder=16,
            cmap='Reds')

        for marker in ["$0x25B2$",
                       "$\u2B2E$"]:
            ax.scatter(
                means4['longitude'].dropna().iloc[-1], # maxs4['longitude'][-1],
                means4['latitude'].dropna().iloc[-1], # maxs4['latitude'][-1],
                #c='tab:blue',
                transform=ccrs.PlateCarree(),
                s=300,
                zorder=101,
                marker=marker,
                color='yellow',
                edgecolors='k')

def plot_column_averaged_methane_concentrations(fig, ax, time, showmodel=True):
    import dictionaries
    from download_glider_data import utils
    import methane_utils as mutils
    from xhistogram.xarray import histogram
    from matplotlib import pyplot as plt
    import xarray
    from scipy.signal import argrelextrema
    import numpy as np
    import string
    import os
    import matplotlib.dates as mdates
    os.environ['HDF5_USE_FILE_LOCKING']='FALSE'

    dataset_ids = dictionaries.timelines['NS_methane_sensor']
    dataset_ids = ['nrt_'+dataset_id for dataset_id in dataset_ids]
    datasets_dict = utils.download_glider_dataset(dataset_ids)
    datasets_list = [gt.load.voto_seaexplorer_dataset(ds) for ds in [*datasets_dict.values()]]
    ds = gt.load.voto_concat_datasets(datasets_list)
    groups = gt.utils.group_by_profiles(ds)
    meantimes = groups.mean(numeric_only=False)['time']
    mld = gt.physics.mixed_layer_depth(ds, 'potential_density', thresh=0.03, verbose=False)
    dsg = ds

    datagaps = []
    for index in range(0,len(datasets_list))[:-1]:
        time1 = datasets_list[index].isel(time=-1).time.values
        time2 = datasets_list[index+1].isel(time=0).time.values
        datagaps.append((time1, time2))

    groups = gt.physics.group_by_profiles(dsg)
    means = groups.mean(numeric_only=False)

    dsg['background'] = mutils.compute_methane_equilibrium(dsg, mode='surface')*1e-3
    dsg['Methane_tempcorr'] = mutils.methane_temperature_correction(
        dsg.mets_raw_temperature, dsg.methane_raw_concentration)*1e-3

    groups = dsg[['profile_num', 'temperature', 'time', 'fdom','latitude', 'longitude',
                    'depth', 'methane_concentration', 'profile_direction', 'background',
                    'Methane_tempcorr']].to_pandas().reset_index().groupby('profile_num')
    means = groups.mean(numeric_only=False)

    order = 50
    maxs = argrelextrema(means.latitude.values, np.greater, order=order)
    mins = argrelextrema(means.latitude.values, np.less, order=order)

    maxs = list(maxs[0])
    mins = list(mins[0])

    maxs.insert(0, 5)

    glons = means['longitude'].to_xarray()
    glats = means['latitude'].to_xarray()
    gtims = means['time'].to_xarray()

    outfile = '/media/coffee/T7/Nordstream/data/output/opendrift_alternative_hourly_100.nc'
    bathy = xarray.open_dataset("../data/input/BAL-MFC_003_006_mask_bathy.nc")
    bathy = bathy
    coordinates = xarray.open_dataset("../data/input/BAL-MFC_003_006_coordinates.nc")
    coordinates = coordinates
    coordinates['volumes'] = coordinates.e3t*coordinates.e2t*coordinates.e1t
    mask = bathy['mask'].values.astype(int)
    np.shape(coordinates['volumes'].values)
    coordinates['volumes'] = coordinates['volumes'].where(mask)

    latbins = coordinates.latitude[:-1:]-(0.5*np.diff(coordinates.latitude))
    latbins = np.append(latbins, coordinates.latitude[-2:].values+0.5*np.diff(coordinates.latitude)[-1])
    lonbins = coordinates.longitude[:-1:]-0.5*np.diff(coordinates.longitude)
    lonbins = np.append(lonbins, coordinates.longitude[-2:].values+0.5*np.diff(coordinates.longitude)[-1])
    zbins   = np.array([-100, -60, 0])

    # decrease amount of data to Glider sampling area
    lonbins = lonbins[(lonbins>14.9) & (lonbins<16.6)]
    latbins = latbins[(latbins>54.9) & (latbins<56.1)]

    dsm = xarray.open_mfdataset(outfile, chunks='auto')
    dsm = dsm.isel(time=slice(0,-1,1))

    glider_max_depth = -60
    dsm = dsm.where((dsm.z>glider_max_depth) & (dsm.z<0))
    weights = dsm.mass.where((dsm.z>glider_max_depth) & (dsm.z<0), other=0.0)
    h = histogram(dsm.lon, dsm.lat,  bins=[lonbins, latbins, ], dim=['trajectory'], weights=weights)

    h2 = h/(coordinates['volumes'].sum(dim='depth')
        ).rename(latitude='lat_bin', longitude='lon_bin') # [ug m-3]
    h2 = h2.sel(lon_bin=slice(14.9,16.6), lat_bin=slice(54.9,56.1))
    h2 = h2*1e-6/(16.04*1e3) # [mol/l]
    h2 = h2.rolling(lon_bin=3, lat_bin=3, center=True).mean()

    result = h2.interp(
        time=gtims,
        lon_bin=glons,
        lat_bin=glats,
        method='linear')

    # compute distance from the leak sites:
    def great_circle(lon1, lat1, lon2, lat2):
        lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
        return 6371 * (
            np.arccos(np.sin(lat1) * np.sin(lat2) + np.cos(lat1)
                    * np.cos(lat2) * np.cos(lon1 - lon2))
        )
    distance = great_circle(dictionaries.leak3['lon'], dictionaries.leak3['lat'], means.longitude, means.latitude)
    import datetime
    means = means.where(means.time<time)
    ax.plot(means['time'],
                means['Methane_tempcorr'],
                c='k', lw=2,
                alpha=0.7, label='glider obs.')

    endpoint = means.last_valid_index()
    if endpoint:
        # if there is already some data at the specified time,
        # e.g. time != 0
        ax.scatter(means.loc[endpoint]['time'],
               means.loc[endpoint]['Methane_tempcorr'],
               c='k', zorder=10, lw=2,
               alpha=1, )

    ax.set_ylim(2e-10,2e-5)
    ax.set_xlim(datetime.datetime(2022,9,26), datetime.datetime(2023,1,3))

    ax.set_yscale('log')
    ax.set_ylabel('[mol l⁻¹]', fontsize=18)
    ax.set_yticks([1e-9, 1e-8, 1e-7, 1e-6, 1e-5])
    ax.tick_params(axis='both', which='major', labelsize=18)

    if showmodel:
        result = result.swap_dims({'profile_num':'time'})
        endpoint = result.sel(time=time, method='nearest')#.values
        if endpoint:
            ax.scatter(endpoint.time,
               endpoint.values,
               c='darkred', zorder=10, lw=2,
               alpha=1, )
        result = result.where(result.time<time)
        result = result
        ax.plot(result.time,
                result.where(result.values>means.background.values, other=means.background.values),
                color='tab:red',
                lw=2, label='num. model')


    for datagap in datagaps:
        ax.axvspan(datagap[0], datagap[1], color='lightgrey', zorder=99, alpha=1)

    ax.grid(zorder=100)
    ax.plot(means.time, means.background, color='grey', alpha=0.8, label='back. conc.')

    ax.legend(fontsize=14, loc='upper right', framealpha=1., handlelength=1).set_zorder(100)

    locator = mdates.AutoDateLocator(minticks=7, maxticks=14)
    formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)

    ax.text(x=0.01, y=0.08, s='column averaged methane concentration', bbox=dict(facecolor='white', alpha=1),
        fontsize=16,
        transform=ax.transAxes,
        zorder=100)

    ax.set_axisbelow(False)
    return fig, ax


def load_ferry_data():
    df_ferry = pandas.read_csv(
    '../data/input/IOW_Ferry_box/IOW_SOOP_Finnmaid_CH4_data_Rehder_Bittig_Glockzin_for_VOTO/'
    'IOW_SOOP_Finnmaid_CH4_data_Rehder_Bittig_Glockzin_for_DLR.txt')
    df_ferry['time'] = pandas.to_datetime(df_ferry['MatlabTime']-719529, unit='D')
    df_ferry = df_ferry.dropna()
    df_ferry = df_ferry.drop(['MatlabTime', 'ExcelTime', 'Day', 'Month', 'Year', 'Hour', 'Minute'], axis=1)
    return df_ferry

def load_glider_data():
    from download_glider_data import utils
    import methane_utils as mutils
    import gsw

    dataset_ids = dictionaries.timelines['NS_methane_sensor']
    dataset_ids = ['nrt_'+dataset_id for dataset_id in dataset_ids]
    datasets_dict = utils.download_glider_dataset(dataset_ids,
        variables=['profile_num', 'temperature', 'salinity', 'mets_raw_temperature',
        'oxygen_concentration', 'mets_raw_concentration', 'longitude', 'latitude'])
    datasets_list = [gt.load.voto_seaexplorer_dataset(ds) for ds in [*datasets_dict.values()]]

    ds = gt.load.voto_concat_datasets(datasets_list)
    ds['Methane_tempcorr'] = mutils.methane_temperature_correction_post_calibration(
        ds.mets_raw_temperature, ds.methane_raw_concentration) # in mol/m3
    ds['oxy_conc'] = ds['oxygen_concentration']/gsw.O2sol(
        ds['salinity'], ds['temperature'],
        ds['pressure'], ds['longitude'],
        ds['latitude'])
    ds['Methane_corr'] = mutils.methane_oxygen_correction_post_calibration(
         ds['Methane_tempcorr'],
         ds['oxy_conc'])
    groups = gt.utils.group_by_profiles(ds)
    dsg = groups.min(numeric_only=False)
    dsg = dsg.to_xarray()
    return dsg