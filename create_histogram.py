from xhistogram.xarray import histogram
import xarray
import os
import sys
import argparse, sys

parser=argparse.ArgumentParser()

parser.add_argument("--runnumber", help="specifiy runnumber, e.g. 100")
parser.add_argument("--calibrationrun", action='store_true', default=False,
    help="if calibrationrun, only output first few days")

args=parser.parse_args()
runnumber = int(args.runnumber)
calibrationrun = args.calibrationrun

if calibrationrun:
    print('running script with flag -calibrationrun=True')

print(f"running script for runnumber {runnumber}")
outfile = f"../data/output/opendrift_alternative_hourly_{runnumber}.nc"
bathy = xarray.open_dataset("../data/input/BAL-MFC_003_006_mask_bathy.nc")
bathy = bathy
coordinates = xarray.open_dataset("../data/input/BAL-MFC_003_006_coordinates.nc")
coordinates = coordinates
coordinates['volumes'] = coordinates.e3t*coordinates.e2t*coordinates.e1t
mask = bathy['mask'].values.astype(int)
coordinates['volumes'] = coordinates['volumes'].where(mask)

latbins = coordinates.latitude[:-1:]-(0.5*np.diff(coordinates.latitude))
latbins = np.append(latbins, coordinates.latitude[-2:].values+0.5*np.diff(coordinates.latitude)[-1])
lonbins = coordinates.longitude[:-1:]-0.5*np.diff(coordinates.longitude)
lonbins = np.append(lonbins, coordinates.longitude[-2:].values+0.5*np.diff(coordinates.longitude)[-1])
zbins = coordinates.depth[:-1:]-0.5*np.diff(coordinates.depth)
zbins = np.append(zbins, coordinates.depth[-2:].values+0.5*np.diff(coordinates.depth)[-1])

print(len(coordinates['depth']))
print(len(zbins))

dsm = xarray.open_mfdataset(outfile, chunks=dict(time=24))
print('reading in outfile...')
if calibrationrun:
    dsm = dsm.isel(time=slice(40,94,2), drop=True)
else:
    # the subsampling here is not strictly necessary, but prevents low-ram machines to
    # run out of memory...
	dsm = dsm.isel(time=slice(0,-1,24), drop=True)

weights = dsm.mass

coordinates = coordinates.rename(latitude='lat_bin', longitude='lon_bin', depth='z_bin')
coordinates = coordinates.transpose('lon_bin', 'lat_bin', 'z_bin')

h = histogram(dsm.lon, dsm.lat, -dsm.z, bins=[lonbins, latbins, zbins], dim=['trajectory'], weights=weights)
h = h.to_dataset(name='mass') # [mol/l]

print('creating histogram...')
h['mass_volatilized'] = histogram(dsm.lon, dsm.lat, -dsm.z, bins=[lonbins, latbins, zbins], dim=['trajectory'], weights=dsm.mass_volatilized)#*1e-6/(16.04*1e3))
h['mass_zsum'] = h.mass.sum(dim='z_bin')
h['mass_volatilized_zsum'] = h.mass_volatilized.sum(dim='z_bin')

print('writing netcdf')
if calibrationrun:
    h.to_netcdf(f'../data/output/{runnumber}_histogram_calibration.nc')
else:
    h.to_netcdf(f'../data/output/{runnumber}_histogram.nc',)
