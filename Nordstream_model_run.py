# instructions how to run:
# change inputparametersXX to run number
# change file inputparametersXX.py to include right outfile_base number
# run python Nordstream_model_run.py

import xarray
import matplotlib.pyplot as plt
import datetime
import matplotlib
import opendrift
import os
import pandas
import numpy as np
import inputparameters_tryout as ip
os.environ['HDF5_USE_FILE_LOCKING']='FALSE'


leak1 = dict(lat=54.876667, lon=15.41, starttime=datetime.datetime(2022,9,26,2)) # NS2A
NS2A = pandas.read_excel('../data/input/NS2 Comparison PBREAK CATHARE_reprocessed.xlsx')
leak1['data'] = NS2A

leak2 = dict(lat=55.535000, lon=15.698333, starttime=datetime.datetime(2022,9,26,19)) # NS1A
NS1A = pandas.read_excel('../data/input/NS1 A Comparison PBREAK CATHARE_reprocessed.xlsx')
leak2['data'] = NS1A

leak3 = dict(lat=55.556667, lon=15.788333, starttime=datetime.datetime(2022,9,26,19)) # NS1B
NS1B = pandas.read_excel('../data/input/NS1 B Comparison PBREAK CATHARE_reprocessed.xlsx')
leak3['data'] = NS1B


# monkeypatch for the chemicaldrift library
def volatilization_custom(self):
    if self.get_config('chemical:transformations:volatilization') is True:
        volatilized_now = np.zeros(self.num_elements_active())
        MolWtCO2 = 44
        MolWtH2O = 18
        MolWt = self.get_config('chemical:transformations:MolWt')
        Undiss_n = 1

        R = 8.206e-05  # (atm m3)/(mol K)
        mixedlayerdepth = self.environment.ocean_mixed_layer_thickness

        # mask of dissolved elements within mixed layer
        W = (self.elements.specie == self.num_lmm) \
            * (-self.elements.z <= mixedlayerdepth)

        mixedlayerdepth = mixedlayerdepth[W]

        T = self.environment.sea_water_temperature[W]
        # temporary fix for missing values
        T[T == 0] = np.median(T)
        S = self.environment.sea_water_salinity[W]
        wind = (self.environment.x_wind[W]
                ** 2 + self.environment.y_wind[W]**2)**.5

        def gas_transfer_velocity(T, S, U):
            A1=1909.4; B1=-120.78; C1=4.1555; D1=-0.080578; E1=0.00065777 # for freshwater ~0PSU
            A2=2101.2; B2=-131.54; C2=4.4931; D2=-0.08676; E2=0.00070663  # for seawater ~35PSU
            # t=means.temperature # Replace this with the mixed layer temperature eventually
            Sc1 = A1 + B1*T + C1*T**2 + D1*T**3 + E1*T**4 #(t in °C). # freshwater
            Sc2 = A2 + B2*T + C2*T**2 + D2*T**3 + E2*T**4 #(t in °C). # seawater
            Sc = Sc1 + (Sc1-Sc2)*S/35
            gtv = 0.251 * wind**2 * (Sc/660)**(-0.5)
            return gtv # [cm h-1]

        def compute_solubility_H(T, units='original'):
            """Environmental Organic Chemistry, John Wiley & Sons, 2016
            Computes Henry’s law solubility constant H
            ------------------------------------------

            Parameters:
            T (float), temperature in K
            units (str), eiter 'original' or 'converted'
            returns the result either in [mol m-3 Pa-1] or in [mg m-3 Pa-1]

            Source: https://acp.copernicus.org/articles/15/4399/2015/acp-15-4399-2015.pdf
            """
            H_cp = 1.4e-5
            HR = 1900
            H = H_cp*np.exp(HR*(1/T-1/298.15)) # [mol m-3 Pa-1]
            return H**-1
            # See sanders table (2015) in papers for
            # explanation, in particular
            # 2.5.1 The Henry volatility defined via concentration pc (KH )

        Henry = compute_solubility_H(
            T+273.15, units='original')  # ((Vp * self.tempcorr("Arrhenius", DH_Vp, T, Tref_Vp)))   \
        #/ (Slb * self.tempcorr("Arrhenius", DH_Slb, T, Tref_Slb))  \
        #* MolWt / 101325.    # atm m3 mol-1

        #k_S_fin = k_S_tot * self.tempcorr("Arrhenius",DH_kSt,TS,Tref_kSt)
        # Calculate mass transfer coefficient water side
        # Schwarzenbach et al., 2016 Eq.(19-20)
        MTCw = ((9e-4)+(7.2e-6*wind**3)) * \
            (MolWtCO2/MolWt)**0.25 / Undiss_n

        # Calculate mass transfer coefficient air side
        # Schwarzenbach et al., 2016 Eq.(19-17)(19-18)(19-19)

        # More complex
        Sca_H2O = 0.62
        MTCaH2O = 0.1 + wind*(6.1+0.63*wind)**0.5 \
            / (13.3*(Sca_H2O)**0.5 + (6.1e-4+(6.3e-5)*wind)**-0.5 - 5 + 1.25*np.log(Sca_H2O))

        MTCa = MTCaH2O * (MolWtH2O/MolWt)**(1/3)

        # Calculate overall volatilization mass tansfer coefficient

        HenryLaw = Henry * (1 + 0.01143 * S) / (R * (T+273.15))
        MTCvol = 1 / (1/MTCw + 1/(MTCa * HenryLaw))

        mixedlayerdepth = self.environment.ocean_mixed_layer_thickness
        # mask of dissolved elements within mixed layer
        W = (self.elements.specie == self.num_lmm) \
            * (-self.elements.z <= mixedlayerdepth)

        mixedlayerdepth = mixedlayerdepth[W]
        K_volatilization = 0.01 * MTCvol / mixedlayerdepth  # (1/s)
        volatilized_now[W] = self.elements.mass[W] * \
            (1-np.exp(-K_volatilization * self.time_step.seconds))

        self.elements.mass_volatilized = self.elements.mass_volatilized + volatilized_now
        self.elements.mass = self.elements.mass - volatilized_now

    else:
        pass

def degradation_custom(self):
    """degradation."""
    if ip.degradation:
        # logger.debug('Calculating overall degradation using overall rate constants')
        degraded_now = np.zeros(self.num_elements_active())
        # Degradation in the water
        k_W_tot = ip.degradation_rate
        W =   (self.elements.specie == self.num_lmm)

        k_W_fin = k_W_tot
        degraded_now[W] = self.elements.mass[W] * (1-np.exp(-k_W_fin * self.time_step.total_seconds()))
        self.elements.mass_degraded_water[W] = self.elements.mass_degraded_water[W] + degraded_now[W]
        self.elements.mass_degraded = self.elements.mass_degraded + degraded_now
        self.elements.mass = self.elements.mass - degraded_now
    else:
        pass


def horizontal_diffusion(self):
    """Move elements with random walk according to given horizontal diffuivity."""
    D_mld = self.get_config('drift:horizontal_diffusivity')
    D_deep = 0.00 # both mld and deep values from https://doi.org/10.5194/gmd-2021-101
    if (D_mld == 0) & (D_deep == 0):
        logger.debug('Horizontal diffusivity is 0, no random walk.')
        return
    dt = np.abs(self.time_step.total_seconds())
    mixedlayerdepth = self.environment.ocean_mixed_layer_thickness
    # mask of dissolved elements within mixed layer
    W = (self.elements.specie == self.num_lmm) \
        * (-self.elements.z <= mixedlayerdepth)

    D = np.where(W, 5.0, 0.01)
    x_vel = self.elements.moving * np.sqrt(2*D/dt) * np.random.normal(
        scale=1, size=len(self.elements.moving))
    y_vel = self.elements.moving * np.sqrt(2*D/dt) * np.random.normal(
        scale=1, size=len(self.elements.moving))
    speed = np.sqrt(x_vel * x_vel + y_vel * y_vel)
    self.update_positions(x_vel, y_vel)
    if ip.degradation:
        degradation_custom(self)

# load bathymetry
bathy = xarray.open_dataset("../data/input/BAL-MFC_003_006_mask_bathy.nc")
bathy = bathy.sel(longitude=slice(10,20), latitude=slice(53,60.02), drop=True)

# depth of the leaks:
for leak in [leak1, leak2, leak3]:
    print(bathy.sel(latitude=leak['lat'], longitude=leak['lon'], method='nearest').deptho.values)

resolution = ip.resolution
n_total_seeds = ip.n_total_seeds

time_step = ip.time_step
mode = ip.mode

if mode == 'dataset-bal-analysis-forecast-phy-hourly':
    model_currents_file='../data/input/BAL_currents_temperature_salinity/out_hourly_alternative*.nc'
    model_ts_file = '../data/input/BAL_currents_temperature_salinity/alt_out_so_thetao*.nc'
elif mode == 'cmems_mod_bal_phy_anfc_PT1h-i':
    model_currents_file = '../data/input/NEMO_alternative/out_hourly_alternative_*'
    model_ts_file = '../data/input/NEMO_alternative/out_hourly_alternative_*'

outfile_base = ip.outfile_base
outfile = outfile_base+'.nc'

Winds = ip.Winds
print(model_currents_file, n_total_seeds, outfile)

speedrun = ip.speedrun # False
tot_mass_max = ip.total_mass_max
# print(f'release in total {tot_mass_max/1000} tons of material into the water AT EACH leak site')
# print(f'each methane seed has a mass of {tot_mass_max/(70*n_total_seeds)} kg') # (four leaks)

# Physical modeling part
from opendrift.models.oceandrift import OceanDrift # used this before
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.chemicaldrift import ChemicalDrift # now using this one
ChemicalDrift.volatilization_custom = volatilization_custom
ChemicalDrift.degradatio_custom = degradation_custom
ChemicalDrift.horizontal_diffusion = horizontal_diffusion # three monkey patches
from opendrift.models.physics_methods import verticaldiffusivity_Sundby1983, verticaldiffusivity_Large1994, verticaldiffusivity_stepfunction

o = ChemicalDrift(loglevel=20)
readBornholm = reader_netCDF_CF_generic.Reader(
    model_currents_file, standard_name_mapping={
            'w': 'upward_sea_water_velocity',
            'ocean_mixed_layer_thickness': 'mlotst'
    })
readBornholm.Dataset = readBornholm.Dataset.drop_duplicates(dim='time')

readTSBornholm = reader_netCDF_CF_generic.Reader(
    model_ts_file, standard_name_mapping={
        'so': 'sea_water_salinity',
        'thetao': 'sea_water_temperature'
    })

readERA5Wind = reader_netCDF_CF_generic.Reader(
    Winds, standard_name_mapping={'u10': 'x_wind',
                                  'v10': 'y_wind',})

readTSBornholm.Dataset['time'] = pandas.to_datetime(readTSBornholm.Dataset.time, unit='D', origin=np.datetime64('1900-01-01'))
readBornholm.Dataset['time'] = pandas.to_datetime(readBornholm.Dataset.time, unit='D', origin=np.datetime64('1900-01-01'))
readERA5Wind.Dataset['time'] = pandas.to_datetime(readERA5Wind.Dataset.time, unit='h', origin=np.datetime64('1900-01-01'))

bathy = reader_netCDF_CF_generic.Reader('../data/input/BAL-MFC_003_006_mask_bathy.nc')
end_time = ip.end_time

o.set_config('environment:fallback:y_wind', 5) # could be usefull if ERA5 is incomplete
o.set_config('environment:fallback:x_wind', 5)

if speedrun:
    print('using simpliedified model: speedrun (neglect diffusion, vertical mixing and more)')
    pass
else:
    print('using complete model: normal mode')
    o.set_config('drift:horizontal_diffusivity', 5) # I believe this was specified somewhere in the NEMO documentation. Search for source.
    o.set_config('drift:vertical_mixing', ip.vertical_drift_mixing)
    o.set_config('vertical_mixing:diffusivitymodel', 'windspeed_Large1994')
    o.set_config('chemical:transformations:degradation', False) # ip.degradation) # replaced by degration_custom
    o.set_config('chemical:transformations:volatilization', ip.volatilization)

    o._add_config({'chemical:transformations:MolWt': {'type': 'float', 'default': 128.1705,
                    'min': 10, 'max': 1000, 'units': 'amu',
                    'level': 3,
                    'description': 'molecular weight'}})

    o.set_config('chemical:transformations:MolWt', 16.04)
    o.set_config('seed:wind_drift_factor', False)
    # See Compilation of Henry’s law constants (version 4.0) for water as solvent, Chapt 2.5.1 for eplanation
    o.set_config('chemical:transformations:Henry', 1/(1.4e-5)*0.00000986923) # factor is pa -> atm conversion

o.add_reader([readBornholm, readTSBornholm, bathy, readERA5Wind])

o.set_config('general:coastline_action', 'previous')
print('latest date of model currents is {}'.format(readBornholm.end_time.isoformat()))
print('model end time is set to {}'.format(end_time))

# running
o.set_config('seed:particle_fraction',0.)
o.set_config('general:coastline_action', 'previous')
o.set_config('seed:LMM_fraction',1.)

#partitioning
o.set_config('chemical:transfer_setup','organics')
o.set_config('chemical:transformations:dissociation','nondiss')
o.set_config('chemical:transformations:LogKOW',0.)
o.set_config('chemical:transformations:TrefKOW',0.)
o.set_config('chemical:transformations:DeltaH_KOC_Sed',0.)
o.set_config('chemical:transformations:DeltaH_KOC_DOM',0.)
o.set_config('chemical:transformations:Setchenow', 0.)

# initial release impementations (some antiexponential release patterns...)
# We went with the CATHARE in the end
timedependence = ip.timedependence
leaks = [leak1, leak2, leak3]
for leakindex,leak in enumerate(leaks):
    print(f'seeding for leak {leakindex}')
    if timedependence == '1_over_t':
        seedendtime = datetime.datetime(2022,12,31)
        t = pandas.date_range(start=leak['starttime'], end=seedendtime,
                periods=(seedendtime-leak['starttime'])/datetime.timedelta(seconds=3600)+1)
        xt = np.linspace(1,len(t), len(t)) # unit hours
        seedfunction = 1/xt # 1/x

    elif timedependence == 'linear_decrease':
        seedendtime = datetime.datetime(2022,10,1,23)
        t = pandas.date_range(start=leak['starttime'], end=seedendtime,
                periods=(seedendtime-leak['starttime'])/datetime.timedelta(seconds=3600)+1)
        xt = np.linspace(1,len(t), len(t))
        seedfunction = (max(xt)-xt) # linear decrease

    elif timedependence == 'PBREAK':
        data = leak['data'][['Time (hr)', 'Time (s)', 'Total Flow Methane (t/hr)']].dropna()
        ds = xarray.Dataset(
            data_vars=dict(
                methane=(["time"], data['Total Flow Methane (t/hr)'])
            ),
            coords=dict(
                time=data['Time (hr)']
            )
        )
        seedendtime = leak['starttime'] + datetime.timedelta(hours=float(ds['time'].max()))
        seedendtime = pandas.to_datetime(seedendtime).round('h')
        t = pandas.date_range(start=leak['starttime'], end=seedendtime, periods=(seedendtime-leak['starttime'])/datetime.timedelta(seconds=3600))
        xt = np.linspace(1,len(t), len(t))
        seedfunction = ds['methane'].sel(time=xt, method='nearest')

    elif timedependence == 'CATHARE':
        data = leak['data'][['Time (hr).1', 'Time (s).1', 'Total Flow Methane (t/hr).1']].dropna()
        ds = xarray.Dataset(
            data_vars=dict(
                methane=(["time"], data['Total Flow Methane (t/hr).1'])
            ),
            coords=dict(
                time=data['Time (hr).1']
            )
        )
        seedendtime = leak['starttime'] + datetime.timedelta(hours=float(ds['time'].max()))
        seedendtime = pandas.to_datetime(seedendtime).round('h')
        t = pandas.date_range(start=leak['starttime'], end=seedendtime, periods=(seedendtime-leak['starttime'])/datetime.timedelta(seconds=3600)+1)
        xt = np.linspace(1,len(t), len(t))
        seedfunction = ds['methane'].sel(time=xt, method='nearest')

    # Function defined, proceed with normalisation
    seedfunction = ip.extranormalization(seedfunction)
    integral = np.trapz(y=seedfunction, x=xt)
    seedfunction = (seedfunction/integral) # normalized funtion over time
    for tindex,seedtime in enumerate(t[:-1]):
        # compute amount of seeds to release in this timestep as a fraction of n_total_seeds
        number = n_total_seeds/len(leaks)*seedfunction[tindex]
        number = int(np.ceil(number))
        mass = tot_mass_max/n_total_seeds
        if ip.dissayanake_depth_distribution:
            o.seed_elements(lat=leak['lat'], lon=leak['lon'],
                            number=int(np.ceil(0.985*number)), radius=ip.radius, radius_type='gaussian',
                            time=seedtime,#[seedtime, seedtime+np.timedelta64(1,'h')],
                            z=np.random.randint(-30,0,int(np.ceil(0.985*number))),
                            mass=mass*1e9,
                            origin_marker=leakindex)
            o.seed_elements(lat=leak['lat'], lon=leak['lon'],
                            number=int(np.ceil(0.015*number)),
                            radius=ip.radius,
                            radius_type='gaussian',
                            time=seedtime,#[seedtime, seedtime+np.timedelta64(1,'h')],
                            z=np.random.randint(-70,-30,int(np.ceil(0.015*number))),
                            mass=mass*1e9,
                            origin_marker=leakindex)
        else:
            o.seed_elements(lat=leak['lat'], lon=leak['lon'], number=number, radius=ip.radius, radius_type='gaussian',
                time=seedtime,#[seedtime, seedtime+np.timedelta64(1,'h')],
                z=np.random.randint(-70,0,number),
                mass=mass*1e9,
                origin_marker=leakindex)

if os.path.isfile(outfile):
    os.remove(outfile)
o.run(time_step=time_step,
      time_step_output=ip.time_skip*time_step,
      export_variables=['status', 'moving', 'age_seconds', 'origin_marker',
                       'lon', 'lat', 'z', 'mass', 'mass_volatilized', 'ocean_mixed_layer_thickness',
                       'mass_degraded'],
      export_buffer_length=10,
      end_time=end_time,
      outfile=outfile)

