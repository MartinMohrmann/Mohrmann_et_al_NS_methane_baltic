import datetime
import numpy as np

outfile_base='../data/output/opendrift_alternative_hourly_194'
speedrun = False

resolution = 'hourly' # 'hourly'
n_total_seeds = 100000
time_step = 3600
time_skip = 4 # data is written into outputfile every nth timestep. Default is 1

# mode = 'cmems_mod_bal_phy_anfc_PT1h-i'
mode = 'dataset-bal-analysis-forecast-phy-hourly'
Winds='../data/input/ERA5_winds_expversum.nc'
timedependence = 'CATHARE' # PBREAK, linear_decrease, 1_over_t, CATHARE
total_mass_max = 10.4e6 # [kg]
end_time = datetime.datetime(2023,1,2)
radius = 5000
vertical_drift_mixing = False
degradation = False
volatilization = True
degradation_rate = 0.1*4.0509259259259263e-07 # from Dissayanake et al 2023, in 1/s
dissayanake_depth_distribution = False # if True, place 98.5 % of methane in surface layer

def extranormalization(seedfunction):
    return seedfunction

