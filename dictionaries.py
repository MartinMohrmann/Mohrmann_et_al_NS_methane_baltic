from cmocean import cm as cmo
import datetime

timelines = dict(
    # Nordstream timelines
    NS1=['SEA077_M11',
         'SEA069_M11', 'SEA077_M13', 'SEA045_M69',
         'SEA077_M15', 'SEA045_M71', 'SEA077_M17', 'SEA045_M73'],  # This is the SAMBA trajectory
    NS2=['SEA045_M67', 'SEA077_M13'],
    NS3=['SEA077_M12'],
    NS4=['SEA076_M8', 'SEA076_M9'],
    NS5=['SEA070_M13', 'SEA070_M14', 'SEA070_M15'],
    NS_methane_sensor=['SEA070_M13', 'SEA070_M14',
                       'SEA070_M15', 'SEA056_M54',
                       'SEA056_M55', 'SEA056_M56', 'SEA056_M57'],)


cmap_dict = dict(
    conservative_temperature=cmo.thermal,
    potential_density=cmo.dense,
    temperature=cmo.thermal,
    salinity=cmo.haline,
    backscatter=cmo.turbid,
    cdom=cmo.matter,
    fdom=cmo.haline,
    chlorophyll=cmo.algae,
    oxygen_concentration=cmo.amp,
    N2=cmo.balance,  # cmo.amp,
    spice=cmo.matter,
    temperature_oxygen=cmo.thermal,
    turbidity=cmo.turbid,
    profile_num=cmo.haline,
    methane_concentration=cmo.thermal,
    methane_raw_concentration=cmo.thermal
    )

limits_dict_below_mld = dict(
    temperature=[3, 11],
    salinity=[8, 15],
    N2=[-0.010, 0.010],
    potential_density=[1006, 1011],
    oxygen_concentration=[0, 200],
    chlorophyll=[0, 0.5],
    cdom=[0, 10],
    turbidity=[0.1, 0.2],
    backscatter=[0, 0.0006])

limits_dict_mld = dict(
    temperature=[15, 19],
    salinity=[7, 10],
    N2=[-0.010, 0.010],
    potential_density=[1003, 1005],
    oxygen_concentration=[250, 350],
    chlorophyll=[0.4, 0.7],
    cdom=[0, 10],
    turbidity=[0.16, 0.23],
    backscatter=[0, 0.0006])

leak1 = dict(lat=54.876667, lon=15.41, starttime=datetime.datetime(2022,9,26,2))
leak2 = dict(lat=55.535000, lon=15.698333, starttime=datetime.datetime(2022,9,26,19))
leak3 = dict(lat=55.556667, lon=15.788333, starttime=datetime.datetime(2022,9,26,19))
leak4 = dict(lat=55.557500, lon=15.779000, starttime=datetime.datetime(2022,9,26,19))

levels = [1e-8, 1e-7, 1e-6, 1e-5, 1e-3]
hatches = ['''\\''', '''////''', '''xx''', '''xxxx''']
spacing = [9, 9, 9, 99] # this has to be adjusted to the levels distance
hatchesdicts = [dict(vmin=levels[i],
                     vmax=levels[i+1],
                     hatch=hatches[i],
                     spacing=spacing[i]) for i in range(0,len(hatches))]

S_lim = 53.5
N_lim = 60
W_lim = 10
E_lim = 20
