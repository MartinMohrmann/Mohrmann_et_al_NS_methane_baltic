import numpy as np

# molar mass of methane
m_CH4 = 16.05  # [g/l]


def compute_solubility_H(T, units='original'):
    """
    Computes Henry’s law solubility constant H
    ------------------------------------------

    Parameters:
    T (float), temperature in K
    units (str), eiter 'original' or 'converted'
    returns the result either in [mol m-3 Pa-1] or in [mg m-3 Pa-1]

    Source: https://acp.copernicus.org/articles/15/4399/2015/acp-15-4399-2015.pdf
    """
    H_cp = 1.4e-5  # mol m-3 Pa-1
    # 14nmol Pa-1 kg-1 = 14*(e-9*e3) = 14e-6
    HR = 1900
    H = H_cp*np.exp(HR*(1/T-1/298.15))
    # if units == 'converted':
    #     H = H*16.05*(1000/0.001) # Mol weight, milligram, m3
    print('compute_solubility_H returning solubillity in [mol m-3 Pa-1]')
    return H


def nnmol_to_glider_units(value):
    """ conversion from nmol to mg m-3"""
    return value*1e-9*m_CH4*(1000/0.001)  # mol*e-9 l-1 * 16.05 g/mol


def glider_units_to_nnmol(value):
    """ conversion from nmol to mg m-3"""
    return value*1e9*(0.001/1000)/m_CH4


def compute_methane_equilibrium(ds, mode):
    """"
    Parameter:
    ----------
    T (float), temperature in K

    Return:
    -------
    molar concentration [mol m-3] of Methane in the water in equilibrium with
    typical atmospheric partial methane pressure.
    """
    T = ds.temperature+273.15
    p_atm = 101325  # [Pa]
    if mode == 'surface':
        x_i = 1900/1e9  # atmospheric Methane mole fraction ca 1900 ppb
        p_i = x_i*p_atm  # Partial pressure of the methane in air
    elif mode == 'in-situ':
        P = ds.pressure*10000  # [Pa]
        p_i = P+p_atm
    # [mol m-3 Pa-1] * [Pa] = [mol m-3]
    conc_mol = compute_solubility_H(T)*p_i
    print("""compute_methane_equilibrium returning molar
        concentration in [mol m-3]""")
    return conc_mol


def methane_temperature_correction(meth_tv, meth_mv):
    """"
    Parameters:
    -----------
    conc_CH4 (numpy.array), methane concentration in [umol/l]
    conc_O2 (nump.array), oxygen concentration in [%]
    """
    conc_CH4 = np.exp(1.667 * np.log(
        (0.083+1.212*np.exp(-meth_tv/0.592)) * (
            (1/meth_mv) - (1/(-20.356+26.263*np.exp(-meth_tv/9.241))))
        ))  # methane concentration in [umol/l]
    return conc_CH4*(1e-6/1e-3)


def methane_oxygen_correction(conc_CH4, conc_O2):
    """"
    corrects the methane values for oxygen concentrations,
    numerical values in formula are sensor-specific and not
    adjusted yet.

    Parameters:
    -----------
    conc_CH4 (numpy.array), methane concentration in [umol/l]
    conc_O2 (nump.array), oxygen concentration in [%]

    Return
    ------
    conc_CH4 (numpy.array), methane concentration in [mol/m-3]
    """
    conc_CH4 = (-0.326+0.360*np.exp(conc_O2/77.779))*conc_CH4  # in [umol/l]
    return conc_CH4*(1e-6/1e-3)


def methane_temperature_correction_post_calibration(meth_tv, meth_mv):
    """"
    Parameters:
    -----------
    conc_CH4 (numpy.array), methane concentration in [umol/l]
    conc_O2 (nump.array), oxygen concentration in [%]
    """
    conc_CH4 = np.exp(1.709 * np.log(
        (0.138+1.645*np.exp(-meth_tv/0.546)) * (
            (1/meth_mv) - (1/(-8.549+15.274*np.exp(-meth_tv/4.589))))
        ))  # methane concentration in [umol/l]
    return conc_CH4


def mets_temperature_from_voltage_post_calibration(meth_tv):
    """
    Parameters:
    -----------
    meth_tv (numpy.array), methane temperature voltage

    Return:
    -------
    gas temperature of methane sensor in Degree Celsius
    """
    t = (meth_tv*21.28)-4.77
    return t


def methane_oxygen_correction_post_calibration(conc_CH4, conc_O2):
    """"
    corrects the methane values for oxygen concentrations,
    numerical values in formula are sensor-specific and not
    adjusted yet.

    Parameters:
    -----------
    conc_CH4 (numpy.array), methane concentration in [umol/l]
    conc_O2 (nump.array), oxygen concentration in [%]

    Return
    ------
    conc_CH4 (numpy.array), methane concentration in [mol/m-3]
    """
    conc_CH4 = (-0.426+0.450*np.exp(conc_O2/0.872))*conc_CH4  # in [umol/l]
    return conc_CH4


def compute_oxygen_from_oxygen_frequency(dsg):
    import gsw

    # Nominal parameters
    # IDENTICAL FOR ALL SENSORS
    A_1=-173.4292
    A_2=249.6339
    A_3=143.3483
    A_4=-21.8492
    B_1=-0.033096
    B_2=0.014259
    B_3=-0.0017

    # SENSOR SPECIFIC CALIBRATION FACTORS
    # Parameters that can be adjusted in post-processing
    Soc_SEA070_M13 = 3.3343e-004
    Foffset_SEA070_M13 = -855.82
    E_SEA070_M13 = 0.036

    # fixed parameters
    A_SEA070_M13 = -4.6474e-003
    B_SEA070_M13 = 1.6144e-004
    C_SEA070_M13 = -2.2489e-006

    # Oxsol_SEA070_M13= np.exp(A_1+A_2.*(100./(data.data['temperature'].values+273.15)) + A_3 * np.log((data.data['temperature']+273.15)./100)+A_4.*(SEA070_M13_PLD_GPCTD_TEMPERATURE+273.15)./100+SEA070_M13_PLD_PSAL.*(B_1+B_2.*(SEA070_M13_PLD_GPCTD_TEMPERATURE+273.15)./100+B_3.*((SEA070_M13_PLD_GPCTD_TEMPERATURE+273.15)./100).^2));
    o2_sol = gsw.O2sol(
        dsg.salinity,
        dsg.temperature,
        dsg.pressure,
        dsg.longitude,
        dsg.latitude)

    PHI_SEA070_M13 = o2_sol * (
        1 + A_SEA070_M13*dsg.temperature + B_SEA070_M13*dsg.temperature **2 +
        C_SEA070_M13*dsg.temperature**3) * np.exp(E_SEA070_M13 * dsg.pressure /
        (dsg.temperature+273.15))

    o2 = Soc_SEA070_M13 * ( dsg.oxygen_frequency + Foffset_SEA070_M13) * PHI_SEA070_M13  #*44660./ (SEA070_M13_PLD_DENS);
    o2_sat = o2*100/o2_sol

    bad = (o2_sat < -150) | (o2_sat > 150)

    o2[bad] = np.NaN
    o2_sat[bad] = np.NaN

    dsg['o2'] = o2
    #dsg['o2_sat'] = o2_sat
    if 'oxygen_concentration' in dsg.variables:
        dsg['oxygen_concentration'] = (['time'], np.where(~np.isnan(dsg.o2), dsg.o2, dsg.oxygen_concentration))
    else:
        # e.g. datasets of the first three glider deployments with oxygen_frequency only
        dsg['oxygen_concentration'] = (['time'], dsg.o2.data)
    dsg['oxygen_saturation'] = (['time'], (100*dsg.oxygen_concentration/o2_sol).data)
    return dsg
