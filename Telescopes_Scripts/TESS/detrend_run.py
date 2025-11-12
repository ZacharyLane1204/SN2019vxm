import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from glob import glob
from copy import deepcopy
from tqdm import tqdm
import warnings

from scipy.optimize import curve_fit

warnings.filterwarnings('ignore')

init_path = '/Users/zgl12/Modules/SN2019vxm/Telescopes_Scripts/TESS/Control_Curves/'
savepath = '/Users/zgl12/Modules/SN2019vxm/Telescopes_Scripts/TESS/Control_Curves/Detrended_Data/'

def linear_model(x, a, c):
    return a * x + c

def quadratic_model(x, a, b, c):
    return a*x**2 + b*x + c

def magnitude(flux, zp):
    m = -2.5 * np.log10(flux) + zp
    return m

def magnitude_error(f, df, dzp):
    return np.sqrt((-2.5 / np.log(10) * df/ (f))**2 + dzp**2)

def flux_to_jansky(flux, zp):
    c = (zp - 8.9) / -2.5
    return flux * 1e6 * 10**c

def jansky_error(mag, mag_error):
    return np.sqrt((-0.4 * np.log(10) * 10**((8.9 - mag)/2.5))**2 * mag_error**2) * 1e6

def jansky_error_flux(flux, d_flux, zp, d_zp):
    
    a = 10**((8.9-zp)/2.5)
    term1 = d_flux * a
    term2 = flux * np.log(10) * a * 0.4 * d_zp
    return np.sqrt(term1**2 + term2**2) * 1e6

def file_proc(file):
    
    df = pd.read_csv(file)

    df = df[(df['flux'] != 0) & ((df['m'] != 15))].reset_index(drop=True)
    
    zp = df['ZP'].values[0]
    zp_e = df['d_ZP'].values[0]

    diffs = np.abs(np.diff(df.MJD))
    indices = np.where(diffs > 0.5)[0]

    ind1 = 80
    ind2 = indices[0] - 30
    ind3 = indices[0]
    ind4 = df.index.values[-120]
    
    return df, zp, zp_e, ind1, ind2, ind3, ind4
    
def indexing_nan(df, ind1, ind2, ind3, ind4):
    mjds = deepcopy(df.MJD.values)
    fluxes = deepcopy(df.flux.values)
    d_fluxes = deepcopy(df.dflux.values)
    
    
    mjds[:ind1] = np.nan
    mjds[ind2:ind3+1] = np.nan
    mjds[ind4:] = np.nan

    fluxes[:ind1] = np.nan
    fluxes[ind2:ind3+1] = np.nan
    fluxes[ind4:] = np.nan

    d_fluxes[:ind1] = np.nan
    d_fluxes[ind2:ind3+1] = np.nan
    d_fluxes[ind4:] = np.nan
    
    return mjds, fluxes, d_fluxes

def fitting_procedure_lin(mjds, fluxes, d_fluxes, ind1, ind2, ind3, ind4):

    guess_lin = [1, 2.1e+04]

    popt_lin_1, _ = curve_fit(linear_model, mjds[ind1:ind2], fluxes[ind1:ind2], sigma = d_fluxes[ind1:ind2], 
                                    p0 = guess_lin, absolute_sigma=True)
    popt_lin_2, _ = curve_fit(linear_model, mjds[ind3+1:ind4-1], fluxes[ind3+1:ind4-1], sigma = d_fluxes[ind3+1:ind4-1], 
                                    p0 = guess_lin, absolute_sigma=True)
    
    # mjds[ind1:ind2], fluxes[ind1:ind2] - linear_model(mjds[ind1:ind2], *popt_lin_1)
    # mjds[ind3+1:ind4-1], fluxes[ind3+1:ind4-1] - linear_model(mjds[ind3+1:ind4-1], *popt_lin_2)
    
    mjd_linear = np.array(mjds[ind1:ind2].tolist() + mjds[ind3+1:ind4-1].tolist())
    flux_linear = np.array((fluxes[ind1:ind2] - linear_model(mjds[ind1:ind2], *popt_lin_1) ).tolist() + 
                           (fluxes[ind3+1:ind4-1] - linear_model(mjds[ind3+1:ind4-1], *popt_lin_2)).tolist())
    d_fluxes_linear = np.array(d_fluxes[ind1:ind2].tolist() + d_fluxes[ind3+1:ind4-1].tolist())
    
    return mjd_linear, flux_linear, d_fluxes_linear

def fitting_procedure_qua(mjds, fluxes, d_fluxes, ind1, ind2, ind3, ind4):

    guess_qua = [1, 1, 2.1e+04]

    popt_qua_1, _ = curve_fit(quadratic_model, mjds[ind1:ind2], fluxes[ind1:ind2], sigma = d_fluxes[ind1:ind2], 
                                    p0 = guess_qua, absolute_sigma=True)
    popt_qua_2, _ = curve_fit(quadratic_model, mjds[ind3+1:ind4-1], fluxes[ind3+1:ind4-1], sigma = d_fluxes[ind3+1:ind4-1], 
                                    p0 = guess_qua, absolute_sigma=True)
    
    # mjds[ind1:ind2], fluxes[ind1:ind2] - quadratic_model(mjds[ind1:ind2], *popt_qua_1)
    # mjds[ind3+1:ind4-1], fluxes[ind3+1:ind4-1] - quadratic_model(mjds[ind3+1:ind4-1], *popt_qua_2)
    
    mjd_quad = np.array(mjds[ind1:ind2].tolist() + mjds[ind3+1:ind4-1].tolist())
    flux_quad = np.array((fluxes[ind1:ind2] - quadratic_model(mjds[ind1:ind2], *popt_qua_1) ).tolist() + 
                           (fluxes[ind3+1:ind4-1] - quadratic_model(mjds[ind3+1:ind4-1], *popt_qua_2)).tolist())
    d_fluxes_quad = np.array(d_fluxes[ind1:ind2].tolist() + d_fluxes[ind3+1:ind4-1].tolist())
    
    return mjd_quad, flux_quad, d_fluxes_quad

def file_info(file):
    number = file.split('_')[-1].split('.')[0]
    sector = file.split('_sector_')[1].split('_')[0]
    method = file.split('_sector_')[1].split('_')[1]
    
    return number, sector, method

def creating_df(file, mjd, flux, d_flux, zp, d_zp, detrend = 'linear'):
    
    number, sector, method = file_info(file)

    saving_df = pd.DataFrame(columns=['MJD', 'flux', 'd_flux', 'uJy', 'd_uJy', 'm', 'd_m', 'ZP', 'd_ZP', 'Method', 'Sector'])
    
    uJy = flux_to_jansky(flux, zp)
    d_uJy = jansky_error_flux(flux, d_flux, zp, d_zp)
    m = magnitude(flux, zp)
    d_m = magnitude_error(flux, d_flux, d_zp)
    
    saving_df['MJD'] = mjd
    saving_df['flux'] = flux
    saving_df['d_flux'] = d_flux
    saving_df['uJy'] = uJy
    saving_df['d_uJy'] = d_uJy
    saving_df['m'] = m
    saving_df['d_m'] = d_m
    saving_df['ZP'] = zp
    saving_df['d_ZP'] = d_zp
    saving_df['Method'] = method
    saving_df['Sector'] = sector
    
    saving_df.to_csv(f'{savepath}tess_sector_{sector}_{method}_{detrend}_{number}.csv', index=False)
        
def medianing(mjds, fluxes, d_fluxes, ind1, ind2, ind3, ind4):
    
    mjd_linear = np.array(mjds[ind1:ind2].tolist() + mjds[ind3+1:ind4-1].tolist())
    flux_linear = np.array((fluxes[ind1:ind2] -  
                            weighted_value_and_uncertainty(fluxes[ind1:ind2], d_fluxes[ind1:ind2])).tolist() + 
                           (fluxes[ind3+1:ind4-1] - 
                            weighted_value_and_uncertainty(fluxes[ind3+1:ind4-1], d_fluxes[ind3+1:ind4-1])).tolist())
    d_fluxes_linear = np.array(d_fluxes[ind1:ind2].tolist() + d_fluxes[ind3+1:ind4-1].tolist())
    
    return mjd_linear, flux_linear, d_fluxes_linear

def weighted_value_and_uncertainty(data, uncertainties):

    # Calculate weights as the inverse of the uncertainties squared
    weights = 1 / uncertainties**2

    weighted_mean = np.nansum(weights * data) / np.nansum(weights)
    weighted_uncertainty = np.sqrt(1 / np.nansum(weights))

    return weighted_mean

def main(file):
    
    df, zp, zp_e, ind1, ind2, ind3, ind4 = file_proc(file)
    mjds, fluxes, d_fluxes = indexing_nan(df, ind1, ind2, ind3, ind4)
    
    d_fluxes[d_fluxes < 1] = 1
    
    # if 'psf' in file:
    #     # mjd_linear, flux_linear, d_fluxes_linear = medianing(mjds, fluxes, d_fluxes, ind1, ind2, ind3, ind4)
    #     mjd_linear, flux_linear, d_fluxes_linear = fitting_procedure_lin(mjds, fluxes, d_fluxes, ind1, ind2, ind3, ind4)
    #     creating_df(file, mjd_linear, flux_linear, d_fluxes_linear, zp, zp_e, detrend = 'std')
    # else:
    
    mjd_linear, flux_linear, d_fluxes_linear = fitting_procedure_lin(mjds, fluxes, d_fluxes, ind1, ind2, ind3, ind4)
    mjd_quad, flux_quad, d_fluxes_quad = fitting_procedure_qua(mjds, fluxes, d_fluxes, ind1, ind2, ind3, ind4)
    
    creating_df(file, mjd_linear, flux_linear, d_fluxes_linear, zp, zp_e, detrend = 'linear')
    creating_df(file, mjd_quad, flux_quad, d_fluxes_quad, zp, zp_e, detrend = 'quad')
    
    
files = sorted(glob(init_path + '*17*.csv'))

# files = ['/Users/zgl12/Modules/SN2019vxm/Data/tess_sector_17_aperture_calib_vxm.csv', '/Users/zgl12/Modules/SN2019vxm/Data/tess_sector_17_psf_calib_vxm.csv']

for file in tqdm(files, desc = 'Detrending...'):
    main(file)
    # break
    
    
    
