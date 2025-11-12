import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from copy import deepcopy
from datetime import datetime
from astropy.time import Time
import emcee
import corner
from matplotlib.ticker import MultipleLocator, StrMethodFormatter
from scipy.optimize import curve_fit, minimize
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d

import warnings 

warnings.filterwarnings("ignore")

og_lc = np.load('/Users/zgl12/Modules/SN2019vxm/Data/tess_sector_18_psf.npy')
ind1 = 80
ind2 = 525


# df = pd.read_csv('/Users/zgl12/Modules/SN2019vxm/Data/tess_sector_18_psf_calib_vxm.csv')

# og_lc = np.row_stack((df.MJD.values, df.flux.values, df.dflux.values))

def residuals(theta, mjd, flux, flux_err):
    
    mjds = mjd[mjd >= theta[3]]
    fluxes = flux[mjd >= theta[3]]
    flux_errs = flux_err[mjd >= theta[3]]
    
    model = power_law(mjds, *theta)
    
    resid = np.nansum((abs(model - fluxes))**2 / flux_errs**2 ) / (len(fluxes) - len(theta))
    if (resid == 0) | (np.isnan(resid)):
        return np.nan
    else:
        return resid

def power_law(x, a, n, c, d): #, f):
    # d = x[0]
    # return a * (x-d)**(n*(x-d) + f) + c
    return a * (x-d)**(n) + c

def ln_likelihood(theta, mjd, flux):
    model = power_law(mjd, *theta)
    chi2 = np.nansum((model - flux)**2)# / flux_err**2)
    return -chi2/2.0

def check_lc_significance(flux, input_frame, flux_sign):

    buffer = 5
    
    frame_start = input_frame - buffer
    frame_end = input_frame + buffer
        
    frames = np.arange(len(flux))
    ind = (frames < frame_start) | (frames > frame_end)
    med = np.nanmedian(flux[ind])
    std = np.nanstd(flux[ind], ddof =1) 
    lcevent = flux[int(frame_start):int(frame_end)]
    
    # Light curve significance
    lc_sig = (lcevent - med) / std

    if flux_sign >= 0:
        sig_max = np.nanmax(lc_sig)
        sig_med = np.nanmean(lc_sig)
        
    else:
        sig_max = abs(np.nanmin(lc_sig))
        sig_med = abs(np.nanmean(lc_sig))
    
    lc_sig = (flux - med) / std
    return sig_max, sig_med, lc_sig * flux_sign

def weighted_value_and_uncertainty(data, uncertainties, mode = 'lc'):

    # Calculate weights as the inverse of the uncertainties squared
    weights = 1 / uncertainties**2

    # Calculate weighted mean
    if mode == 'lc':
        weighted_mean = np.nansum(weights * data) / np.nansum(weights)
        weighted_uncertainty = np.sqrt(1 / np.nansum(weights))
    elif mode == 'field':
        weighted_mean = np.nansum(weights * data, axis = 0) / np.nansum(weights, axis = 0)
        weighted_uncertainty = np.sqrt(1 / np.nansum(weights, axis = 0))

    return weighted_mean, weighted_uncertainty

def binned_averages(lc, exptime = '30', bin_size = 6, verbose = True, mode = 'lc'):

    if mode == 'lc':
        num_data_points = len(lc[0])
    elif mode == 'field':
        num_data_points = lc.shape[0]

    if verbose:
        print(f'Binning data in {bin_size} hour slots')
        if exptime == '200':
            print(f"Exposure time is {3.33} mins")
        else:
            print(f"Exposure time is {exptime} mins")
          
    if exptime == '30':
        points_per_bin = bin_size * 2
    elif exptime == '10':
        points_per_bin = bin_size * 6
    elif exptime == '200':
        points_per_bin = bin_size * 18
    
    num_bins = num_data_points // points_per_bin
    remainder = num_data_points % points_per_bin
    
    if remainder > 0:
        num_bins += 1

    blc, bins_median, bins_std = errors(lc, num_bins, points_per_bin, mode = mode)

    return blc, bins_median, bins_std

def errors(lc, num_bins, points_per_bin, mode = 'lc'):
    bins_median = []
    bins_std = []
    blc = []
    print('Number of bins:', num_bins)
    for i in range(num_bins):
        start = i * points_per_bin
        end = (i + 1) * points_per_bin

        if mode == 'lc':
            bin_lc = lc[:, start:end]
            blc.append(np.nanmean(bin_lc[0]))
            med, std = weighted_value_and_uncertainty(bin_lc[1], bin_lc[2], mode = mode)
            bins_median.append(med)
            bins_std.append(std)
        elif mode == 'field':
            bin_lc = lc[start:end]
            blc.append(None)
            med = np.nanmean(bin_lc, axis = 0)
            bins_median.append(med)
            bins_std.append(None)

    blc = np.array(blc)
    bins_median = np.array(bins_median)
    bins_std = np.array(bins_std)
    
    return blc, bins_median, bins_std

def gather_data(og_lc, ind1, ind2, ind3):    
    min_error = 1

    lc = deepcopy(og_lc)

    lc[2][lc[2] < min_error] = min_error

    lc[:,:ind1] = np.nan
    lc[:,ind2:ind3] = np.nan

    lc[1][ind1:ind2] -= np.nanmedian(lc[1][ind1:ind2])

    model = deepcopy(lc[:,ind3:])

    guess = [2.7, 1.4, -0.8, 58804]

    res = minimize(residuals, guess, args = (model[0], model[1], model[2]))
    
    mjds = model[0][model[0] >= res.x[3]]
    fluxes = model[1][model[0] >= res.x[3]]
    flux_errs = model[2][model[0] >= res.x[3]]
    
    popt, pcov = curve_fit(power_law, mjds, fluxes, p0 = res.x, sigma = flux_errs, absolute_sigma = True, maxfev = 1000)
    
    resid = residuals(popt, model[0], model[1], model[2])
    
    return popt, pcov, resid

def running_optimise():

    final_resids_df = pd.DataFrame(columns = ['index', 'residuals', 'a', 'n', 'c', 'd', 'd_a', 'd_n', 'd_c', 'd_d'])

    for ind in np.arange(555, 650):

        try:
            popt, pcov, resid = gather_data(og_lc, ind1, ind2, ind)
        except:
            continue
        
        final_resids_list = [ind, resid, *popt, *np.sqrt(np.diag(pcov))]
        
        final_resids_df.loc[len(final_resids_df)] = final_resids_list
        
    opt_val, add_uncert, min_ind = weighting(final_resids_df)
        
    return final_resids_df, opt_val, add_uncert, min_ind

def weighting(final_resids_df):
    ds= final_resids_df[abs(final_resids_df['residuals'].idxmin() - final_resids_df.index.values) < 10]['d'].values
    d_ds= final_resids_df[abs(final_resids_df['residuals'].idxmin() - final_resids_df.index.values) < 10]['d_d'].values

    d_ds = d_ds[np.isfinite(d_ds)]
    ds = ds[np.isfinite(ds)]

    opt_val, add_uncert = weighted_value_and_uncertainty(ds, d_ds, mode = 'lc')
    
    return opt_val, add_uncert, final_resids_df['residuals'].idxmin()

def recreating_lc(min_ind, min_error = 1):

    ind3 = 552

    popt, pcov, resid = gather_data(og_lc, ind1, ind2, ind3 = min_ind)

    lc = deepcopy(og_lc)

    model_mjds = np.linspace(og_lc[0][ind3], 58815, 1001)

    lc[2][lc[2] < min_error] = min_error

    lc[1:,:ind1] = np.nan
    lc[1:,ind2:ind3] = np.nan
    lc[1:,ind2:ind3] = np.nan

    lc[1][ind1:ind2] -= np.nanmedian(lc[1][ind1:ind2])

    lc[1][ind3:] -= popt[2]

    blc = binned_averages(lc, exptime = '30', bin_size = 6, verbose = True)
    
    np.save('../../Data/tess_binned_lc_18.npy', blc)
    
    return lc, blc, popt, pcov, 
    
def main(min_error = 1):
    final_resids_df, opt_val, add_uncert, min_ind = running_optimise()
    
    lc, blc, popt, pcov = recreating_lc(min_ind, min_error = min_error)
    
    return final_resids_df, opt_val, add_uncert, min_ind, lc, blc, popt, pcov