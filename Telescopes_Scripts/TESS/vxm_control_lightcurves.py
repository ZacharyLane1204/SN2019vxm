import tessreduce as tr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from tqdm import tqdm
import os

from hidden_prints.hidden_prints import HiddenPrints

savepath = '/Users/zgl12/Modules/SN2019vxm/Telescopes_Scripts/TESS/Control_Curves/'

def magnitude(flux, zp = 20.023960740994795):
    m = -2.5 * np.log10(flux) + zp
    return m

def magnitude_error(f, df, dzp):
    return np.sqrt((-2.5 / np.log(10) * df/ (f))**2 + dzp**2)

def flux_to_jansky(flux, zp = 20.44):
    c = (zp - 8.9) / -2.5
    return flux * 1e6 * 10**c

def jansky_error(mag, mag_error):
    return np.sqrt((-0.4 * np.log(10) * 10**((8.9 - mag)/2.5))**2 * mag_error**2) * 1e6

def jansky_error_flux(flux, d_flux, zp, d_zp):
    
    a = 10**((8.9-zp)/2.5)
    term1 = d_flux * a
    term2 = flux * np.log(10) * a * 0.4 * d_zp
    return np.sqrt(term1**2 + term2**2)

def tess_reducing(ras, decs, sector):
    
    for i in tqdm(range(len(ras)), desc='Controls'):
        ra = ras[i]
        dec = decs[i]
        # try:
        # with HiddenPrints():
        tess_minor_reducing(ra, dec, sector, i)
        # except:
        #     print(f'Error with RA: {ra}, Dec: {dec}, Sector: {sector}')
        #     continue
    
def tess_minor_reducing(ra, dec, sector, i):
    """
    Function to reduce TESS data for a given RA and Dec.
    Parameters:
    ra (float): Right Ascension of the target.
    dec (float): Declination of the target.
    sector (int): TESS sector number.
    """   

    obs = tr.spacetime_lookup(ra, dec, time=58804, print_table=False)

    methods = ['psf']#, 'psf']

    obs_df = pd.DataFrame(obs, columns=['RA', 'Dec', 'Sector', 'Camera', 'CCD',  'Covers'])

    index = int(obs_df[obs_df['Sector'] == sector].index.values[0])
    
    for method in methods:

        saving_df = pd.DataFrame(columns=['MJD', 'flux', 'dflux', 'uJy', 'duJy', 'm', 'dm', 'ZP', 'd_ZP', 'Method', 'Sector'])
        
        if os.path.exists(f'{savepath}tess_sector_{sector}_{method}_calib_{i+1}.csv'):
            print('Existing File, Skipping...')
        elif os.path.exists(f'{savepath}tess_sector_{sector}_{method}_assume_{i+1}.csv'):
            print('Existing File, Skipping...')
        else:
        
            try:
                tess = None
                with HiddenPrints():
                    tess = tr.tessreduce(obs_list=obs[index], plot=False, phot_method=method, 
                                        sourcehunt=False, calibrate = True, verbose = False, num_cores=3)
                
                mag = magnitude(tess.lc[1], zp = tess.zp)
                mag_err = magnitude_error(tess.lc[1], tess.lc[2], tess.zp_e)
                
                zp = tess.zp
                zp_e = tess.zp_e
                calib = 'calib'
                
                mag[~np.isfinite(mag)] = 25
                
            except:
                tess = None
                with HiddenPrints():
                    tess = tr.tessreduce(obs_list=obs[index], plot=False, phot_method=method, 
                                sourcehunt=False, calibrate = False, verbose = False, num_cores=3)
        
                temp_df = pd.read_csv(savepath + f'tess_sector_18_aperture_calib_{i+1}.csv')
                zp = temp_df['ZP'].values[0]
                zp_e = temp_df['d_ZP'].values[0]
                
                mag = magnitude(tess.lc[1], zp = zp)
                mag_err = magnitude_error(tess.lc[1], tess.lc[2], zp_e)
                calib = 'assume'
                
                mag[~np.isfinite(mag)] = 25
            
            
            saving_df['MJD'] = tess.lc[0]
            saving_df['flux'] = tess.lc[1]
            saving_df['dflux'] = tess.lc[2]
            saving_df['uJy'] = flux_to_jansky(tess.lc[1], zp = zp)
            saving_df['duJy'] = jansky_error_flux(tess.lc[1], tess.lc[2], zp, zp_e) * 1e6
            saving_df['m'] = mag
            saving_df['dm'] = mag_err
            saving_df['ZP'] = zp
            saving_df['d_ZP'] = zp_e
            saving_df['Method'] = method
            saving_df['Sector'] = sector
                
            saving_df.to_csv(f'{savepath}tess_sector_{sector}_{method}_{calib}_{i+1}.csv', index=False)