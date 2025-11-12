
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from glob import glob
import os
import requests
from bs4 import BeautifulSoup
from tqdm import tqdm
from copy import deepcopy
import warnings

import healpy as hp

from dynesty import NestedSampler
from dynesty.utils import resample_equal
from dynesty import plotting as dyplot

from matplotlib.ticker import MultipleLocator, StrMethodFormatter
from matplotlib.lines import Line2D

from scipy.interpolate import interp1d

from astropy.io import fits
from astropy.table import Table

from gdt.core import data_path
from gdt.missions.fermi.gbm.localization import GbmHealPix
from gdt.missions.fermi.gbm.phaii import GbmPhaii
from gdt.core.plot.lightcurve import Lightcurve
from gdt.core.binning.binned import rebin_by_time
from gdt.core.binning.unbinned import bin_by_time
from gdt.core.plot.spectrum import Spectrum
from gdt.missions.fermi.gbm.tte import GbmTte
from gdt.core.background.primitives import BackgroundRates
from gdt.core.background.fitter import BackgroundFitter
from gdt.core.background.binned import Polynomial

warnings.filterwarnings('ignore')


def download_fermi_files(download_folder, url):

    # url = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2019/bn191117006/current/"

    response = requests.get(url)

    if response.status_code == 200:
        soup = BeautifulSoup(response.text, 'html.parser')

        links = soup.find_all('a', href=True)
        fit_files = [link['href'] for link in links if link['href'].endswith('.fit')]

        for file in tqdm(fit_files, desc="Downloading files"):
            file_url = url + file
            file_path = download_folder + file
            # print(f"Downloading {file_url} to {file_path}")

            file_response = requests.get(file_url)

            if file_response.status_code == 200:
                with open(file_path, 'wb') as f:
                    f.write(file_response.content)
                # print(f"Downloaded: {file}")
            else:
                print(f"Failed to download {file}")
    else:
        print("Failed to retrieve the webpage")

class GDT_interactions():
    def __init__(self, folder, download_files = True, url = None, detector = None,
                 find_spatial_probability = True, ra = None, dec = None, 
                 bin_by_time = 1, time_ref = 0, bkgd_times = [(-100.0, -20.0), (40.0, 100.0)], 
                 poly_order_bkg = 1, time_interp = np.linspace(-100, 100, 101), erange = (8, 900.0),
                 src_time = (-2, 6.0), init_plot = False):
        
        self.bin_by_time = bin_by_time
        self.time_ref = time_ref
        self.bkgd_times = bkgd_times
        self.poly_order_bkg = poly_order_bkg
        self.time_interp = time_interp
        self.erange = erange
        self.src_time = src_time
        self.init_plot = init_plot
        
        if detector is None:
            raise ValueError("Detector needs to be given")
        else:
            self.detector = detector
        
        if folder[-1] != '/':
            folder += '/'
        
        # self.folder = folder
        
        if download_files:
            if url:
                download_fermi_files(folder, url)
            else:
                raise ValueError('Requires a URL for downloading files')
        
        
        healpix_file = glob(folder + '*healpix*')[0]
        self.midsection = healpix_file.split('.fit')[0].split('_')[3]
        
        self.healpix_file = healpix_file
        self.folder = folder
        
        if find_spatial_probability:
            if (ra is not None) & (dec is not None):
                self.src_prob = self.spatial_probability(ra, dec)
            else:
                raise ValueError("If location probabilty is desired, RA and DEC of target needs to be given")
        
    # def spatial_probability(self, ra, dec):
    #     loc = GbmHealPix.open(self.healpix_file)
    #     # src_prob = loc.source_probability(ra, dec, prior =0.5) # probability that the source is associated with that detector
        
    #     tol = 1e-3
    #     prior = 0.5
        
    #     prev = 3
        
    #     src_prob = None
        
    #     src_prob = loc.source_probability(ra, dec, prior = prior)
        
    #     print('Spatial Probability:', src_prob)
    #     return src_prob
    
    def spatial_probability(self, ra, dec):
        loc = GbmHealPix.open(self.healpix_file)

        # Get per-pixel probability at source location
        p = loc.probability(ra, dec, per_pixel=True)

        # Compute the null hypothesis likelihood (uniform over the sphere)
        u = 1.0 / (4.0 * np.pi)
        u *= hp.nside2resol(loc.nside)**2

        # Log-likelihood: P(ℐ | prior) = p*prior + u*(1 - prior)
        def log_likelihood(theta):
            prior = theta[0]
            if prior <= 0.0 or prior >= 1.0:
                return -np.inf  # Invalid prior, zero likelihood
            likelihood = p * prior + u * (1.0 - prior)
            if likelihood <= 0.0:
                return -np.inf  # Log safety
            return np.log(likelihood)

        # Uniform prior on prior ∈ [0, 1]
        def prior_transform(uv):
            return [uv[0]]  # uv ∈ [0,1] already

        # Run nested sampling
        sampler = NestedSampler(log_likelihood, prior_transform, ndim=1)
        sampler.run_nested(dlogz=1e-5)
        results = sampler.results

        # Resample posterior
        samples = resample_equal(results.samples, results.logwt - results.logz[-1])
        pi_samples = samples[:, 0]

        # Compute marginalized spatial probability
        probs = [loc.source_probability(ra, dec, prior=pi) for pi in pi_samples]
        posterior_mean = np.mean(probs)
        posterior_std = np.nanpercentile(probs, 16), np.nanpercentile(probs, 84)

        print('Spatial Probability (Bayesian posterior mean):', posterior_mean)
        print('Spatial Probability (Bayesian posterior median):', np.median(probs))
        print('Posterior stddev:', posterior_std)
        print('Posterior lower:', np.nanstd(probs, ddof = 1))

        return posterior_mean
    
    def tte_phaii_analysis(self):
        
        tte = GbmTte.open(self.folder + f'glg_tte_n{str(self.detector)}_{self.midsection}_v00.fit')
        # tte = 
        phaii = tte.to_phaii(bin_by_time, self.bin_by_time, time_ref=self.time_ref)
        
        fitter = BackgroundFitter.from_phaii(phaii, Polynomial, time_ranges=self.bkgd_times)

        fitter.fit(order=self.poly_order_bkg)

        back_rates = fitter.interpolate_bins(self.time_interp[:-1], self.time_interp[1:])
        
        if self.init_plot:
            lcplot = Lightcurve(data=phaii.to_lightcurve(energy_range=self.erange), 
                                background=back_rates.integrate_energy(*self.erange))
        

        src_lc = phaii.to_lightcurve(time_range=self.src_time, energy_range=self.erange)
        
        lc_bin_centers = phaii.to_lightcurve(energy_range=self.erange).centroids
        lc_bin_widths = phaii.to_lightcurve(energy_range=self.erange).widths
        lc_flux = phaii.to_lightcurve(energy_range=self.erange).rates
        lc_uncert = phaii.to_lightcurve(energy_range=self.erange).rate_uncertainty

        back_time = back_rates.integrate_energy(*self.erange).time_centroids
        back_ratee = back_rates.integrate_energy(*self.erange).rates

        m = (back_ratee[1] - back_ratee[0]) / (back_time[1] - back_time[0])

        f = interp1d(back_time, back_ratee, kind='linear', fill_value='extrapolate')

        spectra_bin_centers = phaii.to_spectrum(time_range=self.src_time, energy_range=self.erange).centroids
        spectra_bin_widths = phaii.to_spectrum(time_range=self.src_time, energy_range=self.erange).widths
        spectra_flux = phaii.to_spectrum(time_range=self.src_time, energy_range=self.erange).rates_per_kev
        spectra_uncert = phaii.to_spectrum(time_range=self.src_time, energy_range=self.erange).rate_uncertainty_per_kev

        spec_bkgd = back_rates.integrate_time(*self.src_time)

        spectra_back_time = spec_bkgd.centroids
        spectra_back_ratee = spec_bkgd.rates_per_kev
        
        corrected_lc_flux = lc_flux - m*lc_bin_centers - f(0)
        
        indices = np.where(np.isin(spectra_back_time, spectra_bin_centers))[0]
        
        corrected_spectra_flux = spectra_flux - spectra_back_ratee[indices]
        
        return lc_bin_centers, corrected_lc_flux, lc_uncert, spectra_bin_centers, corrected_spectra_flux, spectra_uncert

