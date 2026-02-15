 
import matplotlib.pyplot as plt 
import numpy as np
import os

from gdt.core.background.binned import Polynomial
from gdt.core.background.fitter import BackgroundFitter
from gdt.core.background.primitives import BackgroundSpectrum
from gdt.core.binning.binned import rebin_by_time
from gdt.core.plot.lightcurve import Lightcurve
from gdt.core.plot.model import ModelFit
from gdt.core.plot.spectrum import Spectrum
from gdt.core.spectra.fitting import SpectralFitterPstat

from gdt.missions.fermi.gbm.collection import GbmDetectorCollection
from gdt.missions.fermi.gbm.finders import TriggerFinder
from gdt.missions.fermi.gbm.response import GbmRsp, GbmRsp2
from gdt.missions.fermi.gbm.tte import GbmTte


def load_data(grb, path=None, dets=None, download=False):

    # connect to data repos and download data
    if download is not False:
        trigfinder = connect_to_protocol(grb)
        trigfinder.get_tte(path, dets=dets, verbose=True)
        trigfinder.get_rsp2(path, cspec=True, ctime=False, 
            dets=dets, verbose=True)
    
    # load the data
    ttes = load_tte(path)
    rsp2s = load_rsp2(path)
    return ttes, rsp2s

def connect_to_protocol(grb_name):
    # connect to the FTP, HTTPS, or AWS protocals in order to download GBM data
    for proto in ['FTP', 'HTTPS', 'AWS']:
        try:
            trigger_ftp = TriggerFinder(grb_name, protocol=proto, timeout=5.)
            break
        except Exception as e:
            print(f"[warn] {proto} failed: {e}")
    if 'trigger_ftp' not in locals():
        raise RuntimeError(f"All protocols failed for {grb_name}")
    return trigger_ftp

def load_tte(path):
    # load and open the data
    tte_files = os.listdir(path)
    tte_files = [os.path.join(path, d) 
                for d in tte_files if 'glg_tte' in d]

    # open the tte data as GDT objects
    ttes = [GbmTte.open(tte) for tte in tte_files] 
    # get the detector names
    dets = [tte.detector for tte in ttes]          
    # save the tte files in a collection
    return GbmDetectorCollection.from_list(ttes, dets=dets, names=dets)

def load_rsp2(path):
    data_files = sorted(os.listdir(path))
    rsp_files = [os.path.join(path, d) for d in data_files \
                if 'cspec' in d and '.rsp2' in d]
    
    rsps = [GbmRsp2.open(f) for f in rsp_files]
    dets = [tte.split('glg_cspec_')[1][:2] for tte in rsp_files]

    return GbmDetectorCollection.from_list(rsps, dets=dets, names=dets)

def fit_background(phaiis, bkgd_ints, polynomials):

    # Rebin background at 1 second
    rebinned_phaiis = phaiis.rebin_time(rebin_by_time, 1.024)
    try:
        dets = phaiis.items
    except:
        dets = [phaiis.detector]

    if isinstance(rebinned_phaiis, list) is False:
        rebinned_phaiis = [rebinned_phaiis]
    if isinstance(phaiis, list) is False:
        phaiis = [phaiis]
    if isinstance(polynomials, list) is False:
        polynomials = [polynomials]
    
    # inialize background fitters with selections and perform the fit
    fitters = [BackgroundFitter.from_phaii(
                phaii, Polynomial, time_ranges=bkgd_int)
                for phaii, bkgd_int in zip(rebinned_phaiis, bkgd_ints)]
    _ = [fitter.fit(order=polynomial) \
                for fitter, polynomial in zip(fitters, polynomials)]

    # extrapolate fit through PHAII files at original timescales
    return [fitter.interpolate_bins(phaii.data.tstart, phaii.data.tstop)
                for fitter, phaii in zip(fitters, phaiis)]

def make_pha(phaii, src, nai_energy=(8., 1000.), bgo_energy=(300., 40000.)):

    # setup
    phas = []
    phaii_src = []
    phaii = phaii.to_list()
    dets = [p.detector for p in phaii]

    # go through each detector
    for d in range(len(dets)):

        # get exact energy limits
        if dets[d][0] == 'n' or dets[d][0] == 'N':
            erange = nai_energy
        else:
            erange = bgo_energy
        ebounds = get_limits(phaii[d], erange)

        # extract the spectra
        pha = phaii[d].to_pha(time_ranges=src, energy_range=ebounds,
            filename=phaii[d].detector)

        # check time binning (there's a floating point error in GDT)
        new_src = [src[0], src[1]]
        if pha.time_range[0] != src[0]:
            new_src[0] += 0.0001
        if pha.time_range[1] != src[1]:
            new_src[1] -= 0.0001
        phas.append(phaii[d].to_pha(time_ranges=new_src, 
            energy_range=ebounds, filename=phaii[d].detector))

        # save the sliced phaii for plotting
        phaii_src.append(phaii[d].slice_time(new_src))        

    phaii_src = GbmDetectorCollection.from_list(phaii_src, dets=dets, names=dets)
    return phas, phaii_src

def get_limits(phaii, erange):
    """Finds the exact detector energy limits from the TTE / CSPEC
    file based on the user-selected range

    Args:
    ----------
    phaii: gdt.missions.fermi.gbm.phaii.Cspec or 
        gdt.missions.fermi.gbm.tte.GbmTte object
        detector data file object
    energy_range: 

    Returns:
    ----------
    low_ebounds[0]: float
        low energy bin edge of the lowest energy bin 
        selected by the user
    high_ebounds[-1]: float
        high energy bin edge of the highest energy bin
        selected by the user 
    """    
    low_ebounds = np.array(phaii.ebounds.low_edges())
    low_ebounds = low_ebounds[low_ebounds >= erange[0]]
    high_ebounds = np.array(phaii.ebounds.high_edges())
    high_ebounds = high_ebounds[high_ebounds <= erange[1]]
    return low_ebounds[0], high_ebounds[-1]

def spectral_analysis(src_phas, bkgds_phas, rsps, models):

    fit_results = []
    for model in models:

        specfitter = SpectralFitterPstat(
            src_phas, bkgds_phas, rsps, method='Nelder-Mead')
        specfitter.fit(model, options={'maxiter': 10000})
        print('{}: {}'.format(model.name, specfitter.message))

        # check if any element in the list is NaN
        if np.isnan(specfitter.parameters).any() == True:
            print('{} has nan parameter values!!!'.format(model.name))
        fit_results.append(specfitter)

        print ('Parameters', specfitter.parameters)
        print ('Uncertainties', specfitter.asymmetric_errors())
        print ('Chi^2/d.o.f.', specfitter.statistic/specfitter.dof, '\n')
        # print ('DOF', specfitter.dof, '\n')

    return fit_results

def plot_lightcurves(phaiis, source=None, bkgd=None, bkgd_ints=None,
    dets=None, nai_energy=(10.,1000.), bgo_energy=(300.,40000.), 
    xlim=None, ylim=None, show=True):

    if dets is None:
        dets = [p.detector for p in phaiis]

    # set the energy ranges
    eranges = [nai_energy if p.detector[0]=='n' 
                else bgo_energy for p in phaiis]

    # make the lightcurves
    lcs = [p.to_lightcurve(energy_range=e) 
                    for p, e in zip(phaiis, eranges)]
    if bkgd is not None:
        bkgd_lcs = [bkgd.integrate_energy(
                    emin=erange[0], emax=erange[1]) 
                    for bkgd, erange in zip(bkgd, eranges)]
    if source is not None:
        src_lcs = [phaii.to_lightcurve(
                    time_range=source, energy_range=erange) 
                    for phaii, erange in zip(phaiis, eranges)]

    # plot data and background fit
    if bkgd is not None:
        lc_plots = [Lightcurve(data=lc, background=blc, label=det) \
                    for lc, blc, det in zip(lcs, bkgd_lcs, dets)]
    else:
        lc_plots = [Lightcurve(data=lc, label=det) 
                    for lc, det in zip(lcs, dets)]
    if source is not None:
        _ = [lc_plot.add_selection(src_lc) 
            for lc_plot, src_lc in zip(lc_plots, src_lcs)]
        _ = [lc_plot.selections[0]._artists[-1].set_label(
            'Source selection') for lc_plot in lc_plots]

    # draw the background selections
    if bkgd_ints is not None:
        for lc_plot, bkgd_int in zip(lc_plots, bkgd_ints):
            lc_plot.ax.axvline(bkgd_int[0][0],
                color=lc_plot.lightcurve.color, linewidth=1., linestyle='--')
            lc_plot.ax.axvline(bkgd_int[0][1], 
                color=lc_plot.lightcurve.color, linewidth=1., linestyle='--')
            lc_plot.ax.axvline(bkgd_int[1][0], 
                color=lc_plot.lightcurve.color, linewidth=1., linestyle='--')
            lc_plot.ax.axvline(bkgd_int[1][1], 
                color=lc_plot.lightcurve.color, linewidth=1., linestyle='--',
                label='Background selection:\n ({}, {}), ({},{})'.format(
                    bkgd_int[0][0], bkgd_int[0][1], 
                    bkgd_int[1][0], bkgd_int[1][1]))
            lc_plot._bkgd._artists[0].set_label('Background fit')

    for plot, det, erange in zip(lc_plots, dets, eranges):
        if xlim is not None:
            plot.xlim = xlim 
        else:
            plot.xlim = (-10, 50)
        if ylim is not None and det[0] != 'b':
            plot.ylim = ylim

        plot.ax.set_title(
            '{} Lightcurve ({:.0f}-{:.0f} keV)'.format(det, *erange))
        plot.ax.legend(loc='upper right')

    if show is not False:
        plt.show()
    else:
        plt.savefig('detector_lightcurve_{}.png'.format(dets[p]), dpi=300)
    plt.close()
    return

def plot_count_spectra(phaiis, back_specs=None, energy_phas=None, show=True):

    phaiis = phaiis.to_list()
    dets = [p.detector for p in phaiis]

    for p in range(len(phaiis)):

        # plot the phaii spectrum over the source selection
        specplot = Spectrum(data=phaiis[p].to_spectrum(), label=dets[p])

        if back_specs is not None:
            _ = specplot.set_background(back_specs[p])
            specplot._bkgd._artists[0].set_label('Background')

        # add the energy selection
        if energy_phas is not None:
            _ = specplot.add_selection(energy_phas[p].data)

            specplot.ax.set_title('{} Spectrum ({:.3f} s - {:.3f} s)'.format(
                specplot.ax.get_label().upper(), 
                energy_phas[p].time_range[0], 
                energy_phas[p].time_range[1]))

            specplot._spec._artists[0].set_label('Source selection')
            specplot.selections[0]._artists[-1].set_label(
                'Source selection ({0:.0f}-{1:.0f} keV)'.format(*energy_phas[p].energy_range))
            specplot.ax.legend(loc="upper right")
        
        if show is not False:
            plt.show()
        else:
            plt.savefig('detector_spectrum_{}.png'.format(dets[p]), dpi=300)

    plt.close()
    return

def plot_spectral_fits(fit_results):

    for fit in fit_results:
        modelplot = ModelFit(fitter=fit)
        modelplot.count_spectrum()
        
        # make it look pretty by setting ymin and ymax
        _, _, data_cnts, data_cnts_err, _ = modelplot._fitter.data_count_spectrum()
        ydata = [np.max(d) for d in data_cnts]
        ymax = np.max(ydata) + 2 * np.max(ydata)
        yerr = [np.min(d) for d in data_cnts_err]
        ymin = np.min(yerr) - np.min(yerr) / 2
        modelplot.ax.set_ylim(bottom=ymin, top=ymax)

        modelplot.ax.grid(which='both')
        modelplot.ax.set_title('{} Model'.format(fit.function_name))
        plt.tight_layout()
        plt.show()
    plt.close()

    for fit in fit_results:
        try:
            modelplot = ModelFit(fitter=fit, interactive=False)
            modelplot.nufnu_spectrum(num_samples=10000)
            modelplot.ylim = (0.1, 10000.0)
            modelplot.ax.grid(which='both')
            modelplot.ax.set_title('{} Model'.format(fit.function_name))
            plt.tight_layout()
            plt.show()
        except:
            print ('Model {} fit plotting failed.'.format(fit.function_name))
    plt.close()
    return