# SN 2019vxm Data analysis and 

Description of the scripts and what analysis is contained in the paper.

## Data contents

Not every data script is completely necessary

- 2019vxm_control_coords.txt -- ATClean control coordinates output (Figure 2)
- 2019vxm-combined-20200213_ap1.flm -- Spectra used for TESS zeropoint calculation (Used in Figure 1)
- all_magnitude_thresholds_gaussian.txt -- The precursor search results for the early rise (Figure 6 & Table 2)
- best_magnitude_thresholds_gaussian.txt -- The precursor search results for the early rise (Figure 6 & Table 2)
- collated_vxm_data_with_tess.csv -- All vxm photometry data (Figure 4 & Figure 5)
- ps1_galaxy.csv -- PanSTARRS1 galaxy information
- schulze_2021.csv -- The properties of Type IIn hosts from Schulze et al. (2021) (Figure 8)
- SN2019vxm_2020-03-13_KAST.txt -- Extinction correction spectra
- tess_binned_lc_18.npy -- TESS binned lightcurve produced from TESSReduce (Figure 1 & Figure 4 & Figure 5)
- tess_sector_18_psf.npy -- TESSReduce rise data used for the main analysis (Figure 1)
- x_nufnu.npy -- Fermi forward modelled spectra samples energies for plotting and analysis (Figure 3)
- y_nufnu.npy -- Fermi forward modelled spectra samples rates for plotting and analysis (Figure 3)

## General_Scripts contents

Not every general script is completely necessary

- extinction_correction.ipynb -- Correcting the extinction from MW, host, and more.
- extinction_V.ipynb -- Correcting the extinction from MW, host, and more in V-band.
- localisation_script.ipynb -- Contains the contours and known X-ray sources (Figure 7)
- mangle.py -- Spectral analysis for bandpass calculations
- plot_vxm.ipynb -- Plotting scripts (Figure 4 & Figure 5)
- precursor_plot.ipynb -- Precursor plotting script (Figure 6)
- run_mangle.ipynb -- Spectral analysis for bandpass calculations running script

## Modelling contents

Not every modelling script is completely necessary

- mosfit_all_lcs_single.py -- MOSFit code for fitting
- mosfit_plotter_single.py -- MOSFit code for fitting
- tensions_multimodal.ipynb -- Multimodal distribution sigma tension calculator

## Probability contents

- probability_brute_force.ipynb -- Brute force probability calculator
- sigma_level.ipynb -- Confidence levels and calculations
- temporal_one_two.ipynb -- Temporal methods one and two for probabilities (Table 3)

## SED_Analysis/Anya contents

- SBIPP_Blast_Master_Plotting.ipynb -- FrankenBLAST plotter

## Telescopes_Scripts contents

### Telescopes_Scripts/ATLAS contents

- plot_atlas_lc_position.ipynb -- Plotting (Figure 2)
- ps1_table_obj.py -- PanSTARRS1 information gathering

### Telescopes_Scripts/Fermi contents

- fermi_script.py -- Fermi light curve script for display (Figure 3)
- fermi_to_plot_analyse.ipynb -- Fermi plotting script (Figure 3)
- gbm_spectral_analysis.ipynb -- Gathering Fermi specta and saving it for plotting (Figure 3)
- gbm_tools2.py -- Spectral analysis and modelling background script (Figure 3)

### Telescopes_Scripts/TESS contents

- TESS_Outlines/* -- The TESS CCD outlines for Sector 18 (Figure 7)
- Control_Curves/Detrended_Data/* -- detrending the control light curve data (Figure 6 & Table 2)
- Control_Curves/Original_Lightcurves/* -- original control light curves data for presursor significance breakout (Figure 6 & Table 2)
- control_run.py -- Gather contol light curves (Figure 6 & Table 2)
- detrend_run.py -- Detrend the control light curves (Figure 6 & Table 2)
- model_vxm_detrend.ipynb -- Detrend vxm light curve (Figure 1)
- plot_tess_rise.ipynb -- Plot the TESS rise (Figure 1)
- plotting_trends.ipynb -- Plotting control curve trends
- tess_lc_rise.py -- Main script for the TESS rise (Figure 1)
- vxm_control_create.ipynb -- Create control curves running script
- vxm_control_lightcurves.py -- Major script that creates the control curves

### Telescopes_Scripts/ZTF contents

- ztf_data.ipynb -- ZTF data extraction, unnecessary as in the data scripts