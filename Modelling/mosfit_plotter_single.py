#%matplotlib inline
#%config InlineBackend.figure_format='retina'

import corner
import json
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from tqdm import tqdm_notebook
from collections import OrderedDict
from mosfit.plotting import bandcolorf
import glob

snname = 'SN2019vxm'
sns.reset_orig()

labels = [
    r"log($M_{\mathrm{CSM}}/\mathrm{M_\odot})$",  # Ejecta mass in solar masses
    r"$M_{\mathrm{ej}}/\mathrm{M_\odot}$",  # Ejecta velocity in km/s
    r"n",  # Kinetic energy in ergs
    r"log($n_{\mathrm{h}}/\mathrm{cm^{-3}})$",
    r"log($R_{\mathrm{0}}/\mathrm{AU})$",
    r"log($\rho_{\mathrm{0}}/\mathrm{g cm^{-3}})$",
    r"s",
    r"log($T_{\mathrm{min}}/\mathrm{K})$",
    r"t$_{\mathrm{exp}}$ / days",
    r"log($\sigma)$",
    r"log($v_{\mathrm{ej}}/\mathrm{km s^{-1}})$",  # Peak time in days
]

plt.rcParams["font.family"] = "serif"
plt.rcParams.update({'font.size': 12})

#folders = '/Users/conorransome/SN_IIn/mosfit_jsons/new_s/walkers/'
folders = '/Users/conorransome/Research/MOSFiT/SN2019vxm_new_converged_withu.json'
#folders = '/Users/conorransome/SN_IIn/SN2023jgf_project/products/'
#folders = '/Users/conorransome/SN_IIn/mosfit_jsons/param_tests/nologmej/fixing/walkers/SN2020tan_plot_plot_ext.json' 
#folders = '/Users/conorransome/SN_IIn/mosfit_jsons/extended_lcs/'  #param_tests/nologmej/fixing/walkers/ SN2020tan_plot_plot_ext
#folders = '/Users/conorransome/SN_IIn/Rubin_DR1/mosfit_jsons/'
#folders = folders+snname+'__all_data'+'.json'

suffix = ''
#for i in folders:


#with open('../products/walkers.json', 'r', encoding = 'utf-8') as f:
#with open('/storage/work/cbr5597/IIn_LC/IIn_LC/ZTFBTS_data/ZTF_SNe/products/SN2018fdt.json', 'r', encoding = 'utf-8') as f:
#with open('/Users/cbr5597/SN_IIn/MOSFITJSON_Rfree/r1/SN2000eo.json', 'r', encoding = 'utf-8') as f:
with open(folders, 'r', encoding = 'utf-8') as f:
    data = json.loads(f.read())
    
    if 'name' not in data:
        data = data[list(data.keys())[0]]

photo = data['photometry']
model = data['models'][0]
snname = data['name']

real_data = len([x for x in photo if 'band' in x and 'magnitude' in x and (
    'realization' not in x or 'simulated' in x)]) > 0

band_attr = ['band', 'instrument', 'telescope', 'system', 'bandset']
band_list = list(set([tuple(x.get(y, '')
                            for y in band_attr) for x in photo
                            if 'band' in x and 'magnitude' in x]))
real_band_list = list(set([tuple(x.get(y, '')
                                 for y in band_attr) for x in photo
                                 if 'band' in x and 'magnitude' in x and (
                                     'realization' not in x or 'simulated' in x)]))
xray_instrument_attr = ['instrument', 'telescope']
xray_instrument_list = list(set([tuple(x.get(y, '')
                            for y in xray_instrument_attr) for x in photo
                            if 'instrument' in x and 'countrate' in x])) 
real_band_list = list(set([tuple(x.get(y, '')
                                 for y in band_attr) for x in photo
                                 if 'band' in x and 'magnitude' in x and (
                                     'realization' not in x or 'simulated' in x)]))
real_xray_instrument_list = list(set([tuple(x.get(y, '')
                                 for y in xray_instrument_attr) for x in photo
                                 if 'instrument' in x and 'countrate' in x and (
                                     'realization' not in x or 'simulated' in x)]))


fig = plt.figure(figsize=(12,10))
plt.gca().invert_yaxis()
plt.gca().set_xlabel('MJD')
plt.gca().set_ylabel('Apparent Magnitude')
used_bands = []
for full_band in tqdm_notebook(band_list, desc='Photo', leave=False):
    (band, inst, tele, syst, bset) = full_band
    try:
        inst_exclusive_list
    except:
        pass
    else:
        if inst not in inst_exclusive_list:
            continue
    extra_nice = ', '.join(list(filter(None, OrderedDict.fromkeys((inst, syst, bset)).keys())))
    nice_name = band + ((' [' + extra_nice + ']') if extra_nice else '')
    
    realizations = [[] for x in range(len(model['realizations']))]
    for ph in photo:
        rn = ph.get('realization', None)
        si = ph.get('simulated', False)
        if rn and not si:
            if tuple(ph.get(y, '') for y in band_attr) == full_band:
                realizations[int(rn) - 1].append((
                    float(ph['time']), float(ph['magnitude']), [
                        float(ph.get('e_lower_magnitude', ph.get('e_magnitude', 0.0))),
                        float(ph.get('e_upper_magnitude', ph.get('e_magnitude', 0.0)))],
                ph.get('upperlimit')))
    numrz = np.sum([1 for x in realizations if len(x)])
    for rz in realizations:
        if not len(rz):
            continue
        xs, ys, vs, us = zip(*rz)
        label = '' if full_band in used_bands or full_band in real_band_list else nice_name
        if max(vs) == 0.0:
            plt.plot(xs, ys, color='green',
                             label=label, linewidth=0.5) #bandcolorf(band)
        else:
            xs = np.array(xs)
            ymi = np.array(ys) - np.array([np.inf if u else v[0] for v, u in zip(vs, us)])
            yma = np.array(ys) + np.array([v[1] for v in vs])
            #plt.fill_between(xs, ymi, yma, color=bandcolorf(band), edgecolor=None,
            #                 label=label, alpha=1.0/numrz, linewidth=0.0)
            #plt.plot(xs, ys, 
            #                 label=label, alpha=0.7, linewidth=0.5)#color='red'bandcolorf(band
            if band == 'g':
                plt.plot(xs, ys, color = 'green', 
                             label=label, alpha=0.5, linewidth=0.5)
            if band == 'r':
                plt.plot(xs, ys, color = 'red',
                             label=label, alpha=0.5, linewidth=0.5)
        if label:
            used_bands = list(set(used_bands + [full_band]))
    if real_data:
        for s in range(2):
            if s == 0:
                cond = False
                symb = 'o'
            else:
                cond = True
                symb = 'v'
            vec = [(float(x['time']), float(x['magnitude']),
                    0.0 if 'upperlimit' in x else float(x.get('e_lower_magnitude', x.get('e_magnitude', 0.0))),
                    float(x.get('e_upper_magnitude', x.get('e_magnitude', 0.0)))) for x in photo
                   if 'magnitude' in x and ('realization' not in x or 'simulated' in x) and
                   'host' not in x and 'includeshost' not in x and
                   x.get('upperlimit', False) == cond and
                   tuple(x.get(y, '') for y in band_attr) == full_band]
            if not len(vec):
                continue
            xs, ys, yls, yus = zip(*vec)
            label = nice_name if full_band not in used_bands else ''
            if band == 'g':
                plt.errorbar(xs, ys, yerr=(yus, yls), color='green', fmt=symb,
                         label=label,
                         markeredgecolor='black', markeredgewidth=1, capsize=1,
                         elinewidth=1.5, capthick=2, zorder=10)
                plt.errorbar(xs, ys, yerr=(yus, yls), color='k', fmt=symb, capsize=2,
                         elinewidth=2.5, capthick=3, zorder=5)
            if band == 'r':
                plt.errorbar(xs, ys, yerr=(yus, yls), color='red', fmt=symb,
                         label=label,
                         markeredgecolor='black', markeredgewidth=1, capsize=1,
                         elinewidth=1.5, capthick=2, zorder=10)
                plt.errorbar(xs, ys, yerr=(yus, yls), color='k', fmt=symb, capsize=2,
                         elinewidth=2.5, capthick=3, zorder=5)
            if label:
                used_bands = list(set(used_bands + [full_band]))
plt.margins(0.02, 0.1)
plt.legend()
#plt.ylim(bottom=23, top = 18)
#plt.xlim(59050, 59275)
plt.show()
fig.savefig('/Users/conorransome/Research/' + snname + suffix + '_new_lc.pdf')
plt.close()

import logging
logging.disable(logging.WARNING)

# Construct walker arrays for corner
corner_input = []
pars = [x for x in model['setup'] if model['setup'][x].get('kind') == 'parameter' and
        'min_value' in model['setup'][x] and 'max_value' in model['setup'][x]]
weights = []
for realization in model['realizations']:
    par_vals = realization['parameters']
    if 'weight' in realization:
        weights.append(float(realization['weight']))
    var_names = ['$' + ('\\log\\, ' if par_vals[x].get('log') else '') +
                 par_vals[x]['latex'] + '$' for x in par_vals if x in pars and 'fraction' in par_vals[x]]
    corner_input.append([np.log10(par_vals[x]['value']) if
                         par_vals[x].get('log') else par_vals[x]['value'] for x in par_vals
                         if x in pars and 'fraction' in par_vals[x]])

    #print(corner_input)

weights = weights if len(weights) else None
ranges = [0.999 for x in range(len(corner_input[0]))]
cfig = corner.corner(np.array(corner_input), labels=labels, quantiles=[0.16, 0.5, 0.84],
                    show_titles=True, weights=weights, range=ranges,  title_kwargs ={"fontsize": 15,
    "bbox": dict(facecolor='white', edgecolor='black', boxstyle='round')}, label_kwargs={"fontsize": 18})
plt.tick_params(axis='both', which='major', labelsize=16)
#cfig.subplots_adjust(hspace=0, wspace=0)
#cfig.tight_layout()
cfig.savefig('/Users/conorransome/Research/' + snname + suffix + '_new_corner.pdf', bbox_inches="tight")#, bbox_inches="tight")


#print(snname)