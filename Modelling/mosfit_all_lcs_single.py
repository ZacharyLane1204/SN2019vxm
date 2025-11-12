import json
import numpy as np
import matplotlib.pyplot as plt
import glob
from tqdm import tqdm
import os

folder = '/Users/conorransome/SN_IIn/mosfit_jsons/extended_lcs/'
folder = '/Users/conorransome/SN_IIn/mosfit_jsons/newrho0/'
folder = '/Users/conorransome/SN_IIn/SN2023jgf_project/mosfit/'
#folder = '/Users/conorransome/SN_IIn/Rubin_DR1/mosfit_jsons/default/'
plot_folder = '/Users/conorransome/SN_IIn/mosfit_jsons/param_tests/nologmej/plots/'
plot_folder = '/Users/conorransome/SN_IIn/SN2023jgf_project/plots/'
#plot_folder = '/Users/conorransome/SN_IIn/Rubin_DR1/plots/default/'

colors = {
    "r": '#D73027', "ZTF-r": '#D73027', "R": '#FC8D59', "rp": '#B2182B', "r'": '#B2182B',
    "g": '#1A9850', "ZTF-g": '#1A9850', "G": '#66BD63', "gp": '#006837', "g'": '#006837',
    "i": '#8E0152', "I": '#4D004B', "ip": '#C51B7D', "i'": '#C51B7D', "z": '#3288BD',
    "Z": '#5E4FA2', "zp": '#2C7BB6', "y": '#A6611A', "Y": '#DFD3B1', "yp": '#8073AC',
    "V": '#8C510A', "B": '#4575B4', "U": '#E31A1C', "u": '#D01C8B', "up": '#F1B6DA',
    "J": '#74ADD1', "H": '#542788', "K": '#000000'
}

# Function to generate and save the plot for each supernova
def plot_supernova(snname, photo, model, plot_folder):
    fig, ax = plt.subplots(figsize=(5, 5))

    constant = 0
    ylimit = []
    band_attr = ['band', 'instrument', 'telescope', 'system', 'bandset']
    real_data = [x for x in photo if 'band' in x and 'magnitude' in x and ('realization' not in x or 'simulated' in x)]
    band_list = list(set([tuple(x.get(y, '') for y in band_attr) for x in photo if 'band' in x and 'magnitude' in x]))
    
    reference_peak_time = None  # Reset for each supernova
    bands_in_legend = set()  # Track bands added to the legend

    for full_band in band_list:
        (band, inst, tele, syst, bset) = full_band

        try:
            realizations = [[] for _ in range(len(model['realizations']))]

            for ph in photo:
                rn = ph.get('realization', None)
                if rn:
                    if tuple(ph.get(y, '') for y in band_attr) == full_band:
                        realizations[int(rn) - 1].append((
                            float(ph['time']), float(ph['magnitude']), [
                                float(ph.get('e_lower_magnitude', ph.get('e_magnitude', 0.0))),
                                float(ph.get('e_upper_magnitude', ph.get('e_magnitude', 0.0)))
                            ],
                            ph.get('upperlimit')
                        ))

            if any(realizations):
                for rz in realizations:
                    if band in colors:
                        if rz:
                            xs_r, ys_r, _, _ = zip(*rz)
                            peak_r_ind = np.argmin(ys_r)
                            peak_time_r = xs_r[peak_r_ind]

                            if reference_peak_time is None:
                                reference_peak_time = peak_time_r

                            xs_r = np.array(xs_r) - reference_peak_time
                            # Plot without the label (for realizations)
                            ax.plot(xs_r, np.array(ys_r) - constant, c=colors[band], alpha=0.05)

                for ph in real_data:
                    if tuple(ph.get(y, '') for y in band_attr) == full_band:
                        x_s = float(ph['time']) - reference_peak_time
                        y_s = float(ph['magnitude'])
                        lower_err = float(ph.get('e_lower_magnitude', ph.get('e_magnitude', 0.0)))
                        upper_err = float(ph.get('e_upper_magnitude', ph.get('e_magnitude', 0.0)))
                        ylimit.append(np.min(y_s))

                        if ph.get('upperlimit', False):
                            # Add to legend if band not already added
                            if band not in bands_in_legend:
                                ax.scatter(x_s, y_s - constant, s=50, marker='v', facecolor='none', edgecolor=colors[band], alpha=1, zorder=1000, label=band)
                                bands_in_legend.add(band)
                            else:
                                ax.scatter(x_s, y_s - constant, s=50, marker='v', facecolor='none', edgecolor=colors[band], alpha=1, zorder=1000)
                        else:
                            # Add to legend if band not already added
                            if band not in bands_in_legend:
                                ax.scatter(x_s, y_s - constant, s=150, marker='o', c=colors[band], edgecolor='black', alpha=1, zorder=1000, label=band)
                                ax.errorbar(x_s, y_s - constant, yerr=[[lower_err], [upper_err]], c=colors[band], alpha=1, ls='none')
                                bands_in_legend.add(band)
                            else:
                                ax.scatter(x_s, y_s - constant, s=150, marker='o', c=colors[band], edgecolor='black', alpha=1, zorder=1000)
                                ax.errorbar(x_s, y_s - constant, yerr=[[lower_err], [upper_err]], c=colors[band], alpha=1, ls='none')

                constant -= 1.5

        except Exception as e:
            print(f"Error processing {snname} with band {band}: {e}")
            continue

    ax.set_xlim(left=-20, right=125)
    ax.invert_yaxis()
    #ax.set_ylim(27, 16.5)

    # if ylimit:
    #     ax.set_ylim(top=np.min(ylimit) - 3, bottom=np.min(ylimit) - constant + 3)

   # ax.text(0.01, 0.95, snname, transform=ax.transAxes, fontsize=12, va='top', ha='left',
   #         bbox=dict(facecolor='white', alpha=0.75, edgecolor='none'))

    #fig.text(0.04, 0.55, 'Magnitude + offset', va='center', rotation=90, fontsize=14)
    #fig.text(0.45, 0.03, 'Time / days from peak', va='center', fontsize=14)
    plt.xlabel('Phase from peak / days', fontsize = 14)
    plt.ylabel('Apparent Magnitude + offset', fontsize = 14)

    # Add legend outside of the plot
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10, title="Band")

    # Save the figure
    plt.savefig(os.path.join(plot_folder, f'{snname}_light_curve.pdf'), format='pdf', dpi=150, bbox_inches='tight')
    plt.close()

# Get a list of all JSON files in the directory
light_curve_files = glob.glob(os.path.join(folder, '*.json'))

# Loop through each file and generate the plot for each supernova
for filename in tqdm(light_curve_files, desc='Plotting Light Curves'):
    with open(filename, 'r', encoding='utf-8') as f:
        data = json.loads(f.read())
        if 'name' not in data:
            data = data[list(data.keys())[0]]

        photo = data['photometry']
        model = data['models'][0]
        snname = data['name']

        plot_supernova(snname, photo, model, plot_folder)
