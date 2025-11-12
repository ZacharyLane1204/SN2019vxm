import tessreduce as tr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from tqdm import tqdm

from vxm_control_lightcurves import * 

df = pd.read_csv('/Users/zgl12/Python_Scripts/SN2019VXM/TESS_Outlines/s18_4_1_control_curves.csv')

ras = df['RA'].values
decs = df['DEC'].values

# print(ras)

tess_reducing(ras, decs, sector = 18)