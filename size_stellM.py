import numpy as np
import math as m
import pandas as pd
import statistics as stats
import matplotlib.pyplot as plt
#import seaborn as sns
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
from scipy import stats
from scipy.optimize import curve_fit
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
plt.rcParams["figure.figsize"] = [10, 8]
plt.rcParams.update({'font.size': 18})

#%%
##-------------------
## read in sample-- GSWLC, Salim 2016 x SGA,
##-------------------
file0 = Table.read('/Users/hannahchristie/Documents/GitHub/size_stellMass_rel/GSWLCxSGA.fits')
file0.keep_columns(['Z', 'LOGMSTAR', 'LOGMSTARERR', 'LOGSFRSED', 'LOGSFRSEDERR', 'D25_LEDA', 'MAG_LEDA', 'SB_D25_LEDA', 'gmag', 'rmag', 'imag'])
df = file0.to_pandas()
df = df[df.Z < 0.05]
df = df[df.LOGMSTAR > -99]
df = df[df.LOGSFRSED > - 99]
#df = df[(df.gmag - df.rmag) > 0.64]
print(df.shape)


#%%
## information on 2885
#sfr_2885 = 2.47
#d25_2885 - ????

#%%

plt.scatter(df.D25_LEDA, df.LOGMSTAR, c = df.LOGSFRSED, cmap='viridis', alpha =0.5)
plt.xlabel('D_25 (arcmin)')
plt.ylabel('LOGMSTAR ($M_\odot$)')
cbar = plt.colorbar()
cbar.set_label('LOGSFR', fontsize = 15)
plt.legend()
plt.show()

#%%
plt.scatter(df.LOGMSTAR, df.D25_LEDA, c = df.gmag - df.rmag, cmap='RdBu', vmin = 0, vmax = 1.2, alpha =0.5)
plt.ylabel('D_25 (arcmin)')
plt.xlabel('LOGMSTAR ($M_\odot$)')
cbar = plt.colorbar()
cbar.set_label('g-r colour', fontsize = 15)
#plt.legend()
plt.show()