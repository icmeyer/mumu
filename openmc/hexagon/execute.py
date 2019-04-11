from hexfunction import hexfunction
import openmc
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pandas as pd
import os, shutil

if not os.path.isdir('work'):
    os.mkdir('work')
os.chdir('work')
if not os.path.isdir('figures'):
    os.mkdir('figures')


fractions = np.linspace(0.1,0.5,10)
pitch = 2 # cm
all_fractions = []
ks = []
k_sds = []
for pack_frac in fractions:
    all_tallies = hexfunction(pitch, pack_frac) #cm
    all_fractions.append(all_tallies)
    ks.append(all_tallies['kinf mean'])
    k_sds.append(all_tallies['kinf sd'])

plt.errorbar(fractions, ks, yerr=k_sds, fmt='.',markersize=12,capsize=3)

plt.ylabel('k inf')
plt.xlabel('packing fraction')
plt.show()

#Write results to csv
os.chdir('..')
all_fractions.to_csv('results_packing_fraction.csv')
