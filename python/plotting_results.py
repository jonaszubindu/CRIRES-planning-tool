#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 11:25:22 2020

@author: jonaszbinden
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import cm

cc = cm.viridis(np.linspace(0,1,10))
# CSV files with eclipse data
filename100 = 'csv_files/Eclipse_events_processed_2020-11-17_365d_SN100edit.csv'
filename150 = 'csv_files/Eclipse_events_processed_2020-11-29_365d_SN150edit.csv'
filename200 = 'csv_files/Eclipse_events_processed_2020-11-17_365d_SN200edit.csv'

def snr_total(n_transits, snr_single):
    return np.sqrt(n_transits * 20) * snr_single

# Loading eclipse data into dataframe
df100 = pd.read_csv(filename100)
df150 = pd.read_csv(filename150)
df200 = pd.read_csv(filename200)

# Loading planet list
planets = pd.read_csv('csv_files/PlanetList.csv')

pl_bmassj_list100 = []
pl_bmassj_list150 = []
pl_bmassj_list200 = []

pl_radj_list100 = []
pl_radj_list150 = []
pl_radj_list200 = []

for name in df100['Name']:
    planet = planets[(planets['pl_name'] == name)]
    pl_radj_list100.append(planet.pl_radj.values[0])
    pl_bmassj_list100.append(planet.pl_bmassj.values[0])
df100['pl_bmassj'] = pl_bmassj_list100
df100['pl_radj'] = pl_radj_list100

for name in df150['Name']:
    # print(name)
    planet = planets[(planets['pl_name'] == name)]
    pl_radj_list150.append(planet.pl_radj.values[0])
    pl_bmassj_list150.append(planet.pl_bmassj.values[0])
df150['pl_bmassj'] = pl_bmassj_list150
df150['pl_radj'] = pl_radj_list150


for name in df200['Name']:
    # print(name)
    planet = planets[(planets['pl_name'] == name)]
    pl_radj_list200.append(planet.pl_radj.values[0])
    pl_bmassj_list200.append(planet.pl_bmassj.values[0])
df200['pl_bmassj'] = pl_bmassj_list200
df200['pl_radj'] = pl_radj_list200

Names100 = []
for name in df100['Name']:
    Names100.append(name)

Names100 = list(dict.fromkeys(Names100))
df_plot100 = pd.DataFrame(Names100)
df_plot100.columns = ['Name']
df_plot100['pl_bmassj'] = np.zeros(len(df_plot100))
df_plot100['pl_radE'] = np.zeros(len(df_plot100))
df_plot100['Number of transits'] = np.zeros(len(df_plot100))
df_plot100['Max Number of Exp'] = np.zeros(len(df_plot100))
df_plot100['Effective Temperature'] = np.zeros(len(df_plot100))

for name in Names100:
    pl_bmassj = df100[df100['Name'] == name]['pl_bmassj'].values[0]
    pl_radE = df100[df100['Name'] == name]['pl_radj'].values[0]*11.2
    star_Teff = df100[df100['Name'] == name]['Effective Temperature'].values[0]
    num_of_trans = df100[df100['Name'] == name]['Number of transits'].values[0]
    max_exp_poss = df100[df100['Name'] == name]['Number of exposures possible'].max()
    i = df_plot100.index[df_plot100['Name'] == name]
    df_plot100.loc[:, ('pl_bmassj')][i] = pl_bmassj
    df_plot100.loc[:, ('pl_radE')][i] = pl_radE
    df_plot100.loc[:, ('Number of transits')][i] = num_of_trans
    df_plot100.loc[:, ('Max Number of Exp')][i] = max_exp_poss
    df_plot100.loc[:, ('Effective Temperature')][i] = star_Teff
    
Names150 = []
for name in df150['Name']:
    Names150.append(name)
Names150 = list(dict.fromkeys(Names150))
df_plot150 = pd.DataFrame(Names150)
df_plot150.columns = ['Name']
df_plot150['pl_bmassj'] = np.zeros(len(df_plot150))
df_plot150['pl_radE'] = np.zeros(len(df_plot150))
df_plot150['Number of transits'] = np.zeros(len(df_plot150))
df_plot150['Max Number of Exp'] = np.zeros(len(df_plot150))
df_plot150['Effective Temperature'] = np.zeros(len(df_plot150))


for name in Names150:
    pl_bmassj = df150[df150['Name'] == name]['pl_bmassj'].values[0]
    pl_radE = df150[df150['Name'] == name]['pl_radj'].values[0]*11.2
    star_Teff = df150[df150['Name'] == name]['Effective Temperature'].values[0]
    num_of_trans = df150[df150['Name'] == name]['Number of transits'].values[0]
    max_exp_poss = df150[df150['Name'] == name]['Number of exposures possible'].max()
    i = df_plot150.index[df_plot150['Name'] == name]
    df_plot150.loc[:, ('pl_bmassj')][i] = pl_bmassj
    df_plot150.loc[:, ('pl_radE')][i] = pl_radE
    df_plot150.loc[:, ('Number of transits')][i] = num_of_trans
    df_plot150.loc[:, ('Max Number of Exp')][i] = max_exp_poss
    df_plot150.loc[:, ('Effective Temperature')][i] = star_Teff
    
Names200 = []
for name in df200['Name']:
    Names200.append(name)
Names200 = list(dict.fromkeys(Names200))
df_plot200 = pd.DataFrame(Names200)
df_plot200.columns = ['Name']
df_plot200['pl_bmassj'] = np.zeros(len(df_plot200))
df_plot200['pl_radE'] = np.zeros(len(df_plot200))
df_plot200['Number of transits'] = np.zeros(len(df_plot200))
df_plot200['Max Number of Exp'] = np.zeros(len(df_plot200))
df_plot200['Effective Temperature'] = np.zeros(len(df_plot200))

for name in Names200:
    pl_bmassj = df200[df200['Name'] == name]['pl_bmassj'].values[0]
    pl_radE = df200[df200['Name'] == name]['pl_radj'].values[0]*11.2
    star_Teff = df200[df200['Name'] == name]['Effective Temperature'].values[0]
    num_of_trans = df200[df200['Name'] == name]['Number of transits'].values[0]
    max_exp_poss = df200[df200['Name'] == name]['Number of exposures possible'].max()
    i = df_plot200.index[df_plot200['Name'] == name]
    df_plot200.loc[:, ('pl_bmassj')][i] = pl_bmassj
    df_plot200.loc[:, ('pl_radE')][i] = pl_radE
    df_plot200.loc[:, ('Number of transits')][i] = num_of_trans
    df_plot200.loc[:, ('Max Number of Exp')][i] = max_exp_poss
    df_plot200.loc[:, ('Effective Temperature')][i] = star_Teff
    
    
df_plot100 = df_plot100[df_plot100['Number of transits'] != 0]
df_plot150 = df_plot150[df_plot150['Number of transits'] != 0]
df_plot200 = df_plot200[df_plot200['Number of transits'] != 0]

snr100 = snr_total(1, 100)
snr150 = snr_total(1, 150)
snr200 = snr_total(1, 200)
    
plt.clf()

fig1 = plt.figure(figsize=(6,5))
ax1 = fig1.add_subplot(111)
binwidth = 1
rw = 0.9
hist_bins = np.arange(int(min(df_plot100['Number of transits'])), int(max(df_plot100['Number of transits'])) + binwidth, binwidth)
ax1.hist([df_plot100['Number of transits'].values, df_plot150['Number of transits'].values, df_plot200['Number of transits'].values], bins=hist_bins, rwidth=rw, color=[cc[1],cc[7],cc[9]], stacked=False, label = [f'minimum SN = {snr100:.1f}', f'minimum SN = {snr150:.1f}', f'minimum SN = {snr200:.1f}'])
# ax1.hist(df_plot150['Number of transits'].values, bins=np.arange(min(df_plot150['Number of transits']), max(df_plot150['Number of transits']) + binwidth, binwidth), label = f'minimum SN = {snr150:.1f}')
# ax1.hist(df_plot200['Number of transits'].values, bins=np.arange(min(df_plot200['Number of transits']), max(df_plot200['Number of transits']) + binwidth, binwidth), label = f'minimum SN = {snr200:.1f}')

ax1.set_xlabel('Number of observable transits')
ax1.set_ylabel('Number of Planets')
fig1.legend(loc='upper right')

fig2 = plt.figure(figsize=(6,5))
ax2 = fig2.add_subplot(111)

binwidth = 1
rw = 0.6 # this resizes the width of the bins on the plot to rw * binwidth (it's just aesthetic), the bins values are still computed for binwidth
hist_bins = np.arange(int(min(df_plot100['pl_radE'])), int(max(df_plot100['pl_radE'])) + binwidth, binwidth)
ax2.hist([df_plot100['pl_radE'].values, df_plot150['pl_radE'].values, df_plot200['pl_radE'].values], bins=hist_bins, rwidth=rw, color=[cc[1],cc[7],cc[9]], stacked=False, label = [f'minimum SN = {snr100:.1f}', f'minimum SN = {snr150:.1f}', f'minimum SN = {snr200:.1f}'])
ax2.set_xlabel('Planet Radius [Earth radii]')
ax2.set_ylabel('Number of Planets')
ax2.set_xticks(hist_bins)
fig2.legend(loc='upper right')

# fig3 = plt.figure(figsize=(6,5))
# ax3 = fig3.add_subplot(111)
# binwidth = 10
# rw = 1.2
# hist_bins = np.arange(int(min(df_plot100['Max Number of Exp'])), int(max(df_plot100['Max Number of Exp'])) + binwidth, binwidth)
# ax3.hist([df_plot100['Max Number of Exp'].values, df_plot150['Max Number of Exp'].values, df_plot200['Max Number of Exp'].values], bins=hist_bins, rwidth=rw, color=[cc[1],cc[7],cc[9]], stacked=False, label = [f'minimum SN = {snr100:.1f}', f'minimum SN = {snr150:.1f}', f'minimum SN = {snr200:.1f}'])

# ax3.set_xlabel('Maximum Number of Exposures possible')
# ax3.set_ylabel('Number of Planets')
# fig3.legend(loc='upper right')

fig4 = plt.figure(figsize=(6,5))
ax4 = fig4.add_subplot(111)
binwidth = 100
rw = 0.9
hist_bins = np.arange(int(min(df_plot100['Effective Temperature'])), int(max(df_plot100['Effective Temperature'])) + binwidth, binwidth)
ax4.hist([df_plot100['Effective Temperature'].values, df_plot150['Effective Temperature'].values, df_plot200['Effective Temperature'].values], bins=hist_bins, rwidth=rw, color=[cc[1],cc[7],cc[9]], stacked=False, label = [f'minimum SN = {snr100:.1f}', f'minimum SN = {snr150:.1f}', f'minimum SN = {snr200:.1f}'])

ax4.set_xlabel('Stellar Effective Temperature')
ax4.set_ylabel('Number of Planets')
fig4.legend(loc='upper right')

fig1.savefig('Plots/Num_of_transHist.eps')
fig2.savefig('Plots/RadiusHist.eps')
# fig3.savefig('Plots/Max_Num_of_exp_possHist.eps')
fig4.savefig('Plots/TeffHist.eps')









