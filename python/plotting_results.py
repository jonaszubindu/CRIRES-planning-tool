#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 11:25:22 2020

@author: jonaszbinden
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

filename100 = 'csv_files/Eclipse_events_processed_2020-11-17_365d_SN100edit.csv'
filename150 = 'csv_files/Eclipse_events_processed_2020-11-29_365d_SN150edit.csv'
filename200 = 'csv_files/Eclipse_events_processed_2020-11-17_365d_SN200edit.csv'
# filename800 = 'csv_files/Eclipse_events_processed_2020-11-17_365d_SN800edit.csv'


df100 = pd.read_csv(filename100)
df150 = pd.read_csv(filename150)
df200 = pd.read_csv(filename200)
# df800 = pd.read_csv(filename800)

planets = pd.read_csv('csv_files/PlanetList.csv')

pl_bmassj_list100 = []
pl_bmassj_list150 = []
pl_bmassj_list200 = []
# pl_bmassj_list800 = []

pl_radj_list100 = []
pl_radj_list150 = []
pl_radj_list200 = []
# pl_radj_list800 = []

for name in df100['Name']:
    # print(name)
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


# for name in df800['Name']:
#     # print(name)
#     planet = planets[(planets['pl_name'] == name)]
#     pl_radj_list800.append(planet.pl_radj.values[0])
#     pl_bmassj_list800.append(planet.pl_bmassj.values[0])
# df800['pl_bmassj'] = pl_bmassj_list800
# df800['pl_radj'] = pl_radj_list800

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

for name in Names100:
    pl_bmassj = df100[df100['Name'] == name]['pl_bmassj'].values[0]
    pl_radE = df100[df100['Name'] == name]['pl_radj'].values[0]*11.2
    num_of_trans = df100[df100['Name'] == name]['Number of transits'].values[0]
    max_exp_poss = df100[df100['Name'] == name]['Number of exposures possible'].max()
    i = df_plot100.index[df_plot100['Name'] == name]
    df_plot100.loc[:, ('pl_bmassj')][i] = pl_bmassj
    df_plot100.loc[:, ('pl_radE')][i] = pl_radE
    df_plot100.loc[:, ('Number of transits')][i] = num_of_trans
    df_plot100.loc[:, ('Max Number of Exp')][i] = max_exp_poss
    
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

for name in Names150:
    pl_bmassj = df150[df150['Name'] == name]['pl_bmassj'].values[0]
    pl_radE = df150[df150['Name'] == name]['pl_radj'].values[0]*11.2
    num_of_trans = df150[df150['Name'] == name]['Number of transits'].values[0]
    max_exp_poss = df150[df150['Name'] == name]['Number of exposures possible'].max()
    i = df_plot150.index[df_plot150['Name'] == name]
    df_plot150.loc[:, ('pl_bmassj')][i] = pl_bmassj
    df_plot150.loc[:, ('pl_radE')][i] = pl_radE
    df_plot150.loc[:, ('Number of transits')][i] = num_of_trans
    df_plot150.loc[:, ('Max Number of Exp')][i] = max_exp_poss
    
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

for name in Names200:
    pl_bmassj = df200[df200['Name'] == name]['pl_bmassj'].values[0]
    pl_radE = df200[df200['Name'] == name]['pl_radj'].values[0]*11.2
    num_of_trans = df200[df200['Name'] == name]['Number of transits'].values[0]
    max_exp_poss = df200[df200['Name'] == name]['Number of exposures possible'].max()
    i = df_plot200.index[df_plot200['Name'] == name]
    df_plot200.loc[:, ('pl_bmassj')][i] = pl_bmassj
    df_plot200.loc[:, ('pl_radE')][i] = pl_radE
    df_plot200.loc[:, ('Number of transits')][i] = num_of_trans
    df_plot200.loc[:, ('Max Number of Exp')][i] = max_exp_poss
    
    
# Names800 = []
# for name in df800['Name']:
#     Names800.append(name)
# Names2 = list(dict.fromkeys(Names800))
# df_plot800 = pd.DataFrame(Names800)
# df_plot800.columns = ['Name']
# df_plot800['pl_bmassj'] = np.zeros(len(df_plot800))
# df_plot800['pl_radE'] = np.zeros(len(df_plot800))
# df_plot800['Number of transits'] = np.zeros(len(df_plot800))
# df_plot800['Max Number of Exp'] = np.zeros(len(df_plot800))

# for name in Names800:
#     pl_bmassj = df800[df800['Name'] == name]['pl_bmassj'].values[0]
#     pl_radE = df800[df800['Name'] == name]['pl_radj'].values[0]*11.2
#     num_of_trans = df800[df800['Name'] == name]['Number of transits'].values[0]
#     max_exp_poss = df800[df800['Name'] == name]['Number of exposures possible'].max()
#     i = df_plot800.index[df_plot800['Name'] == name]
#     df_plot200.loc[:, ('pl_bmassj')][i] = pl_bmassj
#     df_plot200.loc[:, ('pl_radE')][i] = pl_radE
#     df_plot200.loc[:, ('Number of transits')][i] = num_of_trans
#     df_plot200.loc[:, ('Max Number of Exp')][i] = max_exp_poss
    
    
# for name in Names800:
#     pl_bmassj = df800[df800['Name'] == name]['pl_bmassj'].values[0] 
#     pl_radE = df800[df800['Name'] == name]['pl_radj'].values[0]*11.2
#     num_of_trans = df800[df800['Name'] == name]['Number of transits'].values[0]
#     max_exp_poss = df800[df800['Name'] == name]['Number of exposures possible'].max()
#     i = df_plot800.index[df_plot200['Name'] == name]
#     df_plot800.loc[:, ('pl_bmassj')][i] = pl_bmassj
#     df_plot800.loc[:, ('pl_radE')][i] = pl_radE
#     df_plot800.loc[:, ('Number of transits')][i] = num_of_trans
#     df_plot800.loc[:, ('Max Number of Exp')][i] = max_exp_poss


df_plot100 = df_plot100[df_plot100['Number of transits'] != 0]
df_plot150 = df_plot150[df_plot150['Number of transits'] != 0]
df_plot200 = df_plot200[df_plot200['Number of transits'] != 0]
# df_plot800 = df_plot800[df_plot800['Number of transits'] != 0]

    
plt.clf()

fig1 = plt.figure(figsize=(6,5))
ax1 = fig1.add_subplot(111)

ax1.hist(df_plot100['pl_radE'].values, bins=16, label = 'minimum SN = 100')
ax1.hist(df_plot150['pl_radE'].values, bins=16, label = 'minimum SN = 150')
ax1.hist(df_plot200['pl_radE'].values, bins=16, label = 'minimum SN = 200')
# ax1.scatter(df_plot800['pl_radE'].values, df_plot800['Max Number of Exp'].values, c='black', marker = 's')
ax1.set_xlabel('Planet Radius [Earth radii]')
# ax1.set_xticks(np.arange(17))
ax1.set_ylabel('Number of observable transits')
fig1.legend(loc='upper right')

fig2 = plt.figure(figsize=(6,5))
ax2 = fig2.add_subplot(111)
binwidth = 1
ax2.hist(df_plot100['pl_radE'].values, bins=np.arange(min(df_plot100['pl_radE']), max(df_plot100['pl_radE']) + binwidth, binwidth), label = 'minimum SN = 100')
ax2.hist(df_plot150['pl_radE'].values, bins=np.arange(min(df_plot150['pl_radE']), max(df_plot150['pl_radE']) + binwidth, binwidth), label = 'minimum SN = 150')
ax2.hist(df_plot200['pl_radE'].values, bins=np.arange(min(df_plot200['pl_radE']), max(df_plot200['pl_radE']) + binwidth, binwidth), label = 'minimum SN = 200')
# ax2.scatter(df_plot800['pl_radE'].values, df_plot800['Max Number of Exp'].values, c='black', marker = 's')
ax2.set_xlabel('Planet Radius [Earth radii]')
# ax2.set_xticks(np.arange(17))
# ax2.set_ylabel('Maximum Number of Possible Exposures')
fig2.legend(loc='upper right')

fig1.savefig('Plots/Final_plot_num_of_transits.eps')
fig2.savefig('Plots/Max_Num_of_exp_poss.eps')