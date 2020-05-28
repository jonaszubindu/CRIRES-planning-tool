#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 14:02:08 2020

Solemnly use to manipulate results, sort and visualize. 
Standard way of usage is sorting candidates according to some ranking and plotting the observable transits during the 
timespan of the data. Change manually which pickled file should get processed. This file has no implementation of
visualization of other targets and was only made to analyze exoplanet transits.


@author: jonaszbinden
GitHub: jonaszubindu
"""
import time
import pandas as pd
import sys
#import numpy as np
import datetime
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import dates as mpl_dates
from classes import Eclipses, Nights, load_Eclipses_from_file
import numpy as np
from astropy import units as u
from astroplan import Observer
import copy


paranal = Observer.at_site('paranal', timezone='Chile/Continental')

dt = datetime.timedelta(days=1)
# d = datetime.datetime(2020, 4, 1, 0, 0, 0) # choose start day manually
d = datetime.date.today()

Max_Delta_days = 7 # choose manually for how many days you want to compute the transits


delta_midnight = np.linspace(-12, 12, 1000) * u.hour # defines number of timesteps per 24 hours

try:
    filename = sys.argv[1]
except IndexError:
    # filename = 'Eclipse_events_processed_2020-05-28-365.pkl'
    filename = 'Eclipse_events_processed_2020-05-28-7.pkl'
    # filename = default_file = 'Eclipse_events_processed_{}-{}.pkl'.format(d, Max_Delta_days)

Eclipses_List2 = load_Eclipses_from_file(filename, Max_Delta_days)

frame1 = []

general1 = []
ranking = []
j = 0
for planet in Eclipses_List2:
    if planet.eclipse_observable != []:
        
        for eclipse in planet.eclipse_observable:
            if eclipse['Number of exposures possible'] == 'Target does not reach 20 exposures':
                pass
            else:
                eclipse1 = copy.deepcopy(eclipse)
    
                ecl_data = []
                ecl_data
                for key in ['Eclipse Begin', 'Eclipse Mid','Eclipse End']:
                    eclipse1[key]['az'] = np.float32(eclipse1[key]['az'])
                    eclipse1[key]['alt'] = np.float32(eclipse1[key]['alt'])
                    ecl_data.append(eclipse1[key])
                    eclipse1.pop(key)
                try :
                    eclipse1.pop('List of Exposure Times')
                except KeyError:
                    pass
                general1.append(eclipse1)
                ecl_data = pd.DataFrame(ecl_data)
                ecl_data.rename(index={0: f"Eclipse Begin : {eclipse['Name']}", 1: "Eclipse Mid", 2: "Eclipse End"}, inplace=True)
                frame1.append(ecl_data)

general1 = pd.DataFrame(general1)

if len(general1)>1:
    df_frame1 = pd.concat(frame1, axis=0)
else:
    df_frame1 = frame1[0]
df_gen1 = general1
df_gen1.sort_values('Number of exposures possible', inplace=True)

for planet in Eclipses_List2:
    df_per_plan = df_gen1.loc[df_gen1['Name'] == planet.name]
    if len(df_per_plan) == 0:
        pass
    else:
        num_exp_mean = df_per_plan['Number of exposures possible'].sum(axis=0)/len(df_per_plan)
        ranking.append(((len(df_per_plan)*num_exp_mean**2),planet.name))

ranking.sort()
df_gen1 = df_gen1.reindex(index=df_gen1.index[::-1])
df_gen1.reset_index(drop=True, inplace=True)

file = filename.split('.')[0]
filename = file + '.csv'
with open(filename, 'w') as f:
    df_gen1.to_csv(f, index=False)
    df_frame1.to_csv(f)
print(f"Data written to {filename}")
    
Nights_paranal = Nights(d, Max_Delta_days)

# Nights_paranal.Calculate_nights_paranal(delta_midnight)


""" Plotting Environment """

if Max_Delta_days > 90:
    for n in range(int(np.floor(Max_Delta_days/90))):
        plt.clf()
        planet_names = []
        fig = plt.figure(figsize=(120,1.2*len(ranking)))
        ax = fig.add_subplot(111)
        plt.style.use('seaborn-notebook')
        mpl.rc('lines', linewidth=6)
        y_range = range(len(ranking))
        j = 0
        for elem in ranking:
            y_planet = [y_range[j],y_range[j]]
            for planet in Eclipses_List2:
                if planet.name == elem[1]:
                    tran_dur = np.float16(planet.transit_duration.to(u.hour))
                    for ecl in planet.eclipse_observable:
                        x_planet = [ecl['Eclipse Begin']['time'].value, ecl['Eclipse End']['time'].value]
                        x_planet_dates = mpl_dates.date2num(x_planet)
                        ax.plot(x_planet, y_planet)
            planet_names.append("{} : {:.3}".format(elem[1], tran_dur))
            j += 1
        
        d = Nights_paranal.date[n*90]
        d_end = Nights_paranal.date[(n+1)*90]
        lims = [d,d_end]
        
        plt.xlim(lims)
        plt.xticks(Nights_paranal.date[n*90:(n+1)*90], fontsize=22)
        plt.xticks(rotation=70)
        plt.yticks(y_range, planet_names, fontsize=22)
        plt.xlabel('Date', fontsize=24)
        plt.ylabel('Planet name : \nTransit duration [h]', fontsize=24)
        
        # plt.tight_layout()
        plt.show()
        fig.savefig(f"{d}-{d_end}-results.eps")
        
        
    """ plotting the last part of unfull months """
    
    plt.clf()
    planet_names = []
    fig = plt.figure(figsize=(120,1.2*len(ranking)))
    ax = fig.add_subplot(111)
    plt.style.use('seaborn-notebook')
    mpl.rc('lines', linewidth=6)
    y_range = range(len(ranking))
    j = 0
    for elem in ranking:
        y_planet = [y_range[j],y_range[j]]
        for planet in Eclipses_List2:
            if planet.name == elem[1]:
                tran_dur = np.float16(planet.transit_duration.to(u.hour))
                for ecl in planet.eclipse_observable:
                    x_planet = [ecl['Eclipse Begin']['time'].value, ecl['Eclipse End']['time'].value]
                    x_planet_dates = mpl_dates.date2num(x_planet)
                    ax.plot(x_planet, y_planet)
        planet_names.append("{} : {:.3}".format(elem[1], tran_dur))
        j += 1
    
    d = Nights_paranal.date[(n+1)*90]
    d_end = Nights_paranal.date[-1]
    lims = [d,d_end]
    
    plt.xlim(lims)
    plt.xticks(Nights_paranal.date[n*90:], fontsize=22)
    plt.xticks(rotation=70)
    plt.yticks(y_range, planet_names, fontsize=22)
    plt.xlabel('Date', fontsize=24)
    plt.ylabel('Planet name : \nTransit duration [h]', fontsize=24)
    
    # plt.tight_layout()
    plt.show()
    fig.savefig(f"{d}-{d_end}-results.eps")
else:
    planet_names = []
    plt.clf()
    fig = plt.figure(figsize=(1.5*len(Nights_paranal.date),1.2*len(ranking)))
    ax = fig.add_subplot(111)
    plt.style.use('seaborn-notebook')
    mpl.rc('lines', linewidth=6)
    y_range = range(len(ranking))
    j = 0
    for elem in ranking:
        y_planet = [y_range[j],y_range[j]]
        for planet in Eclipses_List2:
            if planet.name == elem[1]:
                tran_dur = np.float16(planet.transit_duration.to(u.hour))
                for ecl in planet.eclipse_observable:
                    x_planet = [ecl['Eclipse Begin']['time'].value, ecl['Eclipse End']['time'].value]
                    x_planet_dates = mpl_dates.date2num(x_planet)
                    ax.plot(x_planet, y_planet)
        planet_names.append("{} : {:.3}".format(elem[1], tran_dur))
        j += 1

    d = Nights_paranal.date[0]
    d_end = Nights_paranal.date[-1]
    lims = [d,d_end]
    
    plt.xlim(lims)
    plt.xticks(Nights_paranal.date, fontsize=22)
    plt.xticks(rotation=70)
    plt.yticks(y_range, planet_names, fontsize=22)
    plt.xlabel('Date', fontsize=24)
    plt.ylabel('Planet name : \nTransit duration [h]', fontsize=24)
    
    # plt.tight_layout()
    plt.show()
    fig.savefig(f"{d}-{d_end}-results.eps")


ranking.reverse()
# # Store final Data
# filename = 'Eclipse_events_processed_{}-{}.pkl'.format(d, Max_Delta_days)
# fun.pickle_dumper_objects(filename, Eclipses_List)

