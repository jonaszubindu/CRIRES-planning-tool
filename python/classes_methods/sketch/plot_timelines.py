#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 14:02:08 2020

This file contains the last part of Transit_List.py and is to play around with the data by hand. Solemnly use to manipulate results, sort and visualize. 
Standard way of usage is sorting candidates according to some ranking and plotting the observable transits during the 
timespan of the computed data. Change manually which pickled file should get processed. 


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

Max_Delta_days = 365 # choose manually for how many days you want to compute the transits

try:
    filename = sys.argv[1]
except IndexError:
    # filename = 'Eclipse_events_processed_2020-05-28_365d.pkl'
    filename = 'Eclipse_events_processed_2020-05-29_365d.pkl'
    # filename = default_file = 'Eclipse_events_processed_{}_{}d.pkl'.format(d, Max_Delta_days)

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
if frame1 == []:
    raise Warning('There are no eclipses with at least 20 exposures to observe in the selected time frame.')
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


""" Plotting Environment """

if Max_Delta_days > 90:
    for n in range(int(np.floor(Max_Delta_days/90)-1)):
        plt.clf()
        planet_names = []
        fig = plt.figure(figsize=(120,1.2*len(ranking)))
        ax = fig.add_subplot(111)
        plt.style.use('seaborn-notebook')
        mpl.rc('lines', linewidth=8)
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
                        ax.plot(x_planet, y_planet, color='blue')
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
    mpl.rc('lines', linewidth=8)
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
                    ax.plot(x_planet, y_planet, color='blue')
        planet_names.append("{} : {:.3}".format(elem[1], tran_dur))
        j += 1
    
    d = Nights_paranal.date[(n+1)*90]
    d_end = Nights_paranal.date[-1]
    lims = [d,d_end]
    
    plt.xlim(lims)
    plt.xticks(Nights_paranal.date[(n+1)*90:], fontsize=22)
    plt.xticks(rotation=70)
    plt.yticks(y_range, planet_names, fontsize=22)
    plt.xlabel('Date', fontsize=24)
    plt.ylabel('Planet name : \nTransit duration [h]', fontsize=24)
    
    # plt.tight_layout()
    plt.show()
    fig.savefig(f"{d}-{d_end}-results.eps")
else:
    
    """ Plots data that do not reach 90 days """
    
    planet_names = []
    plt.clf()
    fig = plt.figure(figsize=(1.5*len(Nights_paranal.date),1.2*len(ranking)))
    ax = fig.add_subplot(111)
    plt.style.use('seaborn-notebook')
    mpl.rc('lines', linewidth=8)
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
                    ax.plot(x_planet, y_planet, color='blue')
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
    fig.savefig("plots/" + f"{d}-{d_end}-results.eps")


ranking.reverse()


#######################################################################################################

"""
Plots the targets of a single night, depending where on the night sky they appear at which time.

Parameters
----------
date : datetime or astropy.time.Time
    date for which the night shall be plottet.
location : astroplan.Observer.location or astropy.coordinates.EarthLocation
    Location of observatory.
obs_obj : class object or list of class objects
    Class object containing information about coordinates, observation times.
mix_types : int (obtional)
    set to zero if you want to only compare mutual transits in the same night.

Returns
-------
None.

"""


fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

delta_midnight = np.linspace(-12, 12, 1000)*u.hour
date = Time(date)
night = date + delta_midnight
frame_obs = AltAz(obstime=night, location=location)
if type(obs_obj) == list:
    for obs_obj in obs_obj:
        if hasattr(obs_obj, 'eclipse_observable'): # If object is related to eclipses
            for eclipse in obs_obj.eclipse_observable:
                if eclipse['Number of exposures possible'] == 'Target does not reach 20 exposures':
                    pass
                else:
                    eclipse1 = copy.deepcopy(eclipse)
                    if eclipse1['Eclipse Mid']['time'].datetime.date() == date:
                        obs_time = eclipse1['Eclipse Mid']['time']
                        t = obs_time.datetime.time()
                        h = (t.hour + t.minute/60 + t.second/3600) * u.hour
                        delta_eclipse = np.linspace(h - planet.transit_duration/2, h * u.hour - planet.transit_duration/2, 100)
                        
                        frame_ecl = frame_obs = AltAz(obstime=obs_time, location=location)
                        obs_ecl = obs_obj.transform_to(frame_ecl)
                        ax1.scatter(delta_eclipse, obs_ecl.alt, color='red', marker=10)
                        obs_altazs = obs_obj.transform_to(frame_obs)
                        ax1.scatter(delta_midnight, obs_altazs.alt,
                                c=obs_altazs.secz, label=obs_obj.name, lw=0, s=8,
                                cmap='viridis') # Plot candidate
        elif mix_types == 1:
            obs_altazs = obs_obj.transform_to(frame_obs)
            ax1.scatter(delta_midnight, obs_altazs.alt,
                    c=obs_altazs.secz, label=obs_obj.name, lw=0, s=8,
                    cmap='viridis') # Plot candidate
        
elif type(obs_obj) != list:
    if hasattr(obs_obj, 'eclipse_observable'): # If object is related to eclipses
        for eclipse in obs_obj.eclipse_observable:
            if eclipse['Number of exposures possible'] == 'Target does not reach 20 exposures':
                pass
            else:
                eclipse1 = copy.deepcopy(eclipse)
                if eclipse1['Eclipse Mid']['time'].datetime.date() == date:
                    obs_time = eclipse1['Eclipse Mid']['time']
                    t = obs_time.datetime.time()
                    h = (t.hour + t.minute/60 + t.second/3600) * u.hour
                    delta_eclipse = np.linspace(h - planet.transit_duration/2, h * u.hour - planet.transit_duration/2, 100)
                    
                    frame_ecl = frame_obs = AltAz(obstime=obs_time, location=location)
                    obs_ecl = obs_obj.transform_to(frame_ecl)
                    ax1.scatter(delta_eclipse, obs_ecl.alt, color='red', marker=10)
                    obs_altazs = obs_obj.transform_to(frame_obs)
                    ax1.scatter(delta_midnight, obs_altazs.alt,
                            c=obs_altazs.secz, label=obs_obj.name, lw=0, s=8,
                            cmap='viridis') # Plot candidate
                        
    else:                
        obs_altazs = obs_obj.transform_to(frame_obs)
        ax1.scatter(delta_midnight, obs_altazs.alt,
                c=obs_altazs.secz, label=obs_obj.name, lw=0, s=8,
                cmap='viridis') # Plot candidate
    second_xticks = obs_altazs.az[0::50]
    ax2.xticks(second_xticks)
    
moon = get_moon(night)
sun = get_sun(night)
moon_altazs = moon.transform_to(frame_obs)
sun_altazs = sun.transform_to(frame_obs)


ax1.plot(delta_midnight, sun_altazs.alt, color='r', label='Sun')
ax1.plot(delta_midnight, moon_altazs.alt, color=[0.75]*3, ls='--', label='Moon')

plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg,
                 sun_altazs.alt < -0*u.deg, color='0.5', zorder=0)
plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg,
                 sun_altazs.alt < -18*u.deg, color='k', zorder=0)

plt.colorbar().set_label('Airmass')
plt.legend(loc='upper left')
ax1.xlim(-12*u.hour, 12*u.hour)
ax1.xticks((np.arange(13)*2-12)*u.hour)

ax1.ylim(0*u.deg, 90*u.deg)
ax1.xlabel('Hours from EDT Midnight')
ax1.ylabel('Altitude [deg]')
plt.show()


# # Store final Data
# filename = 'Eclipse_events_processed_{}_{}d.pkl'.format(d, Max_Delta_days)
# fun.pickle_dumper_objects(filename, Eclipses_List)

