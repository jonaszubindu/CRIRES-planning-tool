#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 14:02:08 2020

@author: jonaszbinden
GitHub: jonaszubindu
"""

import pandas as pd
import sys
#import numpy as np
import datetime
from matplotlib import pyplot as plt
from matplotlib import dates as mpl_dates
from classes import Eclipses, load_Eclipses_from_file, Nights
import numpy as np
from astropy import units as u
from astroplan import Observer

""" REDEFINITION:Location and UTC offset Paranal """
# paranal_loc = EarthLocation(lat=-24.627 * u.deg, lon=-70.405 * u.deg, height=2635.43 * u.m)
utcoffset = -4 * u.hour
paranal = Observer.at_site('paranal', timezone='Chile/Continental')

dt = datetime.timedelta(days=1)
# d = datetime.datetime(2020, 4, 1, 0, 0, 0) # choose start day manually
today = datetime.date.today()-datetime.timedelta(days=2)
d = datetime.datetime(today.year, today.month, today.day, 0, 0, 0)


# Max_Delta_days = int(max(Per_max) * 2)
Max_Delta_days = 30 # choose manually for how many days you want to compute the transits


delta_midnight = np.linspace(-12, 12, 1000) * u.hour # defines number of timesteps per 24 hours
d_end = d + dt*Max_Delta_days


Nights_paranal = Nights(d, d_end, LoadFromPickle = 1)
d = d.date()
d = d.isoformat() # start date from which the nights in paranal are calculated
if Nights_paranal.loaded == 0:
    raise Warning('Data for observation range not available: {}-{}'.format(d, Max_Delta_days))



try:
    filename = sys.argv[1]
except IndexError:
    filename = default_file = 'Eclipse_events_processed3_2020-05-20-30.pkl'
    # filename = default_file = 'Eclipse_events_processed3_{}-{}.pkl'.format(d, Max_Delta_days)

Eclipses_List = load_Eclipses_from_file(filename, Max_Delta_days)

Eclipse_Observable = []
Eclipse_Mid_Observable = []

for planet in Eclipses_List:
    if planet.eclipse_observable != []:
        Eclipse_Observable.append(planet)
    elif planet.eclipse_mid_observable != []:
        Eclipse_Mid_Observable.append(planet)
    else:
        Eclipses_List.remove(planet)
        


# Filter data like this:
# df_plot[df_plot['Name'] == 'GJ 9827 b']
# or
# filt = (df_plot['Name'] == 'GJ 9827 b') Generator like!
# df_plot[filt], filt can also be used in .loc[filt, 'obs_time']
# or through changing the indexing
# df_plot.set_index('Name', inplace=True)
# and
# df_plot.loc['GJ 9827 b']
# write df_plot.reset_index(inplace=True) to reset the indexing
# or
# group_df = df_plot.groupby('Name')
# group_df.get_group('GJ 9827 b')









# # sorting data

# group_df = df_plot.groupby('Name')
# planet_names = [iterator[0] for iterator in group_df]
# planet_elem = [[elem] for elem in planet_names]
# for i in planet_elem:
#     i.extend(list(group_df.get_group(i[0])['obs_time'][:]))


# Planets = plot_elem()
# [Planets.add_planet(planet = plot_elem(elem)) for elem in planet_elem]
# for planet in Planets.planets:
#     m = [datetime.datetime.fromisoformat(i) for i in planet.obs_time]
#     planet.obs_time = m

# Nights_paranal_dates['Date'] = pd.to_datetime(Nights_paranal_dates['Date'])

# # grab the obs_time column and split it into time and date
# def split_date_time(Column):
#     new_col_date = []
#     new_col_time = []
#     times = [str.split(i, sep=' ') for i in Column]
#     for obs_time in times:
#         new_col_date.append(obs_time[0])
#         new_col_time.append(obs_time[1])
#     return new_col_date, new_col_time

# df_plot['Date'], df_plot['Eclipse Mid'] = split_date_time(df_plot['obs_time'])
# _, df_plot['Eclipse Begin'] = split_date_time(df_plot['Eclipse Begin'])
# _, df_plot['Eclipse End'] = split_date_time(df_plot['Eclipse End'])

# df_plot.drop(labels='obs_time',axis=1, inplace=True) # drop/remove obs_time column

# Night_groups = df_plot.groupby('Date')
# Night_groups = [Night_groups.get_group(i) for i in df_plot['Date']] # List with all groups of planetary transits that could be observed together!

# Night_groups




# """Plotting Environment"""

# plt.figure(figsize=(30,10))
# y_range = range(len(Planets.planets))
# j = 0
# for planet in Planets.planets:
#     y_planet = [y_range[j] for time in planet.obs_time]
#     j += 1
#     plt.scatter(planet.obs_time, y_planet, marker='|' , s=100, )


# plt.style.use('seaborn')
# plt.tight_layout()

# lims = [Nights_paranal_dates['Date'][i] for i in range(len(Nights_paranal_dates))]
# #lims = lims[0::5]

# plt.xticks(ticks=lims)
# plt.xticks(rotation=70)
# plt.yticks(y_range,planet_names)
# plt.show()



