#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 14:02:08 2020

@author: jonaszbinden
"""

import pandas as pd
import sys
#import numpy as np
import datetime
from matplotlib import pyplot as plt
from matplotlib import dates as mpl_dates

try:
    file_name = sys.argv[1]
except IndexError:
    file_name = default_file = 'test_obs_file3'

try:
    Nights_paranal_table = pd.read_csv('Nights_time_at_paranal.csv', index_col = 0)
    Nights_paranal_dates = pd.read_csv('dates_timespan.csv', index_col = 0)
except Exception:
    print('Error: files for the observation timerange not found')


try:
    df_plot = pd.read_csv(file_name + '.csv', index_col = 0)
except FileNotFoundError:
    df_plot = pd.read_csv(default_file, index_col = 0)
    print("can't parse {}.csv, taking default file: {}".format(file_name, default_file))

df_plot = df_plot[df_plot['Primary eclipse observable?']==True]


df_plot = df_plot.reset_index(drop=True)

class Obs_times:

    def __init__(self, obs_time = None):
        # self.name = obs_time[0]
        # self.obs_time = obs_time[1:]
        # or equivalently for pedagogic purposes:
        if obs_time is None:
            self.name = []
            self.obs_time = []
        else:
            try:
                self.__setattr__('name', obs_time[0])
                self.__setattr__('obs_time', obs_time[1:])
            except:
                raise Warning('Cannot create attribute {} for object {}'.format(obs_time[0],self))


class plot_elem(Obs_times):

    def __init__(self, info = ['Planet name','Transit Observable?'], planet_list = None):
        super().__init__(info)
        if planet_list is None:
            self.planets = []
        else:
            self.planets = planet_list

    def add_planet(self, planet):
        if planet not in self.planets:
            self.planets.append(planet)

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


# sorting data

group_df = df_plot.groupby('Name')
planet_names = [iterator[0] for iterator in group_df]
planet_elem = [[elem] for elem in planet_names]
for i in planet_elem:
    i.extend(list(group_df.get_group(i[0])['obs_time'][:]))


Planets = plot_elem()
[Planets.add_planet(planet = plot_elem(elem)) for elem in planet_elem]
for planet in Planets.planets:
    m = [datetime.datetime.fromisoformat(i) for i in planet.obs_time]
    planet.obs_time = m

Nights_paranal_dates['Date'] = pd.to_datetime(Nights_paranal_dates['Date'])

# grab the obs_time column and split it into time and date
def split_date_time(Column):
    new_col_date = []
    new_col_time = []
    times = [str.split(i, sep=' ') for i in Column]
    for obs_time in times:
        new_col_date.append(obs_time[0])
        new_col_time.append(obs_time[1])
    return new_col_date, new_col_time

df_plot['Date'], df_plot['Eclipse Mid'] = split_date_time(df_plot['obs_time'])
_, df_plot['Eclipse Begin'] = split_date_time(df_plot['Eclipse Begin'])
_, df_plot['Eclipse End'] = split_date_time(df_plot['Eclipse End'])

df_plot.drop(labels='obs_time',axis=1, inplace=True) # drop/remove obs_time column

Night_groups = df_plot.groupby('Date')
Night_groups = [Night_groups.get_group(i) for i in df_plot['Date']] # List with all groups of planetary transits that could be observed together!

Night_groups




"""Plotting Environment"""

plt.figure(figsize=(30,10))
y_range = range(len(Planets.planets))
j = 0
for planet in Planets.planets:
    y_planet = [y_range[j] for time in planet.obs_time]
    j += 1
    plt.scatter(planet.obs_time, y_planet, marker='|' , s=100, )


plt.style.use('seaborn')
plt.tight_layout()

lims = [Nights_paranal_dates['Date'][i] for i in range(len(Nights_paranal_dates))]
#lims = lims[0::5]

plt.xticks(ticks=lims)
plt.xticks(rotation=70)
plt.yticks(y_range,planet_names)
plt.show()



