#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 16:18:01 2021

@author: jonaszbinden
"""
import astroplan
import astropy
import astropy.units as u
# from astropy import table
from astropy.time import Time, TimeDelta
import math as m

#from astroplan import EclipsingSystem
#import astropy.coordinates
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_moon
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style, quantity_support
plt.style.use(astropy_mpl_style)
quantity_support()
from astropy.coordinates import get_sun
from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive
# from astroquery.exoplanet_orbit_database import ExoplanetOrbitDatabase
# import astroquery.open_exoplanet_catalogue as oec
import datetime
import time
from astroplan import Observer, FixedTarget
# from threading import Thread


from astroplan import download_IERS_A, get_IERS_A_or_workaround
from classes_methods.Helper_fun import help_fun_logger


import pickle
import classes_methods.Helper_fun as fun
from classes_methods.Helper_fun import help_fun_logger
import logging
import pandas as pd
from classes_methods.classes import load_Eclipses_from_file, Nights 

filename_pickles_1 = 'Eclipse_events_processed_2021-10-01_185d_1.pkl'
filename_pickles_2 = 'Eclipse_events_processed_2021-10-01_185d_2.pkl'
filenames_pickles = [filename_pickles_1, filename_pickles_2]

filename_csv = '/Users/jonaszbinden/Desktop/Target Selection Paper/Results/P108_metric_2D/P108_targets_sorted.xlsx'

df = pd.read_excel(filename_csv)
Max_Delta_days = 185 #int((filename_csv.split('_')[-1].split('.')[0]).split('d')[0])
d = datetime.date.fromisoformat(filename_pickles_1.split('_')[-3])
name_list = df['Name']
name_list = name_list.drop_duplicates()
Eclipses_List_new = []
for filename in filenames_pickles:
    Eclipses_List = load_Eclipses_from_file(filename, Max_Delta_days, path = '/Users/jonaszbinden/Desktop/Target Selection Paper/Results/P108_metric_2D/original_files/')
    time.sleep(3)
    for name in name_list:
        for planet in Eclipses_List:
            if planet.name == name:
                print(planet.name, 'here')
                Eclipses_List_new.append(planet)
Eclipses_List_new.pop(-2)
delta_midnight = np.linspace(-12, 12, 1000) * u.hour
Nights_paranal = Nights(d, Max_Delta_days, LoadFromPickle=0)
Nights_paranal.Calculate_nights_paranal(delta_midnight)
beg_end_night = []
for night in Nights_paranal.night:
    beg_end_night.append({'beg': night[0], 'end': night[-1]})

ranking, df_gen, df_frame, _ = fun.data_sorting_and_storing(Eclipses_List_new, write_to_csv=1)
ranked_events, Obs_events = fun.postprocessing_events(d, Max_Delta_days, Nights, Eclipses_List_new)
fun.xlsx_writer(filename, df_gen, df_frame, Obs_events)
ranking = fun.plotting_transit_data(d, Max_Delta_days, ranking, Eclipses_List_new, Nights, ranked_events)


            
            
    
    