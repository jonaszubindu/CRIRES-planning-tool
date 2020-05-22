#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 13:55:10 2020

Original File Transit List.

@author: jonaszbinden
GitHub: jonaszubindu
"""


import requests
import os
import json
from json import JSONDecodeError
import logging
import pickle

logging.basicConfig(filename = 'Transit_List.log', filemode='w', level=logging.DEBUG, format='%(asctime)s-%(levelname)s-%(message)s')

def connect(host='http://exoplanetarchive.ipac.caltech.edu/'): # Nasa Exoplanet Archive
    """Check Internet Connection to Nasa Exoplanet Archive"""
    req = requests.get(host)  # Python 3.x
    if req.ok: 
        print('Connected to {}'.format(host))
    else:
        raise Warning('Check connection to {}, response code:{}'.format(host,req.status_code))
    return req

response = connect()

import Helper_fun as fun
import astroplan
import astropy
import astropy.units as u
# from astropy import table
from astropy.time import Time

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
from astroplan import Observer
# from threading import Thread

import csv_file_import
from astroplan import download_IERS_A, get_IERS_A_or_workaround
from Helper_fun import help_fun_logger
from classes import Exoplanets, Nights, Eclipses

""" Update most recent IERS data """
get_IERS_data = 'yes'

try:
    if get_IERS_data == 'yes':
        download_IERS_A(show_progress=True)
        print('IERS data successfully downloaded')
    else:
        try:
            get_IERS_A_or_workaround()  # For now, in future always download the most recent ones
            print('IERS data successfully retrieved')
        except Exception:
            download_IERS_A(show_progress=True)
            print('IERS data successfully downloaded')
except Exception:
    print('No input given, downloading IERS data...')
    download_IERS_A(show_progress=True)
    print('IERS data successfully downloaded')


try:
    """ Name_list includes all the exoplanet's names downloaded via Request_Table_NasaExoplanetArchive. """
    Name_list = csv_file_import.main()
except:
    raise Warning('csv file is corrupted or not available')
    logging.error('csv file is corrupted or not available')
    
for name in Name_list: 
    """ Check for any non string-type elements in the list of planet names to parse. """
    if type(name) is not str:
        Warning('Name is not a string, cannot parse {}'.format(name))
        logging.error('Name is not a string, cannot parse {}'.format(name))
    

Exoplanets = Exoplanets()

for name in Name_list[0:5]: #ALL PLANETS
    try:
        """Load Planets from Nasa Exoplanet Archive"""
        Exoplanets.Planet_finder(name)
    except Exception:
        Exoplanets.Fail.append(name)
        print('Planet ' + name + ' added to list "Fail"\n')
        Exoplanets.Fail.append(Name_list.index(name))

""" Check all Planets if they have the necessary properties in the database to process for transit observations """
Exoplanets.hasproperties() # work with Parse_list_transits in the same way as in the Transit_Example to retrieve all necessary transit information
print('Found {} planets in Nasa Archive with Transit data'.format(len(Exoplanets.Parse_planets_Nasa)))

# ----------------------------------------------------------------------------------------------------------------------------
# Constructions for Functions using other ways to add Exoplanets

# def add_manually(self): #define a funcion to add manyally a planet
    #    if #input = 'skip'
    #        pass
    #    else:
    #        #ask the user to input the different columns.
    #        pass

# def add_exoplanets(self, name):  # maybe nicer with tables
    #     """name of planet must be in form strings such that ExoplanetOrbitDatabase can read them"""
    #     try:
    #         Exo = ExoplanetOrbitDatabase.query_planet(
    #             name)  # General Exoplanet Data Source
    #         if not Exo:
    #             Warning('Exoplanet' + name + 'could not be parsed')
    #             self.ExoError.append(name)
    #     except KeyError:
    #         pass


# for planet in cat.findall('.//planet'):
            #     try:
            #         if oec.findvalue(planet, 'name') == name:
            #             print(oec.findvalue(planet, 'name'))
            #             pass  # include adding planet data from the third Archive here
        # else:
        #             Warning('Planet ' + name +
        #                     ' is in no available Database, add manually\n')
        #             # self.add_manually()


# To access for instance the name of a current list element write Exoplanets.Exoplanets_List[1]['NAME']
#-------------------- sky coordinates for the current time write Exoplanets.Exoplanets_List[1]['sky_coord']
# When calculating the sky_coordinates for a specific observation, the DATABASE for the particular EXOPLANET must be UPDATED!

# -----------------------------------------------------------------------------------------------------------------------------

d = datetime.date.today()

""" Definition of maximum number of days to plan observations into the future """
Per_max = [] 
for planet in Exoplanets.Parse_planets_Nasa:
    Per_max.append(np.float64(planet['pl_orbper'] / u.day)) # Maximum orbital period of exoplanets

# Max_Delta_days = int(max(Per_max) * 2)
Max_Delta_days = 30 # choose manually for how many days you want to compute the transits
midnight = datetime.time(0,0,0)
delta_midnight = np.linspace(-12, 12, 1000) * u.hour # defines number of timesteps per 24 hours

Nights_paranal = Nights(d, Max_Delta_days, LoadFromPickle=0)
"""Generates the class object Nights and calculates the nights for paranal between d and d_end"""
# Nights_paranal.Calculate_nights_paranal(delta_midnight, WriteToPickle = 1) 


"""Altitude constraints definition"""
Altcons = astroplan.AltitudeConstraint(min=+30 * u.deg, max=None)  

"""Airmass constraints definition"""
Airmasscons = astroplan.AirmassConstraint(min=None, max=1.7)  

"""Astronomical Nighttime constraints definition: begin and end of each night at paranal as AtNightConstraint.twilight_astronomical"""
Night_cons_per_night = astroplan.AtNightConstraint.twilight_astronomical()

"""
Some docs here
"""
obs_time = Time(datetime.datetime.combine(Nights_paranal.date[0], midnight))
Eclipses_List = []
for planet in Exoplanets.Parse_planets_Nasa:
    Planet = Eclipses(planet, Max_Delta_days)
    Eclipses_List.append(Planet)
    Planet.Observability(obs_time, Nights_paranal, Night_cons_per_night, Altcons, Airmasscons, check_eclipse=1)

"""
For each observable planet:

Do stuff with the input file 'etc-form.json' here:
use: ETC.update_etc_form(**kwargs) from Etc_form_class

Then write the whole file again as a json file: etc-form.json
with ETC.write_etc_format_file()
and run it with Etc_form_class.run_etc_calculator

"""

"""
Calculates the median of the signal to noise ratio achievable in transits that allow around or more than 20 single exposures
during the transit. The Number of exposures possible is added to the list eclipse_observable and eclipse_mid_observable 
for comparison. Each exposure is optimised to have NDIT between 16 and 32 with a minimum S/N = 100. The resulting S/N ratios 
of each exposure are used to compute the overall median. More values like DIT, NDIT, SN of each exposure for each transit 
could be stored as well, but this has not been implemented yet. If one gets stuck in this loop due to internet connection or 
unexpected errors, one may rerun the code from here, instead of rerunning everything again. 
"""
for planet in Eclipses_List:
    
    for eclipse in planet.eclipse_observable:
        try :
            fun.SN_Transit_Observation_Optimization(eclipse, planet)
        except Warning as w:
            print(w)
            print('Something went wrong in:{}:{}, taking next observation...'.format(planet.name, eclipse['obs_time']))
            logging.exception(w)
            logging.error('Something went wrong in:{}:{}, taking next observation...'.format(planet.name, eclipse['obs_time']))
            break
        except Exception as e:
            #go back to menu
            print(e)
            logging.exception(e)
            print('Catched random exception')
            
            print('Shall we save what has been computed so far to a picklefile? You may load that pickle file anytime later and continue from there. Just use the function pickled_items to load manually and investigate the output with next(output) or use load_planets_from_pickle to generate a list of Eclipses instances again like Eclipses_List')
            save = input('Do you want to save? y/n')
            if save == 'y':
                d = d.date()
                d = d.isoformat() # start date from which the nights in paranal are calculated
                filename = 'Eclipse_events-{}-{}.pkl'.format(d, Max_Delta_days)
                fun.pickle_dumper_objects(filename, Eclipses_List)
            elif save == 'n':
                break
            else:
                break
        
# Store final Data
filename = 'Eclipse_events_processed3_{}-{}.pkl'.format(d, Max_Delta_days)
fun.pickle_dumper_objects(filename, Eclipses_List)

