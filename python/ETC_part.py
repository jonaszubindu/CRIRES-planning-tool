#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 13:39:37 2020

This is to only run the ETC calculator part with objects loaded from pickle file.

@author: jonaszbinden
"""
import numpy as np
import Helper_fun as fun
# import Etc_form_class
import pickle
# import numpy as np
# import Etc_form_class
from astropy import units as u
from astroplan import Observer
import datetime
# import time
# import json
# from json import JSONDecodeError
# import requests
import logging
from classes import Eclipses, load_Eclipses_from_file

logging.basicConfig(filename = 'ETC_call_Transits.log', filemode='w', level=logging.DEBUG, format='%(asctime)s-%(levelname)s-%(message)s')


""" REDEFINITION: UTC offset Paranal """

dt = datetime.timedelta(days=1)
# d = datetime.datetime(2020, 4, 1, 0, 0, 0) # choose start day manually
today = datetime.date.today()-datetime.timedelta(days=1)
d = datetime.datetime(today.year, today.month, today.day, 0, 0, 0)


# Max_Delta_days = int(max(Per_max) * 2)
Max_Delta_days = 30 # choose manually for how many days you want to compute the transits


delta_midnight = np.linspace(-12, 12, 1000) * u.hour # defines number of timesteps per 24 hours
d_end = d + dt*Max_Delta_days

# """
# For each observable planet:

# Do stuff with the input file 'etc-form.json' here:
# use: ETC.update_etc_form(**kwargs) from Etc_form_class

# Then write the whole file again as a json file: etc-form.json
# with ETC.write_etc_format_file()
# and run it with Etc_form_class.run_etc_calculator



# Load planets from filename.pkl
d = d.date()
d = d.isoformat() # start date from which the nights in paranal are calculated
filename = 'Eclipse_events_processed_{}-{}.pkl'.format(d, Max_Delta_days)
Eclipses_List = load_Eclipses_from_file(filename, Max_Delta_days)

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

    

