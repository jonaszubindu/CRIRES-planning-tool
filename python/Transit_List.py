#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 13:55:10 2020

Original File Transit List. This file contains all the routines to use the functions, classes and class methods to compute observability of targets,
transits of exoplanets and to call the Exposure Time Calculator (ETC) from ESO to calculate S/N signal to noise ratio for observations with CRIRES+.

The different functionalities can be accessed via a menu popping up when running this file.

More functionalities can of course be added and the ETC part can be extended to other instruments used at the VLT
in Paranal, Chile.

Part of this tool could also be extracted and used for other observatories like the first part about observability
under certain constraints.

This tool was created as part of my master thesis: "Planning observations of terrestrial Exoplanets around M type Stars
with CRIRES+" during the peculiar times of Covid-19 in 2020.

Documentation can be found in the README file of this bundle or on GitHub or in my master thesis which can be found here: "COMING SOON".

Examples on how to use:
-----------------------

call Transit_List.py - The tool will welcome you and check if you have internet connection. If your internet connection is working,
you should see now a menu like:

    Choose one of the following options:
 1: Run full transit calculation
 2: run call ETC part for a list of transits
 3: run single transit planning
 4: run single target planning
 5: Plotting data of some result file
Enter number:


The different obtions will start the following procedures
1: Run full transit calculation -

2: run call ETC part for a list of transits -

3: run single transit planning -

4: run single target planning -

5: Plotting data of some result file -




@author: jonaszbinden
GitHub: jonaszubindu
"""


import numpy as np
import pandas as pd
import requests
import matplotlib.pyplot as plt
import os
import json
from json import JSONDecodeError
import logging
import pickle
import sys
import copy
import datetime
import time
import astroplan
from astroplan import download_IERS_A, get_IERS_A_or_workaround, Observer
import astropy
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_moon
from astropy.visualization import astropy_mpl_style, quantity_support
from astropy.coordinates import get_sun
from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive
import classes_methods.Helper_fun as fun
from classes_methods.Helper_fun import help_fun_logger
from classes_methods.classes import Exoplanets, Nights, Eclipses, load_Eclipses_from_file
from classes_methods import csv_file_import
from classes_methods.misc import misc
from astropy.utils import iers

plt.style.use(astropy_mpl_style)
quantity_support()

# from astroquery.exoplanet_orbit_database import ExoplanetOrbitDatabase
# import astroquery.open_exoplanet_catalogue as oec
# from threading import Thread

logging.basicConfig(filename='Transit_List.log', filemode='w',
                    level=logging.DEBUG, format='%(asctime)s-%(levelname)s-%(message)s')

def connect(host='http://exoplanetarchive.ipac.caltech.edu/'):  # Nasa Exoplanet Archive
    """Check Internet Connection to Nasa Exoplanet Archive"""
    req = requests.get(host)  # Python 3.x
    if req.ok:
        print('Connected to {}'.format(host))
    else:
        raise Warning('Check connection to {}, response code:{}'.format(
            host, req.status_code))
    return req


response = connect()






""" Ask for menu input """
k = misc.user_menu(menu=('Run full transit calculation', 'Run call ETC part for a list of transits',
                         'Run single transit planning', 'Run single target planning', 'Plotting data of some result file'))


""" Location and UTC offset Paranal """
paranal = Observer.at_site('paranal', timezone='Chile/Continental')

midnight = datetime.time(0, 0, 0)


""" Altitude constraints definition """
Altcons = astroplan.AltitudeConstraint(min=+30 * u.deg, max=None)

""" Airmass constraints definition """
Airmasscons = astroplan.AirmassConstraint(min=None, max=1.7)

""" Astronomical Nighttime constraints definition: begin and end of each night at paranal as AtNightConstraint.twilight_astronomical """
Night_cons_per_night = astroplan.AtNightConstraint.twilight_astronomical()

""" Moon Constraint """
Mooncons = astroplan.MoonSeparationConstraint(min=+45 * u.deg, max=None)

constraints = [Night_cons_per_night, Altcons, Airmasscons, Mooncons]


##########################################################################################################

if k == 1:

    """ Computes the observability of a list of Candidates """

    d = misc.ask_for_value(
        msg='Enter start date like 2020-06-01 or press enter to use todays date ')
    Max_Delta_days = misc.ask_for_value(
        msg='Enter number of days to compute transits for ')
    Max_Delta_days = int(Max_Delta_days)
    if d == '':
        d = datetime.date.today()
    else:
        d = datetime.date.fromisoformat(d)

    dt = datetime.timedelta(days=1)
    d_end = Max_Delta_days * dt + d
    ETC_calculator = misc.ask_for_value(
        msg='Do you want to call the ETC calculator to process the results S/N ratio? (WARNING : Only works with stable internet connection!) y/n ')

    print(
        f"*** Running full transit analysis for transits between {d} and {d_end} ***")

    """ Update most recent IERS data """
    get_IERS_data = 'yes' # not working at the moment, problem seams to be on IERS side.
    timeoutcount = 0
    
    iers.Conf.iers_auto_url.set('https://datacenter.iers.org/data/9/finals2000A.all') # temporary fix for iers data
    success = 0
    while timeoutcount < 5 and success == 0:
        timeoutcount += 1
        try:
            download_IERS_A(show_progress=True)
            print('IERS data successfully downloaded')
            success = True
        except Exception as e:
            print(e, timeoutcount)
            
    if success == 0:    
        get_IERS_A_or_workaround()
        print('IERS data successfully retrieved')
    
    
#    try:
#        if get_IERS_data == 'yes':
#            download_IERS_A(show_progress=True)
#            print('IERS data successfully downloaded')
#        else:
#            try:
#                get_IERS_A_or_workaround()  # For now, in future always download the most recent ones
#                print('IERS data successfully retrieved')
#            except Exception:
#                download_IERS_A(show_progress=True)
#                print('IERS data successfully downloaded')
#    except Exception:
#        print('No input given, downloading IERS data...')
#        download_IERS_A(show_progress=True)
#        print('IERS data successfully downloaded')
    
    


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

    # initialize list to sort exoplanet candidates and check if data are available to calculate observability.
    Exoplanets = Exoplanets()

    for name in Name_list:
        try:
            """Load Planets from Nasa Exoplanet Archive"""
            Exoplanets.Planet_finder(name)
        except Exception:
            Exoplanets.Fail.append(name)
            print('Planet ' + name + ' added to list "Fail"\n')
            Exoplanets.Fail.append(name)

    """ Check all Planets if they have the necessary properties in the database to process for transit observations """
    Exoplanets.hasproperties()  # work with Parse_list_transits in the same way as in the Transit_Example to retrieve all necessary transit information
    print('Found {} planets in Nasa Archive with Transit data'.format(
        len(Exoplanets.Parse_planets_Nasa)))

    """ Generates the class object Nights and calculates the nights for paranal between d and d_end """
    Nights_paranal = Nights(d, Max_Delta_days, LoadFromPickle=0)

    """ Running methods to compute observability of single transits during the given timespan. """

    obs_time = Time(datetime.datetime.combine(
        Nights_paranal.date[0], midnight))
    Eclipses_List = []
    for planet in Exoplanets.Parse_planets_Nasa:
        Planet = Eclipses(Max_Delta_days, planet)
        Eclipses_List.append(Planet)
        Planet.Observability(obs_time, Nights_paranal,
                             constraints=constraints, check_eclipse=1)

##########################################################################################################

if k == 1 or k == 2:

    """ ETC part to process list of observable candidates """

    if k == 2:
        d = misc.ask_for_value(
            msg='Enter date like 2020-05-28 of the file of transit data you want to use ')
        Max_Delta_days = misc.ask_for_value(
            msg='Enter timespan in days to the appropriate file ')
        filename = 'Eclipse_events_processed_{}_{}d.pkl'.format(
            d, Max_Delta_days)
        misc.wait_for_enter(
            msg=f"Do you want to load file : {filename} to feed to ETC?")
        Eclipses_List = load_Eclipses_from_file(filename, Max_Delta_days)
        Max_Delta_days = int(Max_Delta_days)
        d = datetime.date.fromisoformat(d)



    if k == 2 or ETC_calculator == 'y':

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
            if planet.name == 'WASP-163 b': #Work around for planets which are missing data but did not get filtered. Get the name of the planet that was last called and enter it here. Can contain several names.
                Eclipses_List.remove(planet)
                
        for planet in Eclipses_List:    
            for eclipse in planet.eclipse_observable:
                try:
                    fun.SN_estimate_num_of_exp(eclipse, planet)
                except Warning as w:
                    print(w)
                    print('Something went wrong in:{}:{}, taking next observation...'.format(
                        planet.name, eclipse['obs time']))
                    logging.exception(w)
                    logging.error('Something went wrong in:{}:{}, taking next observation...'.format(
                        planet.name, eclipse['obs time']))
                    break
                except Exception as e:
                    # go back to menu
                    print(e)
                    logging.exception(e)
                    print('\a')
                    print('Catched random exception')

                    print('Shall we save what has been computed so far to a picklefile? You may load that pickle file anytime later and continue from there. Just use the function pickled_items to load manually and investigate the output with next(output) or use load_planets_from_pickle to generate a list of Eclipses instances again like Eclipses_List')
                    save = input('Do you want to save? y/n ')
                    if save == 'y':
                        d = d.isoformat()  # start date from which the nights in paranal are calculated
                        filename = 'Eclipse_events_processed_{}_{}d.pkl'.format(
                            d, Max_Delta_days)
                        fun.pickle_dumper_objects(filename, Eclipses_List)
                        sys.exit()
                    elif save == 'n':
                        sys.exit()
                    else:
                        sys.exit()

    # Store final Data
    try:
        if save == 'n':
            pass
    except Exception:
        filename = 'Eclipse_events_processed_{}_{}d.pkl'.format(
            d, Max_Delta_days)
        fun.pickle_dumper_objects(filename, Eclipses_List)

##########################################################################################################

if k == 3:

    """ Only run single Planet candidate for observability and compute exact number of possible exposures thorugh the ETC """

    d = misc.ask_for_value(
        msg='Enter start date like 2020-06-01 or press enter to use todays date ')
    Max_Delta_days = misc.ask_for_value(
        msg='Enter number of days to compute transits for ')
    Max_Delta_days = int(Max_Delta_days)
    if d == '':
        d = datetime.date.today()
    else:
        d = datetime.date.fromisoformat(d)

    dt = datetime.timedelta(days=1)
    d_end = Max_Delta_days * dt + d
    name = misc.ask_for_value(
        msg='Enter name like KELT-10 b of planet you want check observability ')
    print(
        f"*** Running single transit analysis for transits between {d} and {d_end} ***")
    Planet = NasaExoplanetArchive.query_planet(name, all_columns=True)
    if len(Planet['pl_name']) == 0:
        print('\a')
        name = input('Planet name does not exist, try again ')
        Planet = NasaExoplanetArchive.query_planet(name, all_columns=True)
        if len(Planet['pl_name']) == 0:
            raise ValueError('planet could not be found')

    Exoplanet = Exoplanets()
    Exoplanet.Exoplanets_List_Nasa.append(Planet)
    Exoplanet.hasproperties()
    Planet = Exoplanet.Parse_planets_Nasa[0]

    """ Generates the class object Nights and calculates the nights for paranal between d and d_end """
    Nights_paranal = Nights(d, Max_Delta_days, LoadFromPickle=0)

    """ ETC part for single candidate exposure time calculation and number of possible exposures. """

    obs_time = Time(datetime.datetime.combine(
        Nights_paranal.date[0], midnight))
    Planet = Eclipses(Max_Delta_days, Planet)
    Planet.Observability(obs_time, Nights_paranal,
                         constraints=constraints, check_eclipse=1)

    for eclipse in Planet.eclipse_observable:
        try:
            fun.SN_Transit_Observation_Optimization(eclipse, Planet)
        except Warning as w:
            print(w)
            print('Something went wrong in:{}:{}, taking next observation...'.format(
                Planet.name, eclipse['obs time']))
            logging.exception(w)
            logging.error('Something went wrong in:{}:{}, taking next observation...'.format(
                Planet.name, eclipse['obs time']))
            break
        except Exception as e:
            # go back to menu
            print(e)
            logging.exception(e)
            print('\a')
            print('Catched random exception')

            print('Shall we save what has been computed so far to a picklefile? You may load that pickle file anytime later and continue from there. Just use the function pickled_items to load manually and investigate the output with next(output) or use load_planets_from_pickle to generate a list of Eclipses instances again like Eclipses_List')
            save = input('Do you want to save? y/n ')
            if save == 'y':
                d = d.isoformat()  # start date from which the nights in paranal are calculated
                filename = 'Eclipse_events_processed_{}_{}d.pkl'.format(
                    d, Max_Delta_days)
                fun.pickle_dumper_objects(filename, Eclipses_List)
                sys.exit()
            elif save == 'n':
                sys.exit()
            else:
                sys.exit()

    # Store final Data
    try:
        if save == 'n':
            pass
    except Exception:
        name = Planet.name.split(' ')
        name = name[0] + name[1]
        filename = '{}_events_processed_{}_{}d.pkl'.format(
            name, d, Max_Delta_days)
        fun.pickle_dumper_objects(filename, Planet)

##########################################################################################################

if k == 4:

    """ Run analysis of single target for its observability, if this should also be available for target lists, I need some more information about where these lists would come from. """

    print('Coming soon')
    sys.exit()

    # print('Run through target properties manually or load from name...')
    # name = misc.ask_for_value(msg='name of target ')
    # target = NasaExoplanetArchive.query_star(name)

    # d = misc.ask_for_value(msg='name of target ')
    # Max_Delta_days = misc.ask_for_value(msg='name of target ')

    # st_Teff = misc.ask_for_value(msg='name of target ')
    # st_jmag = misc.ask_for_value(msg='name of target ')
    # number_of_time_steps_per_night = misc.ask_for_value(msg='how many time steps per night, press enter for default value=1000')
    # DIT = misc.ask_for_value(msg='name of target ')
    # NDIT = misc.ask_for_value(msg='name of target ')

    # """ Generates the class object Nights and calculates the nights for paranal between d and d_end """
    # Nights_paranal = Nights(d, Max_Delta_days, LoadFromPickle=0)

    # if delta_midnight == '':
    #     delta_midnight = np.linspace(-12, 12, 1000) * u.hour # defines number of timesteps per 24 hours

    # """
    # Some docs here
    # """

    # target = Target(name, star_Teff, star_jmag)
    # target.target_observability(Nights_paranal, constraints=constraints, delta_midnight=delta_midnight)

    # obs_time = misc.ask_for_value(msg='Choose a observation time for which you want the signal to noisu ratio S/N')
    # try :
    #     fun.SN_Ratio_Target(obs_time, target)
    # except Warning as w:
    #     print(w)
    #     print('Something went wrong in:{}:{}, taking next observation...'.format(Planet.name, eclipse['obs time']))
    #     logging.exception(w)
    #     logging.error('Something went wrong in:{}:{}, taking next observation...'.format(Planet.name, eclipse['obs time']))
    #     break
    # except Exception as e:
    #     #go back to menu
    #     print(e)
    #     logging.exception(e)
    #     print('\a')
    #     print('Catched random exception')

    #     print('Shall we save what has been computed so far to a picklefile? You may load that pickle file anytime later and continue from there. Just use the function pickled_items to load manually and investigate the output with next(output) or use load_planets_from_pickle to generate a list of Eclipses instances again like Eclipses_List')
    #     save = input('Do you want to save? y/n ')
    #     if save == 'y':
    #         d = d.isoformat() # start date from which the nights in paranal are calculated
    #         filename = 'Eclipse_events_processed_{}_{}d.pkl'.format(d, Max_Delta_days)
    #         fun.pickle_dumper_objects(filename, Eclipses_List)
    #         sys.exit()
    #     elif save == 'n':
    #         sys.exit()
    #     else:
    #         sys.exit()

##########################################################################################################


if k == 1 and ETC_calculator == 'n':
    sys.exit()

if k == 1 or k == 2 or k == 3 or k == 4:

    """ Storing data and plotting data """
    if k == 3:
        Eclipses_List = Planet

    ranking, df_gen, df_frame = fun.data_sorting_and_storing(Eclipses_List, filename, write_to_csv=1)
    ranked_events, Obs_events = fun.postprocessing_events(d, Max_Delta_days, Nights, Eclipses_List)
    fun.xlsx_writer(filename, df_gen, df_frame, Obs_events)
    k2 = misc.user_menu(menu=('Plot candidates over full period', 'Plot single night of (mutual) target(s)', 'Get target finder image '))

if k == 5:
    """ Plotting data of some result file """

    k2 = misc.user_menu(menu=('Plot candidates over full period', 'Plot single night of (mutual) target(s)', 'Get target finder image '))

    if k2 == 1 or k2 == 2:
        filename = misc.ask_for_value(
            msg='Enter filename with data to plot:  ')
        if filename.split('.')[-1] == 'pkl':
            pass
        else:
            filename = filename.split('.')[0] + '.pkl'
        
        d = datetime.date.fromisoformat(filename.split('_')[-2])
        Max_Delta_days = int((filename.split('_')[-1].split('.')[0]).split('d')[0])
        Eclipses_List = load_Eclipses_from_file(filename, Max_Delta_days)


if k2 == 1:
    """ Plotting candidates over full period """
    ranking, df_gen, df_frame = fun.data_sorting_and_storing(Eclipses_List, write_to_csv=0)
    ranked_events, Obs_events = fun.postprocessing_events(d, Max_Delta_days, Nights, Eclipses_List)
    fun.xlsx_writer(filename, df_gen, df_frame, Obs_events)
    ranking = fun.plotting_transit_data(d, Max_Delta_days, ranking, Eclipses_List, Nights, ranked_events)

if k2 == 2:
    """ Plot single night of (mutual) target(s) """
    d = misc.ask_for_value(
        msg='Enter date like 2020-05-28 of the night you want to investigate, CAUTION: dates are regarded with respect to UTC ')
    d = datetime.date.fromisoformat(d)
    midn = datetime.time(0, 0, 0)
    d = datetime.datetime.combine(d, midn)
    name = misc.ask_for_value(
        msg=f"Do you want to plot all mutually observable targets for the night of the {d}? Press enter, otherwise write the name of the target you want to plot ")
    found = 0
    if name != '':
        # if name == target.name:
        #     fun.plot_night(d, location = paranal.location, obs_obj = target)
        # else:
        for planet in Eclipses_List:
            if planet.name == name:
                fun.plot_night(d, location=paranal.location, obs_obj=planet)
                found = 1
    if found == 0:
        print('Did not get valid name, plotting all candidates...')
        fun.plot_night(d, location=paranal.location, obs_obj=Eclipses_List)

if k2 == 3:
    """ Get target finder image """
    name = misc.ask_for_value(
        msg="Write the name of the target you want to plot ")
    try:
        fun.find_target_image(name)
    except NameError:
        name = misc.ask_for_value(msg=f"{name} does not exist, try again ")
        fun.find_target_image(name)
