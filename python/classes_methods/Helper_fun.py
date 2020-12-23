#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 11:53:09 2020

This file contains helper Functions used in Transit_List.py

@author: jonaszbinden
GitHub: jonaszubindu
"""
import pickle
import numpy as np
from classes_methods import Etc_form_class
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_moon, get_sun
from astroplan import moon_phase_angle, Observer, FixedTarget
from astropy.time import Time
import datetime
import time
import json
from json import JSONDecodeError
import requests
import logging
from classes_methods import misc
import astroplan.constraints
import astropy

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import dates as mpl_dates
import copy
import pandas as pd
from astroplan.plots import plot_finder_image
from astropy.visualization import astropy_mpl_style, quantity_support
plt.style.use(astropy_mpl_style)
quantity_support()
import os
import xlsxwriter
from datetimerange import DateTimeRange

""" Location and UTC offset Paranal """
paranal = Observer.at_site('paranal', timezone='Chile/Continental')

##########################################################################################################

def help_fun_logger(orig_fun):
    """ Function to log execution of other functions """
    logging.info('succesfully ran function:{}'.format(orig_fun.__name__))

    return orig_fun

##########################################################################################################

@help_fun_logger
def pickled_items(filename):
    """ Unpickle a file of pickled data. """
    path = os.getcwd() + '/picklefiles/'

    with open(path + filename, "rb") as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError as e:
                logging.exception(e)
                break

##########################################################################################################

@help_fun_logger
def Etc_calculator_Texp(obs_obj, obs_time, snr=500):
    """
        Optimizes NDIT for the S/N minimum defined by ''snr'' for a given DIT for a certain
        observation target ''obs_obj'' at a certain observation time ''obs_time''.

        Parameters
        ----------
        obs_obj : class object
            class object of a observation target.

        obs_time : astropy.time.core.Time
            Time in UTC of observation.

        snr : float
            Minimum S/N ratio that should be reached in a complete exposure

        Returns
        -------
        Exposure_time : float
            Time to compute a full set of exposures with NDIT single exposures of duration DIT.

        DIT : float
            Exposure time for single exposure.

        NDIT : float
            number of single exposures to create one complete set of exposure.

        output : namespace object
            Object containing the full output data from the ETC.

        ETC : etc_form object
            etc_form class instance with input data for the ETC.

    """
    NDIT_opt = 24 # NDIT should optimally be between 16-32
    ETC = Etc_form_class.etc_form(inputtype = "snr-Templ")
    gsmag = obs_obj.star_jmag
    if gsmag < 9.3:
        gsmag = 9.3 # old crires, brighter targets were dimmed by a filter down to 9.3
    moon_target_sep, moon_phase, airmass, _ = airmass_moon_sep_obj_altaz(obs_obj, obs_time) #add moon_target_sep
    ETC.update_etc_form(snr = snr, temperature = obs_obj.star_Teff, brightness = obs_obj.star_jmag, airmass = airmass, moon_target_sep = moon_target_sep, moon_phase = moon_phase, gsmag = gsmag)

    ETC.write_etc_format_file()
    try:
        NDIT, output = ETC.run_etc_calculator(obs_obj.name, obs_time)
    except Exception as e:
        print(type(e))
        if type(e) ==  json.decoder.JSONDecodeError: # catches errors in the etc-form.json input file for the ETC calculator
            # Routine to fix the JSONDecodeError, tries to find the false input.
            ETC.etc_debugger("snr-Templ", obs_obj.name, obs_time, snr = snr, temperature = obs_obj.star_Teff, brightness = obs_obj.star_jmag, airmass = airmass, moon_target_sep = moon_target_sep, moon_phase = moon_phase, gsmag = gsmag)
    # Routine to change ndit to 16-32 and change dit accordingly:
    cycles = 0
    while NDIT < 16 or NDIT > 32:
        Exposure_time = NDIT*ETC.input.timesnr.dit
        DIT_new = Exposure_time/NDIT_opt # determine DIT for NDIT=24
        ETC.input.timesnr.dit = DIT_new
        ETC.write_etc_format_file() # write new DIT into 'etc-form.json'
        logging.info('executed cycle:{}, DIT:{}, NDIT{}'.format(cycles, DIT_new, NDIT))
        try:
            NDIT, output = ETC.run_etc_calculator(obs_obj.name, obs_time) # recalculate the new NDIT
        except Warning:
            raise Warning('DIT seams not feasable input value: {}'.format(DIT_new))
        if NDIT == 0:
            raise Warning('NDIT not available from etc-calculator')
        if cycles > 5:
            raise Warning('too many tries to bring NDIT between 16-32')
        cycles += 1

    DIT = ETC.input.timesnr.dit
    Exposure_time = NDIT*DIT # seconds
    logging.info('Final values: Exposure time:{}, DIT:, NDIT:{}'.format(Exposure_time, DIT, NDIT))

    return Exposure_time, DIT, NDIT, output, ETC

##########################################################################################################

@help_fun_logger
def Etc_calculator_SN(obs_obj, obs_time, ndit, dit):
    """
        Calculates solely the S/N ratio for a given ''dit'' and ''ndit'' for a certain observation
        target ''obs_obj'' at a certain observation time ''obs_time''. CAUTION: This function has not been used
        or tested much yet.

        Parameters
        ----------
        obs_obj : class object
            class object of an observation target.

        obs_time : astropy.time.core.Time
            Time in UTC of observation.

        ndit : int
            Number of frames taken during a full single exposure.
        dit : float
            Exposure time for each frame.

        Returns
        -------
        output : namespace object
            Object containing the full output data from the ETC.

        ETC : etc_form object
            etc_form class instance with input data for the ETC.
    """
    ETC = Etc_form_class.etc_form(inputtype = "ndit-Templ")
    gsmag = obs_obj.star_jmag
    if gsmag < 9.3:
        gsmag = 9.3 # Check again why that one is
    moon_target_sep, moon_phase, airmass, _ = airmass_moon_sep_obj_altaz(obs_obj, obs_time) #add moon_target_sep
    ETC.update_etc_form(temperature = obs_obj.star_Teff, brightness = obs_obj.star_jmag, airmass = airmass, dit = dit, ndit = ndit, moon_target_sep = moon_target_sep, moon_phase=moon_phase, gsmag = gsmag)

    ETC.write_etc_format_file()
    try:
        NDIT, output = ETC.run_etc_calculator(obs_obj.name, obs_time)
    except Exception as e:
        print(type(e))
        if type(e) ==  json.decoder.JSONDecodeError:
            # Routine to fix the JSONDecodeError
            ETC.etc_debugger("ndit-Templ", obs_obj.name, obs_time, temperature = obs_obj.star_Teff, brightness = obs_obj.star_jmag, airmass = airmass, moon_target_sep = moon_target_sep, moon_phase = moon_phase, gsmag = gsmag)

    return output, ETC

##########################################################################################################

@help_fun_logger
def calculate_SN_ratio(sn_data):
    """
        Calculates the median of the signal to noise S/N ratio data ''sn_data''.

        Parameters
        ----------
        sn_data : list
            Containing the S/N ratio data of which the median should be calculated.

        Returns
        -------
        median_SN : float
            Median of the S/N ratio data.
        min_SN : float
            minimum S/N.
        max_SN : float
            maximum S/N.
    """
    #Find the median.
    SN_sorted = np.sort(sn_data)
    length = len(SN_sorted)
    if (length % 2 == 0):
        median_SN = (SN_sorted[(length)//2] + SN_sorted[(length)//2-1]) / 2
    else:
        median_SN = SN_sorted[(length-1)//2]
    min_SN = min(SN_sorted)
    max_SN = max(SN_sorted)

    return median_SN, min_SN, max_SN

##########################################################################################################


@help_fun_logger
def extract_out_data(outputs):
    """
        Function to extract the S/N ratio data from the ''output'' file generated by the ETC.

        Parameters
        ----------
        outputs : namespace object or list
            Object or list of objects containing the full output data from the ETC.

        Returns
        -------
        SN_data : list
            Contains a list of all data from the output(s) of the ETC.
    """
    SN_data = []
    if type(outputs) == list:
        for output in outputs:
            for data in output.data.orders:
                for det in data.detectors:
                    SN_data.extend(det.data.snr.snr.data)
    else:
        output = outputs
        for data in output.data.orders:
            for det in data.detectors:
                SN_data.extend(det.data.snr.snr.data)
    return SN_data

##########################################################################################################

@help_fun_logger
def airmass_moon_sep_obj_altaz(obs_obj, obs_time, location=paranal.location):
    """
        This function calculates the moon target separation, moon phase (moon sun separation), airmass factor and local coordinates to observe
        the object ''obs_obj'' at ''obs_time'' at the location given in ''location'', which is normally paranal.

        Parameters
        ----------
        obs_obj : class object
            instance of class Eclipses, Can also be any other .

        obs_time : astropy.time.core.Time
            Time in UTC of observation.

        location : astropy.coordinates.EarthLocation, optional
            location of the observatory. The default is paranal.location.

        Returns
        -------
        moon_target_sep : float
            angular seperation between moon and target in degrees.

        moon_phase : float
            angular seperation between moon and sun in degrees.

        airmass : float
            Airmass factor at observation time and location.

        obs_altazs : astropy.coordinates.AltAz object
            Azimuth and Altitude in deg at which the object can be observed at the chosen time and location.

    """
    if type(obs_obj) == FixedTarget:
        obs_coor = obs_obj.coord
    else:
        obs_coor = obs_obj.Coordinates.coord
    frame_obs = AltAz(obstime=obs_time, location=location)
    obs_altazs = obs_coor.transform_to(frame_obs)
    airmass = obs_altazs.secz
    moon = get_moon(obs_time).transform_to(frame_obs)
    sun = get_sun(obs_time).transform_to(frame_obs)
    moon_target_sep = moon.separation(obs_altazs) # calculates the moon target separation
    moon_target_sep = moon_target_sep.deg * u.deg
    moon_phase = sun.separation(moon)
    # moon_phase = moon_phase_angle(time=obs_time)
    # moon_phase = moon_phase.to(u.deg)
    z = 90 * u.deg - obs_altazs.alt
    zmoon = 90 * u.deg - moon.alt
    sep_min = np.abs(z-zmoon)
    sep_max = np.abs(z+zmoon)
    if moon_target_sep > sep_max:
        moon_target_sep = sep_max
    elif moon_target_sep < sep_min:
        moon_target_sep = sep_min

    moon_target_sep = moon_target_sep.value
    moon_phase = moon_phase.value
    airmass = airmass.value
    moon_alt = moon.alt.value
    moon_tar_sep = (moon_target_sep, moon_alt)

    return moon_tar_sep, moon_phase, airmass, obs_altazs

##########################################################################################################

@help_fun_logger
def pickle_dumper_objects(filename, Objects):
    """
        Simple function to store class objects or list of class objects ''Objects'' as .pkl file under ''filename''.

        Parameters
        ----------
        filename : str
            filename under which the data should be stored.

        Objects : class object or list of class objects
            class object to store.

    """
    path = os.getcwd() + '/picklefiles/'

    with open(path + filename, 'wb') as out:
        if type(Objects) == list:
            for elem in Objects:
                for key in list(elem.__dict__.keys()):
                    object_attribute = elem.__getattribute__(key)
                    pickle.dump(object_attribute, out, -1)
        else:
            elem = Objects
            for key in list(elem.__dict__.keys()):
                object_attribute = elem.__getattribute__(key)
                pickle.dump(object_attribute, out, -1)
    print('Successfully pickled file {}'.format(filename))


##########################################################################################################

@help_fun_logger
def SN_Transit_Observation_Optimization(eclipse, planet):
    """
        Calculates exactly how many exposures are possible to take during a single transit and adds the data
        to the object Eclipses.eclipse_observable. This function gets only called for single targets cause
        it calls the ETC calculator several times. This is to make sure that changes in the exposure time
        during the transit can be accounted for.

        Parameters
        ----------
        eclipse : object from list Eclipses.eclipse_observable.
            contains all the data about a single transit observation.

        planet : instance of class Eclipses.
            Eclipses instance with all the data about the planet in question.

    """

    """ Checking if eclipse has already been processed """
    try:
        test_empty = eclipse['Number of exposures possible']
        if test_empty != None:
            print('{} has already been processed, skipping...'.format(planet.name))
    except KeyError:
        print('{} gets fed to ETC calculator for best observations'.format(planet.name))
        logging.info('{} gets fed to ETC calculator for best observations'.format(planet.name))
        obs_time = eclipse['Eclipse Mid']['time']

        Transit_dur = planet.transit_duration.to(u.second).value # in seconds

        """ Initial Calculation of Exposure time for Transit: eclipse """
        # for different S/N ratio add argument 'snr'= , to every Etc_calculator_Texp function
        Exposure_time, _, _, output, _ = Etc_calculator_Texp(planet, obs_time) #obtimising NDIT for each single exposure with S/N min = 100 in seconds

        """ Find how many single exposures are possible to take """

        Exposure_times = []
        Exposure_times.append(Exposure_time)
        SN_data_overall = []
        num_exp = 1
        range_obs_times = Exposure_time
        time_between_exposures = 10 # buffer between two exposures in seconds
        n = 1
        median_SN_single_exp = []

        delta = datetime.timedelta(seconds = int(np.ceil(Exposure_time/2 + time_between_exposures)))
        obs_time_up = obs_time + delta
        obs_time_down = obs_time - delta

        while Transit_dur > range_obs_times:
            print('number of exposures: {}'.format(num_exp))
            Exposure_time_up, _, _, output, _ = Etc_calculator_Texp(planet, obs_time_up) # obtimising NDIT for each single exposure with S/N min = 100 in seconds
            Exposure_time_down, _, _, output, _ = Etc_calculator_Texp(planet, obs_time_down) # obtimising NDIT for each single exposure with S/N min = 100 in seconds
            Exposure_times.append(Exposure_time_up)
            Exposure_times.append(Exposure_time_down)
            delta_up = datetime.timedelta(seconds = int(np.ceil(Exposure_time_up + time_between_exposures)))
            delta_down = datetime.timedelta(seconds = int(np.ceil(Exposure_time_down + time_between_exposures)))
            obs_time_up = obs_time_up + delta_up
            obs_time_down = obs_time_down - delta_down
            num_exp = 2*n
            range_obs_times = (obs_time_up-obs_time_down).sec
            n += 1

            # Observations.extend(obs_times)
            SN_data = extract_out_data(output)
            SN_data_overall.extend(SN_data)
            median_SN, _, _ = calculate_SN_ratio(SN_data)
            median_SN_single_exp.append(median_SN)
        Exposure_times.sort()

        if num_exp < 20 and range_obs_times > Transit_dur:
            print('Time to reach 20 exposure exceeded, number of possible exposures: {}'.format(num_exp))
            eclipse['Number of exposures possible'] = num_exp
            eclipse['Comment'] = 'Reaching 20 exposures with S/N = 100 exceeds Transit Duration'
        elif num_exp >= 20: #and range_obs_times <= Transit_dur:
            print('Reached 20 Exposures in: {} seconds during a {} seconds long Transit'.format(np.ceil(range_obs_times), Transit_dur))
            Median_SN, _, _ = calculate_SN_ratio(SN_data_overall)
            median_SN_single_exp.sort()
            eclipse['Number of exposures possible'] = num_exp # estimates the number of exposures possible according to the transit duration and the maximum exposure time calculated reaching 20 exposures
            eclipse['Time necessary to reach 20 exposures [s]'] = np.ceil(range_obs_times)
            eclipse['S/N overall median'] = Median_SN
            eclipse['Minimum S/N'] = median_SN_single_exp[0]
            eclipse['Maximum S/N'] = median_SN_single_exp[-1]
            eclipse['Minimum Exposure Time'] = Exposure_times[0]
            eclipse['Maximum Exposure Time'] = Exposure_times[-1]
            eclipse['List of Exposure Times'] = Exposure_times

##########################################################################################################

def SN_estimate_num_of_exp(eclipse, planet):
    """
        Calculates the exposure time to reach 16 < NDIT < 32 for Transit mid, begin and end
        and from the maximum exposure time estimates the number of possible exposure during the whole transit.
        Also stores data about minimum, maximum and medium S/N ratio and adds the data
        to the object Eclipses.eclipse_observable.

        Parameters
        ----------
        eclipse : object from list Eclipses.eclipse_observable.
            contains all the data about a single transit observation.

        planet : instance of class Eclipses.
            Eclipses instance with all the data about the planet in question.

    """

    """ Checking if eclipse has already been processed """
    try:
        test_empty = eclipse['Number of exposures possible']
        if test_empty != None:
            print('{} has already been processed, skipping...'.format(planet.name))

    except KeyError:
        print('Eclipse {} {} gets fed to ETC calculator for best observations'.format(planet.name, eclipse['obs time']))
        logging.info('{} gets fed to ETC calculator for best observations'.format(planet.name))

        obs_time = eclipse['Eclipse Mid']['time']
        # obs_time_begin = eclipse['Eclipse Begin']['time']
        # obs_time_end = eclipse['Eclipse End']['time']
        Transit_dur = planet.transit_duration.to(u.second).value # in seconds
        # for different S/N ratio add argument 'snr'= , to every Etc_calculator_Texp function
        Exposure_time, _, _, output, _ = Etc_calculator_Texp(planet, obs_time) #obtimising NDIT for each single exposure with S/N min = 100 in seconds
        # Exposure_time_begin, _, _, output, _ = Etc_calculator_Texp(planet, obs_time_begin) #obtimising NDIT for each single exposure with S/N min = 100 in seconds
        # Exposure_time_end, _, _, output, _ = Etc_calculator_Texp(planet, obs_time_end) #obtimising NDIT for each single exposure with S/N min = 100 in seconds
        # Exposure_times = [Exposure_time_begin, Exposure_time_mid, Exposure_time_mid] # get max exposure time
        # Exposure_times.sort()

        SN_data = extract_out_data(output)
        median_SN, min_SN, max_SN = calculate_SN_ratio(SN_data)
        num_exp_possible = int(np.floor(Transit_dur/Exposure_time))
        if num_exp_possible >= 20:
            eclipse['Number of exposures possible'] = num_exp_possible # estimates the number of exposures possible according to the transit duration and the maximum exposure time
            eclipse['S/N median'] = median_SN
            eclipse['Minimum S/N'] = min_SN
            eclipse['Maximum S/N'] = max_SN
            eclipse['Average Exposure Time [s]'] = Exposure_time
            # eclipse['Minimum Exposure Time'] = Exposure_times[0]
            # eclipse['Maximum Exposure Time'] = Exposure_times[-1]
        else:
            eclipse['Number of exposures possible'] = num_exp_possible
            # eclipse['Estimated number of exposures'] = num_exp_possible
            eclipse['S/N median'] = median_SN
            eclipse['Minimum S/N'] = min_SN
            eclipse['Maximum S/N'] = max_SN
            eclipse['Average Exposure Time [s]'] = Exposure_time
            # eclipse['Maximum Exposure Time'] = Exposure_times[-1]

##########################################################################################################

@help_fun_logger
def data_sorting_and_storing(Eclipses_List, filename=None, write_to_csv=1):
    """
        Sorting and storing final data from ''Eclipses_List'' as csv files to ''filename'',
        For now this only works with Eclipses, will include later functionality
        to sort and store more general targets and maybe different functions to plot different kinds of data.
        Might contain more types of output than it has now.

        Parameters
        ----------
        Eclipses_List : list
            Contains all the objects from class Eclipses that have been loaded.

        filename : str, obtional
            If this function is called independent of the observability runs, include filename
            from which Eclipse data should get loaded for sorting and processing.

        Returns
        -------
        ranking : list
            Ranking of the Transits according to (Number of exposures possible)**2 * (number of Transit in computed timespan).

    """
    frame1 = []

    general1 = []
    ranking = []
    num_trans = []

    if type(Eclipses_List) != list:
        Eclipses_List = [Eclipses_List]

    for planet in Eclipses_List:
        if planet.eclipse_observable != []:

            for eclipse in planet.eclipse_observable:
                try:
                    # if eclipse['Number of exposures possible'] == 'Target does not reach 20 exposures':
                    #     pass
                    # else:
                    eclipse1 = copy.deepcopy(eclipse)

                    ecl_data = []
                    ecl_data
                    for key in ['Eclipse Begin', 'Eclipse Mid', 'Eclipse End']:
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
                    ecl_data.rename(index={0: f"Eclipse Begin : {eclipse['Name']}", 1: f"Eclipse Mid : {eclipse['Name']}", 2: f"Eclipse End : {eclipse['Name']}"}, inplace=True)
                    frame1.append(ecl_data)
                except Exception:
                    print(f"{planet.name} could not be processed by the ETC and is excluded from the final processing")


    if frame1 == []:
        raise Warning('There are no eclipses with at least 20 exposures to observe in the selected time frame.')
    general1 = pd.DataFrame(general1)

    if len(general1)>1:
        df_frame1 = pd.concat(frame1, axis=0)
    else:
        df_frame1 = frame1[0]
    df_gen1 = general1

    for planet in Eclipses_List:
        df_per_plan = df_gen1.loc[(df_gen1['Name'] == planet.name) & (df_gen1['Number of exposures possible'] >= 20)]
        if len(df_per_plan) == 0:
            num_exp_mean = 0
        else:
            num_exp_mean = df_per_plan['Number of exposures possible'].sum(axis=0)/len(df_per_plan)
        ranking.append(((len(df_per_plan)*num_exp_mean**2), planet.name))
        num_trans.append((len(df_per_plan), planet.name))

    df_frame1_sorted = df_frame1.sort_values(by='time')
    df_gen_ranking = []
    df_gen_num_trans = []
    for elem, elem1 in zip(ranking, num_trans):
        for name in df_gen1['Name']:
            if elem[1] == name:
                df_gen_ranking.append(elem[0])
                df_gen_num_trans.append(elem1[0])
    df_gen1['rank'] = df_gen_ranking
    df_gen1['Number of transits'] = df_gen_num_trans
    df_gen1.sort_values('Number of exposures possible', inplace=True)

    df_gen1 = df_gen1.reindex(index=df_gen1.index[::-1])
    df_gen1.reset_index(drop=True, inplace=True)
    ranking.sort()
    df_gen1.drop(columns = ['Primary eclipse observable?'], inplace=True)
    if write_to_csv == 1:
        path = os.getcwd() + '/csv_files/'
        if filename == None:
            d = datetime.date.today()
            Max_Delta_days = 'unknown'
            filename = f"Eclipse_events_processed_{d}_{Max_Delta_days}d.csv"
        else:
            file = filename.split('.')[0]
            filename = file + '.csv'
            
        with open(path + filename, 'w') as f:
            df_gen1.to_csv(f, index=False)
            df_frame1_sorted.to_csv(f)
        print(f"Data written to {filename}")

    return ranking, df_gen1, df_frame1, num_trans

##########################################################################################################

@help_fun_logger
def plotting_transit_data(d, Max_Delta_days, ranking, Eclipses_List, Nights, ranked_events=None):
    """
        Plotting final data in ''Eclipses_List'' for the time span given in ''Nights'' or from date ''d'' for ''Max_Delta_days'' days.
        For now this only works for Eclipses, will include later functionality to plot general
        targets and maybe different functions to plot different kinds of data. Might contain more types of output
        than it has now.

        Parameters
        ----------
        Max_Delta_days : int
            Date span to plot the data for.

        ranking : list
            Ranking of the Transits according to (Number of exposures possible)**2 * (number of Transit in computed timespan).

        Eclipses_List : list
            contains Eclipses class objects, which should be plotted.

        Nights : class object
            Class object containing night data for the range of days to plot the results.

        Returns
        -------
        ranking : list
            ranking in reversed order.
    """

    """ Generates the class object Nights and calculates the nights for paranal between d and d_end """
    
    if type(Eclipses_List) != list:
        Eclipses_List = [Eclipses_List]
    
    Nights = Nights(d, Max_Delta_days, LoadFromPickle=0)
    d_orig = d
    # delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    
    if Max_Delta_days > 90:
        for n in range(int(np.floor(Max_Delta_days/90))):
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
                for planet in Eclipses_List:
                    if planet.name == elem[1]:
                        tran_dur = np.float16(planet.transit_duration.to(u.hour))
                        for ecl in planet.eclipse_observable:
                            x_planet = [ecl['Eclipse Begin']['time'].value, ecl['Eclipse End']['time'].value]
                            if ecl['obs time error'] > 1 / 24:
                                ax.plot(x_planet, y_planet, color='red')
                            elif ecl['Number of exposures possible'] < 20:
                                ax.plot(x_planet, y_planet, color='orange')
                            else:
                                ax.plot(x_planet, y_planet, color='blue')
                planet_names.append("{} : {:.3}".format(elem[1], tran_dur))
                j += 1
            
            d = Nights.date[n*90]
            d_end = Nights.date[(n+1)*90]
            lims = [d,d_end]

            plt.xlim(lims)
            plt.xticks(Nights.date[n*90:(n+1)*90], fontsize=22)
            plt.xticks(rotation=70)
            plt.yticks(y_range, planet_names, fontsize=22)
            plt.xlabel('Date', fontsize=24)
            plt.ylabel('Planet name : \nTransit duration [h]', fontsize=24)

            plt.tight_layout()
            plt.show()

            path = os.getcwd() + '/Plots'
            fig.savefig(f"{path}/{d_orig}_{Max_Delta_days}d_{d}-{d_end}-results.eps", bbox_inches='tight', dpi=100)


        """ plotting the last part of unfull months """
        n =int(np.floor(Max_Delta_days/90))
        d = Nights.date[n*90]
        d_end = Nights.date[-1]
        lims = [d,d_end]
        
        plt.clf()
        planet_names = []
        fig = plt.figure(figsize=(1.5*len(Nights.date),1.2*len(ranking))) #figsize=(3*len(lims),1.4*len(ranking))
        ax = fig.add_subplot(111)
        plt.style.use('seaborn-notebook')
        mpl.rc('lines', linewidth=8)
        y_range = range(len(ranking))
        j = 0
        for elem in ranking:
            y_planet = [y_range[j],y_range[j]]
            for planet in Eclipses_List:
                if planet.name == elem[1]:
                    tran_dur = np.float16(planet.transit_duration.to(u.hour))
                    for ecl in planet.eclipse_observable:
                        x_planet = [ecl['Eclipse Begin']['time'].value, ecl['Eclipse End']['time'].value]
                        if ecl['obs time error'] > 1 / 24:
                            ax.plot(x_planet, y_planet, color='orange')
                        elif ecl['Number of exposures possible'] < 20:
                            ax.plot(x_planet, y_planet, color='red')
                        else:
                            ax.plot(x_planet, y_planet, color='blue')
            planet_names.append("{} : {:.3}".format(elem[1], tran_dur))
            j += 1

        plt.xlim(lims)
        plt.xticks(Nights.date[n*90:], fontsize=22)
        plt.xticks(rotation=70)
        plt.yticks(y_range, planet_names, fontsize=22)
        plt.xlabel('Date', fontsize=24)
        plt.ylabel('Planet name : \nTransit duration [h]', fontsize=24)

        plt.tight_layout()
        plt.show()

        path = os.getcwd() + '/Plots'
        fig.savefig(f"{path}/{d_orig}_{Max_Delta_days}d_{d}-{d_end}-results.eps", bbox_inches='tight', dpi=100)
    else:
        planet_names = []
        # plt.clf()
        fig = plt.figure(figsize=(1.5*len(Nights.date),1.2*len(ranking)))
        ax = fig.add_subplot(111)
        plt.style.use('seaborn-notebook')
        mpl.rc('lines', linewidth=8)
        y_range = range(len(ranking))
        j = 0
        for elem in ranking:
            y_planet = [y_range[j],y_range[j]]
            for planet in Eclipses_List:
                if planet.name == elem[1]:
                    tran_dur = np.float16(planet.transit_duration.to(u.hour))
                    for ecl in planet.eclipse_observable:
                        x_planet = [ecl['Eclipse Begin']['time'].value, ecl['Eclipse End']['time'].value]
                        if ecl['obs time error'] > 1 / 24:
                            ax.plot(x_planet, y_planet, color='red')
                        elif ecl['Number of exposures possible'] < 20:
                            ax.plot(x_planet, y_planet, color='orange')
                        else:
                            ax.plot(x_planet, y_planet, color='blue')
            planet_names.append("{} : {:.3}".format(elem[1], tran_dur))
            j += 1
        
        d = Nights.date[0]
        d_end = Nights.date[-1] + datetime.timedelta(days = 1)
        lims = [d,d_end]
                
        # fig.legend(loc='upper left')
        plt.xlim(lims)
        plt.xticks(Nights.date, fontsize=18)
        plt.xticks(rotation=70)
        plt.yticks(y_range, planet_names, fontsize=18)
        plt.xlabel('Date', fontsize=24)
        plt.ylabel('Planet name : \nTransit duration [h]', fontsize=18)
        plt.tight_layout()
        plt.show()

        path = os.getcwd() + '/Plots'
        fig.savefig(f"{path}/{d_orig}_{Max_Delta_days}d_{d}-{d_end}-results.eps", bbox_inches='tight', dpi=100)


    ranking.reverse()
    return ranking

##########################################################################################################

def find_target_image(name):
    """
        Plotting function to get a quick find image for the target.

        Parameters
        ----------
        name : str
            Name of target.

        Returns
        -------
        Figure showing the target relative on the night sky.

    """
    messier1 = FixedTarget.from_name(name)
    ax, hdu = plot_finder_image(messier1)
    plt.show()

##########################################################################################################

def plot_night(date, location, obs_obj, mix_types = 1):
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

    """
    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    # ax2 = ax1.twiny()

    delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    d = Time(date)
    night = d + delta_midnight
    frame_obs = AltAz(obstime=night, location=location)

    moon = get_moon(night)
    sun = get_sun(night)
    moon_altazs = moon.transform_to(frame_obs)
    sun_altazs = sun.transform_to(frame_obs)

    ax1.fill_between(delta_midnight, 0*u.deg, 90*u.deg,
                     sun_altazs.alt < -0*u.deg, color='0.5', zorder=0)
    ax1.fill_between(delta_midnight, 0*u.deg, 90*u.deg,
                     sun_altazs.alt < -18*u.deg, color='k', zorder=0)

    ax1.plot(delta_midnight, sun_altazs.alt, color=[0.75]*3, label='Sun')
    ax1.plot(delta_midnight, moon_altazs.alt, color=[0.75]*3, ls='--', label='Moon')
    
    no_list = 0
    warning = 0

    if type(obs_obj) == list and len(obs_obj) > 1:
        """ plotting for list of objects """
        no_ecl_observable = 0
        for obs_obj in obs_obj:
            if hasattr(obs_obj, 'eclipse_observable'): # If object is related to eclipses
                for eclipse in obs_obj.eclipse_observable:
                    eclipse1 = copy.deepcopy(eclipse)
                    if eclipse1['Eclipse Mid']['time'].datetime.date() == date.date():
                        obs_time = eclipse1['Eclipse Mid']['time']
                        t = obs_time.datetime.time()
                        h = (t.hour + t.minute/60 + t.second/3600) * u.hour
                        delta_eclipse = np.linspace(h - obs_obj.transit_duration.to(u.hour)/2, h + obs_obj.transit_duration.to(u.hour)/2, 100)
                        delta_eclipse_frame = np.linspace(- obs_obj.transit_duration/2, + obs_obj.transit_duration/2, 100)
                        transit = Time(obs_time) + delta_eclipse_frame
                        frame_ecl = AltAz(obstime=transit, location=location)
                        obs_ecl = obs_obj.Coordinates.coord.transform_to(frame_ecl)

                        obs_altazs = obs_obj.Coordinates.coord.transform_to(frame_obs)
                        if any(obs_altazs[obs_altazs.alt > 85 * u.deg]): #checking if target reaches zenith angle
                            warning = 1
                        im = ax1.scatter(delta_midnight, obs_altazs.alt, label=obs_obj.name, lw=0, s=8,
                                cmap='viridis', vmin=-10, vmax=10) # plot candidate
                        if eclipse['Number of exposures possible'] >= 20:
                            ax1.scatter(delta_eclipse, obs_ecl.alt, color='red', lw=3, s=8)
                        no_ecl_observable = 0
                    else:
                        no_ecl_observable = 1

            elif mix_types == 1:
                obs_altazs = obs_obj.Coordinates.coord.transform_to(frame_obs)
                im = ax1.scatter(delta_midnight, obs_altazs.alt, label=obs_obj.name, lw=0, s=8,
                        cmap='viridis', vmin=-10, vmax=10) # Plot candidate

    elif type(obs_obj) == list and len(obs_obj) == 1:
        obs_obj = obs_obj[0]
        no_list = 1
    else:
        pass
        no_list = 1

    if no_list == 1:
        no_ecl_observable = 0
        """ Plotting for single objects """
        if hasattr(obs_obj, 'eclipse_observable'): # If object is related to eclipses
            for eclipse in obs_obj.eclipse_observable:
                eclipse1 = copy.deepcopy(eclipse)
                if eclipse1['Eclipse Mid']['time'].datetime.date() == date.date():
                    obs_time = eclipse1['Eclipse Mid']['time']
                    t = obs_time.datetime.time()
                    h = (t.hour + t.minute/60 + t.second/3600) * u.hour
                    delta_eclipse = np.linspace(h - obs_obj.transit_duration/2, h + obs_obj.transit_duration/2, 100)
                    delta_eclipse_frame = np.linspace(- obs_obj.transit_duration/2, + obs_obj.transit_duration/2, 100)
                    transit = Time(obs_time) + delta_eclipse_frame
                    frame_ecl = AltAz(obstime=transit, location=location)
                    obs_ecl = obs_obj.Coordinates.coord.transform_to(frame_ecl)

                    obs_altazs = obs_obj.Coordinates.coord.transform_to(frame_obs)
                    if any(obs_altazs[obs_altazs.alt > 85 * u.deg]): #checking if target reaches zenith angle
                        warning = 1
                    im = ax1.scatter(delta_midnight, obs_altazs.alt,
                            c=obs_altazs.secz.value, label=obs_obj.name, lw=0, s=8,
                            cmap='viridis', vmin=-10, vmax=10) # plot candidate
                    if eclipse['Number of exposures possible'] >= 20:
                        ax1.scatter(delta_eclipse, obs_ecl.alt, color='red', lw=3, s=8) # plot transit
                else:
                    no_ecl_observable = 1

    if no_ecl_observable == 1:
        obs_altazs = obs_obj.Coordinates.coord.transform_to(frame_obs)
        im = ax1.scatter(delta_midnight, obs_altazs.alt, label=obs_obj.name, lw=0, s=8,
                cmap='viridis', vmin=-10, vmax=10) # Plot candidate
        # phi = np.linspace(0, np.pi, 20)
        # second_xticks = obs_altazs.az[np.int16(np.floor(500*(1+np.tanh(phi/2))))]
        # ax2.set_xlim(obs_altazs.az[0], obs_altazs.az[-1])
        # ax2.set_xticks(second_xticks)
        # ax2.set_xlabel('Azimuth Target [deg]')
    if no_list == 1:
        fig.colorbar(im).set_label('Airmass')
        
    fig.legend(loc='upper left')
    ax1.set_xlim(-12*u.hour, 12*u.hour)
    ax1.set_xticks((np.arange(13)*2-12)*u.hour)
    if warning == 1:
        ax1.set_title(f"{date.date()}, WARNING: Target close to Zenith!")
    else:
        ax1.set_title(f"{date.date()}")
    ax1.set_ylim(0*u.deg, 90*u.deg)
    ax1.set_xlabel('Hours from EDT Midnight')
    ax1.set_ylabel('Altitude [deg]')
    plt.show()

    path = os.getcwd() + '/Plots'
    fig.savefig(f"{path}/{d}-single_night.eps")
    
##########################################################################################################

def xlsx_writer(filename, df_gen, df_frame, ranked_obs_events = None):
    """
        Function to call for customized creation of excel files to store the Exoplanet candidate data.
        This function can be changed in any suitable way to highlight or modify prefered cell formats.
    
        Parameters
        ----------
        filename : str
            filename under which the xlsx file should be stored.
    
        df_gen : pandas DataFrame
            dataframe containing the candidate data to store.
    
        df_frame : pandas DataFrame
            dataframe containing the observation times to store.
    
        ranked_obs_events : pandas DataFrame
            dataframe containing the ranked observation time events to store.
    
        Returns
        -------
        Stores xlsx file to csv_file folder.
    
    """
    path = os.getcwd() + '/csv_files/'
    
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(path + filename.split('.')[0] + '.xlsx', engine='xlsxwriter')
    
    df_frame.reset_index(inplace=True)
    
    workbook = writer.book
    # Set up a format
    # book_format = workbook.add_format(properties={'bold': True, 'font_color': 'red'})
    cell_format = workbook.add_format()
    
    cell_format.set_pattern(1)  # This is optional when using a solid fill.
    cell_format.set_bg_color('red')  # Highlights the background of the cell
    
    # Create a sheet
    worksheet1 = workbook.add_worksheet('Candidates')
    worksheet1.set_column(0, 1, 25)
    worksheet1.set_column(2, 12, 20)
    # worksheet2 = workbook.add_worksheet('Observations')
    # worksheet2.set_column(0, 12, 25)
    if type(ranked_obs_events) == pd.core.frame.DataFrame:
        worksheet3 = workbook.add_worksheet('Ranked Observations')
        worksheet3.set_column(0, 12, 25)
        for col_num, header in enumerate(ranked_obs_events.keys()):
            worksheet3.write(0, col_num, header)
    else:
        pass
    
    
    # Write the headers
    for col_num, header in enumerate(df_gen.keys()):
        worksheet1.write(0, col_num, header)
        
    
    # for col_num, header in enumerate(df_frame.keys()):
    #     worksheet2.write(0, col_num, header)
    
        
    obs_time = []
    # Save the data from the OrderedDict into the excel sheet
    for row_num, row_data in enumerate(df_gen.values):
        for col_num, cell_data in enumerate(row_data):
            if (col_num == 2 and cell_data > 1 / 24) or (col_num == 6 and cell_data < 20):
                obs_time.append(df_gen['obs time'][row_num])
                try:
                    worksheet1.write(row_num + 1, col_num, cell_data, cell_format)
                except TypeError:
                    if type(cell_data) == astropy.time.Time:
                        cell_data = cell_data.value.isoformat()
                    else:
                        cell_data = cell_data.value
                    worksheet1.write(row_num + 1, col_num, cell_data, cell_format)
            else:
                try:
                    worksheet1.write(row_num + 1, col_num, cell_data)
                except TypeError:
                    if type(cell_data) == astropy.time.Time:
                        cell_data = cell_data.value.isoformat()
                    else:
                        cell_data = cell_data.value
                    worksheet1.write(row_num + 1, col_num, cell_data)
    
    # # Save the data from the OrderedDict into the excel sheet
    # for row_num in range(int(len(df_frame.values) / 3)):
    #     row_num = row_num * 3
    #     obs_t = copy.deepcopy(df_frame['time'][row_num+1]) # to compare to transit mid time
    #     for col_num, _ in enumerate(df_frame.values[row_num]):
    #         tim_delta = [np.abs((obs_t - obs).value) == 0.0 for obs in obs_time]
    #         if any(tim_delta):
    #             try:
    #                 worksheet2.write(row_num + 1, col_num, df_frame.values[row_num][col_num], cell_format)
    #                 worksheet2.write(row_num + 2, col_num, df_frame.values[row_num + 1][col_num], cell_format)
    #                 worksheet2.write(row_num + 3, col_num, df_frame.values[row_num + 2][col_num], cell_format)
    #             except TypeError:
    #                 if type(df_frame.values[row_num][1]) == astropy.time.core.Time:
    #                     df_frame.loc[row_num,'time'] = df_frame.values[row_num][col_num].value.isoformat()
    #                     df_frame.loc[row_num + 1, 'time'] = df_frame.values[row_num +1][col_num].value.isoformat()
    #                     df_frame.loc[row_num + 2, 'time'] = df_frame.values[row_num +2][col_num].value.isoformat()
    #                 else:
    #                     df_frame.loc[row_num, 'time'] = df_frame.values[row_num][col_num].value
    #                     df_frame.loc[row_num + 1, 'time'] = df_frame.values[row_num + 1][col_num].value
    #                     df_frame.loc[row_num + 2, 'time'] = df_frame.values[row_num + 2][col_num].value
        
    #                 worksheet2.write(row_num + 1, col_num, df_frame.values[row_num][col_num], cell_format)
    #                 worksheet2.write(row_num + 2, col_num, df_frame.values[row_num + 1][col_num], cell_format)
    #                 worksheet2.write(row_num + 3, col_num, df_frame.values[row_num + 2][col_num], cell_format)
        
    #         else:
    #             try:
    #                 worksheet2.write(row_num + 1, col_num, df_frame.values[row_num][col_num])
    #                 worksheet2.write(row_num + 2, col_num, df_frame.values[row_num + 1][col_num])
    #                 worksheet2.write(row_num + 3, col_num, df_frame.values[row_num + 2][col_num])
    #             except TypeError:
    #                 if type(df_frame.values[row_num][1]) == astropy.time.core.Time:
    #                     df_frame.loc[row_num, 'time'] = df_frame.values[row_num][col_num].value.isoformat()
    #                     df_frame.loc[row_num + 1, 'time'] = df_frame.values[row_num + 1][col_num].value.isoformat()
    #                     df_frame.loc[row_num + 2, 'time'] = df_frame.values[row_num + 2][col_num].value.isoformat()
    #                 else:
    #                     df_frame.loc[row_num, 'time'] = df_frame.values[row_num][col_num].value
    #                     df_frame.loc[row_num + 1, 'time'] = df_frame.values[row_num + 1][col_num].value
    #                     df_frame.loc[row_num + 2, 'time'] = df_frame.values[row_num + 2][col_num].value
    #                 worksheet2.write(row_num + 1, col_num, df_frame.values[row_num][col_num])
    #                 worksheet2.write(row_num + 2, col_num, df_frame.values[row_num + 1][col_num])
    #                 worksheet2.write(row_num + 3, col_num, df_frame.values[row_num + 2][col_num])
    
    if type(ranked_obs_events) == pd.core.frame.DataFrame:
        ranked_obs_events.reset_index(inplace=True)
        ranked_obs_events.drop(columns='level_0', inplace=True)
        for row_num in range(int(len(ranked_obs_events.values) / 3)):
            row_num = row_num * 3
            obs_dat = copy.deepcopy(ranked_obs_events['date'][row_num+1])
            obs_tim = copy.deepcopy(ranked_obs_events['time'][row_num+1]) # to compare to transit mid time
            obs_t = Time(datetime.datetime.combine(obs_dat, obs_tim))
            # print(obs_t)
            for col_num, _ in enumerate(ranked_obs_events.values[row_num]):
                tim_delta = [np.abs((obs_t - obs).value) < 1E-10 for obs in obs_time]
                # print(tim_delta)
                if any(tim_delta):
                    
                    if type(ranked_obs_events.values[row_num][col_num]) != datetime.date and type(ranked_obs_events.values[row_num][col_num]) != datetime.time:
                        worksheet3.write(row_num + 1, col_num, ranked_obs_events.values[row_num][col_num], cell_format)
                        worksheet3.write(row_num + 2, col_num, ranked_obs_events.values[row_num + 1][col_num], cell_format)
                        worksheet3.write(row_num + 3, col_num, ranked_obs_events.values[row_num + 2][col_num], cell_format)
                    else:
                        
                        ranked_obs_events.loc[row_num, 'date'] = ranked_obs_events.loc[row_num, 'date'].isoformat()
                        ranked_obs_events.loc[row_num + 1, 'date'] = ranked_obs_events.loc[row_num + 1, 'date'].isoformat()
                        ranked_obs_events.loc[row_num + 2, 'date'] = ranked_obs_events.loc[row_num + 2, 'date'].isoformat()
                    
                        ranked_obs_events.loc[row_num, 'time'] = ranked_obs_events.loc[row_num, 'time'].isoformat()
                        ranked_obs_events.loc[row_num + 1, 'time'] = ranked_obs_events.loc[row_num + 1, 'time'].isoformat()
                        ranked_obs_events.loc[row_num + 2, 'time'] = ranked_obs_events.loc[row_num + 2, 'time'].isoformat()
        
                        worksheet3.write(row_num + 1, col_num, ranked_obs_events.values[row_num][col_num], cell_format)
                        worksheet3.write(row_num + 2, col_num, ranked_obs_events.values[row_num + 1][col_num], cell_format)
                        worksheet3.write(row_num + 3, col_num, ranked_obs_events.values[row_num + 2][col_num], cell_format)
        
                else:
                    if type(ranked_obs_events.values[row_num][col_num]) != datetime.date and type(ranked_obs_events.values[row_num][col_num]) != datetime.time:
                        worksheet3.write(row_num + 1, col_num, ranked_obs_events.values[row_num][col_num])
                        worksheet3.write(row_num + 2, col_num, ranked_obs_events.values[row_num + 1][col_num])
                        worksheet3.write(row_num + 3, col_num, ranked_obs_events.values[row_num + 2][col_num])
                    else:
                        
                        ranked_obs_events.loc[row_num, 'date'] = ranked_obs_events.loc[row_num, 'date'].isoformat()
                        ranked_obs_events.loc[row_num + 1, 'date'] = ranked_obs_events.loc[row_num + 1, 'date'].isoformat()
                        ranked_obs_events.loc[row_num + 2, 'date'] = ranked_obs_events.loc[row_num + 2, 'date'].isoformat()
                    
                        ranked_obs_events.loc[row_num, 'time'] = ranked_obs_events.loc[row_num, 'time'].isoformat()
                        ranked_obs_events.loc[row_num + 1, 'time'] = ranked_obs_events.loc[row_num + 1, 'time'].isoformat()
                        ranked_obs_events.loc[row_num + 2, 'time'] = ranked_obs_events.loc[row_num + 2, 'time'].isoformat()
        
                        worksheet3.write(row_num + 1, col_num, ranked_obs_events.values[row_num][col_num])
                        worksheet3.write(row_num + 2, col_num, ranked_obs_events.values[row_num + 1][col_num])
                        worksheet3.write(row_num + 3, col_num, ranked_obs_events.values[row_num + 2][col_num])
    else:
        pass
    # print(obs_time)
    # Close the workbook
    workbook.close()


##########################################################################################################

def postprocessing_events(d, Max_Delta_days, Nights, Eclipses_List):
    """
    

    Parameters
    ----------
    filename : string
        filename of pickled file from which the data to process should get retrieved.

    Raises
    ------
    Warning
        raised if the comparison between nights and the number of mutual targets is illogical.

    Returns
    -------
    ranking_dates : list
        contains the ranked observation events, collections of nights in sequence with good targets.
    
    Obs_events : pandas DataFrame
        dataframe containing the ranked events.

    """
    d = datetime.datetime.combine(d, datetime.time(0,0,0))
    # Max_Delta_days = int(filename.split('_')[4][:-5])
    Nights = Nights(d, Max_Delta_days)
    
    ranking, df_gen, df_frame, _ = data_sorting_and_storing(
        Eclipses_List, write_to_csv=0)
    
    
    df_frame_date = []
    df_frame_time = []
    for elem in df_frame['time']:
        df_frame_date.append(elem.value.date())
        df_frame_time.append(elem.value.time())
    
    # df_frame.reset_index(inplace=True)
    for n in range(int(len(df_frame['time'])/3-1)):
        start = df_frame['time'][3*n]
        end = df_frame['time'][(n+1)*3-1]
        if end < start:
            # print(f"we got a problem with {df_frame.loc[3*n]}-{df_frame.loc[(n+1)*3-1]}")
            raise Warning('smth went wrong in handling the times')
            
    df_frame.drop(columns=['time'], inplace=True)
    df_frame.insert(0, 'date', df_frame_date)
    df_frame.insert(1, 'time', df_frame_time)
    # df_frame.sort_values(by=df_frame['date'])  # sorts eclipses in time and date.
    
    
    ranking_dates = []
    date_section = []
    
    
    for n in range(int(len(df_frame['date']) / 3)):
        date_section.append(df_frame[n * 3:(n + 1) * 3])
    # for date_sec in date_section:
    #     ranking_dates.append((len(date_sec) / 3, date_sec))
    
    z = 0
    for n, date in enumerate(Nights.date):
        date_sec = []
        Nights.date[n] = [date,0]
        for date_obj in date_section:
            # date_obj.reset_index(inplace=True)
            if date_obj['date'][0] == Nights.date[n][0].date():
                Nights.date[n][1] = 1
                if n > 0 and Nights.date[n-1][1] == 1:
                    ranking_dates[-1][0] = ranking_dates[-1][0] + 1
                    ranking_dates[-1][1].append(date_obj)
                else:    
                    date_sec.append(date_obj)
        if date_sec == []:
            pass
        else:
            ranking_dates.append([len(date_sec), date_sec])
            z += 1
    
    for date_obJ in ranking_dates:  # for each date
        ran = []
        date_obJ = list(date_obJ)
        for ecl in date_obJ[1]:
            ecl.reset_index(inplace=True)
            start = ecl['date'][0].strftime("%Y-%m-%dT") + ecl['time'][0].strftime("%H:%M:%S-0000")  # begin of eclipse
            end = ecl['date'][2].strftime("%Y-%m-%dT") + ecl['time'][2].strftime("%H:%M:%S-0000")  # end of eclipse
            ran.append((ecl, DateTimeRange(start, end)))
    
        # for range1 in ran:
        #     for range2 in ran:
        #         if range1[1] == range2[1]:
        #             pass
        #         else:
        #             if range1[1].start_datetime < range1[1].end_datetime and range2[1].start_datetime < range2[1].end_datetime:
        #                 inter_sect = range1[1].intersection(range2[1])
        #                 if inter_sect != None:
        #                     date_obJ[0] += -1 / 2
        #             elif range1[1].start_datetime < range1[1].end_datetime:
        #                 print(f"{range1[1].start_datetime} > {range1[1].end_datetime} in {range1[0][['index','date','time']]}")
                        
        #             elif range2[1].start_datetime < range2[1].end_datetime:
        #                 print(f"{range2[1].start_datetime} > {range2[1].end_datetime} in {range2[0][['index','date','time']]}")

    for obs in (ranking_dates):
        final_ranking_per_night = 0
        for single_obs in obs[1]:
            for n in range(len(df_gen)):
                if single_obs.loc[0]['index'].split(':')[1][1:] == df_gen.loc[n]['Name']:
                    if df_gen.loc[n]['Number of exposures possible'] < 20:
                        pass
                    else:
                        final_ranking_per_night += df_gen.loc[n]['Number of exposures possible']
        obs[0] = final_ranking_per_night
    
                
    ranking_dates.sort(key=lambda lis: lis[0], reverse=True)
            
    obs_events = []
    for obs_event in ranking_dates:
        obs_events.extend(obs_event[1])
    Obs_events = pd.concat(obs_events)
    
    return ranking_dates, Obs_events

