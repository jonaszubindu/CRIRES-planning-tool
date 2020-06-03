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
import Etc_form_class
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_moon, get_sun
from astroplan import moon_phase_angle, Observer, FixedTarget
import datetime
import time
import json
from json import JSONDecodeError
import requests
import logging
from misc import misc
import astroplan.constraints

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import dates as mpl_dates
import copy
import pandas as pd

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
    with open(filename, "rb") as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError as e:
                logging.exception(e)
                break

##########################################################################################################

@help_fun_logger  
def Etc_calculator_Texp(obs_obj, obs_time, snr=100):
    """
    Optimizes NDIT for the S/N minimum defined by snr for a given dit for a certain 
    observation target at a certain observation time.

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
    output : namespace object
        Object containing the full output data from the ETC calculator.

    """
    NDIT_opt = 24 # NDIT should optimally be between 16-32
    ETC = Etc_form_class.etc_form(inputtype = "snr")
    if snr != 100:
        ETC.update_etc_form(snr=snr)
    
    moon_target_sep, moon_phase, airmass, _ = airmass_moon_sep_obj_altaz(obs_obj, obs_time) #add moon_target_sep
    ETC.update_etc_form(temperature = float(obs_obj.star_Teff/u.K), brightness = obs_obj.star_jmag, airmass = airmass, moon_target_sep = moon_target_sep, moon_phase = moon_phase)
 
    ETC.write_etc_format_file()
    try:
        NDIT, output = ETC.run_etc_calculator(obs_obj.name, obs_time)
    except Exception as e:
        print(type(e))
        if type(e) ==  json.decoder.JSONDecodeError: # catches errors in the etc-form.json input file for the ETC calculator
            # Routine to fix the JSONDecodeError, tries to find the false input.
            ETC.etc_debugger(obs_obj.name, obs_time, temperature = float(obs_obj.star_Teff/u.K), brightness = obs_obj.star_jmag, airmass = airmass, moon_target_sep = moon_target_sep, moon_phase = moon_phase)
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
    Calculates solely the S/N ratio for a given dit and ndit for a certain observation target at a certain observation time.

    Parameters
    ----------
    obs_obj : class object 
        class object of a observation target.
        
    obs_time : astropy.time.core.Time
        Time in UTC of observation.
        
    ndit : int
        Number of frames taken during a full single exposure.
    dit : float
        Exposure time for each frame.

    Returns
    -------
    output : namespace object
        Object containing the full output data from the ETC calculator.

    """
    ETC = Etc_form_class.etc_form(inputtype = "ndit")
    moon_target_sep, moon_phase, airmass, _ = airmass_moon_sep_obj_altaz(obs_obj, obs_time) #add moon_target_sep
    ETC.update_etc_form(temperature = float(obs_obj.star_Teff/u.K), brightness = obs_obj.star_jmag, airmass = airmass, dit = dit, ndit = ndit, moon_target_sep = moon_target_sep, moon_phase=moon_phase)
 
    ETC.write_etc_format_file()
    try:
        NDIT, output = ETC.run_etc_calculator(obs_obj.name, obs_time)
    except Exception as e:
        print(type(e))
        if type(e) ==  json.decoder.JSONDecodeError:
            # Routine to fix the JSONDecodeError
            ETC.etc_debugger(obs_obj.name, obs_time, temperature = float(obs_obj.star_Teff/u.K), brightness = obs_obj.star_jmag, airmass = airmass, moon_target_sep = moon_target_sep, moon_phase = moon_phase)
            
    return output, ETC

##########################################################################################################

@help_fun_logger  
def calculate_SN_ratio(sn_data):
    """
    Calculates the median of the signal to noise S/N ratio data

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

@help_fun_logger  
def extract_out_data(outputs):
    """
    Function to extract the S/N ratio data from the output file generated by the ETC calculator

    Parameters
    ----------
    outputs : namespace object or list
        Object or list of objects containing the full output data from the ETC calculator.

    Returns
    -------
    SN_data : list
        Contains a list of all data from the output(s) of the ETC calculator.

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
    This function calculates the moon distance, moon phase angle, airmass factor and local coordinates to observe 
    the object obs_obj.

    Parameters
    ----------
    obs_obj : class object
        instance of class Eclipses, Can also be any other .
        
    obs_time : astropy.time.core.Time
        Time in UTC of observation.
        
    location : TYPE, optional
        location of the observatory. The default is paranal.location.

    Returns
    -------
    moon_target_sep : TYPE
        DESCRIPTION.
    moon_phase : TYPE
        DESCRIPTION.
    airmass : np.float64
        Airmass factor at observation time and location.
        
    obs_altazs : TYPE
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


    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    Objects : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    with open(filename, 'wb') as out:
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
    

    Parameters
    ----------
    eclipse : TYPE
        DESCRIPTION.
    planet : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

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
    

    Parameters
    ----------
    eclipse : TYPE
        DESCRIPTION.
    planet : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

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
        obs_time_begin = eclipse['Eclipse Begin']['time']
        obs_time_end = eclipse['Eclipse End']['time']
        Transit_dur = planet.transit_duration.to(u.second).value # in seconds
        # for different S/N ratio add argument 'snr'= , to every Etc_calculator_Texp function
        Exposure_time_mid, _, _, output, _ = Etc_calculator_Texp(planet, obs_time) #obtimising NDIT for each single exposure with S/N min = 100 in seconds
        Exposure_time_begin, _, _, output, _ = Etc_calculator_Texp(planet, obs_time_begin) #obtimising NDIT for each single exposure with S/N min = 100 in seconds
        Exposure_time_end, _, _, output, _ = Etc_calculator_Texp(planet, obs_time_end) #obtimising NDIT for each single exposure with S/N min = 100 in seconds
        Exposure_times = [Exposure_time_begin, Exposure_time_mid, Exposure_time_mid] # get max exposure time
        Exposure_times.sort()
        
        SN_data = extract_out_data(output)
        median_SN, min_SN, max_SN = calculate_SN_ratio(SN_data)
        num_exp_possible = int(np.floor(Transit_dur/Exposure_times[-1]))
        if num_exp_possible >= 20:
            eclipse['Number of exposures possible'] = num_exp_possible # estimates the number of exposures possible according to the transit duration and the maximum exposure time
            eclipse['S/N median'] = median_SN
            eclipse['Minimum S/N'] = min_SN
            eclipse['Maximum S/N'] = max_SN
            eclipse['Minimum Exposure Time'] = Exposure_times[0]
            eclipse['Maximum Exposure Time'] = Exposure_times[-1]
        else:
            eclipse['Number of exposures possible'] = 'Target does not reach 20 exposures'
            eclipse['Estimated number of exposures'] = num_exp_possible
    
##########################################################################################################    

@help_fun_logger  
def data_sorting_and_storing(Eclipses_List, filename):
    """
    Sorting and storing final data as csv files, For now this only works with Eclipses, will include later functionality 
    to sort and store more general targets and maybe different functions to plot different kinds of data. 
    Might contain more types of output than it has now.

    Parameters
    ----------
    Eclipses_List : TYPE
        DESCRIPTION.
    filename : TYPE
        DESCRIPTION.

    Returns
    -------
    ranking : list
        Ranking of the Transits according to (Number of exposures possible)**2 * (number of Transit in computed timespan).

    """
    frame1 = []
    
    general1 = []
    ranking = []
    
    for planet in Eclipses_List:
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
    
    for planet in Eclipses_List:
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
    
    return ranking

##########################################################################################################  

@help_fun_logger  
def plotting_transit_data(Max_Delta_days, ranking, Eclipses_List, Nights):
    """
    Plotting final data, For now this only works for Eclipses, will include later functionality to plot general 
    targets and maybe different functions to plot different kinds of data. Might contain more types of output 
    than it has now.

    Parameters
    ----------
    Max_Delta_days : int
        Date span to plot the data for.
    ranking : list
        Ranking of the Transits according to (Number of exposures possible)**2 * (number of Transit in computed timespan).
    Eclipses_List : TYPE
        DESCRIPTION.
    Nights : TYPE, optional
        DESCRIPTION. The default is Nights_paranal.

    Returns
    -------
    ranking : list
        ranking in reversed order.

    """  
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
                            x_planet_dates = mpl_dates.date2num(x_planet)
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
            for planet in Eclipses_List:
                if planet.name == elem[1]:
                    tran_dur = np.float16(planet.transit_duration.to(u.hour))
                    for ecl in planet.eclipse_observable:
                        x_planet = [ecl['Eclipse Begin']['time'].value, ecl['Eclipse End']['time'].value]
                        x_planet_dates = mpl_dates.date2num(x_planet)
                        ax.plot(x_planet, y_planet, color='blue')
            planet_names.append("{} : {:.3}".format(elem[1], tran_dur))
            j += 1
        
        d = Nights.date[(n+1)*90]
        d_end = Nights.date[-1]
        lims = [d,d_end]
        
        plt.xlim(lims)
        plt.xticks(Nights.date[n*90:], fontsize=22)
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
                        x_planet_dates = mpl_dates.date2num(x_planet)
                        ax.plot(x_planet, y_planet, color='blue')
            planet_names.append("{} : {:.3}".format(elem[1], tran_dur))
            j += 1
    
        d = Nights.date[0]
        d_end = Nights.date[-1]
        lims = [d,d_end]
        
        plt.xlim(lims)
        plt.xticks(Nights.date, fontsize=22)
        plt.xticks(rotation=70)
        plt.yticks(y_range, planet_names, fontsize=22)
        plt.xlabel('Date', fontsize=24)
        plt.ylabel('Planet name : \nTransit duration [h]', fontsize=24)
        
        # plt.tight_layout()
        plt.show()
        fig.savefig(f"{d}-{d_end}-results.eps")
    
    
    ranking.reverse()
    return ranking

##########################################################################################################