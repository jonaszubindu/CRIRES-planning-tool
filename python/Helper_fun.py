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

""" Location and UTC offset Paranal """
paranal = Observer.at_site('paranal', timezone='Chile/Continental')


def help_fun_logger(orig_fun):
    """ Function to log execution of other functions """
    logging.info('succesfully ran function:{}'.format(orig_fun.__name__))
    
    return orig_fun
    
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
    min_SN = SN_sorted[0]
    max_SN = SN_sorted[1]

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
    
    
    