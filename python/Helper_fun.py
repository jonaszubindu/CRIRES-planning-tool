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
from astroplan import moon_phase_angle, Observer
import datetime
import time
import json
from json import JSONDecodeError
import requests
import logging

""" Location and UTC offset Paranal """
paranal_loc = EarthLocation(lat=-24.627 * u.deg, lon=-70.405 * u.deg, height=2635.43 * u.m)
utcoffset = -4 * u.hour
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
    
    _, _, airmass, _ = airmass_moon_sep_obj_altaz(obs_obj, obs_time) #add moon_target_sep
    ETC.update_etc_form(temperature = float(obs_obj.star_Teff/u.K), brightness = obs_obj.star_jmag, airmass = airmass) #, moon_target_sep = moon_target_sep)
 
    ETC.write_etc_format_file()
    try:
        NDIT, output = ETC.run_etc_calculator(obs_obj.name, obs_time)
    except Exception as e:
        print(type(e))
        if type(e) ==  json.decoder.JSONDecodeError: # catches errors in the etc-form.json input file for the ETC calculator
            # Routine to fix the JSONDecodeError, tries to find the false input.
            ETC.etc_debugger(obs_obj.name, obs_time, temperature = float(obs_obj.star_Teff/u.K), brightness = obs_obj.star_jmag, airmass = airmass)

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
    _, _, airmass, _ = airmass_moon_sep_obj_altaz(obs_obj, obs_time) #add moon_target_sep
    ETC.update_etc_form(temperature = float(obs_obj.star_Teff/u.K), brightness = obs_obj.star_jmag, airmass = airmass, dit = dit, ndit = ndit) #, moon_target_sep = moon_target_sep)
 
    ETC.write_etc_format_file()
    try:
        NDIT, output = ETC.run_etc_calculator(obs_obj.name, obs_time)
    except Exception as e:
        print(type(e))
        if type(e) ==  json.decoder.JSONDecodeError:
            # Routine to fix the JSONDecodeError
            ETC.etc_debugger(obs_obj.name, obs_time, temperature = float(obs_obj.star_Teff/u.K), brightness = obs_obj.star_jmag, airmass = airmass)
        elif type(e) == ConnectionError:
            #STORE OBS DATA#
            input('Connection Error: Check Internet connection and press enter when problem resolved:')
        else:
            raise e
            
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
    if type(obs_obj) == SkyCoord:
        obs_coor = obs_obj
    else:
        obs_coor = obs_obj.Coordinates
    frame_obs = AltAz(obstime=obs_time, location=location)
    obs_altazs = obs_coor.transform_to(frame_obs)
    # moon_altazs = get_moon(obs_time).transform_to(frame_obs)
    moon_phase = moon_phase_angle(time=obs_time)
    # moon_target_sep = obs_altazs - moon_altazs
    moon_target_sep = 0
    airmass = np.float64(obs_altazs.secz)
    
    return moon_target_sep, moon_phase, airmass, obs_altazs


    
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
           
            
@help_fun_logger  
def SN_Transit_Observation_Optimization(eclipse, planet):
    """ Checking if eclipse has already been processed """
    try:
        test_empty = eclipse['Total Exposure Time']
        if test_empty != None:
            print('{} has already been processed, skipping...'.format(planet.name))
    except KeyError:
        print('{} gets fed to ETC calculator for best observations'.format(planet.name))
        logging.info('{} gets fed to ETC calculator for best observations'.format(planet.name))
        obs_time = eclipse['Eclipse Mid']['time']
        Transit_dur = (np.float128(planet.transit_duration/u.day))*24*3600 # in seconds
        
        """ Initial Calculation of Exposure time for Transit: eclipse """
        
        Exposure_time, _, _, output, _ = Etc_calculator_Texp(planet, obs_time) #obtimising NDIT for each single exposure with S/N min = 100 in seconds
        
        """ Find how many single exposures are possible to take """
        
        # Observations = []
        SN_data_overall = []
        num_exp = 1
        range_obs_times = Exposure_time
        time_between_exposures = 10 # buffer between two exposures in seconds
        n = 1
        median_SN_single_exp = []
        
        delta = datetime.timedelta(seconds = int(np.ceil(Exposure_time/2 + time_between_exposures)))
        obs_time_up = obs_time + delta
        obs_time_down = obs_time - delta
        
        while num_exp <= 20 and Transit_dur > range_obs_times:
            print('number of exposures: {}'.format(num_exp))
            Exposure_time_up, _, _, output, _ = Etc_calculator_Texp(planet, obs_time_up) # obtimising NDIT for each single exposure with S/N min = 100 in seconds
            Exposure_time_down, _, _, output, _ = Etc_calculator_Texp(planet, obs_time_down) # obtimising NDIT for each single exposure with S/N min = 100 in seconds
            
            delta_up = datetime.timedelta(seconds = int(np.ceil(Exposure_time_up + time_between_exposures)))
            delta_down = datetime.timedelta(seconds = int(np.ceil(Exposure_time_down + time_between_exposures)))
            obs_time_up = obs_time + delta_up
            obs_time_down = obs_time - delta_down
            # obs_times = [obs_time_down, obs_time_up]
            num_exp = 2*n
            range_obs_times = (obs_time_up-obs_time_down).sec
            n += 1
            
            # Observations.extend(obs_times)
            SN_data = extract_out_data(output)
            SN_data_overall.extend(SN_data)
            median_SN, _, _ = calculate_SN_ratio(SN_data)
            median_SN_single_exp.append(median_SN) 
        # Observations.sort()
                            
        if range_obs_times > Transit_dur:
            print('Time to reach 20 exposure exceeded, number of possible exposures: {}'.format(num_exp))
            eclipse['Number of exposures possible'] = num_exp
            eclipse['Comment'] = 'Reaching 20 exposures with S/N = 100 exceeds Transit Duration'
        elif num_exp >= 20: #and range_obs_times <= Transit_dur:
            print('Reached 20 Exposures in: {} during a {} Transit'.format(range_obs_times, Transit_dur))
            Median_SN, _, _ = calculate_SN_ratio(SN_data_overall)
            median_SN_single_exp.sort()
            eclipse['Number of exposures possible'] = num_exp
            eclipse['Time necessary to reach 20 exposures'] = range_obs_times
            eclipse['S/N overall median'] = Median_SN
            eclipse['minimum S/N'] = median_SN_single_exp[0]
            eclipse['maximum S/N'] = median_SN_single_exp[-1]
            
    