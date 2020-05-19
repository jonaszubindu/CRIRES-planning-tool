#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 11:53:09 2020

This file contains helper Functions used in Transit_List.py

@author: jonaszbinden
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
paranal = Observer.at_site('paranal', timezone='Etc/GMT-4')


# def error_logger(orig_fun):
#     """ Function to log all occurring errors """
    
#     @wraps(orig_fun)
#     def wrapper(*args, **kwargs):
        
#         return orig_fun(*args, **kwargs)
    
#     return wrapper

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
        if type(e) ==  json.decoder.JSONDecodeError:
            # Routine to fix the JSONDecodeError
            ETC.etc_debugger(temperature = float(obs_obj.star_Teff/u.K), brightness = obs_obj.star_jmag, airmass = airmass)
        elif type(e) == ConnectionError:
            #STORE OBS DATA#
            input('Connection Error: Check Internet connection and press enter when problem resolved:')
        else:
            raise e

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
            ETC.etc_debugger(temperature = float(obs_obj.star_Teff/u.K), brightness = obs_obj.star_jmag, airmass = airmass)
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
    