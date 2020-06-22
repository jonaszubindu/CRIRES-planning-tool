#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 10:45:38 2020

Test file to test the computation of the changing variables that the ETC calculator is fed with.

@author: jonaszbinden
"""
from astropy.coordinates import SkyCoord, AltAz, get_moon
from astropy.time import Time
from astropy import units as u
from astroplan import Observer, moon_phase_angle
import numpy as np
import datetime
import pprint

""" Location and UTC offset Paranal """
paranal = Observer.at_site('paranal', timezone='Chile/Continental')

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
        obs_coor = obs_obj.Coordinates.coord
    frame_obs = AltAz(obstime=obs_time, location=location)
    obs_altazs = obs_coor.transform_to(frame_obs)
    airmass = obs_altazs.secz
    moon = get_moon(obs_time).transform_to(frame_obs)
    moon_target_sep = moon.separation(obs_altazs) # calculates the moon target separation
    moon_target_sep = moon_target_sep.deg * u.deg
    moon_phase = moon_phase_angle(time=obs_time)
    moon_phase = moon_phase.to(u.deg)
    z = 90 * u.deg - obs_altazs.alt
    zmoon = 90 * u.deg - moon.alt
    sep_min = np.abs(z-zmoon)
    sep_max = np.abs(z+zmoon)
    if moon_target_sep > sep_max: 
        moon_target_sep = sep_max
    elif moon_target_sep < sep_min:
        moon_target_sep = sep_min
        
    moon_target_sep = np.float64(moon_target_sep.value)
    moon_phase = np.float64(moon_phase.value)
    airmass = np.float64(airmass.value)
    moon_alt = moon.alt.value
    moon_tar_sep = (moon_target_sep, moon_alt)
    
    return moon_tar_sep, moon_phase, airmass, obs_altazs
    


obs_obj = SkyCoord.from_name('KELT-10 b')
obs_time = Time(datetime.datetime.fromisoformat('2020-06-20 04:52:09.346232'))

moon_target_sep, moon_phase, airmass, _ = airmass_moon_sep_obj_altaz(obs_obj, obs_time) #add moon_target_sep
# pprint.pprint({'moon_alt': moon.alt, 'moon_az' : moon.az, 'moon_target_sep' : moon_target_sep, 'moon_phase' : moon_phase, 'airmass' : airmass, 'Az' : obs_altazs.az, 'Alt' : obs_altazs.alt, 'sep_min' : sep_min, 'sep_max' : sep_max})







