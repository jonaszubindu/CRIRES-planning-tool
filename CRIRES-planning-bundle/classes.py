#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 13:49:49 2020

Classes for Transit_List observability

@author: jonaszbinden
"""

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
from astroplan import Observer, FixedTarget
# from threading import Thread

import csv_file_import
from astroplan import download_IERS_A, get_IERS_A_or_workaround
from Helper_fun import help_fun_logger



import pickle
import Helper_fun as fun
from Helper_fun import help_fun_logger
import logging


""" Location and UTC offset Paranal """
paranal = Observer.at_site('paranal', timezone='Chile/Continental')

##########################################################################################################

class Exoplanets:
    """
    some docs here
    """
    @help_fun_logger
    def __init__(self):
        """
        Initialize lists to sort the planets retrieved from Nasa Exoplanet Archive according
        to the available data to compute Transits.

        """
        self.Exoplanets_planet_objects = []
        self.Exoplanets_List_Nasa = []
        self.Exoplanet_not_Nasa = []
        
        self.Transit_data_missing = []
        self.Transit_data_avail = []
        self.Parse_planets_Nasa = []

    @help_fun_logger
    def Planet_finder(self, name):
        """ Checking if Planet can be found in Nasa Exoplanet Archive """
        Planet_try = NasaExoplanetArchive.query_planet(name, all_columns=True)
        print('Planet ' + name + ' found in Nasa Exoplanet Archive\n')
        self.Exoplanets_List_Nasa.append(Planet_try)
        if not Planet_try:
            print('Planet not in Nasa Exoplanet Archive\n')
            Exoplanets.Exoplanet_not_Nasa.append(name)
            
    ##########################################################################################################
    
    @help_fun_logger
    def hasproperties(self):
        """
        some docs here

        Returns
        -------
        None.

        """
        for planet in self.Exoplanets_List_Nasa:
            flag = None
            print('Checking Planet ' + planet['pl_name'][0] + ' for Transit Data')
            if np.ma.is_masked(planet['pl_tranmid'][0]) is True:
                print('Planet ' + planet['pl_name'][0] +
                      ' has no data about transit mid\n')
                flag = True
            if np.ma.is_masked(planet['pl_orbper'][0]) is True:
                print('Planet ' + planet['pl_name'][0] +
                      ' has no data about orbit period\n')
                flag = True
            if np.ma.is_masked(planet['pl_trandur'][0]) is True:
                print('Planet ' + planet['pl_name'][0] +
                      ' has no data about transit duration\n')
                flag = True
            try:
                sky_coords = SkyCoord.from_name(planet['pl_name'][0])
            except Exception:
                try:
                    sky_coords = planet['sky_coord']
                    print('no Sky coordinates found for ' + planet['pl_name'][0])
                except:
                    flag = True
            if not sky_coords:
                flag = True
            if not flag:
                self.Transit_data_avail.append(planet['pl_name'][0])
                self.Parse_planets_Nasa.append(planet)
                print('Planet ' + planet['pl_name'][0] +
                      ' added to Transit_data_avail\n')
            else:
                self.Transit_data_missing.append(planet['pl_name'][0])
                print('Planet ' + planet['pl_name'][0] + ' added to Transit_data_missing\n')

##########################################################################################################

class Nights(object):
    """
        Calculates the nights at paranal for a certain start date and end date. Retrieves the sun coordinates
        for each night from astroplan.
    
        Parameters
        ----------
        d : datetime
            Start date from which the nights at paranal are computed.
            
        d_end : datetime
            End date until which the nights at paranal are computed.
    """
    @help_fun_logger
    def __init__(self, d, Max_Delta_days, LoadFromPickle = 0):
        dt = datetime.timedelta(days=1)
        d_end = d + dt*Max_Delta_days
        self.Max_Delta_days = Max_Delta_days
        if LoadFromPickle == 1:
            """ Check if there exist pkl files with night data for the preferred range, if not: initializes empty lists for Nights instance """
            try:
                d_str = d.isoformat() # start date from which the nights in paranal are calculated
                filename = 'Nights_paranal_{}-{}.pkl'.format(d_str, self.Max_Delta_days)
                nights = fun.pickled_items(filename)
                self.start = nights.__next__()
                self.end = nights.__next__()
                self.Max_Delta_days = nights.__next__()
                self.date = nights.__next__()
                self.night = nights.__next__()
                self.loaded = 1
            except Exception as e:
                print(e)
                self.loaded = 0
                d_str = d.isoformat() # start date from which the nights in paranal are calculated
                filename = 'Nights_paranal_{}-{}.pkl'.format(d_str, self.Max_Delta_days)
                print('No Night data found for {}, computing nights...'.format(filename))
        
        self.start = d
        self.end = d_end
        self.Max_Delta_days = (d_end - d).days
        self.date = []
        # self.night = []
        self.loaded = 0
            
        
        for k in range(self.Max_Delta_days):
            date = self.start + dt * k # list with datetime objects for the midnights in question
            self.date.append(date)

    ##########################################################################################################

    @help_fun_logger
    def Calculate_nights_paranal(self, delta_midnight, observatory = paranal, WriteToPickle = 0):
        """
        Calculates the nights at observatory, default=paranal for a certain start date and end date. Retrieves the sun coordinates
        for each night from astroplan.
    
        Parameters
        ----------
        
        WriteToPickle : int
            Object Nights gets written into a pickle file.
        
        Returns
        -------
        self.dates : list
            Contains datetime.date objects for each night between d and d_end.      
        
        self.coords : list
            Contains dict with sun coordinates for each time in Nights.night.
            
        self.night : list
            Contains lists with nighttimes for each night. The number timesteps for each nights is defined in 
            delta_midnight.
    
        """
        
        if self.loaded == 1:
            """ If the nights at paranal with startdate d could be loaded from file, yield: """
            print('Nights loaded from file, continueing with processing planets for observability')
        else:
            
            """ All times are in UTC respectively """
            print('Calculating the nights of paranal from the {} until the {}'.format(self.start,self.end))
            self.night = []
            midnight = datetime.time(0,0,0)
            for date in self.date:
                
                midnight_datetime = datetime.datetime.combine(date,midnight) # list with datetime objects for the midnights in question
                midnight_at_site_UTC = observatory.datetime_to_astropy_time(midnight_datetime) # Time object for each midnight gets created in UTC, midnight in UTC.
                
                Night_paranal = midnight_at_site_UTC + delta_midnight # in UTC  
                # compute frame AltAz for get_sun
                frame_24_h_paranal = AltAz(obstime=Night_paranal, location=observatory.location)
                sunaltazs_24_h = get_sun(Night_paranal).transform_to(frame_24_h_paranal)
        
                night = []
                for n, _ in enumerate(delta_midnight):
                    """ Calculates the night times for each night """
                    if sunaltazs_24_h[n].alt < -18 * u.deg:
                        night.append(str(sunaltazs_24_h[n].obstime.value))
                self.night.append(Time(night))
                
            
            if WriteToPickle == 1:
                """Write Nights_paranal_table to file"""
                d = self.start
                d = d.isoformat() # start date from which the nights in paranal are calculated
                filename = 'Nights_paranal_{}-{}.pkl'.format(d, self.Max_Delta_days)
                fun.pickle_dumper_objects(filename, self)

##########################################################################################################                    

class Eclipses:
    """
    Initialization of Eclipse class. For a planet the necessary data for eclipse observation get initialized here. 
    
    Parameters:
    -------------------
    name : string
        'pl_name': name of the planet
        
    epoch : astropy.time.core.Time
        'pl_tranmid': the mid time of the next transit
    
    period : astropy.units.quantity.Quantity
        'pl_orbper': orbital period of the planet around its host star in u.day
    
    transit_duration : astropy.units.quantity.Quantity
        'pl_trandur': duration of a transit in u.day
        
    Coordinates : astropy.coordinates.sky_coordinate.SkyCoord (ICRS)
        'sky_coord': right ascension and azimuth of host star in degrees
    
    eccentricity : float
        'pl_eccen': eccentricity of the orbit of the planet around its host star
    
    star_Teff : astropy.units.quantity.Quantity
        'st_Teff': Effective temperature of the host star in u.K (Kelvin)
    
    star_jmag : float
        'st_j': Magnitude of the host star in the J-band
    
    Max_Delta_days : int
        Days for which the eclipses get computed
    
    -------------------
    
    Other parameters get initialized as empty lists and get assigned later.
    More parameters can be added manually. The parameters all come from the 
    'NasaExoplanetArchive.query_planet(name, all_columns=True)' function from
    astroquery. For the filtering of the Exoplanet candidates refer to 
    'Request_Table_NasaExoplanetArchive.py' and the file used to filter the Archive:
    'Nasa_Archive_Selection.txt', in this file you may find information to look up keywords
    that can be used to add additional information to a instance of Eclipses.
    
    -------------------
    
    Furthermore initializing an instance of this class calls astroplan.EclipsingSystem and creates an instance of
    EclipsingSystem, using the already initialized data. The instance is stored under self.Planets_eclipse.
    Additionally the number of eclipses possible in the evaluated timespan is computed and stored in self.num_eclipses.
    
    """
    @help_fun_logger
    def __init__(self, Max_Delta_days, planet = None):
        if planet == None:
            """ Create Empty class to load data from file """
            self.name = None
            self.epoch = None
            self.period = None
            self.transit_duration = None
            self.eccentricity = None
            self.star_Teff = None
            self.star_jmag = None
            self.Planets_eclipse = None
            self.num_eclipses = None
            
        else:
            """ Initialize Eclipse instance from Nasa query_planet object """
            self.name = planet['pl_name'][0]
            self.epoch = Time(planet['pl_tranmid'][0], format='jd')
            self.period = planet['pl_orbper'][0]
            self.transit_duration = planet['pl_trandur'][0] * u.day
            self.eccentricity = planet['pl_orbeccen'][0]
            self.star_Teff = planet['st_teff'][0]
            self.star_jmag = planet['st_j'][0]
            
        self.target_observable = [] 
        self.eclipse_observable = [] 
        
        if planet != None:
            
            Planet_next_eclipse = astroplan.EclipsingSystem(primary_eclipse_time=self.epoch, orbital_period=self.period, duration=self.transit_duration)
            """
            WARNING:
                There are currently two major caveats in the implementation of
                ''EclipsingSystem''. The secondary eclipse time approximation is
                only accurate when the orbital eccentricity is small, and the eclipse
                times are computed without any barycentric corrections. The current
                implementation should only be used forapproximate mid-eclipse times for
                low eccentricity orbits, with event durations longer than the
                barycentric correction error (<=16 minutes).
                
                Shortest Transit duration found in the candidates so far is 20 minutes. 6. May 2020.
                
                From EclipsingSystem.__doc__
                
            """
            self.Planets_eclipse = Planet_next_eclipse # not written to output file
            
            self.num_eclipses = int(np.floor(Max_Delta_days /(self.period / u.day))) 
            
            try:
                self.Coordinates = FixedTarget.from_name(self.name)
            except Exception:
                self.Coordinates = FixedTarget(planet['sky_coord'][0], name=self.name)
        
    ##########################################################################################################  
        
    @help_fun_logger
    def Observability(self, obs_time, Nights, constraints, check_eclipse, check_target=0, delta_midnight=None):
        """
        Some docs here

        Returns
        -------
        None.

        """
        Planet_next_eclipse_Times = self.Planets_eclipse.next_primary_eclipse_time(obs_time, n_eclipses=self.num_eclipses)
        print(self.name + ' is getting processed')
        for date in Nights.date:
            
            if check_target == 1:
                """ Check if Nights object has attribute nights to calculate observability of target """
                if hasattr(Nights,'night'):
                    k = list.index(Nights.date,date)
                    night = Nights.night[k]
                else:
                    Nights.Calculate_nights_paranal(delta_midnight)
                    
            for planet_next_eclipse_by_date in Planet_next_eclipse_Times:
                """ Loop over all eclipses coming up in the given timespan of object planet """
                
                if date == planet_next_eclipse_by_date.datetime.date():  # Check which eclipse can be observed in which night
                    
                    if check_eclipse == 1:
                        
                        Planet_next_eclipse_per_night_MID = planet_next_eclipse_by_date
                        Planet_next_eclipse_per_night_BEGIN = Planet_next_eclipse_per_night_MID - self.transit_duration / 2
                        Planet_next_eclipse_per_night_END = Planet_next_eclipse_per_night_MID + self.transit_duration / 2
                        
                        Planet_Eclipes_NIGHT = [Planet_next_eclipse_per_night_BEGIN, Planet_next_eclipse_per_night_MID, Planet_next_eclipse_per_night_END]  # Begin, midpoint and end of transit
        
                        """ Computes observability of the Transit """
                        ecl_obs = astroplan.is_event_observable(constraints=constraints, observer=paranal, target=self.Coordinates, times=Planet_Eclipes_NIGHT)
                        
                        if all(ecl_obs[0] == True):
                            print('{} total Eclipse is observable'.format(self.name))
                            airmass_moon_sep_obj_altaz_RESULT = [fun.airmass_moon_sep_obj_altaz(self, tim) for tim in Planet_Eclipes_NIGHT]
                            moon_target_sep = [out[0] for out in airmass_moon_sep_obj_altaz_RESULT]
                            moon_phase = [out[1] for out in airmass_moon_sep_obj_altaz_RESULT]
                            airmass = [out[2] for out in airmass_moon_sep_obj_altaz_RESULT]
                            obs_altazs = [out[3] for out in airmass_moon_sep_obj_altaz_RESULT]
                            # print(moon_target_sep, moon_phase, airmass, obs_altazs)
                            self.eclipse_observable.append({
                                'Name': self.name,
                                'obs_time': Planet_next_eclipse_per_night_MID,
                                'Primary eclipse observable?': ecl_obs[0][0],
                                'Transit Length': self.transit_duration.to(u.hour),
                                'Eclipse Begin': {'time' : Planet_next_eclipse_per_night_BEGIN, 
                                                  'airmass' : airmass[0],
                                                  'moon sep' : moon_target_sep[0][0],
                                                  'moon phase' : moon_phase[0],
                                                  'az' : obs_altazs[0].az,
                                                  'alt' : obs_altazs[0].alt
                                                  },
                                'Eclipse Mid' : {'time' : Planet_next_eclipse_per_night_MID,
                                                 'airmass' : airmass[1],
                                                 'moon sep' : moon_target_sep[1][0],
                                                 'moon phase' : moon_phase[1],
                                                 'az' : obs_altazs[1].az,
                                                 'alt' : obs_altazs[1].alt
                                                 },
                                'Eclipse End': {'time' : Planet_next_eclipse_per_night_END,
                                                'airmass' : airmass[2],
                                                'moon sep' : moon_target_sep[2][0],
                                                'moon phase' : moon_phase[2],
                                                'az' : obs_altazs[2].az,
                                                'alt' : obs_altazs[2].alt
                                                }})
                            
                    if check_target == 1:
                        """ Check if target observable independent of Transit, can not be turned on through menu yet. """
                        tar_obs = astroplan.is_event_observable(constraints=constraints, observer=paranal, target=self.Coordinates, times=night)   
                        if any(tar_obs[0] == True):
                            print('{} Target is observable without any primary eclipse'.format(self.name))
                            for n, tar in enumerate(tar_obs[0]):
                                if tar == True:
                                    moon_target_sep, moon_phase, airmass, obs_altazs = fun.airmass_moon_sep_obj_altaz(self, night[n])
                                    self.target_observable.append({
                                        'Name': self.name,
                                        'Object w. o. primary eclipse observable?': tar,
                                        'Obs Data': {'time' : night[n],
                                                         'airmass' : airmass,
                                                         'moon sep' : moon_target_sep[0],
                                                         'moon phase' : moon_phase,
                                                         'az' : obs_altazs.az,
                                                         'alt' : obs_altazs.alt
                                                         }})

                        # Alt_constraints = Altcons.compute_constraint(times=night, observer=paranal, targets=self.Coordinates)
            
                        # Airmass_constraints = Airmasscons.compute_constraint(times=night, observer=paranal, targets=self.Coordinates)
            
                        # Night_constraint = night_cons.compute_constraint(times=night, observer=paranal, targets=self.Coordinates)
            
            # if check_target == 1:
            #     """ Write computed constraints to lists """
            #     self.Airmass_window.append(Alt_constraints)
            #     self.Alt_window.append(Airmass_constraints)
            #     self.Time_window.append(Night_constraint)
            

##########################################################################################################            

@help_fun_logger
def load_Eclipses_from_file(filename, Max_Delta_days):
    """
    

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.

    Returns
    -------
    Eclipses_List : TYPE
        DESCRIPTION.

    """
    Eclipses_List = []
    def att_identifier(att):
        if type(att) == int:
                Planet.num_eclipses = att
        elif type(att) == FixedTarget:
            Planet.Coordinates = att
        elif type(att) == astroplan.periodic.EclipsingSystem:
            Planet.Planets_eclipse = att
        else:
            return att
        
    planet = fun.pickled_items(filename)
    att = None
    while True:
        Planet = Eclipses(Max_Delta_days)
        try:
            if att != None:
                print(att)
                Planet.name = att
            else:
                Planet.name = next(planet)
                
            Planet.epoch = next(planet)
            Planet.period = next(planet)
            Planet.transit_duration = next(planet)
            Planet.eccentricity = next(planet)
            Planet.star_Teff = next(planet)
            Planet.star_jmag = next(planet)
            att = None
            while att == None:
                att = next(planet)
                att = att_identifier(att)
            
            target_observable = att
            
            eclipse_observable = next(planet)
            
            att = None
            while att == None:
                att = next(planet)
                att = att_identifier(att)
            
            Planet.eclipse_observable.extend(eclipse_observable)
            Planet.target_observable.extend(target_observable)
            Eclipses_List.append(Planet)
        except StopIteration:
            print('Eclipses_List has been loaded from {}'.format(filename))
            logging.info('Eclipses_List has bin loaded from {}'.format(filename))
            Planet.eclipse_observable.extend(eclipse_observable)
            Planet.target_observable.extend(target_observable)
            Eclipses_List.append(Planet)
            break
    
    return Eclipses_List


##########################################################################################################


class Targets:
    """
    Initialization of Target class. For a star or other object the necessary data for observations get initialized here. 
    Use NasaExoplanetArchive.query_star to initialize or do manually:
        
    target = NasaExoplanetArchive.query_star('name')
    name = target['st_name'][0]
    star_Teff = target['st_teff'][0]
    star_jmag = target['st_j'][0]
    
    Parameters:
    -------------------
    name : string
        'name': name of the object
        
    Coordinates : astropy.coordinates.sky_coordinate.FixedTarget (ICRS)
        'sky_coord': right ascension and azimuth of host star in degrees
    
    star_Teff : astropy.units.quantity.Quantity
        'st_Teff': Effective temperature of the host star in u.K (Kelvin)
    
    star_jmag : float
        'st_j': Magnitude of the host star in the J-band
    
    -------------------
    
    Other parameters get initialized as empty lists and get assigned later.
    More parameters can be added manually. The parameters all come from the 
    'NasaExoplanetArchive.query_star(name, all_columns=True)' function from
    astroquery.
    
    -------------------
    
    """
    @help_fun_logger
    def __init__(self, name, star_Teff, star_jmag, Coordinates=None):
            
        """ Initialize Eclipse instance from Nasa query_planet object """
        self.name = name
        self.star_Teff = star_Teff
        self.star_jmag = star_jmag
            
        self.target_observable = [] 
        
        try:
            self.Coordinates = FixedTarget.from_name(self.name)
        except Exception:
            self.Coordinates = FixedTarget(Coordinates, name=self.name)
        
    ##########################################################################################################
    
    def target_observable(self, Nights, constraints, delta_midnight=None):
        """
        Calculates for which times during the time span of Nights, the target is observable under the given constraints.
        LATER : Could include plotting of target observability.

        Parameters
        ----------
        Nights : class
            Nights at Paranal for which to compute if the target is observable.
        constraints : list
            list of Astroplan constraints to constrain the observability.
        delta_midnight : numpy.linspace, Obtional
            grid of timesteps within 24 hours for which the observation should be calculated.

        """
        if delta_midnight == None:
            delta_midnight = np.linspace(-12, 12, 1000) * u.hour # defines number of timesteps per 24 hours
            
        print(self.name + ' is getting processed')
        for date in Nights.date:
            """ Check if Nights object has attribute nights to calculate observability of target """
            if hasattr(Nights,'night'):
                k = list.index(Nights.date, date)
                night = Nights.night[k]
            else:
                Nights.Calculate_nights_paranal(delta_midnight)
                night = Nights.night[k]
        
        
            """ Check if target observable """
            tar_obs = astroplan.is_event_observable(constraints=constraints, observer=paranal, target=self.Coordinates, times=night)   
            if any(tar_obs[0] == True):
                print('{} Target is observable without any primary eclipse'.format(self.name))
                for n, tar in enumerate(tar_obs[0]):
                    if tar == True:
                        moon_target_sep, moon_phase, airmass, obs_altazs = fun.airmass_moon_sep_obj_altaz(self, night[n])
                        self.target_observable.append({
                            'Name': self.name,
                            'Object observable?': tar,
                            'Obs Data': {'time' : night[n],
                                             'airmass' : airmass,
                                             'moon sep' : moon_target_sep[0],
                                             'moon phase' : moon_phase,
                                             'az' : obs_altazs.az,
                                             'alt' : obs_altazs.alt
                                             }})

##########################################################################################################
                