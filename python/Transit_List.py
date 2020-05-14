#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 13:55:10 2020

@author: jonaszbinden
"""


import requests
import os
import json


def connect(host='https://exoplanetarchive.ipac.caltech.edu/'): # Nasa Exoplanet Archive
    """Check Internet Connection to Nasa Exoplanet Archive"""
    req = requests.get(host)  # Python 3.x
    if req.ok: 
        print('Connected to {}'.format(host))
    else:
        raise Warning('Check connection to {}, response code:{}'.format(host,req.status_code))
    return req

response = connect()


import astroplan
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
import Etc_form_class

"""Update most recent IERS data"""
get_IERS_data = 'yes'

try:
    if get_IERS_data == 'yes':
        download_IERS_A(show_progress=True)
        print('IERS data successfully downloaded')
    else:
        try:
            get_IERS_A_or_workaround()  # For now, in future always download the most recent ones
            print('IERS data successfully retrieved')
        except:
            download_IERS_A(show_progress=True)
            print('IERS data successfully downloaded')
except:
    print('No input given, downloading IERS data...')
    download_IERS_A(show_progress=True)
    print('IERS data successfully downloaded')



class Exoplanets:
    """
    some docs here
    """
    
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

    
    def Planet_finder(self, name):
        """Checking if Planet can be found in Nasa Exoplanet Archive"""
        Planet_try = NasaExoplanetArchive.query_planet(name, all_columns=True)
        print('Planet ' + name + ' found in Nasa Exoplanet Archive\n')
        self.Exoplanets_List_Nasa.append(Planet_try)
        if not Planet_try:
            print('Planet not in Nasa Exoplanet Archive\n')
            Exoplanets.Exoplanet_not_Nasa.append(name)
            


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



try:
    """Name_list includes all the exoplanet's names downloaded via Request_Table_NasaExoplanetArchive."""
    Name_list = csv_file_import.main()
except:
    raise Warning('csv file is corrupted')

for name in Name_list: # Check where name list comes from to make more efficient
    """Check for any non string-type elements in the list of planet names to parse."""
    if type(name) is not str:
        Warning('Name is not a string, cannot parse {}'.format(name))
    

Exoplanets = Exoplanets()

for name in Name_list[:5]:
    try:
        """Load Planets from Nasa Exoplanet Archive"""
        Exoplanets.Planet_finder(name)
    except Exception:
        Exoplanets.Fail.append(name)
        print('Planet ' + name + ' added to list "Fail"\n')
        Exoplanets.Fail.append(Name_list.index(name))

"""Check all Planets if they have the necessary properties in the database to process for transit observations"""
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

"""Location and UTC offset Paranal"""
paranal_loc = EarthLocation(lat=-24.627 * u.deg, lon=-70.405 * u.deg, height=2635.43 * u.m)
utcoffset = -4 * u.hour
paranal = Observer.at_site('paranal', timezone='Etc/GMT-4')


dt = datetime.timedelta(days=1)
d = datetime.datetime(2020, 4, 1, 0, 0, 0)


# Exoplanets_NO_TRANSITS_DATA = []
# Exoplanets_NO_T0 = []  # List of Exoplanets without time T0 for transit reference

# for n in range(len(Exoplanets.Exoplanets_List)):
#     if type(Exoplanets.Exoplanets_List[n]['T0']) == np.ma.core.MaskedConstant:
#         Exoplanets_NO_T0.append(Exoplanets.Exoplanets_List[n]['NAME'])
#     else:
#         try:
#             Exoplanets.name.append(Exoplanets.Exoplanets_List[n]['NAME'])
#             Exoplanets.T0.append(Exoplanets.Exoplanets_List[n]['T0'])
#             Exoplanets.period.append(Exoplanets.Exoplanets_List[n]['PER'])
#         except:
#             Exoplanets_NO_TRANSITS_DATA.append(
#                 Exoplanets.Exoplanets_List[n]['NAME'])


# For those planets without a T0, check if one can find a T0 in a different DATABASE
# for name in Exoplanets_NO_T0:
#     try:
#         Exoplanets.Exoplanets_altern_T0.altern_Planet_finder(name)
#     except:
#         Warning('No other database entry found for ' + name + '.')

# From here on only Nasa Archive is used

"""Definition of maximum number of days to plan observations into the future"""
Per_max = [] 
for planet in Exoplanets.Parse_planets_Nasa:
    Per_max.append(np.float64(planet['pl_orbper'] / u.day)) # Maximum orbital period of exoplanets

# Max_Delta_days = int(max(Per_max) * 2)
Max_Delta_days = 30 # choose manually for how many days you want to compute the transits


delta_midnight = np.linspace(-12, 12, 1000) * u.hour
d_end = d + dt*Max_Delta_days



class Nights:
    """contains empty lists for dates, coordinates of the sun for the specific nighttimes at Paranal and the nighttimes"""
    def __init__(self):
        self.date = []
        self.coord = []
        self.night = []

# This could be generalized for other observatories as well
    def Calculate_nights_paranal(self, d, d_end, DataFrame = 0, WriteToCSV = 0):
        """
        Calculates the nights at paranal for a certain start date and end date. Retrieves the sun coordinates
        for each night from astroplan.
    
        Parameters
        ----------
        d : datetime
            Start date from which the nights at paranal are computed.
        d_end : datetime
            End date until which the nights at paranal are computed.
    
        Returns
        -------
        Nights_paranal_table : pandas DataFrame
            Contains the sun coordinates form Nights.coord.
        Nights_paranal_dates : pandas DataFrame
            Contains the dates from Nights.date.
    
        """
        midnights_2020_til_2024 = []
        Nights_paranal = []
        frame_2020_til_2024 = []
        sunaltazs_2020_til_2024 = []
        # if LoadFromCSV == 1:
        #     """Check if there exist csv files with night data for the preferred range"""
        #     try:    
        #         Nights_paranal_table = pd.read_csv('Nights_time_at_paranal.csv')
        #         Nights_paranal_dates = pd.read_csv('dates_timespan.csv')
        #     except Exception:
        #         print('No Night data found, computing nights...')
                
        print('Calculating the nights of paranal from the ' +
              str(d.date()) + ' until the ' + str(d_end.date()) + '.\n')
        for k in range(Max_Delta_days):
            midnights_2020_til_2024.append(d + dt * k)
            Nights_paranal.append(Time(str(midnights_2020_til_2024[k])) - utcoffset + delta_midnight)  # introduce
            # compute frame AltAz for get_sun
            frame_2020_til_2024.append(AltAz(obstime=Nights_paranal[k], location=paranal.location))
            Sun_per_day = get_sun(Nights_paranal[k]).transform_to(
                frame_2020_til_2024[k])
            # access sun coord. for a specific date and time via sunaltazs_2020_til_2024[0](day)[0](time)
            sunaltazs_2020_til_2024.append(Sun_per_day)
    
            Nights.date.append(sunaltazs_2020_til_2024[k][0].obstime.datetime.date())
            for n in range(len(delta_midnight)):
                Time_sun = (
                    sunaltazs_2020_til_2024[k][n].obstime.value.split(" ")[1])
                if sunaltazs_2020_til_2024[k][n].alt < -18 * u.deg:
                    Nights.coord.append({
                        'Date': str(Nights.date[k]),
                        'Time': Time_sun,
                        'Az': sunaltazs_2020_til_2024[k][n].az,
                        'Alt': sunaltazs_2020_til_2024[k][n].alt})

        """"Splitting the nights into single nights and create a list with single night times in UTC -> Nights.night"""
        for n in range(len(Nights.date)):
            night = []
            for k in range(len(Nights.coord)):
                T = Nights.coord[k]['Date'] + ' ' + Nights.coord[k]['Time']  # UTC
                # time_list.append(T)
                if Nights.coord[k]['Date'] == str(Nights.date[n]):
                    night.append(T)
            Nights.night.append(Time(night))
        
        if DataFrame == 1:
            """Write the nights to a pandas DataFrame"""
            Nights_paranal_dates = pd.DataFrame(Nights.dates)
            Nights_paranal_table = pd.DataFrame(Nights)  # All nights in Paranal in UTC
            self.Nights_paranal_table.rename(columns={0: 'Date'}, inplace=True)
        if WriteToCSV == 1:
            """Write Nights_paranal_table to csv file"""
            Nights_paranal_table.to_csv('Nights_time_at_paranal.csv')
            Nights_paranal_dates.to_csv('dates_timespan.csv')
        
            
            



Nights = Nights()
"""Generates the class object Nights and calculates the nights for paranal between d and d_end"""
Nights.Calculate_nights_paranal(d, d_end) 


"""Altitude constraints definition"""
Altcons = astroplan.AltitudeConstraint(min=+30 * u.deg, max=None)  

"""Airmass constraints definition"""
Airmasscons = astroplan.AirmassConstraint(min=None, max=1.7)  

"""Astronomical Nighttime constraints definition: begin and end of each night at paranal as LocalTimeConstraint, written to list Night_cons_per_night"""
Night_cons_per_night = []
for m in range(len(Nights.night)):  # 
    dn_min = Nights.night[m][0]
    dn_max = Nights.night[m][-1]
    Night_cons_per_night.append(astroplan.LocalTimeConstraint(
        min=dn_min.datetime.time(), max=dn_max.datetime.time()))


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
    def __init__(self, name, epoch, period, transit_duration, sky_coords, eccentricity, star_Teff, star_jmag):
        self.name = name
        self.epoch = Time(epoch, format='jd')
        self.period = period
        self.transit_duration = transit_duration
        self.eccentricity = eccentricity
        self.star_Teff = star_Teff
        self.star_jmag = star_jmag
        self.num_eclipses = []
        self.target_observable = []
        self.eclipse_observable = []
        self.eclipse_mid_observable = []
        self.Airmass_window = []
        self.Alt_window = []
        self.Time_window = []
        
        
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
        self.Planets_eclipse = Planet_next_eclipse
        
        self.num_eclipses = int(np.floor(Max_Delta_days /(self.period / u.day)))
        
        try:
            self.Coordinates = SkyCoord.from_name(self.name)
        except Exception:
            self.Coordinates = sky_coords
        

    def Observability(self):
        """
        Some docs here

        Returns
        -------
        None.

        """
        
        obs_time = Nights.night[0][0]
        Planet_next_eclipse_Times = self.Planets_eclipse.next_primary_eclipse_time(obs_time, n_eclipses=self.num_eclipses)
        print(self.name + ' is getting processed')
        for date, night, night_cons in zip(Nights.date,Nights.night,Night_cons_per_night):
        
            for planet_next_eclipse_by_date in Planet_next_eclipse_Times:
                """Loop over all eclipses coming up in the given timespan of object planet"""
                
                if date == planet_next_eclipse_by_date.datetime.date():  # Check which eclipse can be observed in which night
                    Planet_next_eclipse_per_night_MID = planet_next_eclipse_by_date
                    Planet_next_eclipse_per_night_BEGIN = Planet_next_eclipse_per_night_MID + self.transit_duration / 2
                    Planet_next_eclipse_per_night_END = Planet_next_eclipse_per_night_MID + self.transit_duration / 2
                    
                    Planet_Eclipes_NIGHT = [Planet_next_eclipse_per_night_BEGIN, Planet_next_eclipse_per_night_MID, Planet_next_eclipse_per_night_END]  # Begin, midpoint and end of transit
    
                    """Careful, this only computes if the complete transit is observable"""
                    ecl_obs = astroplan.is_event_observable(constraints=[Altcons, Airmasscons, night_cons], observer=paranal, target=self.Coordinates, times=Planet_Eclipes_NIGHT)
                    """This computes if mid transit is observable"""
                    ecl_mid_obs = astroplan.is_event_observable(constraints=[Altcons, Airmasscons, night_cons], observer=paranal, target=self.Coordinates, times=Planet_next_eclipse_per_night_MID)
                    
                    tar_obs = astroplan.is_event_observable(constraints=[Altcons, Airmasscons], observer=paranal, target=self.Coordinates, times=night)
                    if all(ecl_obs[0] == True):
                        _, _, airmass, obs_altazs = [airmass_moon_sep_obj_altaz(self, tim) for tim in Planet_Eclipes_NIGHT]
                        print('{} total Eclipse is observable'.format(self.name))
                        self.eclipse_observable.append({
                            'Name': self.name,
                            'obs_time': Planet_next_eclipse_per_night_MID,
                            'Primary eclipse observable?': ecl_obs[0][0],
                            'Transit Length': self.transit_duration,
                            'Eclipse Begin': {'time' : Planet_next_eclipse_per_night_BEGIN, 
                                              'airmass' : airmass[0],
                                              'az' : obs_altazs[0].az,
                                              'alt' : obs_altazs[0].alt
                                              },
                            'Eclipse Mid' : {'time' : Planet_next_eclipse_per_night_MID,
                                             'airmass' : airmass[1],
                                             'az' : obs_altazs[1].az,
                                             'alt' : obs_altazs[1].alt
                                             },
                            'Eclipse End': {'time' : Planet_next_eclipse_per_night_END,
                                            'airmass' : airmass[2],
                                            'az' : obs_altazs[2].az,
                                            'alt' : obs_altazs[2].alt
                                            }})
                    elif ecl_mid_obs[0][0] == True:
                        _, _, airmass, obs_altazs = airmass_moon_sep_obj_altaz(self, Planet_next_eclipse_per_night_MID)
                        print('{} mid Eclipse is observable'.format(self.name))
                        self.eclipse_mid_observable.append({
                            'Name' : self.name,
                            'Primary eclipse observable?' : ecl_mid_obs[0][0],
                            'Eclipse' : {'time' : Planet_next_eclipse_per_night_MID,
                                             'airmass' : airmass,
                                             'az' : obs_altazs.az,
                                             'alt' : obs_altazs.alt
                                             }})
                    elif any(tar_obs[0] == True) and all(ecl_obs[0] == False):
                        print('{} Target is observable without any primary eclipse'.format(self.name))
                        for n, tar in enumerate(tar_obs[0]):
                            if tar == True:
                                _, _, airmass, obs_altazs = airmass_moon_sep_obj_altaz(self, night[n])
                                self.target_observable.append({
                                    'Name': self.name,
                                    'Object w. o. primary eclipse observable?': tar,
                                    'Obs Data': {'time' : night[n],
                                                     'airmass' : airmass,
                                                     'az' : obs_altazs.az,
                                                     'alt' : obs_altazs.alt
                                                     }})

            Alt_constraints = Altcons.compute_constraint(times=night, observer=paranal, targets=self.Coordinates)
            
            Airmass_constraints = Airmasscons.compute_constraint(times=night, observer=paranal, targets=self.Coordinates)
            
            Night_constraint = night_cons.compute_constraint(times=night, observer=paranal, targets=self.Coordinates)
        
            self.Airmass_window.append(Alt_constraints)
            self.Alt_window.append(Airmass_constraints)
            self.Time_window.append(Night_constraint)

def Etc_calculator_Texp(obs_obj, obs_time, snr=100):
    
    NDIT_opt = 24 # NDIT should optimally be between 16-32
    ETC = Etc_form_class.etc_form(inputtype = "snr")
    if snr != 100:
        ETC.update_etc_form(snr=snr)
    
    _, _, airmass, _ = airmass_moon_sep_obj_altaz(obs_obj, obs_time) #add moon_target_sep
    ETC.update_etc_form(temperature = float(obs_obj.star_Teff/u.K), brightness = obs_obj.star_jmag, airmass = airmass) #, moon_target_sep = moon_target_sep)
 
    ETC.write_etc_format_file()
    
    # Routine to change ndit to 16-32 and change dit accordingly:
    NDIT, output = ETC.run_etc_calculator()
    cycles = 0
    while NDIT < 16 or NDIT > 32:
        
        Exposure_time = NDIT*ETC.input.timesnr.dit
        DIT_new = Exposure_time/NDIT_opt # determine DIT for NDIT=24
        ETC.input.timesnr.dit = DIT_new 
        ETC.write_etc_format_file() # write new DIT into 'etc-form.json'
        print(cycles, DIT_new, NDIT)
        NDIT, output = ETC.run_etc_calculator() # recalculate the new NDIT
        if NDIT == 0:
            raise Exception('NDIT not available from etc-calculator')
        if cycles > 5:
            print('too many tries to bring NDIT between 16-32')
            break
        cycles += 1
    print(NDIT, cycles)
    DIT = ETC.input.timesnr.dit
    print(DIT,NDIT)
    Exposure_time = NDIT*DIT # seconds
    return Exposure_time, DIT, NDIT, output

def Etc_calculator_SN(obs_obj, obs_time, ndit, dit):
    
    ETC = Etc_form_class.etc_form(inputtype = "ndit")
    _, _, airmass, _ = airmass_moon_sep_obj_altaz(obs_obj, obs_time) #add moon_target_sep
    ETC.update_etc_form(temperature = float(obs_obj.star_Teff/u.K), brightness = obs_obj.star_jmag, airmass = airmass, dit = dit, ndit = ndit) #, moon_target_sep = moon_target_sep)
 
    ETC.write_etc_format_file()
    _, output = ETC.run_etc_calculator()
    
    return output


def calculate_SN_ratio(output):
    #Find the median.
    SN_data = []
    for data in output.data.orders:
        SN_data.extend(data['detectors'][0]['data']['snr']['snr']['data'])
    SN_sorted = np.sort(SN_data)
    length = len(SN_sorted)
    if (length % 2 == 0):
        median_SN = (SN_sorted[(length)//2] + SN_sorted[(length)//2-1]) / 2
    else:
        median_SN = SN_sorted[(length-1)//2]
    min_SN = SN_sorted[0]
    max_SN = SN_sorted[1]

    return median_SN, min_SN, max_SN


def airmass_moon_sep_obj_altaz(obs_obj, obs_time, location=paranal.location):
    
    obs_coor = obs_obj.Coordinates
    frame_obs = AltAz(obstime=obs_time, location=location)
    obs_altazs = obs_coor.transform_to(frame_obs)
    # moon_altazs = get_moon(obs_time).transform_to(frame_obs)
    moon_phase = astroplan.moon_phase_angle(time=obs_time)
    # moon_target_sep = obs_altazs - moon_altazs
    moon_target_sep = 0
    airmass = np.float64(obs_altazs.secz)
    
    return moon_target_sep, moon_phase, airmass, obs_altazs



"""
Some docs here
"""
obs_time = Nights.night[0][0]
table_eclipse_observable = []
table_object_observable = []
table_eclipse_mid_observable = []
Eclipses_List = []
# for planet in Exoplanets.Parse_planets_Nasa:
planet = Exoplanets.Parse_planets_Nasa[0]
Planet = Eclipses(planet['pl_name'][0], planet['pl_tranmid'][0], planet['pl_orbper'][0], planet['pl_trandur'][0] * u.day, planet['sky_coord'][0], planet['pl_orbeccen'][0], planet['st_teff'][0], planet['st_j'][0])
Planet.Observability()
Eclipses_List.append(Planet)
"""
Calculates the median of the signal to noise ratio achievable in transits that allow around or more than 20 single exposures
during the transit. The Number of exposures possible is added to the list eclipse_observable and eclipse_mid_observable 
for comparison. Each exposure is optimised to have NDIT between 16 and 32 with a minimum S/N = 100. The resulting S/N ratios 
are used to compute the median. More values like DIT, NDIT, SN of each exposure for each transit could be stored as well, 
not implemented yet.
"""
obs_obj = Planet
for eclipse1 in Planet.eclipse_observable:
    if eclipse1['Primary eclipse observable?'] == True:
        print('{} gets fed to ETC calculator for best observations'.format(Planet.name))
        obs_time = eclipse1['Eclipse Mid']
        Exposure_time, DIT, NDIT, output = Etc_calculator_Texp(obs_obj, obs_time) #obtimising NDIT for each single exposure with S/N min = 100 in seconds
        Transit_dur = (np.float128(Planet.transit_duration/u.day))*24*3600 # in seconds
        Number_of_exposures_possible = Transit_dur/Exposure_time
        eclipse1['Total Exposure Time'] = Exposure_time
        eclipse1['DIT'] = DIT
        eclipse1['NDIT'] = NDIT
        eclipse1['Number of Exposures possible'] = int(np.floor(Number_of_exposures_possible))
        if Number_of_exposures_possible >= 10:
            eclipse1['COMMENT'] = "More than 20 exposures possible per transit"
            time_between_exposures = 10 # buffer between two exposures in seconds
            dt = datetime.timedelta(seconds = (Exposure_time + time_between_exposures))
            Exposure_time_single_transit = []
            DIT_single_transit = []
            NDIT_single_transit = []
            num_exp = int(np.floor(Number_of_exposures_possible/2))
            obs_times_up = [obs_time + dt*n for n in range(num_exp)]
            obs_times_down = [obs_time - dt*n for n in range(num_exp)]
            obs_times_down.reverse()
            obs_times = obs_times_down + obs_times_up
            for obs_time in obs_times:
                """
                These values get stored only temporarily, if these data should be stored for later, 
                for example to have data about how each exposure should be taken, then
                write them into some new class or table
                """
                Exposure_time, DIT, NDIT, output = Etc_calculator_Texp(obs_obj, obs_time) #obtimising NDIT for each single exposure with S/N min = 100 in seconds
                Exposure_time_single_transit.append(Exposure_time)
                DIT_single_transit.append(DIT)
                NDIT_single_transit.append(NDIT)
            median_SN, _, _ = calculate_SN_ratio(output)
            eclipse1['Median Signal to Noise ratio'] = median_SN
for eclipse2 in Planet.eclipse_mid_observable:
    if eclipse2['Primary eclipse observable?'] == True:
        print('{} gets fed to ETC calculator for best observations'.format(Planet.name))
        obs_time = eclipse2['Eclipse']['time']
        Exposure_time, DIT, NDIT, output = Etc_calculator_Texp(obs_obj, obs_time) #obtimising NDIT for each single exposure with S/N min = 100 in seconds
        Transit_dur = (np.float128(Planet.transit_duration/u.day))*24*3600 # in seconds
        Number_of_exposures_possible = Transit_dur/Exposure_time
        eclipse2['Total Exposure Time'] = Exposure_time
        eclipse2['DIT'] = DIT
        eclipse2['NDIT'] = NDIT
        # eclipse2['Number of exposures possible'] = Number_of_exposures_possible
        if Number_of_exposures_possible >= 10:
            eclipse2['COMMENT'] = "Not full transit observable, compare with Number of exposures"
            time_between_exposures = 10 # buffer between two exposures in seconds
            dt = datetime.timedelta(seconds = (Exposure_time + time_between_exposures))
            Exposure_time_single_transit = []
            DIT_single_transit = []
            NDIT_single_transit = []
            num_exp = int(np.floor(Number_of_exposures_possible/2))
            obs_times_up = [obs_time + dt*n for n in range(num_exp)]
            obs_times_down = [obs_time - dt*n for n in range(num_exp)]
            obs_times_down.reverse()
            obs_times = obs_times_down + obs_times_up
            for obs_time in obs_times:
                """
                These values get stored only temporarily, if these data should be stored for later,
                for example to have data about how each exposure should be taken, then 
                write them into some new class or table
                """
                Exposure_time, DIT, NDIT, output = Etc_calculator_Texp(obs_obj, obs_time) #obtimising NDIT for each single exposure with S/N min = 100 in seconds
                Exposure_time_single_transit.append(Exposure_time)
                DIT_single_transit.append(DIT)
                NDIT_single_transit.append(NDIT)
            median_SN, _, _ = calculate_SN_ratio(output)
            eclipse2['Median Signal to Noise ratio'] = median_SN
            
pd_ecl_obs = pd.DataFrame(data=Planet.eclipse_observable)
pd_targ_obs = pd.DataFrame(data=Planet.target_observable)
pd_ecl_mid_obs = pd.DataFrame(data=Planet.eclipse_mid_observable)
table_eclipse_mid_observable.append(pd_ecl_mid_obs)
table_eclipse_observable.append(pd_ecl_obs)
table_object_observable.append(pd_targ_obs)



# """
# For each observable planet:

# Do stuff with the input file 'etc-form.json' here:
# use: ETC.update_etc_form(**kwargs) from Etc_form_class

# Then write the whole file again as a json file: etc-form.json
# with ETC.write_etc_format_file()
# and run it with Etc_form_class.run_etc_calculator


# The following section is not yet fully working and first needs some debugging

# """
# for planet in Eclipses_List:
    # """
    # Calculates the median of the signal to noise ratio achievable in transits that allow around or more than 20 single exposures
    # during the transit. The Number of exposures possible is added to the list eclipse_observable and eclipse_mid_observable 
    # for comparison. Each exposure is optimised to have NDIT between 16 and 32 with a minimum S/N = 100. The resulting S/N ratios 
    # are used to compute the median. More values like DIT, NDIT, SN of each exposure for each transit could be stored as well, 
    # not implemented yet.
    # """
    # obs_obj = planet
    # for eclipse1, eclipse2 in zip(planet.eclipse_observable, planet.eclipse_mid_observable):
    #     if eclipse1['Primary eclipse observable?'] == True:
    #         obs_time = eclipse1['Eclipse Mid']
    #         Exposure_time, DIT, NDIT, SN, _ = Etc_calculator_Texp(obs_obj, obs_time)
    #         Number_of_exposures_possible = planet.transit_duration/Exposure_time
    #         eclipse1['Total Exposure Time'] = Exposure_time
    #         eclipse1['DIT'] = DIT
    #         eclipse1['NDIT'] = NDIT
    #         eclipse1['Number of Exposures possible'] = Number_of_exposures_possible
    #         if Number_of_exposures_possible >= 20:
    #             eclipse1['COMMENT'] = "More than 20 exposures possible per transit"
    #             time_between_exposures = 10 # buffer between two exposures in seconds
    #             dt = datetime.timedelta(seconds = (Exposure_time + time_between_exposures))
    #             Exposure_time_single_transit = []
    #             DIT_single_transit = []
    #             NDIT_single_transit = []
    #             num_exp = np.floor(Number_of_exposures_possible/2)
    #             obs_times_up = [obs_time + dt*n for n in num_exp]
    #             obs_times_down = [obs_time - dt*n for n in num_exp]
    #             obs_times = obs_times_down.reverse() + obs_times_up
    #             for obs_time in obs_times:
    #                 """
    #                 These values get stored only temporarily, if these data should be stored for later, 
    #                 for example to have data about how each exposure should be taken, then
    #                 write them into some new class or table
    #                 """
    #                 Exposure_time, DIT, NDIT, SN, _ = Etc_calculator_Texp(obs_obj, obs_time) #obtimising NDIT for each single exposure with S/N min = 100
    #                 Exposure_time_single_transit.append(Exposure_time)
    #                 DIT_single_transit.append(DIT)
    #                 NDIT_single_transit.append(NDIT)
    #             median_SN, _, _ = calculate_SN_ratio(output)
    #             eclipse1['Median Signal to Noise ratio'] = median_SN
            
    #     elif eclipse2['Primary eclipse observable?'] == True:
    #         obs_time = eclipse2['obs_time']
    #         Exposure_time, DIT, NDIT, SN, _ = Etc_calculator_Texp(obs_obj, obs_time) #obtimising NDIT for each single exposure with S/N min = 100
    #         Number_of_exposures_possible = planet.transit_duration/Exposure_time
    #         eclipse2['Total Exposure Time'] = Exposure_time
    #         eclipse2['DIT'] = DIT
    #         eclipse2['NDIT'] = NDIT
    #         # eclipse2['Number of exposures possible'] = Number_of_exposures_possible
    #         if Number_of_exposures_possible >= 20:
    #             eclipse2['COMMENT'] = "Not full transit observable, compare with Number of exposures"
    #             time_between_exposures = 10 # buffer between two exposures in seconds
    #             dt = datetime.timedelta(seconds = (Exposure_time + time_between_exposures))
    #             Exposure_time_single_transit = []
    #             DIT_single_transit = []
    #             NDIT_single_transit = []
    #             num_exp = np.floor(Number_of_exposures_possible/2)
    #             obs_times_up = [obs_time + dt*n for n in num_exp]
    #             obs_times_down = [obs_time - dt*n for n in num_exp]
    #             obs_times = obs_times_down.reverse() + obs_times_up
    #             for obs_time in obs_times:
    #                 """
    #                 These values get stored only temporarily, if these data should be stored for later,
    #                 for example to have data about how each exposure should be taken, then 
    #                 write them into some new class or table
    #                 """
    #                 Exposure_time, DIT, NDIT, SN, output = Etc_calculator_Texp(obs_obj, obs_time)
    #                 Exposure_time_single_transit.append(Exposure_time)
    #                 DIT_single_transit.append(DIT)
    #                 NDIT_single_transit.append(NDIT)
    #             median_SN, _, _ = calculate_SN_ratio(output)
    #             eclipse2['Median Signal to Noise ratio'] = median_SN

# table_eclipse_mid_observable = pd.concat(table_eclipse_mid_observable)
# table_eclipse_observable = pd.concat(table_eclipse_observable)
# table_object_observable = pd.concat(table_object_observable)

# text_file_name = input('Write Filename to save file: ')
# def_text_file_name = 'Observation_Timetable_Eclipse'
# def_text_file_name2 = 'Observation_Timetable_Objects.csv'

# if not text_file_name:
#     text_file_name = def_text_file_name

# text_file_name = text_file_name + '.csv'
# table_eclipse_observable.to_csv(text_file_name)
# direc = os.getcwd()
# print(text_file_name + ' is created in ' + direc)
# table_object_observable.to_csv(def_text_file_name2)


    

