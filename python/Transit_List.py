#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 13:55:10 2020

@author: jonaszbinden
"""


import urllib.request
import os


def connect(host='http://exoplanets.org/table'):
    try:
        urllib.request.urlopen(host)  # Python 3.x
        print('internet connected')
    except:
        print('No internet!')
        return False


if connect() is False:
    raise Warning('No Internet Connection, abort!')

import astroplan
import astropy.units as u
from astropy import table
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
from astroquery.exoplanet_orbit_database import ExoplanetOrbitDatabase
import astroquery.open_exoplanet_catalogue as oec
import datetime
import time
from astroplan import Observer
from threading import Thread

import csv_file_import
from astroplan import download_IERS_A, get_IERS_A_or_workaround

max_wait_time = 5  # seconds

get_IERS_data = 'None'


# def check():
#     time.sleep(5)
#     if get_IERS_data != None:
#         return


# Thread(target=check).start()

get_IERS_data = input(
    'Do you want to download IERS data first or retrieve from cached data?: write yes or press enter ')


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

    def __init__(self):
        self.Exoplanets_List = []
        self.Exoplanets_List_Nasa = []
        self.ExoError = []
        self.Fail = []
        self.Exoplanets_altern_T0 = []
        self.Exoplanet_not_Nasa = []

        self.period = []
        self.name = []
        self.T0 = []

    def add_exoplanets(self, name):  # maybe nicer with tables
        """name of planet must be in form strings such that ExoplanetOrbitDatabase can read them"""
        try:
            Exo = ExoplanetOrbitDatabase.query_planet(
                name)  # General Exoplanet Data Source
            if not Exo:
                Warning('Exoplanet' + name + 'could not be parsed')
                self.ExoError.append(name)
        except KeyError:
            pass

    def altern_Planet_finder(self, name):  # maybe nicer with tables
        # Nasa Archive, primary source for Exoplanet data if possible
        Planet_try = NasaExoplanetArchive.query_planet(name, all_columns=True)
        print('Planet ' + name + ' found in Nasa Exoplanet Archive\n')
        self.Exoplanets_List_Nasa.append(Planet_try)
        if not Planet_try:
            print('Planet not in Nasa Exoplanet Archive, try next Database\n')
            Exoplanets.Exoplanet_not_Nasa.append(name)
            # for planet in cat.findall('.//planet'):
            #     try:
            #         if oec.findvalue(planet, 'name') == name:
            #             print(oec.findvalue(planet, 'name'))
            #             pass  # include adding planet data from the third Archive here
        # else:
        #             Warning('Planet ' + name +
        #                     ' is in no available Database, add manually\n')
        #             # self.add_manually()


def hasproperties(planet_properties):  # maybe nicer with tables
    Transit_data_missing = []
    Transit_data_avail = []
    Parse_planets_Nasa = []
    for n in range(len(planet_properties)):
        planet = planet_properties[n]
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
            Transit_data_avail.append(planet['pl_name'][0])
            Parse_planets_Nasa.append(planet)
            print('Planet ' + planet['pl_name'][0] +
                  ' added to Transit_data_avail\n')
        else:
            Transit_data_missing.append(planet['pl_name'][0])
            print('Planet ' + planet['pl_name'][0] +
                  ' added to Transit_data_missing\n')

    return Transit_data_missing, Transit_data_avail, Parse_planets_Nasa

    # def add_manually(self): #define a funcion to add manyally a planet
    #    if #input = 'skip'
    #        pass
    #    else:
    #        #ask the user to input the different columns.
    #        pass


try:
    # Name_list includes all the exoplanet's names. Ref_list lists the source for the exoplanet's data.
    Name_list = csv_file_import.main()
except:
    raise Warning('csv file is corrupted')


Exoplanets = Exoplanets()

for k in range(len(Name_list)):
    if type(Name_list[k]) is not str:
        Name_list[k] = str(Name_list[k])
        Warning('Name is not a string, cannot parse ' + Name_list[k])

#Transit_data_missing = []
#Transit_data_available = []

for name in Name_list:
    try:
        #        Exoplanets.add_exoplanets(name)
        Exoplanets.altern_Planet_finder(name)
    except Exception:
        Exoplanets.Fail.append(name)
        print('Planet ' + name + ' added to list "Fail"\n')
        Exoplanets.Fail.append(Name_list.index(name))

# work with Parse_list_transits in the same way as in the Transit_Example to retrieve all necessary transit information
Transit_data_missing, Transit_data_available, Parse_planets_Nasa = hasproperties(
    Exoplanets.Exoplanets_List_Nasa)
print('Found {} planets in Nasa Archive with Transit data'.format(
    len(Parse_planets_Nasa)))


# To access for instance the name of a current list element write Exoplanets.Exoplanets_List[1]['NAME']
#-------------------- sky coordinates for the current time write Exoplanets.Exoplanets_List[1]['sky_coord']
# When calculating the sky_coordinates for a specific observation, the DATABASE for the particular EXOPLANET must be UPDATED!


# Location and UTC offset Paranal
paranal_loc = EarthLocation(
    lat=-24.627 * u.deg, lon=-70.405 * u.deg, height=2635.43 * u.m)
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


Per_max = []
for i in range(len(Parse_planets_Nasa)):
    Per_max.append(np.float64(Parse_planets_Nasa[i]['pl_orbper'] / u.day))

# define the maximum of days it should plan observations into the future in one step
Max_Delta_days = int(max(Per_max) * 2)
midnights_2020_til_2024 = []
delta_midnight = np.linspace(-12, 12, 1000) * u.hour

Nights_paranal = []
frame_2020_til_2024 = []
sunaltazs_2020_til_2024 = []


class Nights:
    def __init__(self):
        self.date = []
        self.coord = []
        self.night = []


Nights = Nights()
d_end = d + dt * Max_Delta_days

"Make that maybe faster in some way"
"Calculation of the nights at paranal and the coordinate frame for paranal"


def Calculate_nights_paranal(Nights, d, d_end):
    print('Calculating the nights of paranal from the ' +
          str(d.date()) + ' until the ' + str(d_end.date()) + '.\n')
    for k in range(Max_Delta_days):
        midnights_2020_til_2024.append(d + dt * k)
        Nights_paranal.append(
            Time(str(midnights_2020_til_2024[k])) - utcoffset + delta_midnight)  # introduce
        # compute frame AltAz for get_sun
        frame_2020_til_2024.append(
            AltAz(obstime=Nights_paranal[k], location=paranal_loc))
        Sun_per_day = get_sun(Nights_paranal[k]).transform_to(
            frame_2020_til_2024[k])
        # access sun coord. for a specific date and time via sunaltazs_2020_til_2024[0](day)[0](time)
        sunaltazs_2020_til_2024.append(Sun_per_day)

        Nights.date.append(
            sunaltazs_2020_til_2024[k][0].obstime.datetime.date())
        for n in range(len(delta_midnight)):
            Time_sun = (
                sunaltazs_2020_til_2024[k][n].obstime.value.split(" ")[1])
            if sunaltazs_2020_til_2024[k][n].alt < -18 * u.deg:
                Nights.coord.append({
                    'Date': str(Nights.date[k]),
                    'Time': Time_sun,
                    'Az': sunaltazs_2020_til_2024[k][n].az,
                    'Alt': sunaltazs_2020_til_2024[k][n].alt})

    Nights_paranal_table = pd.DataFrame(
        Nights.coord)  # All nights in Paranal in UTC
    Nights_paranal_dates = pd.DataFrame(Nights.date)
    Nights_paranal_dates.rename(columns={0: 'Date'}, inplace=True)
    return Nights_paranal_table, Nights_paranal_dates


Nights_paranal_table, Nights_paranal_dates = Calculate_nights_paranal(
    Nights, d, d_end)
Nights_paranal_table.to_csv('Nights_time_at_paranal.csv')
Nights_paranal_dates.to_csv('dates_timespan.csv')
time_list = []

"Splitting the nights in to single nights, this can be combined with the paranal night calculation part"
for n in range(len(Nights.date)):
    night = []
    for k in range(len(Nights.coord)):
        T = Nights.coord[k]['Date'] + ' ' + Nights.coord[k]['Time']  # UTC
        time_list.append(T)
        if Nights.coord[k]['Date'] == str(Nights.date[n]):
            night.append(T)
    Nights.night.append(Time(night))

Time_list = Time(time_list)

Target_list = []
Night_cons_per_night = []

# Include moon data here


Altcons = astroplan.AltitudeConstraint(
    min=+30 * u.deg, max=None)  # Altitude constraints definition
Airmasscons = astroplan.AirmassConstraint(
    min=None, max=1.7)  # Airmass constraints definition
for n in range(len(Parse_planets_Nasa)):
    # How fixed resp. timedep are these target coordinates?
    Target_list.append(astroplan.FixedTarget(
        Parse_planets_Nasa[n]['sky_coord'], name=Parse_planets_Nasa[n]['pl_name'][0]))

for m in range(len(Nights.night)):  # Astronomical Nighttime constraints definition
    dn_min = Nights.night[m][0]
    dn_max = Nights.night[m][-1]
    Night_cons_per_night.append(astroplan.LocalTimeConstraint(
        min=dn_min.datetime.time(), max=dn_max.datetime.time()))


class Eclipses:
    def __init__(self):
        self.name = []
        self.epoch = []
        self.period = []
        self.transit_duration = []
        self.Planets_eclipse = []
        self.Alt_window = []
        self.Airmass_window = []
        self.Coordinates = []
        self.Time_window = []


Eclipses = Eclipses()

"Create Eclipse Objects here and check if and when the primary eclipses should be observable. The results"
"are added to the class objects Eclipses"

target_observable = []
eclipse_observable = []
Eclipse_observable_datetime = []
num_eclipses = []

obs_time = Nights.night[0][0]


# Ugly, How do one access the columns of the planet objects like lists?
for k in range(len(Parse_planets_Nasa)):
    Eclipses.name.append(Parse_planets_Nasa[k]['pl_name'][0])
    Eclipses.epoch.append(
        Time(Parse_planets_Nasa[k]['pl_tranmid'], format='jd'))
    Eclipses.period.append(Parse_planets_Nasa[k]['pl_orbper'])
    Eclipses.transit_duration.append(
        Parse_planets_Nasa[k]['pl_trandur'] * u.day)  # unit is days according to Astroplan
    Planet_next_eclipse = astroplan.EclipsingSystem(primary_eclipse_time=Eclipses.epoch[k],
                                                    orbital_period=Eclipses.period[k], duration=Eclipses.transit_duration[k])
    Eclipses.Planets_eclipse.append(Planet_next_eclipse)
    num_eclipses.append(
        int(np.floor(Max_Delta_days / (Parse_planets_Nasa[k]['pl_orbper'][0] / u.day))))
    num_e = int(np.floor(Max_Delta_days /
                         (Parse_planets_Nasa[k]['pl_orbper'][0] / u.day)))
    Planet_next_eclipse_Times = Planet_next_eclipse.next_primary_eclipse_time(
        obs_time, n_eclipses=num_e)
    try:
        Eclipses.Coordinates.append(SkyCoord.from_name(
            Parse_planets_Nasa[k]['pl_name'][0]))
    except Exception:
        Eclipses.Coordinates.append(Parse_planets_Nasa[k]['sky_coord'])

    print(Parse_planets_Nasa[k]['pl_name'][0] + ' is getting processed')

    for n in range(len(Nights.date)):
        obs_date = Nights.date[n]
        obs_night = Nights.night[n]
        # Loop over all eclipses coming up in the given timespan of object [k]
        for m in range(len(Planet_next_eclipse_Times)):
            ple_by_date = Planet_next_eclipse_Times[m]
            if obs_date == ple_by_date.datetime.date():  # Check which eclipse can be observed in which night
                Planet_next_eclipse_per_night_MID = Planet_next_eclipse_Times[m]
                Planet_next_eclipse_per_night_BEGIN = Planet_next_eclipse_per_night_MID + \
                    Parse_planets_Nasa[k]['pl_trandur'] / 2 * u.day
                Planet_next_eclipse_per_night_END = Planet_next_eclipse_per_night_MID + \
                    Parse_planets_Nasa[k]['pl_trandur'] / 2 * u.day
                Planet_Eclipes_NIGHT = [Planet_next_eclipse_per_night_BEGIN, Planet_next_eclipse_per_night_MID,
                                        Planet_next_eclipse_per_night_END]  # Begin, midpoint and end of transit

                """Careful, this only computes if the complete transit is observable"""
                ecl_obs = astroplan.is_event_observable(constraints=[
                                                        Altcons, Airmasscons, Night_cons_per_night[n]], observer=paranal, target=Eclipses.Coordinates[k], times=Planet_Eclipes_NIGHT)
                ecl_mid_obs = astroplan.is_event_observable(constraints=[
                    Altcons, Airmasscons, Night_cons_per_night[n]], observer=paranal, target=Eclipses.Coordinates[k], times=Planet_next_eclipse_per_night_MID)
                tar_obs = astroplan.is_event_observable(constraints=[
                                                        Altcons, Airmasscons], observer=paranal, target=Eclipses.Coordinates[k], times=Nights.night[n])
                eclipse_observable.append({
                    'Name': Parse_planets_Nasa[k]['pl_name'][0],
                    'obs_time': Planet_next_eclipse_per_night_MID,
                    'Primary eclipse observable?': ecl_obs[0][0],
                    'Eclipse Begin': Planet_next_eclipse_per_night_BEGIN[0],
                    'Eclipse End': Planet_next_eclipse_per_night_END[0]})
                for l in range(len(Nights.night[n])):
                    target_observable.append({
                        'Name': Parse_planets_Nasa[k]['pl_name'][0],
                        'Object observable?': tar_obs[0][l],
                        'Obs Night Time': Nights.night[n][l]})

        Alt_constraints = Altcons.compute_constraint(
            times=Nights.night[n], observer=paranal, targets=Eclipses.Coordinates[k])
        Airmass_constraints = Airmasscons.compute_constraint(
            times=Nights.night[n], observer=paranal, targets=Eclipses.Coordinates[k])
        Night_constraint = Night_cons_per_night[n].compute_constraint(
            times=Nights.night[n], observer=paranal, targets=Eclipses.Coordinates[k])
        Eclipses.Airmass_window.append(Alt_constraints)
        Eclipses.Alt_window.append(Airmass_constraints)
        Eclipses.Time_window.append(Night_constraint)


try:
    del table_eclipse_observable
    del table_object_observable
except:
    pass

table_eclipse_observable = pd.DataFrame(data=eclipse_observable)
table_object_observable = pd.DataFrame(data=target_observable)


text_file_name = input('Write Filename to save file: ')
def_text_file_name = 'Observation_Timetable_Eclipse'
def_text_file_name2 = 'Observation_Timetable_Objects.csv'

if not text_file_name:
    text_file_name = def_text_file_name

text_file_name = text_file_name + '.csv'
table_eclipse_observable.to_csv(text_file_name)
direc = os.getcwd()
print(text_file_name + ' is created in ' + direc)
table_object_observable.to_csv(def_text_file_name2)

# print(table_eclipse_observable)
