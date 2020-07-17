#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 15:26:43 2020

File containing example functions from astropy and astroplan on which the whole software tool is built.

@author: jonaszbinden
"""

from astropy.time import Time
import astropy.units as u
import astroplan
#import astropy.coordinates
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_moon
import numpy as np

import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style, quantity_support
plt.style.use(astropy_mpl_style)
quantity_support()
from astropy.coordinates import get_sun
from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive
from astroquery.exoplanet_orbit_database import ExoplanetOrbitDatabase
#primary_eclipse_time = Time(2452826.628514, format='jd')
#orbital_period = 3.52474859 * u.day
#eclipse_duration = 0.1277 * u.day
#
#hd209458 = EclipsingSystem(primary_eclipse_time=primary_eclipse_time,
#                           orbital_period=orbital_period, duration=eclipse_duration,
#                           name='HD 209458 b')



planet_properties = NasaExoplanetArchive.query_planet('TRAPPIST-1 b', all_columns=True)


# get relevant planet properties
epoch = Time(planet_properties['pl_tranmid'], format='jd')
period = planet_properties['pl_orbper']
transit_duration = planet_properties['pl_trandur'] * u.day

# Create an EclipsingSystem object for HD 209458
#from astroplan import EclipsingSystem
trappist1b = astroplan.EclipsingSystem(primary_eclipse_time=epoch, orbital_period=period,
                             duration=transit_duration)

# Calculate next three mid-transit times which occur after ``obs_time``
obs_time = Time('2017-01-01 12:00')
trappist1b.next_primary_eclipse_time(obs_time, n_eclipses=3)


paranal_loc = EarthLocation(lat = -24.627*u.deg, lon = -70.405*u.deg, height = 2635.43*u.m)
utcoffset = -4*u.hour
time = Time('2020-04-01 12:00') - utcoffset


GJ9827b_prop = NasaExoplanetArchive.query_planet('GJ-9827 b', all_columns=True)
GJ9827b = SkyCoord.from_name('GJ9827b')


epoch2 = Time(GJ9827b_prop['pl_tranmid'], format = 'jd')
period2 = GJ9827b_prop['pl_orbper']
transit_duration2 = GJ9827b_prop['pl_trandur'] * u.day

GJ9827b_trans = astroplan.EclipsingSystem(primary_eclipse_time=epoch2, orbital_period=period2,
                          duration=transit_duration2)

midnight = Time('2020-7-13 00:00:00') - utcoffset

delta_midnight = np.linspace(-12, 12, 1000)*u.hour
times_July12_to_13 = midnight + delta_midnight
frame_July12_to_13 = AltAz(obstime=times_July12_to_13, location=paranal_loc)
sunaltazs_July12_to_13 = get_sun(times_July12_to_13).transform_to(frame_July12_to_13)

GJ9827baltaz = GJ9827b.transform_to(AltAz(obstime=time,location=paranal_loc))


delta_midnight = np.linspace(-2, 10, 100)*u.hour
obstime=midnight+delta_midnight

paranal = astroplan.Observer(paranal_loc, timezone='Etc/GMT-4')

Altcons = astroplan.AltitudeConstraint(min = +30*u.deg, max = None)
Airmasscons = astroplan.AirmassConstraint(min = None, max = 3.0)

Alt_constraints = Altcons.compute_constraint(times=obstime, observer=paranal, targets=GJ9827b)
Airmass_constraints = Airmasscons.compute_constraint(times=obstime, observer=paranal, targets=GJ9827b)

observable = astroplan.is_observable(constraints = [Altcons, Airmasscons], observer=paranal, targets=GJ9827b, times=obstime)

midnight1 = Time('2020-7-13 00:00:00') - utcoffset
midnight2 = Time('2021-1-13 00:00:00') - utcoffset

frame_July13night1 = AltAz(obstime=midnight1+delta_midnight, location=paranal_loc)
frame_July13night2 = AltAz(obstime=midnight2+delta_midnight, location=paranal_loc)
GJ9827baltazs_July13night1 = GJ9827b.transform_to(frame_July13night1)
GJ9827baltazs_July13night2 = GJ9827b.transform_to(frame_July13night2)

GJ9827bairmasss_July13night = GJ9827baltazs_July13night1.secz

plt.plot(delta_midnight, GJ9827bairmasss_July13night)
plt.xlim(-2, 10)
plt.ylim(1, 4)
plt.xlabel('Hours from EDT Midnight')
plt.ylabel('Airmass [Sec(z)]')
plt.show()

delta_midnight = np.linspace(-12, 12, 1000)*u.hour
times_July12_to_13 = midnight + delta_midnight
frame_July12_to_13 = AltAz(obstime=times_July12_to_13, location=paranal_loc)

GJ9827baltazs_July12_to_13 = GJ9827b.transform_to(frame_July12_to_13)

moon_July12_to_13 = get_moon(times_July12_to_13)
moonaltazs_July12_to_13 = moon_July12_to_13.transform_to(frame_July12_to_13)


plt.plot(delta_midnight, sunaltazs_July12_to_13.alt, color='r', label='Sun')
plt.plot(delta_midnight, moonaltazs_July12_to_13.alt, color=[0.75]*3, ls='--', label='Moon')
plt.scatter(delta_midnight, GJ9827baltazs_July12_to_13.alt,
            c=GJ9827baltazs_July12_to_13.az, label='GJ-9827 b', lw=0, s=8,
            cmap='viridis')
plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg,
                 sunaltazs_July12_to_13.alt < -0*u.deg, color='0.5', zorder=0)
plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg,
                 sunaltazs_July12_to_13.alt < -18*u.deg, color='k', zorder=0)
plt.colorbar().set_label('Azimuth [deg]')
plt.legend(loc='upper left')
plt.xlim(-12*u.hour, 12*u.hour)
plt.xticks((np.arange(13)*2-12)*u.hour)
plt.ylim(0*u.deg, 90*u.deg)
plt.xlabel('Hours from EDT Midnight')
plt.ylabel('Altitude [deg]')
plt.show()


K2_24_c_transit = ExoplanetOrbitDatabase.query_planet('K2-24 c')
K2_24_c = SkyCoord.from_name('K2-24 c')

from astroplan.plots import plot_finder_image
from astroplan import FixedTarget
import matplotlib.pyplot as plt

messier1 = FixedTarget.from_name("M1")
ax, hdu = plot_finder_image(messier1)
plt.show()