#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 13:00:07 2020

@author: jonaszbinden
"""
import numpy as np


def main():
    print("__name__ is = __main__ now")


import time
from threading import Thread

answer = None


def check():
    time.sleep(5)
    if answer != None:
        return
    print("Too Slow")


Thread(target=check).start()

answer = raw_input("Input something: ")


class Name:
    def __init__(self, name):
        self.name = name

    @classmethod
    def do_something(self, val1, val2):
        val = val1 * val2
        print("I can use this anywhere")
        return val


class Dog(Name):

    tricks = ['looping', 'trick']
    numbers = [1, 2, 3, 4, 5]             # mistaken use of a class variable

    def __init__(self, Name):
        self.name = Name.name
        self.owner = []

    def add_trick(self, trick):
        self.tricks.append(trick)

    def add_owner(self, owner):
        self.owner.append(owner)

    @staticmethod
    # x.num = x.calculate_something() to make this work, can create a new instance object by writing x.num =
    def calculate_something(numbers):
        sum_of_first_and_second = numbers[0] + numbers[1]
        return sum_of_first_and_second


if __name__ == "__main__":
    main()


# getattr(x,'name') is equivalent to x.name
# hasattr...
# setattr...


try:
    from astroplan import download_IERS_A
    download_IERS_A()
    print("IERS data successfully downloaded")
except:
    print("an error occured")

# continue doing some code here:
x = range(5)
print(x)

import datetime
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_moon, get_sun
from astropy.time import Time

dt = datetime.timedelta(days=1)
d = datetime.datetime(2020, 4, 1, 0, 0, 0)

paranal = EarthLocation(lat=-24.627 * u.deg, lon=-
                        70.405 * u.deg, height=2635.43 * u.m)
utcoffset = -4 * u.hour


midnights_2020_til_2024 = []
delta_midnight = np.linspace(-12, 12, 1000) * u.hour

Nights_paranal = []
frame_2020_til_2024 = []
sunaltazs_2020_til_2024 = []


for k in range(20):
    midnights_2020_til_2024.append(d + dt * k)
    Nights_paranal.append(
        Time(str(midnights_2020_til_2024[k])) - utcoffset + delta_midnight)
    frame_2020_til_2024.append(
        AltAz(obstime=Nights_paranal[k], location=paranal))
    Sun_per_day = get_sun(Nights_paranal[k]).transform_to(
        frame_2020_til_2024[k])
    sunaltazs_2020_til_2024.append(Sun_per_day)


paranal = EarthLocation(lat=-24.627 * u.deg, lon=-
                        70.405 * u.deg, height=2635.43 * u.m)
utcoffset = -4 * u.hour
time = Time('2020-04-01 12:00') - utcoffset


GJ9827b_prop = NasaExoplanetArchive.query_planet('GJ-9827 b', all_columns=True)
GJ9827b = SkyCoord.from_name('GJ9827b')


epoch2 = Time(GJ9827b_prop['pl_tranmid'], format='jd')
period2 = GJ9827b_prop['pl_orbper']
transit_duration2 = GJ9827b_prop['pl_trandur'] * u.day

GJ9827b_trans = EclipsingSystem(primary_eclipse_time=epoch2, orbital_period=period2,
                                duration=transit_duration2)

midnight = Time('2020-7-13 00:00:00') - utcoffset

delta_midnight = np.linspace(-12, 12, 1000) * u.hour
times_July12_to_13 = midnight + delta_midnight
frame_July12_to_13 = AltAz(obstime=times_July12_to_13, location=paranal)
sunaltazs_July12_to_13 = get_sun(
    times_July12_to_13).transform_to(frame_July12_to_13)

GJ9827baltaz = GJ9827b.transform_to(AltAz(obstime=time, location=paranal))


nights = []


import requests
r = requests.get('https://xkcd.com/353/', timeout=3)
# print(r.ok)
# print(r.text)

import os
path = os.getcwd()
# for f in os.listdir(): #input default = current working directory
#    print(f)
# find more about name splitting etc. in the video about automated files in python
