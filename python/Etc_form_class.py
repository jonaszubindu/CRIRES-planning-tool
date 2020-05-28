#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  2 10:30:00 2020

This file contains the class etc_form to read in, change and update the input file for the ETC calculator 'etc-form.json'.
IMPORTANT: Do not change the files 'etc-form-default-snr.json' or 'etc-form-default-ndir.json'.


@author: jonaszbinden
GitHub: jonaszubindu
"""
import os
# import pandas as pd
import time
import json
from json import JSONEncoder, JSONDecodeError
import etc_cli
import argparse
import requests
import logging
from functools import wraps
from misc import misc

def Etc_logger(orig_fun):
    """ Function to log execution of other functions """

    @wraps(orig_fun)
    def wrapper(*args, **kwargs):
        logging.info('ran {} with args:{}, and kwargs:{}'.format(orig_fun.__name__, args, kwargs))
        return orig_fun(*args, **kwargs)
        
    return wrapper

try:
    from types import SimpleNamespace as Namespace
except ImportError:
    # Python 2.x fallback
    from argparse import Namespace

"""Warnings"""
JSONDecodeWarning = Warning('Something went wrong processing the etc-form file... I will try to find the error for you')
NDITWarning = Warning('NDIT not available from output file')

def DecodeWarning(key,value):
    DecodeWarning = FutureWarning(f"the error is related to the present {key} input value: {value.__str__()}")
    return DecodeWarning

ErrorNotFoundWarning = DeprecationWarning('Sorry, I cannot find the error, check the file etc-format.json and try to run it manually on the ETC calculator webpage. \n Maybe she can help you find the error...')

ditSTD = 10 # default value for DIT
nditSTD = 1 # default value for NDIT
  
class FormEncoder(JSONEncoder):
    def default(self, o): 
        return o.__dict__
            
           
class etc_form:
    """
    Include ETC constraints here as a different mode to compute
    additional independent constraints
    
    This can be advanced by any method to change some input parameter of
    'etc-form.json' for any type of targets.
    
    WARNING: If the general structure of the file changes due to for instance
             change from inputtype "Spectrum" to "Emission Line", this must be regarded
             when adding methods to alter 'etc-form.json'. Might conflict with other methods!
    """
    @Etc_logger
    def __init__(self, inputtype):
        """
        Initializes 'etc-form-default.json' via pandas to a dataframe object.
        
        Parameters:
        -----------
        inputtype : string
            specify if the ETC-calculator should be run in S/N mode or in 
            NDIT mode.

        """
        try:
            if inputtype == "snr":
                with open('etc-form-default-snr.json') as args: 
                    etc_obj = json.load(args, object_hook=lambda d: Namespace(**d))
            elif inputtype == "ndit":
                with open('etc-form-default-ndit.json') as args: 
                    etc_obj = json.load(args, object_hook=lambda d: Namespace(**d))
            else:
                raise KeyError("wrong inputtype: {}".format(inputtype))
        except FileNotFoundError:
            raise FileNotFoundError("File 'etc-form-default-{}.json' is not existing or not in current directory".format(inputtype))
        self.input = etc_obj
    
    @Etc_logger
    def update_etc_form(self, **kwargs):
        """
        changes constrains in 'etc-form.json'
        
        Parameters: 
        -----------
        Keyword arguments recognized by update_etc_form:
            
        airmass : float

        moon_target_sep : list
            Two values, first value is moon_target_separation, second value is moon_alt, altitude above horizon of the moon.
        
        moon_phase : float
            Illumination of the moon, also known as moon_sun_separation.
        
        snr : int or float
            Minimum signal to noise ratio S/N.
        
        dit : int or float
            DIT exposure time for single exposure.
        
        ndit : int
            NDIT number of single exposures for one single observation.
            
            NDIT*DIT = Texp total exposure time for one single observation.
        
        inputtype : string
            snr or ndit depending on ETC calculator should calculate the NDIT for a certain minimum S/N 
            or S/N for a certain NDIT.
                   
        temperature : float
            Effective temperature of the target object.
            
        brightness : float
            Object brightness, standard is J-band magnitude, system: AB.
        
        others:...
        
        """
        
        if "airmass" in kwargs:
            self.input.sky.airmass = kwargs.get("airmass")
            # self.input.sky.airmass = 12 # Chabis Test
        if "moon_target_sep" in kwargs:
            self.input.sky.moon_target_sep = kwargs.get("moon_target_sep")[0]
            self.input.sky.moon_alt = kwargs.get("moon_target_sep")[1]
        if "moon_phase" in kwargs:
            self.input.sky.moon_sun_sep = kwargs.get("moon_phase")

        if "snr" in kwargs:
            self.input.timesnr.snr.snr = kwargs.get("snr")
        else:
            self.input.timesnr.snr.snr = 100 # default signal to noise ratio: 100
        if "dit" in kwargs:
            self.input.timesnr.dit = kwargs.get("dit")
        if "temperature" in kwargs:
            self.input.target.sed.spectrum.params.temperature = kwargs.get("temperature")
        if "brightness" in kwargs:
            self.input.target.brightness.params.mag = kwargs.get("brightness")
        if "inputtype" in kwargs:
            self.input.timesnr.inputtype = kwargs.get("inputtype")
        if self.input.timesnr.inputtype == "snr":
            if kwargs.get("dit") == None:
                self.input.timesnr.dit = ditSTD
            else:    
                self.input.timesnr.dit = kwargs.get("dit")
        elif self.input.timesnr.inputtype == "ndit":
            if kwargs.get("dit") == None:
                self.input.timesnr.ndit = ditSTD
            else:    
                self.input.timesnr.ndit = kwargs.get("ndit")
            if kwargs.get("ndit") == None:
                self.input.timesnr.ndit = ditSTD
            else:    
                self.input.timesnr.ndit = kwargs.get("ndit")
            
    @Etc_logger     
    def write_etc_format_file(self):
        """
        Writes self.etc to a new JSON file named 'etc-form.json' such
        that it can be interpreted by the ETC online-calculator.

        """
        Etc_write = self.input
        with open('etc-form.json','w') as Dump:
            json.dump(Etc_write, Dump, indent=2, cls=FormEncoder)
    
    @Etc_logger
    def run_etc_calculator(self, name, tim):
        """
        Runs ETC calculator through commandline and asks for output data file
    
        Returns
        -------
        NDIT : int
            Number of single exposures with DIT to reach 
            signal to noise S/N as defined in 'etc-form.json'.
        output : pandas DataFrame
            DataFrame object containing the output data from the
            ETC calculator
    
        """
        success = 0
        while success == 0:
            """ 
            tries to access the ETC calculator, in case of unstable internet connection, the user is asked to solve the problem 
            manually and confirm it by pressing 'enter/return' except a JSONDecodeError gets caught.
            If success = 1 the loop breaks and the output gets processed.
            
            """
            try:
                CallETC(args = ['crires', 'etc-form.json', '-o', 'etc-data.json'])
                print('ETC calculator successfully called for {},{}'.format(name, tim))
                success = 1
            except Exception as e:
                if type(e) == requests.exceptions.ConnectionError:
                    try: 
                        print('could not connect to ETC server, trying again...')
                        time.sleep(5) # waits for better server connection if connectionseams unstable
                        CallETC(args = ['crires', 'etc-form.json', '-o', 'output1.json'])
                    except Exception as e:
                        if type(e) == requests.exceptions.ConnectionError:
                            print(ConnectionError('Could not establish VPN connection to ETC server'))
                            misc.wait_for_enter(msg='Connection Error: Check Internet connection and press enter when problem resolved:')
    
                elif type(e) == json.decoder.JSONDecodeError:
                    raise e
                else:
                    raise e
        
        
        time.sleep(1)
        with open('etc-data.json') as args: 
            output = json.load(args, object_hook=lambda d: Namespace(**d)) # writes the ETC data to a output namespace object
        try:
            NDIT = output.data.time.ndit
        except Exception:
            raise NDITWarning # Warning
            NDIT = 0
         
        if self.input.timesnr.inputtype == "ndit":
            NDIT = self.input.timesnr.ndit
        return NDIT, output
    
    @Etc_logger 
    def etc_debugger(self, name, tim, temperature, brightness, airmass, moon_phase, moon_target_sep):
        """
        This tries to find the error in the etc-format file. As soon as the ETC calculator gets updated with better input error handling
        this function must be updated or replaced by additional error handling in the functions running the ETC calculator.

        Parameters
        ----------
        JSONDecodeError : Exception
            Handle of the JSONDecodeError that occurred while running the ETC calculator.
        temperature : float
            Temperature input parameter that was used.
            
        brightness : float
            Brightness input parameter that was used.
            
        airmass : float
            Airmass input parameter that was used.
            
        moon_phase : float
            Illumination of the moon, also known as moon_sun_separation.
            
        moon_target_sep : list
            Two values, first value is moon_target_separation, second value is moon_alt, altitude above horizon of the moon.

        Raises
        ------
        JSONDecodeError
            If the errornous parameter was found, raises the JSONDecodeError and reviels the faulty parameter.

        Returns
        -------
        None. If no raises occur, the etc_debugger tells the user that it has not found the error and gives the problem 
        back to the user

        """
        print('Something went wrong processing the etc-form file... I will try to fix it for you')
        print(name, tim, temperature, brightness, airmass, moon_phase, moon_target_sep)
        os.system('cp etc-form.json etc-form-copy.json')
        cls = type(self)
        ETC = cls.__new__(cls)
        ETC.__init__('snr')        
        ETC.update_etc_form(temperature = temperature)
        ETC.write_etc_format_file()
        try:
            NDIT, output = ETC.run_etc_calculator(name,tim)
        except JSONDecodeError:
            raise DecodeWarning('temperature',temperature) # Warning
        cls = type(self)
        ETC = cls.__new__(cls)
        ETC.__init__('snr')   
        ETC.update_etc_form(brightness = brightness)
        ETC.write_etc_format_file()
        try:
            NDIT, output = ETC.run_etc_calculator(name,tim)
        except JSONDecodeError:
            raise DecodeWarning('brightness', brightness) # Warning
        cls = type(self)
        ETC = cls.__new__(cls)
        ETC.__init__('snr')   
        ETC.update_etc_form(airmass = airmass)
        ETC.write_etc_format_file()
        try:
            NDIT, output = ETC.run_etc_calculator(name,tim)
        except JSONDecodeError: 
            raise DecodeWarning('airmass', airmass)
        cls = type(self)
        ETC = cls.__new__(cls)
        ETC.__init__('snr')   
        ETC.update_etc_form(moon_phase = moon_phase)
        ETC.write_etc_format_file()
        try:
            NDIT, output = ETC.run_etc_calculator(name,tim)
        except JSONDecodeError: 
            raise DecodeWarning('moon_phase', moon_phase)
            cls = type(self)
        ETC = cls.__new__(cls)
        ETC.__init__('snr')   
        ETC.update_etc_form(moon_target_sep = moon_target_sep)
        ETC.write_etc_format_file()
        try:
            NDIT, output = ETC.run_etc_calculator(name,tim)
        except JSONDecodeError: 
            raise DecodeWarning('moon_target_sep', moon_target_sep)
        print('I will continue with the next planet for now...')
        raise ErrorNotFoundWarning # Warning
        
@Etc_logger
def CallETC(args):
    """This part is extracted from etc-cli.py and is included here to ensure better error handling."""
    #-------------------------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(
    description='Call an ETC with input parameters and optionally an uploaded spectrum.\n'
    + 'Print the resulting JSON on stdout or optionally a file.\n'
    + 'Examples: \n'
    + './etc_cli.py crires etc-form.json -o output1.json\n'
    + './etc_cli.py crires etc-form-uploading.json -u upload.dat -o output2.json',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('etcname',
                        help='Name of instrument ETC to call, e.g. 4most')

    parser.add_argument('postdatafile',
                        help='Name of JSON file with ETC input parameters,\nlike the ETC input form')

    parser.add_argument('-u,', '--upload', dest="uploadfile",
                        help='Name of file with spectrum to upload.\nSee https://etc.eso.org/observing/etc/doc/upload.html')

    parser.add_argument('-c', '--collapse', action='store_true',
                        help='collapse output JSON data arrays to a short indicative strings')

    parser.add_argument('-i', '--indent', type=int, nargs='?', const=4,
                        help='Format the output JSON with indentation (default 4)')

    parser.add_argument('-o', '--outputfile', dest="outputfile",
                        help='Send the output to file')

    
    
    args = parser.parse_args(args=args)
    
    # baseurl = 'http://localhost:8000/observing/etc/etcapi/'
    baseurl = 'https://etctest.hq.eso.org/observing/etc/etcapi/'
    # baseurl = 'https://etc.eso.org/observing/etc/etcapi/'

    etcName = etc_cli.getEtcUrl(args.etcname)

    url = baseurl + etc_cli.getEtcUrl(etcName)

    jsondata = etc_cli.callEtc(args.postdatafile, url, args.uploadfile).json()

    etc_cli.output(jsondata, args)
    # ------------------------------------------------------------------------------------------------

