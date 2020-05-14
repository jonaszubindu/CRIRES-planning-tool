#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  2 10:30:00 2020

This file contains the class etc_form to read in, change and update the input file for the ETC calculator 'etc-form.json'.
IMPORTANT: Do not change the file 'etc-form-default.json'


@author: jonaszbinden
"""
import os
import pandas as pd
import time
import json
from json import JSONEncoder
try:
    from types import SimpleNamespace as Namespace
except ImportError:
    # Python 2.x fallback
    from argparse import Namespace


ditSTD = 10 # default value for DIT
nditSTD = 1 # default value for NDIT
  
class FormEncoder(JSONEncoder):
    def default(self, o): return o.__dict__
           
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
    

    def update_etc_form(self, **kwargs):
        """
        changes constrains in 'etc-form.json'
        
        Parameters: 
        -----------
        Keyword arguments recognized by update_etc_form:
            
        airmass : float

        moon_target_sep : float
        
        moon_sun_sep : float
        
        snr : int or float
            minimum signal to noise ratio S/N 
        
        dit : int or float
            DIT exposure time for single exposure
        
        ndit : int
            NDIT number of single exposures for one single observation
            
            NDIT*DIT = Texp total exposure time for one single observation
        
        inputtype : string
            snr or ndit depending on ETC calculator should calculate the NDIT for a certain minimum S/N 
                    or S/N for a certain NDIT
                   
        temperature : float
            effective temperature of the target object
            
        brightness : float
            object brightness, standard is J-band magnitude, system: AB
        
        others:...
        
        """
        
        if "airmass" in kwargs:
            self.input.sky.airmass = kwargs.get("airmass")
        if "moon_target_sep" in kwargs:
            self.input.sky.moon_target_sep = kwargs.get("moon_target_sep")
        if "moon_sun_sep" in kwargs:
            self.input.sky.moon_sun_sep = kwargs.get("moon_sun_sep")

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
            
            
    def write_etc_format_file(self):
        """
        Writes self.etc to a new JSON file named 'etc-form.json' such
        that it can be interpreted by the ETC online-calculator.

        """
        Etc_write = self.input
        with open('etc-form.json','w') as Dump:
            json.dump(Etc_write, Dump, indent=2, cls=FormEncoder)
    

    def run_etc_calculator(self):
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
        try:    
            return_code = os.system('./etc_cli.py crires etc-form.json -o etc-data.json')
            if return_code == 256:
                raise Warning("VPN connection to ESO server for ETC calculator is not established!")
        except Exception:
            raise Warning('Something went wrong with the ETC calculator')
            
        time.sleep(1)
        with open('etc-data.json') as args: 
            output = json.load(args, object_hook=lambda d: Namespace(**d))
        try:
            NDIT = output.data.time.ndit
        except Exception:
            raise Warning('NDIT not available from output file')
            NDIT = 0
         
        if self.input.timesnr.inputtype == "ndit":
            NDIT = self.input.timesnr.ndit
        return NDIT, output
