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

ditSTD = 10 # default value for DIT
nditSTD = 1 # default value for NDIT

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
        if inputtype == "snr":
            etc_obj = pd.read_json('etc-form-default-snr.json')
        elif inputtype == "ndit":
            etc_obj = pd.read_json('etc-form-default-ndit.json')
        else:
            raise KeyError("wrong inputtype: {}".format(inputtype))
        self.etc = etc_obj

    def update_etc_form(self, **kwargs):
        """
        changes constrains in 'etc-form.json'
        
        Parameters: Keyword arguments recognized by update_etc_form
        -----------
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
        
        others:...
        
        """
        
        if "airmass" in kwargs:
            self.etc.sky.airmass = kwargs.get("airmass")
        if "moon_target_sep" in kwargs:
            self.etc.sky.moon_target_sep = kwargs.get("moon_target_sep")
        if "moon_sun_sep" in kwargs:
            self.etc.sky.moon_sun_sep = kwargs.get("moon_sun_sep")

        if "snr" in kwargs:
            self.etc.timesnr.snr = kwargs.get("snr")
        else:
            self.etc.timesnr.snr = 100 # default signal to noise ratio: 100
        if "dit" in kwargs:
            self.etc.timesnr.dit = kwargs.get("dit")
        if "temperature" in kwargs:
            self.etc.target.sed.loc[:, ('spectrum','params','temperature')] = kwargs.get("temperature")
        if "inputtype" in kwargs:
            self.etc.timesnr.inputtype = kwargs.get("inputtype")
        if self.etc.timesnr.inputtype == "snr":
            try:
                self.etc.timesnr.dit = kwargs.get("dit")
            except Exception:
                raise Warning(
                    'dit not defined, using standard:{}'.format(ditSTD))
                self.etc.timesnr.dit = ditSTD
        elif self.etc.timesnr.inputtype == "ndit":
            try:
                self.etc.timesnr.dit = kwargs.get("dit")
            except Exception:
                raise Warning(
                    'dit not defined, using standard:{}'.format(ditSTD))
                self.etc.timesnr.dit = ditSTD
            try:
                self.etc.timesnr.ndit = kwargs.get("ndit")
            except Exception:
                raise Warning(
                    'dit not defined, using standard:{}'.format(nditSTD))
                self.etc.timesnr.ndit = nditSTD

    def write_etc_format_file(self):
        """
        Writes self.etc to a new JSON file named 'etc-form.json' such
        that it can be interpreted by the ETC online-calculator.

        """
        Etc_write = self.etc
        Etc_write.to_json('etc-form.json', indent=2)

    

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
        return_code = os.system('./etc_cli.py crires etc-form.json -o etc-data.json')
        time.sleep(1)
        output = pd.read_json('etc-data.json')
        try:
            NDIT = output.data.time['ndit']
        except AttributeError:
            raise Warning('NDIT not available from output file')
            NDIT = 0
         
        if self.etc.timesnr.inputtype == "ndit":
            NDIT = self.etc.timesnr.ndit
        return NDIT, output
