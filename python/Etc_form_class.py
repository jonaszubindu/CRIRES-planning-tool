
import json
import os
import astropy
import astroplan
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
    def __init__(self):
        """
        Initializes 'etc-form-default.json' via pandas to a dataframe object.

        Returns
        -------
        etc-form as object:
            
        self.etc.

        """
        etc_obj = pd.read_json('etc-form-default.json')
        self.etc = etc_obj

    def update_etc_form(self, **kwargs):
        """
        
        changes constrains in 'etc-form.json'
        
        Returns
        -------
        None.
        
        """
        
        if "airmass" in kwargs:
            etc_form.Airmass(self, kwargs.get("airmass"))
        if "moon_target_sep" in kwargs:
            etc_form.Moon_target_sep(self, kwargs.get("moon_target_sep"))
        if "moon_sun_sep" in kwargs:
            etc_form.Moon_sun_sep(self, kwargs.get("moon_sun_sep"))

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

        Returns
        -------
        None.

        """
        Etc_write = self.etc
        Etc_write.to_json('etc-form.json', indent=2)

    def Airmass(self):
        pass

    def Moon_target_sep(self):
        pass

    def Moon_sun_sep(self):
        pass



    def run_etc_calculator(self):
        """
        Runs ETC calculator through commandline and asks for output data file
    
        Returns
        -------
        NDIT : INT
            Number of single exposures with DIT to reach 
            signal to noise S/N as defined in 'etc-form.json'.
        output : PANDAS DATAFRAME
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
