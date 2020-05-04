"""
Include ETC constraints here as a different mode to compute
additional independent constraints

This can be advanced by any method to change some input parameter of
'etc-form.json' for any type of targets.

WARNING: If the general structure of the file changes due to for instance
         change from inputtype "Spectrum" to "Emission Line", this must be regarded
         when adding methods to alter 'etc-form.json'. Might conflict with other methods!
"""
import json
import os
import astropy
import astroplan
import pandas as pd
import time

ditSTD = 10

class etc_form:

    def __init__(self):
        etc_obj = pd.read_json('etc-form.json')
        self.etc = etc_obj

    def update_etc_form(self, **kwargs):
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
                    'dit not defined, using standards:{}'.format(ditSTD))
                self.ETCForm.dit = ditSTD

    def write_etc_format_file(self):
        Etc_write = self.etc
        Etc_write.to_json('etc-form.json', indent=2)

    def Airmass(self):
        pass

    def Moon_target_sep(self):
        pass

    def Moon_sun_sep(self):
        pass



def run_etc_calculator():
    return_code = os.system('./etc_cli.py crires etc-form.json -o etc-data.json')
    time.sleep(1)
    output = pd.read_json('etc-data.json')

    NDIT = output.data.time['ndit']
    return NDIT
