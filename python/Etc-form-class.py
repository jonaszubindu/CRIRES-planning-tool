"""
Include ETC constraints here as a different mode to compute additional independent constraints
"""
import json
import etc_cli.py as etc
import os
import astropy
import astroplan

ditSTD = 10

class etc_form:

    def __init__(self, ETCForm):
        with open('etc-form.json', 'r') as etc_form_raw:
            ETCForm = json.load(Input)

        self = ETCForm
        # for key, value in ETCForm.items():
        #     self.key = value

    def update_etc_form(self, **kwargs):

        if "airmass" in kwargs:
            Airmass(self, kwargs.get("airmass"))
        if "moon_target_sep" in kwargs:
            Moon_target_sep(self, kwargs.get("moon_target_sep"))
        if "moon_sun_sep" in kwargs:
            Moon_sun_sep(self, kwargs.get("moon_sun_sep"))

        if "dit" in kwargs:
            self.dit = kwargs("dit")
        if "temperature" in kwargs:
            self.temperature = kwargs.get("temperature")
        if "spectrumtype" in kwargs:
            self.spectrumtype = kwargs.get("spectrumtype")
        if "inputtype" in kwargs:
            self.inputtype = kwargs.get("inputtype")
        if inp_type == "snr":
            try:
                self.dit = kwargs.get("dit")
                self.median = kwargs.get("median")
            except Exception:
                raise Warning(
                    'dit and/or median not defined, using standards:{},{}'.format(ditSTD, medianSTD))
                self.dit = ditSTD
                self.median = medianSTD
        if inp_type == "ndit":
            self.ndit = kwargs.get("ndit")

    def Airmass(self):
        pass

    def Moon_target_sep(self):
        pass

    def Moon_sun_sep(self):
        pass

    def write_etc_format_file(object):
        with open('etc-form.json', 'w') as etc_form_end:
            json.dump(object.__dict__, etc_form_end, indent=2)


def run_etc_calculator():
    os.system('./etc-cli.py crires etc-form.json -o output1.json')
    with open('output1.json', 'r') as output:
        Output = json.load(output)

    NDIT = Output['NDIT']
    return NDIT
