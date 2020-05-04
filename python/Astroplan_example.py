#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 17:19:08 2020

@author: jonaszbinden
"""

from astroplan.plots import plot_finder_image
from astroplan import FixedTarget
import matplotlib.pyplot as plt

messier1 = FixedTarget.from_name("M1")
ax, hdu = plot_finder_image(messier1)
plt.show()