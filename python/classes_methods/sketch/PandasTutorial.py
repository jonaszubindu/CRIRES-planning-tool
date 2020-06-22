#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 12:30:06 2020

@author: jonaszbinden
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import style
style.use('ggplot')

web_stats = {'Day':[1,2,3,4,5,6],
             'Visitors':[43,53,34,45,64,34],
             'Bounce_rate':[65,72,62,64,54,66]}

df = pd.DataFrame(web_stats)

print(df)

print(df.set_index('Day'))

df = df.set_index('Day')
#or
df.set_index('Day', inplace=True)

print(df['Visitors'])
#or
print(df.Visitors)

print(df[['Visitors','Bounce_rate']])

print(df.Visitors.tolist()) #doesn't work with the Visitors and Bounce_rate
