#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 10:00:27 2020

@author: jonaszbinden
"""

import pandas as pd
import sys, os

path = os.getcwd()
default_file = 'PlanetListTest'

try:
    file_name = sys.argv[1]
except IndexError:
    file_name = default_file


if len(file_name) == 0:
    file_name = default_file
    print('Default file ' + default_file + ' is being used')
elif os.path.exists(path + '/' + file_name + '.csv') == False:
    message = "Error:file {} does not exist in current directory".format(sys.argv[1])
    err_msg = NameError(message)
    raise err_msg

def main():
    df = pd.read_csv(path + '/' + file_name + '.csv')
    df_names = df['pl_name']
    return df_names


if __name__ == "__main__":
    main()
