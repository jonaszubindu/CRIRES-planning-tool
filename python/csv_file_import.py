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
#print("Enter a file name with the planets you want to process. If you want to use the default list, press enter")
#file_name = input("enter file name with planetary names: ")
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

#df_dict_sources = df_dict['']
def main():
    df = pd.read_csv(path + '/' + file_name + '.csv')
    #df_dict = pd.DataFrame.to_dict(df, orient='list')
    #df_dict_names = df_dict['NAME']
    #df_json = pd.DataFrame.to_json(df)
    #Ref_list = pd.DataFrame(df, columns=['NAME','ORBREF','ORBURL'])
    df_names = df['pl_name']
    return df_names


if __name__ == "__main__":
    main()
