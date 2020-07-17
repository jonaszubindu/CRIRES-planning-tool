#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 12:14:08 2020

@author: jonaszbinden
"""


# import xlsxwriter
import pandas as pd
from classes_methods.Helper_fun import data_sorting_and_storing
import datetime
from classes_methods.classes import Nights
import sys
from datetimerange import DateTimeRange
from classes_methods.classes import load_Eclipses_from_file
import astropy.units as u
import astropy
ranked_obs_events = None

try:
    filename = sys.argv[1]
except Exception:
    filename = 'Eclipse_events_processed_2020-07-09_30d.pkl'

d = datetime.datetime.fromisoformat(filename.split('_')[3])
Max_Delta_days = int(filename.split('_')[4][:-5])
Nights = Nights(d, Max_Delta_days)

Eclipses_List = load_Eclipses_from_file(filename, Max_Delta_days)
ranking, df_gen, df_frame = data_sorting_and_storing(
    Eclipses_List, write_to_csv=0)


# df_gen.reset_index(inplace=True)

# def xlsx_writer(filename, df_gen, df_frame, ranked_obs_events = None):
# """
#     Function to call for customized creation of excel files to store the Exoplanet candidate data.
#     This function can be changed in any suitable way to highlight or modify prefered cell formats.

#     Parameters
#     ----------
#     filename : str
#         filename under which the xlsx file should be stored.

#     df_gen : pandas DataFrame
#         dataframe containing the candidate data to store.

#     df_frame : pandas DataFrame
#         dataframe containing the observation times to store.

#     ranked_observations : pandas DataFrame
#         dataframe containing the ranked observation time events to store.

#     Returns
#     -------
#     Stores xlsx file to csv_file folder.

# """

# Create a Pandas Excel writer using XlsxWriter as the engine.
writer = pd.ExcelWriter(filename.split('.')[0] + '.xlsx', engine='xlsxwriter')

df_frame.reset_index(inplace=True)

workbook = writer.book
# Set up a format
# book_format = workbook.add_format(properties={'bold': True, 'font_color': 'red'})
cell_format = workbook.add_format()

cell_format.set_pattern(1)  # This is optional when using a solid fill.
cell_format.set_bg_color('red')  # Highlights the background of the cell

# Create a sheet
worksheet1 = workbook.add_worksheet('Candidates')
worksheet1.set_column(0, 1, 25)
worksheet1.set_column(2, 12, 20)
worksheet2 = workbook.add_worksheet('Observations')
worksheet2.set_column(0, 12, 25)
if not ranked_obs_events:
    pass
else:
    worksheet3 = workbook.add_worksheet('Ranked Observations')
    worksheet3.set_column(0, 12, 25)
    for col_num, header in enumerate(ranked_obs_events.keys()):
        worksheet3.write(0, col_num, header)


# Write the headers
for col_num, header in enumerate(df_gen.keys()):
    worksheet1.write(0, col_num, header)
    

for col_num, header in enumerate(df_frame.keys()):
    worksheet2.write(0, col_num, header)

    
obs_time = []
# Save the data from the OrderedDict into the excel sheet
for row_num, row_data in enumerate(df_gen.values):
    for col_num, cell_data in enumerate(row_data):
        if col_num == 2 and cell_data > 1 / 24 * u.day:
            obs_time.append(df_gen['obs time'][row_num])
            try:
                worksheet1.write(row_num + 1, col_num, cell_data, cell_format)
            except TypeError:
                if type(cell_data) == astropy.time.Time:
                    cell_data = cell_data.value.isoformat()
                else:
                    cell_data = cell_data.value
                worksheet1.write(row_num + 1, col_num, cell_data, cell_format)
        else:
            try:
                worksheet1.write(row_num + 1, col_num, cell_data)
            except TypeError:
                if type(cell_data) == astropy.time.Time:
                    cell_data = cell_data.value.isoformat()
                else:
                    cell_data = cell_data.value
                worksheet1.write(row_num + 1, col_num, cell_data)
    
# Save the data from the OrderedDict into the excel sheet
for row_num in range(int(len(df_frame.values) / 3)):
    row_num = row_num * 3
    for col_num, _ in enumerate(df_frame.values[row_num]):
        if df_frame.values[row_num][1] == any(obs_time):
            try:
                worksheet2.write(row_num + 1, col_num, df_frame.values[row_num][col_num], cell_format)
                worksheet2.write(row_num + 2, col_num, df_frame.values[row_num + 1][col_num], cell_format)
                worksheet2.write(row_num + 3, col_num, df_frame.values[row_num + 2][col_num], cell_format)
            except TypeError:
                if type(df_frame.values[row_num][1]) == astropy.time.core.Time:
                    df_frame.loc[row_num,'time'] = df_frame.values[row_num][col_num].value.isoformat()
                    df_frame.loc[row_num + 1, 'time'] = df_frame.values[row_num +1][col_num].value.isoformat()
                    df_frame.loc[row_num + 2, 'time'] = df_frame.values[row_num +2][col_num].value.isoformat()
                else:
                    df_frame.loc[row_num, 'time'] = df_frame.values[row_num][col_num].value
                    df_frame.loc[row_num + 1, 'time'] = df_frame.values[row_num + 1][col_num].value
                    df_frame.loc[row_num + 2, 'time'] = df_frame.values[row_num + 2][col_num].value
    
                worksheet2.write(row_num + 1, col_num, df_frame.values[row_num][col_num], cell_format)
                worksheet2.write(row_num + 2, col_num, df_frame.values[row_num + 1][col_num], cell_format)
                worksheet2.write(row_num + 3, col_num, df_frame.values[row_num + 2][col_num], cell_format)
    
        else:
            try:
                worksheet2.write(row_num + 1, col_num, df_frame.values[row_num][col_num])
                worksheet2.write(row_num + 2, col_num, df_frame.values[row_num + 1][col_num])
                worksheet2.write(row_num + 3, col_num, df_frame.values[row_num + 2][col_num])
            except TypeError:
                if type(df_frame.values[row_num][1]) == astropy.time.core.Time:
                    df_frame.loc[row_num, 'time'] = df_frame.values[row_num][col_num].value.isoformat()
                    df_frame.loc[row_num + 1, 'time'] = df_frame.values[row_num + 1][col_num].value.isoformat()
                    df_frame.loc[row_num + 2, 'time'] = df_frame.values[row_num + 2][col_num].value.isoformat()
                else:
                    df_frame.loc[row_num, 'time'] = df_frame.values[row_num][col_num].value
                    df_frame.loc[row_num + 1, 'time'] = df_frame.values[row_num + 1][col_num].value
                    df_frame.loc[row_num + 2, 'time'] = df_frame.values[row_num + 2][col_num].value
                worksheet2.write(row_num + 1, col_num, df_frame.values[row_num][col_num])
                worksheet2.write(row_num + 2, col_num, df_frame.values[row_num + 1][col_num])
                worksheet2.write(row_num + 3, col_num, df_frame.values[row_num + 2][col_num])
    
if type(ranked_obs_events) == pd.core.frame.DataFrame:
    ranked_obs_events.reset_index(inplace=True)
    ranked_obs_events.drop(columns='level_0', inplace=True)
    for row_num in range(int(len(ranked_obs_events.values) / 3)):
        row_num = row_num * 3
        for col_num, _ in enumerate(ranked_obs_events.values[row_num]):
            if ranked_obs_events.values[row_num][1] == any(obs_time):
                if type(ranked_obs_events.values[row_num][col_num]) != datetime.date and type(ranked_obs_events.values[row_num][col_num]) != datetime.time:
                    worksheet3.write(row_num + 1, col_num, ranked_obs_events.values[row_num][col_num], cell_format)
                    worksheet3.write(row_num + 2, col_num, ranked_obs_events.values[row_num + 1][col_num], cell_format)
                    worksheet3.write(row_num + 3, col_num, ranked_obs_events.values[row_num + 2][col_num], cell_format)
                else:
                    
                    ranked_obs_events.loc[row_num, 'date'] = ranked_obs_events.loc[row_num, 'date'].isoformat()
                    ranked_obs_events.loc[row_num + 1, 'date'] = ranked_obs_events.loc[row_num + 1, 'date'].isoformat()
                    ranked_obs_events.loc[row_num + 2, 'date'] = ranked_obs_events.loc[row_num + 2, 'date'].isoformat()
                
                    ranked_obs_events.loc[row_num, 'time'] = ranked_obs_events.loc[row_num, 'time'].isoformat()
                    ranked_obs_events.loc[row_num + 1, 'time'] = ranked_obs_events.loc[row_num + 1, 'time'].isoformat()
                    ranked_obs_events.loc[row_num + 2, 'time'] = ranked_obs_events.loc[row_num + 2, 'time'].isoformat()
    
                    worksheet3.write(row_num + 1, col_num, ranked_obs_events.values[row_num][col_num], cell_format)
                    worksheet3.write(row_num + 2, col_num, ranked_obs_events.values[row_num + 1][col_num], cell_format)
                    worksheet3.write(row_num + 3, col_num, ranked_obs_events.values[row_num + 2][col_num], cell_format)
    
            else:
                if type(ranked_obs_events.values[row_num][col_num]) != datetime.date and type(ranked_obs_events.values[row_num][col_num]) != datetime.time:
                    worksheet3.write(row_num + 1, col_num, ranked_obs_events.values[row_num][col_num])
                    worksheet3.write(row_num + 2, col_num, ranked_obs_events.values[row_num + 1][col_num])
                    worksheet3.write(row_num + 3, col_num, ranked_obs_events.values[row_num + 2][col_num])
                else:
                    
                    ranked_obs_events.loc[row_num, 'date'] = ranked_obs_events.loc[row_num, 'date'].isoformat()
                    ranked_obs_events.loc[row_num + 1, 'date'] = ranked_obs_events.loc[row_num + 1, 'date'].isoformat()
                    ranked_obs_events.loc[row_num + 2, 'date'] = ranked_obs_events.loc[row_num + 2, 'date'].isoformat()
                
                    ranked_obs_events.loc[row_num, 'time'] = ranked_obs_events.loc[row_num, 'time'].isoformat()
                    ranked_obs_events.loc[row_num + 1, 'time'] = ranked_obs_events.loc[row_num + 1, 'time'].isoformat()
                    ranked_obs_events.loc[row_num + 2, 'time'] = ranked_obs_events.loc[row_num + 2, 'time'].isoformat()
    
                    worksheet3.write(row_num + 1, col_num, ranked_obs_events.values[row_num][col_num])
                    worksheet3.write(row_num + 2, col_num, ranked_obs_events.values[row_num + 1][col_num])
                    worksheet3.write(row_num + 3, col_num, ranked_obs_events.values[row_num + 2][col_num])
else:
    pass

# Close the workbook
workbook.close()

