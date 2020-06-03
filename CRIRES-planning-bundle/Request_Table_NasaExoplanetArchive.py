#!/usr/bin/env python

"""
    Request Confirmed Exoplanets Table from Nasa Exoplanet Archive

    This script opens a file with constraints and columns that should constrain the Nasa exoplanets archive data.
    The script contains two important information:

    which columns do you want to import in your Exoplanet table

    and

    with which contraints should the table be filtered.

    The script looks automatically for constraints and columns in a file called Nasa_Archive_Selection.txt.
    It is important that columns are defined as COLUMN and constraints as CONSTRAINT for the script to find them.
    Please do not add any special characters to a column or constraint. Write the constraint in the format constraint < value
    explicitly with spaces, like this the module will find the details of the constraint. Make sure that the defined constraints are
    also columns of the table you request. Otherwise the constraints are not applicable.
    The script creates a URL to request for the exoplanet table and filters the initial table after the constraints.
    It stores a .csv file of that table that can be imported to Transit_List.py via csv_file_import.py
"""

import pandas as pd
import re

def from_sexagesimal_to_deg(args):
    args = list(filter(None, re.split(r"d|m|s",args)))
    args = [float(arg) for arg in args]
    deg = args[0] + args[1]/60 + args[2]/3600
    return deg

with open('Nasa_Archive_Selection.txt', 'r') as file:

    f_contents = file.readlines()

f_constraints = []
f_column = []
for line in f_contents:
    if 'CONSTRAINT' in line:
        f_constraints.append(line)
    elif 'COLUMN' in line:
        f_column.append(line)

f_req_columns = [line.split("COLUMN")[1] for line in f_column]
f_req_columns = [line[1:line.find(':')] for line in f_req_columns]
f_URL_req = ''
for column in f_req_columns:
    f_URL_req += '{},'.format(column)
f_URL_req = f_contents[0].split(' ')[0] + f_URL_req + f_contents[0].split(' ')[1]

# Remove duplicates
res = []
for i in f_req_columns:
    if i not in res:
        res.append(i)
f_req_columns = res

cons = []
for keyword in f_req_columns:
    for constrain in f_constraints:
        if ' '+keyword+' ' in constrain:
            cons.append(keyword + (constrain.split(keyword,1)[1]).split('\n')[0])

for n in range(len(cons)):
    try:
        cons[n] = cons[n].split('(')[1]
    except IndexError:
        pass
    cons[n] = cons[n].split(')')[0]
keys, conds, values = [[] for _ in range(3)]
for con in cons:
    key, cond, value = con.split(' ')
    keys.append(key), conds.append(cond), values.append(float(value))

df_Nasa_Archive = pd.read_csv(f_URL_req.split('\n')[0])

new_keys = []
degs = []
numb = []
copy_keys = keys.copy()
for key in copy_keys:
    #print(key)
    try:
        if type(df_Nasa_Archive[key][0]) == str and re.findall(r"[dms]", df_Nasa_Archive[key][0]) != []:
            new_key = key + '_deg'
            new_keys.append(new_key)
            for elem in df_Nasa_Archive[key]:
                degs.append(from_sexagesimal_to_deg(elem))
            df_Nasa_Archive[new_key] = degs
        elif type(df_Nasa_Archive[key][0]) == str:
            for elem in df_Nasa_Archive[key]:
                numb.append(float(elem))
            df_Nasa_Archive[key] = numb
    except KeyError:
        print('{} is not in table'.format(key))
        i = keys.index(key)
        no_cond = conds.pop(i)
        no_val = values.pop(i)
        no_key = keys.pop(i)
        print('{}{}{} cannot be constrained'.format(no_key,no_cond,no_val))
    else:
        pass

for new_key in new_keys:
    for i, key in enumerate(keys):
        #print(new_key,key)
        if key in new_key:
            keys[i] = new_key

filt = " & ".join(["(df_Nasa_Archive['{0}'] {1}= {2})".format(key, cond, value) for key, cond, value in zip(keys,conds,values)])

filt = eval(filt)
try:
    df_Nasa_Archive_filtered = df_Nasa_Archive[filt]
except Exception:
    print('could not constrain table, something went wrong')
file_name = input('Write name to store file: ')
file_name = file_name + '.csv'
default_file_name = 'PlanetList.csv'
try:
    df_Nasa_Archive_filtered.to_csv(file_name)
except Exception:
    df_Nasa_Archive_filtered.to_csv(default_file_name)











