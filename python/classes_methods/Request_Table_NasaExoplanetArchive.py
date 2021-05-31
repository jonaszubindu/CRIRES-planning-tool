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
explicitly with spaces similar to the other constraints. The logic symbol < and > are inclusive(>=, <=). 
Like this the module will find the details of the constraint. Make sure that the defined constraints are
also columns of the table you request. Otherwise the constraints are not applicable.
The script creates a URL to request for the exoplanet table and filters the initial table after the constraints.
It stores a .csv file of that table that can be imported to Transit_List.py via csv_file_import.py
    
"""

import pandas as pd
import re
import os
import wget

def from_sexagesimal_to_deg(args):
    args = list(filter(None, re.split(r"d|m|s",args)))
    args = [float(arg) for arg in args]
    deg = args[0] + args[1]/60 + args[2]/3600
    return deg

path = os.getcwd()
path = path[:-15]
# print(path)

with open(path + 'Nasa_Archive_Selection.txt', 'r') as file:
    # Loading file with contstraints for planets
    f_contents = file.readlines()

f_constraints = []
f_column = []

""" 
    Reads in the column names and constraints given from the Nasa_Archive_Selection.txt file
    and processes them 

"""

for line in f_contents:
    
    try:
        if 'CONSTRAINT' in line:
            f_constraints.append(line)
        elif 'COLUMN' in line:
            f_column.append(line)
    except SyntaxError:
        raise Warning('Illegal lines used in Nasa_Archive_Selection - file')


f_req_columns = [line.split("COLUMN")[1] for line in f_column]
f_req_columns = [line[1:line.find(':')] for line in f_req_columns]
# f_URL_req = ''
# for column in f_req_columns:
#     f_URL_req += '{},'.format(column)
# f_URL_req = f_contents[0].split(' ')[0] + f_URL_req + f_contents[0].split(' ')[1]


#####################################################################################################

""" Excerpt from Fabio Lesjak's planning tool """

properties = ['hostname', 'pl_letter', 'pl_name', 'pl_orbper', 'pl_orbsmax', 
              'pl_radj', 'pl_bmassj', 'ra', 'dec', 'pl_orbincl', 'pl_orbeccen', 
              'pl_orbpererr1', 'pl_orbpererr2', 'sy_vmag', 'sy_hmag', 'sy_jmag', 'sy_kmag', 
              'st_teff', 'st_rad', 'st_mass', 'pl_eqt', 'pl_trandep', 'pl_trandur', 
              'pl_tranmid', 'pl_tranmiderr1', 'pl_tranmiderr2']
    
    
urlRoot = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query="
select = "select+"
for p in properties:
    select = ''.join([select, ',', p])
select = select[:7] + select[8:]
table = "+from+pscomppars"
outformat = "&format=csv"

url = ''.join((urlRoot, select, table, outformat))

# f_URL_req = url
path = path + 'csv_files/'
filename = 'PlanetAll.csv'
if os.path.exists(path+filename):
    os.remove(path+filename) # if exist, remove it directly

wget.download(url, out=path+filename)


#####################################################################################################


# Remove duplicates
res = []
for i in f_req_columns:
    if i not in res:
        res.append(i)
f_req_columns = res

cons = []
for keyword in f_req_columns:
    for constrain in f_constraints:
        if keyword+' ' in constrain:
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

"""
    Actually calling the Nasa Exoplanet Archive is done with the pandas function pandas.read_csv(filepath_or_buffer), the usage is as follows:
    
        filepath_or_buffer : str, path object or file-like object
        Any valid string path is acceptable. The string could be a URL. Valid
        URL schemes include http, ftp, s3, and file. For file URLs, a host is
        expected. A local file could be: file://localhost/path/to/table.csv.

"""
# print(os.getcwd())
df_Nasa_Archive = pd.read_csv(path + 'PlanetAll.csv') #(f_URL_req.split('\n')[0])

""" Converting sexagismal values stored as str to floats that can be constrained """

new_keys = []
degs = []
numb = []
copy_keys = keys.copy()
for key in copy_keys:
    # print(key)
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
        # print(new_key,key)
        if key in new_key:
            keys[i] = new_key

""" Filters the data of df_Nasa_Archive """ 

if cons != []:
    
    filt = " & ".join(["(df_Nasa_Archive['{0}'] {1}= {2})".format(key, cond, value) for key, cond, value in zip(keys,conds,values)])
    print(filt)
    filt = eval(filt)

    try:
        df_Nasa_Archive_filtered = df_Nasa_Archive[filt]
    except Exception:
        print('could not constrain table, something went wrong')
else:
    df_Nasa_Archive_filtered = df_Nasa_Archive

""" 
    Stores the data to a csv file. The candidate data can be reviewed here already but the data 
    in the csv file is not important for the next steps, only the names of the candidates 
"""

file_name = input('Write name to store file: [PlanetList.csv]')
if file_name != '':
    file_name = file_name
else:
    file_name = 'PlanetList.csv'

# path = path + 'csv_files/'
df_Nasa_Archive_filtered.to_csv(path + file_name)









