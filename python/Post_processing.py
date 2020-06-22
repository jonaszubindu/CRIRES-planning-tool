import pandas as pd
from classes_methods.Helper_fun import data_sorting_and_storing
import datetime
from classes_methods.classes import Nights
import sys
from datetimerange import DateTimeRange
from classes_methods.classes import load_Eclipses_from_file

try:
    filename = sys.argv[1]
except Exception:
    filename = 'Eclipse_events_processed_2020-06-17_365d.pkl'

d = datetime.datetime.fromisoformat(filename.split('_')[3])
Max_Delta_days = int(filename.split('_')[4][:-5])
Nights = Nights(d, Max_Delta_days)

Eclipses_List = load_Eclipses_from_file(filename, Max_Delta_days)
ranking, df_gen, df_frame = data_sorting_and_storing(
    Eclipses_List, write_to_csv=0)


df_frame_date = []
df_frame_time = []
for elem in df_frame['time']:
    df_frame_date.append(elem.value.date())
    df_frame_time.append(elem.value.time())

# df_frame.reset_index(inplace=True)
for n in range(int(len(df_frame['time'])/3-1)):
    start = df_frame['time'][3*n]
    end = df_frame['time'][(n+1)*3-1]
    if end < start:
        # print(f"we got a problem with {df_frame.loc[3*n]}-{df_frame.loc[(n+1)*3-1]}")
        raise Warning('smth went wrong in handling the times')
        
df_frame.drop(columns=['time'], inplace=True)
df_frame.insert(0, 'date', df_frame_date)
df_frame.insert(1, 'time', df_frame_time)
# df_frame.sort_values(by=df_frame['date'])  # sorts eclipses in time and date.


ranking_dates = []
date_section = []


for n in range(int(len(df_frame['date']) / 3)):
    date_section.append(df_frame[n * 3:(n + 1) * 3])
# for date_sec in date_section:
#     ranking_dates.append((len(date_sec) / 3, date_sec))


for date in Nights.date:
    date_sec = []
    for date_obj in date_section:
        # date_obj.reset_index(inplace=True)
        if date_obj['date'][0] == date.date():
            date_sec.append(date_obj)
    if date_sec == []:
        pass
    else:
        ranking_dates.append((len(date_sec), date_sec))

for date_obj in ranking_dates:  # for each date
    ran = []
    date_obj = list(date_obj)
    for ecl in date_obj[1]:
        ecl.reset_index(inplace=True)
        start = ecl['date'][0].strftime("%Y-%m-%dT") + ecl['time'][0].strftime("%H:%M:%S-0000")  # begin of eclipse
        end = ecl['date'][2].strftime("%Y-%m-%dT") + ecl['time'][2].strftime("%H:%M:%S-0000")  # end of eclipse
        ran.append((ecl, DateTimeRange(start, end)))

    for range1 in ran:
        for range2 in ran:
            if range1[1] == range2[1]:
                pass
            else:
                if range1[1].start_datetime < range1[1].end_datetime and range2[1].start_datetime < range2[1].end_datetime:
                    inter_sect = range1[1].intersection(range2[1])
                    if inter_sect != None:
                        date_obj[0] += -1 / 2
                elif range1[1].start_datetime < range1[1].end_datetime:
                    print(f"{range1[1].start_datetime} > {range1[1].end_datetime} in {range1[0][['index','date','time']]}")
                    
                elif range2[1].start_datetime < range2[1].end_datetime:
                    print(f"{range2[1].start_datetime} > {range2[1].end_datetime} in {range2[0][['index','date','time']]}")

ranking_dates.sort(key=lambda tup: tup[0])
