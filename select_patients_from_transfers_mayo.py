#This code reads in the transfers file created before and subselects specific groups based on a range of criteria
# can subselect based on patient characteristics eg demographics or locations that the specific patient pathway intersects with
# can also select by certain dates to get all the transfers on specific dates.


import numpy as np
import itertools
import functools

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import pandas as pd
#from mpl_toolkits import mplot3d
#from matplotlib import cm
#from matplotlib import colors
import networkx as nx
#from collections import Counter
#from itertools import chain
#from collections import defaultdict
from datetime import datetime
import datetime as dt

#this function allows separation of the time and date stamp from EPIC
def get_separate_date_time(datetimeentry):
    print(datetimeentry)
    if type(datetimeentry) == float:
        return datetime.max
    else:
        #this returns the date in a format where the hours and days can be accessed eg d.year or d.minute
        separate_date_time = datetime.strptime(datetimeentry,"%Y-%m-%d %H:%M:%S")
        return separate_date_time

#change it to a datetime object in python
def make_into_time(dt_entry):
    dt = datetime.strptime(dt_entry, "%Y-%m-%d")
    return dt


# return only the day part of the time stamp and convert to datetime object
def get_transfer_day(date):
    strdate = str(date)
    fmt = "%Y-%m-%d"
    try:
        dd = datetime.strptime(strdate, fmt)
    except ValueError as v:
        ulr = len(v.args[0].partition('unconverted data remains: ')[2])
        if ulr:
            d = datetime.strptime(strdate[:-ulr], fmt)
            date_day = d.date()
            dd = date_day.strftime('%Y-%m-%d')
        else:
            raise v
    return dd

#find the previous day - this is useful for making moving windows in time to then calculate networks for the time period in the window
def get_previous_day(date):
    date_in_date_format = get_transfer_day(date)
    prev_day = datetime.date(date_in_date_format) - dt.timedelta(1)
    prev_day = prev_day.strftime('%Y-%m-%d')
    return prev_day

#find the next day - useful for 3 day windows
def get_next_day(date_n):
    date_in_date_format = get_transfer_day(date_n)
    next_day = datetime.date(date_in_date_format) + dt.timedelta(1)
    next_day = next_day.strftime('%Y-%m-%d')
    return next_day

#read in the transfer file - can be modified (using strain addition)
alltransfers = pd.read_csv("transfer_strain.csv")

#select all adult patients 16 and above if the analysis looks only at adults
adult_transfers= alltransfers.loc[alltransfers['age']>16]
#select the appropriate columns only
adult_transfers = adult_transfers[['ptid','transfer_dt','dt_adm', 'dt_dis', 'spec', 'age', 'asa', 'breach_percentage', 'Date', 'bedsfree', 'strain', 'from', 'to']]
adult_transfers.to_csv('adult_transfers.csv')

#Selection based on strain paraemter: here we select transfers on specific dates with low breach percentage ie days where A&E was very full
transfers_lowed = adult_transfers[adult_transfers['breach_percentage'] < 0.6955]
#transfers_lowed.to_csv('transfers_lowedpercentage.csv')

#The next section finds the transfers that occurred on the day before and the day after a full day in the ER - we assumed that this would reflect the system under strain.
#find the days with low ED percentage first
print(transfers_lowed['transfer_dt'].map(get_transfer_day))
transfers_lowed['day_of_transfer'] = transfers_lowed['transfer_dt'].map(get_transfer_day)
low_ed_perc_dates = transfers_lowed['day_of_transfer'].unique()
low_ed_prev_day = []
low_ed_next_day = []
#find the day before and after
for i in low_ed_perc_dates:
    prev_day = get_previous_day(i)
    next_day = get_next_day(i)
    low_ed_prev_day.append(prev_day)
    low_ed_next_day.append(next_day)

all_dates_low_ed = set(low_ed_prev_day + list(low_ed_perc_dates) + low_ed_next_day)
adult_transfers['day_of_transfer'] = adult_transfers['transfer_dt'].map(get_transfer_day)
transfers_around_low_ed_ind = adult_transfers[adult_transfers['day_of_transfer'].isin(all_dates_low_ed)]
transfers_around_low_ed_ind.to_csv('transfers_around_low_ed_perc.csv') # these are the transfers that occur on the day before, on the day and the day after a full ER.

#The next section does the same for days where the ER was quite empty
#select transfers on specific dates with low breach percentage ie days where A&E was very empty
transfers_highed = adult_transfers[adult_transfers['breach_percentage'] >0.9685]
transfers_lowed.to_csv('transfers_lowedpercentage.csv')

#select patients for the day before, the day of and the day after an empty A&E
#find the days with low ED percentage
transfers_highed['day_of_transfer'] = transfers_highed['transfer_dt'].map(get_transfer_day)
high_ed_perc_dates = transfers_highed['day_of_transfer'].unique()
high_ed_prev_day = []
high_ed_next_day = []
for i in high_ed_perc_dates:
    prev_day = get_previous_day(i)
    next_day = get_next_day(i)
    high_ed_prev_day.append(prev_day)
    high_ed_next_day.append(next_day)

all_dates_high_ed = set(high_ed_prev_day + list(high_ed_perc_dates) + high_ed_next_day)
transfers_around_high_ed_ind = adult_transfers[adult_transfers['day_of_transfer'].isin(all_dates_high_ed)]
transfers_around_high_ed_ind.to_csv('transfers_around_high_ed_perc.csv')



# The next section selects specific subgroups based on demographics:

#Select all the patients who at some point in their stay were in ICU or neurocritical care
#wards = {'ADD GENERAL ICU', 'ADD NEURO ICU', }
#icu_patient_ids = set(adult_transfers.loc[adult_transfers['from'].isin(wards)]['ptid'].unique())
#icu_patient_records = adult_transfers.loc[adult_transfers['ptid'].isin(icu_patient_ids)]
#icu_patient_records.to_csv('transfers_all_icu.csv', header=True, index=False)

#all adult HDU or ICU patients - so including the other high dependency areas
wards = {'ADD GENERAL ICU', 'ADD NEURO ICU', 'ADD D4 IDA UNIT', 'ADD CORONARY CARE UNIT', 'ADD TRANSPLANT HDU'}
icu_patient_ids = set(adult_transfers.loc[adult_transfers['from'].isin(wards)]['ptid'].unique())
icu_patient_records = adult_transfers.loc[adult_transfers['ptid'].isin(icu_patient_ids)]
icu_patient_records.to_csv('transfers_icu.csv', header=True, index=False)


#select specific transfers from the busy ER list
wards = {'ADD GENERAL ICU', 'ADD NEURO ICU', 'ADD D4 IDA UNIT', 'ADD CORONARY CARE UNIT', 'ADD TRANSPLANT HDU'}
icu_patient_ids = set(transfers_around_low_ed_ind.loc[transfers_around_low_ed_ind['from'].isin(wards)]['ptid'].unique())
icu_patient_records = transfers_around_low_ed_ind.loc[transfers_around_low_ed_ind['ptid'].isin(icu_patient_ids)]
icu_patient_records.to_csv('transfers_lowed_hdu_2309.csv', header=True, index=False)

#select specific transfers fromt he calm ER list
wards = {'ADD GENERAL ICU', 'ADD NEURO ICU', 'ADD D4 IDA UNIT', 'ADD CORONARY CARE UNIT', 'ADD TRANSPLANT HDU'}
icu_patient_ids = set(transfers_around_high_ed_ind.loc[transfers_around_high_ed_ind['from'].isin(wards)]['ptid'].unique())
icu_patient_records = transfers_around_high_ed_ind.loc[transfers_around_high_ed_ind['ptid'].isin(icu_patient_ids)]
icu_patient_records.to_csv('transfers_high_hdu_2309.csv', header=True, index=False)



# select patients who were admitted under a specific specialty
specialities = {'Trauma', 'Orthopaedics'}
t_o_patient_ids = set(alltransfers.loc[alltransfers['spec'].isin(specialities)]['ptid'].unique())
t_o_patient_records = alltransfers.loc[alltransfers['ptid'].isin(t_o_patient_ids)]


trauma_spec={'Trauma'}
trauma_ids = set(alltransfers.loc[alltransfers['spec'].isin(trauma_spec)]['ptid'].unique())
trauma_records = alltransfers.loc[alltransfers['ptid'].isin(trauma_ids)]
trauma_adult_records = trauma_records.loc[trauma_records['age'] >16]
trauma_adult_records.to_csv('transfers_trauma_adult.csv', header = True, index = False)

#select all elderly trauma patients
#age_old = {'80','81','82','83','84','85','86','87','88','89','90','91','92','93','94','95'}
#t_o_old_patient_ids = set(t_o_patient_records.loc[t_o_patient_records['age'].isin(age_old)]['ptid'].unique())
#t_o_old_patient_records = t_o_patient_records.loc[t_o_patient_records['ptid'].isin(t_o_old_patient_ids)]
#t_o_old_patient_records.to_csv('transfers_old_tando.csv', header=True, index=False)

#select all ASA 3 and 4 adult patients to find out what happens to higher risk patients
#asacategory= {'3','4'}
#asa34_adult_ids = set(adult_transfers.loc[adult_transfers['asa'].isin(asacategory)]['ptid'].unique())
#asa34_adult_records = adult_transfers.loc[adult_transfers['ptid'].isin(asa34_adult_ids)]
#asa34_adult_records.to_csv('transfers_adult_asa34.csv', header=True, index=False)

#select all paediatric patients
#paeds_transfers = alltransfers.loc[alltransfers['age'] < 16]
#paeds_transfers.to_csv('transfers_paeds_all.csv', header=True, index=False)


