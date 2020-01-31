#Copyright 2020 Katharina Kohler

#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at

#    https://www.apache.org/licenses/LICENSE-2.0

#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

#this file allows the addition of system strain data to the transfers file - eg how busy the ED is or other more advanced features
#this step is not necessary for the intial creation of a network - it is mainly useful for further analysis into how the system behaves under strain or when developing temporal networks - ie a network for each day/week/month

import matplotlib
#from mpl_toolkits import mplot3d
#from matplotlib import cm
#from matplotlib import colors
import networkx as nx
#from collections import Counter
#from itertools import chain
#from collections import defaultdict
from datetime import datetime
import pandas as pd

def get_date_number(dt):
    return dt.year * 10000 + dt.month * 100 + dt.day

#convert a date to the datetime format needed
def get_transfer_day(date):
    strdate = str(date)
    fmt = "%Y-%m-%d"
    try:
        d = datetime.strptime(strdate, fmt)
    except ValueError as v:
        ulr = len(v.args[0].partition('unconverted data remains: ')[2])
        if ulr:
            d = datetime.strptime(strdate[:-ulr], fmt)
        else:
            raise v
    return d

#here the example is to add on bed state ie the percentage fo beds occupied as a measure for system strain.

max_beds = 1155 # maximal number of beds +1 to not end up with infinite strain

def get_free_beds(beds_occupied):
    return max_beds - beds_occupied

#reead int he transfers file created previously
transfer_data = pd.read_csv("transfers.csv")
transfer_data['date'] = pd.to_datetime(transfer_data['transfer_dt'], format="%Y-%m-%d %H:%M")
transfer_data['date_as_number'] = transfer_data['date'].map(get_date_number)

#read in the ER performance data - the details depend on wheich measures are used for the specific example
ed_performance = pd.read_csv("ed_performance_all.csv")
# need transfer date only in a separate column
ed_performance['date'] = pd.to_datetime(ed_performance['day'], format='%d/%m/%Y')
ed_performance['date_number'] = ed_performance['date'].map(get_date_number)
ed_performance.drop(['date'], axis=1, inplace=True)
ed_performance.set_index('date_number', drop=True, inplace=True)
transfer_data_ed = transfer_data.join(ed_performance, on='date_as_number', how='left')


#add on bedstate information to give further details ont he system pressure

bedstate_info = pd.read_csv("all_beds_info.csv")
bedstate_info['date'] = pd.to_datetime(bedstate_info['Date'], format='%Y-%m-%d')
bedstate_info['date_number'] = bedstate_info['date'].map(get_date_number)
bedstate_info.drop(['date'], axis = 1, inplace = True)
bedstate_info.set_index('date_number', drop = True, inplace = True)
transfer_data_beded= transfer_data_ed.join(bedstate_info, on = 'date_as_number', how = 'left')

transfer_data_beded['bedsfree'] = transfer_data_beded['Total Occupied'].map(get_free_beds)
#the following measure was a try at developing a system strain measure - it was not completely successful and is a part of the project still under development.
transfer_data_beded['strain'] = 1.0/transfer_data_beded.bedsfree * transfer_data_beded.breach_percentage *100
transfer_strain = transfer_data_beded.drop(['date_as_number','date', 'day'], axis=1)


transfer_strain.to_csv('transfer_strain.csv', header=True, index=False)