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

#First part of the code:
# reading in EPIC data (locations with time stamps from a range of source csv files.
# collating the locations into a timeline for each patient and then creating the transfers file
# transfers file contains all the transfers for each patient in the data set with demographics attached.

#import the necessary libraries
import pandas as pd
from datetime import datetime
import numpy as np
from collections import deque, namedtuple

#separate out time and date from EPIC time stamp
def get_separate_date_time(datetimeentry):
    print(datetimeentry)
    if type(datetimeentry) == float:
        return datetime.max
    else:
        #this returns the date in a format where the hours and days can be accessed eg d.year or d.minute
        separate_date_time = datetime.strptime(datetimeentry,"%Y-%m-%d %H:%M:%S")
        return separate_date_time


#function to remove "ADD" in our cae representing Addenbrookes hospital. This can be adjusted to remove whichever pieces of location name are not needed
def simplify_theatre_entries(df: pd.DataFrame):
    df_nonan = df[df['adt_department_name'].notna()]
    theatre_rows = df_nonan[df_nonan['adt_department_name'].str.contains('THEATRE')]
    df_nonan.loc[theatre_rows.index, 'adt_department_name'] = 'THEATRE'
    return df_nonan

#reading in data: Here the files are called ADM_POINT_aug.csv, SURGERIES_aug.csv, ENC_POINT_aug.csv, adm_info (demographics)
#each one contains information from a different part of the EPIC database
#admpoint contains the transfers of all the patients between wards
admpoint = pd.read_csv("ADM_POINT_aug.csv")
# only keep the columns actually needed
admpoint = admpoint[['STUDY_SUBJECT_DIGEST', 'in_dttm', 'out_dttm', 'adt_department_name']]
#add on a column indicating the origin of an entry in this case adm point for later when joining on information from other files
s_length = len(admpoint['in_dttm'])
admpoint['data_origin'] = np.repeat('adm', s_length, axis=0)
# Rename the 'STUDY_SUBJECT_DIGEST' column to 'ptid' to signify patient ID. This is needed for later joining the files
admpoint.rename(index=str, columns={'STUDY_SUBJECT_DIGEST': 'ptid'}, inplace=True)

#convert the date and time to a datetime object
admpoint['in_dttm'] = pd.to_datetime(admpoint['in_dttm'])
admpoint['out_dttm'] = pd.to_datetime(admpoint['out_dttm'])
#remove excessive information and nonsense entries
amdpoint = simplify_theatre_entries(admpoint)

#surgeriesinfo contains details about the surgery
surgeriesinfo = pd.read_csv("SURGERIES_aug.csv")
#encinfo contains information about scans, investigations and clinic visits
encinfo = pd.read_csv("ENC_POINT_aug.csv")

#mark the point at which all the files have been read in
print('reading in done')

#change the format of the newly read in files so that they can be joined to admpoint

surg_df = surgeriesinfo[['STUDY_SUBJECT_DIGEST', 'case_start', 'case_end', 'prov_name']]
s_length = len(surg_df['case_start']) #length of series that needs to be added into the new columns

surg_df['data_origin'] = np.repeat ('surg', s_length, axis = 0)
surg_df.rename(index=str, columns={'STUDY_SUBJECT_DIGEST': 'ptid'}, inplace=True)
surg_df.rename(index=str, columns={'case_start': 'in_dttm'}, inplace=True)
surg_df.rename(index=str, columns={'case_end': 'out_dttm'}, inplace=True)
surg_df.rename(index=str, columns={'prov_name': 'adt_department_name'}, inplace = True)
surg_df['in_dttm'] = pd.to_datetime(surg_df['in_dttm'], dayfirst=True)
surg_df['out_dttm'] = pd.to_datetime(surg_df['out_dttm'], dayfirst=True)
surg_df = simplify_theatre_entries(surg_df)

#change the formatting of encinfo to be able to join onto admpoint
enc_df = encinfo[['STUDY_SUBJECT_DIGEST', 'at_time','enctype', 'dep_name']]
#replace the empty entries in the department name with the encounter type eg operation.
empty_indices = enc_df['dep_name'] == ''
enc_df.loc[empty_indices, 'dep_name'] = enc_df.loc[empty_indices, 'enctype']
enc_df = enc_df[['STUDY_SUBJECT_DIGEST', 'at_time', 'dep_name']]
s_length = len(enc_df['at_time'])
enc_df['out_dttm'] = np.repeat(1,s_length, axis =0)
enc_df['data_origin'] = np.repeat ('enc', s_length, axis = 0)
enc_df.rename(index=str, columns={'STUDY_SUBJECT_DIGEST': 'ptid'}, inplace=True)
enc_df.rename(index=str, columns={'at_time': 'in_dttm'}, inplace=True)
enc_df.rename(index=str, columns={'dep_name': 'adt_department_name'}, inplace = True)
enc_df['in_dttm'] = pd.to_datetime(enc_df['in_dttm'])
enc_df['out_dttm'] = enc_df['in_dttm']
enc_df = simplify_theatre_entries(enc_df)

#make dataframe containing all the information by joining the three ddataframes together, then sort by patient id
#convert in_dttm to a datetime object to be able to sort by it later
full_info = pd.concat([admpoint, enc_df, surg_df], ignore_index=True)
# full_info['in_dttm'] = pd.to_datetime(full_info['in_dttm'])
# full_info['out_dttm'] = pd.to_datetime(full_info['out_dttm'])
full_info = full_info.sort_values(by = ['ptid', 'in_dttm'], ascending = [True, True])
full_info.reset_index(drop=True, inplace=True)

#add on the information from the other files that is needed in addition to the transfer data in admpoint
#adminfo contains demographic data for the patients - eg age and admission specialty, ASA grading...
adminfo = pd.read_csv("ADM_INFO_aug.csv")
#pick the columns in the secondary files that are actually needed.
adminfo = adminfo[['adm_hosp', 'dis_hosp', 'specialty', 'admAge', 'STUDY_SUBJECT_DIGEST']]
# Set the index of the adminfo dataframe to the value we want to join to.
adminfo.set_index('STUDY_SUBJECT_DIGEST', drop=True, inplace=True)
# Join the columns in adminfo onto the admpoint dataframe based on patient ID.
full_info = full_info.join(adminfo, on='ptid', how='left')
print('joining on adminfo')


#joining on the information from the surgeries dataframe
surg_extra = surgeriesinfo[['asa_rating_c', 'STUDY_SUBJECT_DIGEST']]
surg_extra.set_index('STUDY_SUBJECT_DIGEST', drop=True, inplace=True)
full_info = full_info.join(surg_extra, on='ptid', how='left')

full_info = full_info[full_info.adt_department_name != ''] # removes all lines where there is no location given
full_info = full_info[full_info.adt_department_name != ' '] # removes all lines where there is no location given

#remove any lines where the location is not at least 3 length - as this is not a real location code
full_info['adt_department_name'] = full_info['adt_department_name'].astype('str')
mask = (full_info['adt_department_name'].str.len() >= 3)
full_info = full_info.loc[mask]


full_info.to_csv('full_info_all.csv', header=True, index=False)

print('full_info created')
## Remove event types we don't care about.
## Event types are: Admission, Transfer Out, Transfer In, Census, Patient Update, Discharge, Hospital Outpatient
admpoint = admpoint[admpoint.evttype != 'Census'] # removes all census lines
admpoint = admpoint[admpoint.evttype != 'Patient Update'] # removes all patient update lines


#at this point the columns of the dataframe should have the following headers: 'ptid', 'transfer_dt', 'from', 'to', 'dt_adm', 'dt_dis', 'spec', 'age', 'asa'
##!!! starting the transfers code

#making a named tuple to store the current transfer tdata
location = namedtuple("location", ["name","dt_in","dt_out","dt_adm","dt_dis","spec","age","asa"])

#now starting the actual transfers code that creates the pathway for each patient

def get_transfers_out(ptid, location_stack, current_dt):
    # In the case of nested locations, this function generates transfers out of those locations to the
    # top-level location.

    # First check if the patient is still in the most recent location at the specified time.
    # If this is true then there are no transfers happening.
    if (len(location_stack) <= 1) or (location_stack[-1].dt_out >= current_dt):
        return []

    transfer_list = []
    #print(len(location_stack))
    #print(location_stack[-1].dt_out)
    #print(current_dt)
    # Pop locations off the stack until we find one that the patient is still in at the current time,
    # or until there are no more locations left.
    # We always leave at least one location on the stack because this function does not handle patient discharge.
    while (len(location_stack) > 1) and (location_stack[-1].dt_out < current_dt):
        loc = location_stack.pop()
        next_loc = location_stack[-1]
        if loc.name != next_loc.name:
            transfer_list.append({'ptid': ptid, 'transfer_dt': loc.dt_out, 'from': loc.name, 'to': next_loc.name, 'dt_adm': loc.dt_adm, 'dt_dis': loc.dt_dis, 'spec': loc.spec, 'age': loc.age, 'asa': loc.asa})

    return transfer_list


def clean_patient_data(patient_data: pd.DataFrame):
    # Remove locations: 'ADD FLUORO', 'nan' - these occur frequently but don't represent transfers as fluoro is done at the patient bedside
    bad_locations = {'ADD FLUORO', 'nan', 'CUH EXT FILM'}
    bad_location_data = patient_data[patient_data['adt_department_name'].isin(bad_locations)]

    if len(bad_location_data) > 0:
        good_data = patient_data.drop(bad_location_data.index, axis=0).reset_index(drop=True)
    else:
        good_data = patient_data.reset_index(drop=True)
    # Replicated locations should be collapsed into a single row with in_dttm as the first in_dttm and out_dttm as
    # the last out_dttm.
    current_loc_index = None
    current_loc = None
    current_loc_dt_out = None
    indices_to_remove = []
    surg_theatre_index = None

    def is_separate_visit(t1: pd.Timestamp, t2: pd.Timestamp):
        # Times different by more one hour or more means separate visits (this was a value judgement, it could have been set at 10 minutes depending on the situation.
        delta_seconds = (t2 - t1).total_seconds()
        return delta_seconds >= 3600

    for i, row in good_data.iterrows():
        loc = row['adt_department_name']
        dt_in = row['in_dttm']
        dt_out = row['out_dttm']

        if (loc != current_loc) or is_separate_visit(current_loc_dt_out, dt_in):
            if (row['data_origin'] == 'surg') and ('THEATRE' in loc):
                # Sometimes a theatre entry from the surgical data is wrong.
                # One way to detect these is if the patient goes back to A&E after the theatre entry.
                # If that happens there is always a correct theatre visit later from the adm source, so remove the bad one.
                # this was detected by looking through the data to find clinically improbable pathways that occur due to incorrect assignment by staff on EPIC.
                surg_theatre_index = i
            elif loc == 'POST-DISCHARGE':
                # If the patient is discharged they are allowed to come to A&E again as a representation

                surg_theatre_index = None
            elif (surg_theatre_index is not None) and (loc == 'ADD EMERGENCY DEPT'):
                # The patient is now in A&E but they were previously in theatre, according to the surgical data.
                # Delete the previous theatre entry.
                indices_to_remove.append(surg_theatre_index)
                surg_theatre_index = None

            current_loc = loc
            current_loc_index = i
            current_loc_dt_out = dt_out
        else:
            # We're already in this location. Remove this duplicate row, and update the out time of the
            # first row for this location to the out time from this row.
            indices_to_remove.append(i)
            current_loc_dt_out = dt_out
            good_data.loc[current_loc_index, 'out_dttm'] = dt_out

    clean_data = good_data.drop(indices_to_remove, axis=0)

    if len(clean_data) == 0:
        ptid = patient_data['ptid'].iloc[0]
        print("Writing pt_%s.csv for patient with no cleaned data." % ptid)
        patient_data.to_csv("pt_%s.csv" % ptid)

    return clean_data

def is_bad_patient(patient_data: pd.DataFrame):
    # If they go from theatre or a ward to A&E, the entry is not good as this is not something that happens clinically.
    # #Occasionally these entries occur due to a patient not being assigned to the correct location
    shouldnt_go_to_ae = False

    for i, row in patient_data.iterrows():
        loc = row['adt_department_name']

        if ('THEATRE' in loc) or ('WARD' in loc):
            shouldnt_go_to_ae = True
        elif loc == 'POST-DISCHARGE':
            shouldnt_go_to_ae = False
        elif (loc == 'ADD EMERGENCY DEPT') and shouldnt_go_to_ae:
            return True

    return False


def get_patient_transfers(ptid, patient_data):
    patient_data = clean_patient_data(patient_data)
    if len(patient_data) == 0:
        return None, None

    if is_bad_patient(patient_data):
        return None, patient_data

    # The stack will contain the previous (location, entry_time, exit_time) tuples.
    location_stack = []
    transfer_list = []
    dt_adm = patient_data['adm_hosp'].iloc[0]

    for i, row in patient_data.iterrows():
        loc = row['adt_department_name']
        dt_in = row['in_dttm']
        dt_out = row['out_dttm']
        dt_dis = row['dis_hosp']
        spec = row['specialty']
        age = row['admAge']
        asa = row['asa_rating_c']

        new_loc = location(loc, dt_in, dt_out, dt_adm, dt_dis, spec, age, asa)

        # If the stack is empty then this is the first location for the patient and we won't
        # be adding any transfers yet.
        if len(location_stack) == 0:
            location_stack.append(new_loc)
        else:
            # If there are nested locations that the patient is no longer in, we need to transfer out of them.
            transfer_list += get_transfers_out(ptid, location_stack, new_loc.dt_in)

            current_loc = location_stack[-1]

            if loc == 'POST-DISCHARGE':
                transfer_list.append(
                    {'ptid': ptid, 'transfer_dt': new_loc.dt_in, 'from': current_loc.name, 'to': 'discharge', 'dt_adm': new_loc.dt_adm, 'dt_dis': new_loc.dt_dis, 'spec': new_loc.spec, 'age': new_loc.age, 'asa': new_loc.asa})
                location_stack = []
            else:
                # Check if the patient entered the next location after leaving the current location.
                # If so, this is a normal transfer, and we need to pop the current location off the stack
                # so it can be replaced by the new location.
                if new_loc.dt_in >= current_loc.dt_out:
                    # This is a normal transfer.
                    # Pop the previous location off the stack because we've left it.
                    location_stack.pop()

                # The patient is now in the new location so put it on top of the stack.
                location_stack.append(new_loc)

                # Add a transfer from the previous location to the new location.
                if current_loc.name != new_loc.name:
                    transfer_list.append(
                        {'ptid': ptid, 'transfer_dt': new_loc.dt_in, 'from': current_loc.name, 'to': new_loc.name, 'dt_adm': new_loc.dt_adm, 'dt_dis': new_loc.dt_dis, 'spec': new_loc.spec, 'age': new_loc.age, 'asa': new_loc.asa})

    if len(location_stack) > 0:
        # In case we are inside a bunch of nested locations, we now need to transfer the patient out of them.
        transfer_list += get_transfers_out(ptid, location_stack, datetime.max)

        # Finally we need to discharge the patient.
        last_loc = location_stack[-1]
        transfer_list.append({'ptid': ptid, 'transfer_dt': last_loc.dt_out, 'from': last_loc.name, 'to': 'discharge','dt_adm': last_loc.dt_adm, 'dt_dis': last_loc.dt_dis, 'spec': last_loc.spec, 'age': last_loc.age, 'asa': last_loc.asa })

    # Return a DataFrame with the transfers.
    transfers = pd.DataFrame(columns=['ptid', 'transfer_dt', 'from', 'to', 'dt_adm', 'dt_dis', 'spec', 'age', 'asa'], data=transfer_list)

    return transfers, patient_data


def get_transfers(location_data: pd.DataFrame):
    num_good_patients = 0
    num_bad_patients = 0
    bad_patient_data = []

    sorted_data = location_data.sort_values(['ptid', 'in_dttm'])

    num_patients = len(sorted_data['ptid'].unique())
    i = 0

    sorted_data.to_csv("sorted_data.csv")
    groups = sorted_data.groupby('ptid')
    all_transfers = []

    for ptid, group in groups:
        patient_transfers, patient_data = get_patient_transfers(ptid, group)
        if patient_transfers is None:
            num_bad_patients += 1
            if patient_data is not None:
                bad_patient_data.append(patient_data)
        else:
            num_good_patients += 1
            all_transfers.append(patient_transfers)

        i += 1
        if (i % 100) == 0:
            print("Finished %s of %s patients. Good patients: %s, bad patients: %s." % (i, num_patients, num_good_patients, num_bad_patients))

    print("Good patients: %s. Bad patients: %s." % (num_good_patients, num_bad_patients))

    if num_bad_patients > 0:
        print("Saving bad patient data.")
        pd.concat(bad_patient_data, ignore_index=True).to_csv('bad_patient_data.csv', header=True, index=False)

    return pd.concat(all_transfers, ignore_index=True).sort_values(by=['ptid', 'transfer_dt'], axis=0).reset_index(drop=True)


all_transfers = get_transfers(full_info)

print("Rows: %s" % len(all_transfers))
all_transfers.to_csv('transfers_with_bad_dates.csv', header=True, index=False)

#set the time period between which the data was collected. All dates that appear outside of this will be disregarded.
first_date = datetime(2015, 1, 1)
last_date = datetime(2018, 6, 1)
after_last_date = all_transfers[all_transfers['transfer_dt'] > last_date]
all_transfers.drop(after_last_date.index, axis=0, inplace=True)

before_first_date = all_transfers[all_transfers['transfer_dt'] < first_date]
if len(before_first_date) > 0:
    all_transfers.drop(before_first_date.index, axis=0, inplace=True)

print("Rows after removing bad dates: %s" % len(all_transfers))

all_transfers.to_csv('transfers.csv', header=True, index=False)
print('transfers file created')
##!!! finish of creating the transfers file
