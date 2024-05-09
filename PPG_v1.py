# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 16:33:10 2024

@author: robbe
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.backends.backend_pdf as pdf
from matplotlib.backends.backend_pdf import PdfPages
import os
from collections import defaultdict
from varname import argname
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
import pickle
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import filedialog
import scipy.stats as stats


"""
constants used in functions
"""
AM_propane = 44.0# atomic weight of propane, used for whsv calculations

Chemicals = { # dictionary of the used molecules, for mass balance calculations
    'Carbon Dioxide': {'C_atoms': 1, 'H_atoms': 0, 'O_atoms': 2},
    'Hydrogen': {'C_atoms': 0, 'H_atoms': 2, 'O_atoms': 0},
    'Argon': {'C_atoms': 0, 'H_atoms': 0, 'O_atoms': 0},
    'Nitrogen': {'C_atoms': 0, 'H_atoms': 0, 'O_atoms': 0},
    'Carbon Monoxide': {'C_atoms': 1, 'H_atoms': 0, 'O_atoms': 1},
    'Methane': {'C_atoms': 1, 'H_atoms': 4, 'O_atoms': 0},
    'Ethene': {'C_atoms': 2, 'H_atoms': 4, 'O_atoms': 0},
    'Ethane': {'C_atoms': 2, 'H_atoms': 6, 'O_atoms': 0},
    'Propylene': {'C_atoms': 3, 'H_atoms': 6, 'O_atoms': 0},
    'Propane': {'C_atoms': 3, 'H_atoms': 8, 'O_atoms': 0},
    'Isobutane': {'C_atoms': 4, 'H_atoms': 10, 'O_atoms': 0},
    'n-Butane': {'C_atoms': 4, 'H_atoms': 10, 'O_atoms': 0},
    'C2-butene': {'C_atoms': 4, 'H_atoms': 8, 'O_atoms': 0}
}
#M_catalyst = 0.20 #gram of catalyst used
"""
functions Module 1
"""
def detect_temperature_ramp(dfr,a,b,c,d):
    """
    takes an excel file and is used to detect and label an temperature ramping
    this is for one ramp experiment, mostly at the beginning of an reactor run
    not for multiple ramp experiments at different flows or pressures

    Parameters
    ----------
    dfr : dataframe
        The dataframe wich contains the reactor data
    a : float
        the FIC-120 used for the ramp experiment
    b : float
        The FIC-130 used for the ramp experiment
    c : float
        The FIC-140 used for the ramp exp
    d : float
        the minimum temperature used for the ramp experiment.

    Returns
    -------
    an updated dataframe, containing a new colum( T_ramp) with either NoRamp or TRamp

    """
    
    df = dfr
    FIC_120_SET_VALUE = a
    FIC_130_SET_VALUE = b
    FIC_140_SET_VALUE = c
    # Find the starting index of the temperature ramp
    start_index = df[(df["TIC-300_PV"] == d) & (df["TIC-300_SP"].diff() > 0)].index.min()

    # Initialize a new column "T_ramp" with default value "NoRamp"
    df["T_ramp"] = "NoRamp"

    if start_index is not None:
        # Extract the portion of the dataframe starting from the detected index
        df_temp_ramp = df.loc[start_index:]

        # Detect the end of the temperature ramp
        end_index = df_temp_ramp[
            (df_temp_ramp["FIC-120"] != FIC_120_SET_VALUE) |
            (df_temp_ramp["FIC-130"] != FIC_130_SET_VALUE) |
            (df_temp_ramp["FIC-140"] != FIC_140_SET_VALUE)
        ].index.min()

        if end_index is not None:
            # Extract the portion of the dataframe up to the detected end index
            df_temp_ramp = df_temp_ramp.loc[:end_index]

            # Set the value of "T_ramp" to "TRamp" for the detected temperature ramp
            df.loc[df_temp_ramp.index, "T_ramp"] = "TRamp"

    return df

def remove_nan_keys(list_of_dicts):
    for dictionary in list_of_dicts:
        keys_to_remove = [key for key, value in dictionary.items() if str(value).lower() == 'nan']
        for key in keys_to_remove:
            del dictionary[key]
    return list_of_dicts

def Conditions_List(file_path):
    """
    generates an list with the correct conditions for each experiment

    Parameters
    ----------
    file_path : excel file
        an excel file with the conditions for each experiment.

    Returns
    list of dictionary, each dictionary has the same keys
    wich includes the name.
    the other keys are defined as the column heads of the excelfile-
    conditions_list : List of dictionaries
        each dictionary has the same keys
        wich includes the name.
        the other keys are defined as the column heads of the excelfile-.

    """
    
    # Read the Excel file into a DataFrame
    df = pd.read_excel(file_path, header=None)

    # Assuming the first row contains the experiment names
    column_names = df.iloc[0]

    # Create a list of dictionaries, where each dictionary represents an experiment with its conditions
    conditions_list = []

    # Iterate through each row and create a dictionary for each experiment
    for index, row in df.iterrows():
        # Skip the first row as it contains column names
        if index == 0:
            continue

        experiment_name = row[0]

        # Create a dictionary for the current experiment
        experiment_dict = {'experiment': experiment_name}

        # Add each condition as a key-value pair to the dictionary
        for col, value in zip(column_names[0:], row[0:]):
            experiment_dict[col] = value

        # Append the dictionary to the list
        conditions_list.append(experiment_dict)

    updated_con_list = remove_nan_keys(conditions_list)
    return updated_con_list


def modify_experiment_names(experiments):
        experiment_counts = defaultdict(int)
        
        for experiment_dict in experiments:
            experiment = experiment_dict['Experiment']
            count = experiment_counts[experiment]
            
            if count > 0:
                experiment_dict['Experiment'] = f"{experiment}_{count}"
            
            experiment_counts[experiment] += 1
        
        return experiments


def general_start_times(dft,con_list):
    """
    gives a list of dictionaries, each dictionary is for one experiment
    keys : experiment, start time (wich is a datetime object of the start date) and end time (datetime objet of end date of experiment)

    Parameters
    ----------
    dft : dataframe
        generated by function detect_temperature_ramp.
    conditions_list : list of dictionaries
        generated by function Conditions_list
        a list of dictionaries, each dictionarie contains the name (key : experiment) of an experiment and the conditions
    Returns
    -------
    experiments : list of dictionaries, each experiment has one dictionarie
        Key:
            End Time (datetime) the end (date) time of each experiment
            Experiment (str) name of experiment
            Start Time (datetime) the start of each experiment (date)

    """
    
   
    dft['Date/Time'] = pd.to_datetime(dft['Date/Time'])

    # Initialize variables to track experiment boundaries
    experiment_start = None
    experiment_type = None

    # Create a list to store the experiment information
    experiments = []

    # Iterate through the rows of the DataFrame
    for index, row in dft.iterrows():
        # Initialize a variable to check if any condition is met
        condition_met = False
        current_experiment = None

        # Iterate through the list of dictionaries
        for conditions in con_list[0:]:  # Skip the first dictionary as it is for 'name'
            # Check if all conditions are met for the current row
            conditions_met = all(
                (pd.isna(value) and pd.isna(row[key])) or (row[key] == value)
                for key, value in conditions.items() if key in row.index
            )
            if conditions_met:
                current_experiment = conditions.get('experiment', None)
                break

        if current_experiment != experiment_type:
            if experiment_type:
                # If there is an ongoing experiment, record it
                experiments.append({
                    'Experiment': experiment_type,
                    'Start Time': experiment_start,
                    'End Time': row['Date/Time']
                })
            experiment_type = current_experiment
            experiment_start = row['Date/Time']

    # Add the last experiment if there is any ongoing experiment
    if experiment_type:
        experiments.append({
            'Experiment': experiment_type,
            'Start Time': experiment_start,
            'End Time': dft['Date/Time'].iloc[-1]
        })

    renamed_experiments = modify_experiment_names(experiments)
    return renamed_experiments

def general_start_hours(dft,con_list):
    """
    same as general_start_times but than with hours instead of dates
    needed for the construction of different plots

    Parameters
    ----------
    dft : dataframe
        dataframe generated by detect_temperature_ramp

    Returns
    -------
    experiments : list of dictionaries
        each experiment has it own dictionary
        keys:
            End hours : float
            Experiment : name of the experiment
            Start hours : float

    """
    # Read the Excel file into a DataFrame
    data_mcb = dft
    data_mcb['Date/Time'] = pd.to_datetime(data_mcb['Date/Time'])
    data_mcb['Hours'] = (data_mcb['Date/Time'] - data_mcb['Date/Time'][0]).dt.total_seconds() / 3600

    # Initialize variables to track experiment boundaries
    experiment_start = None
    experiment_type = None

    # Create a list to store the experiment information
    experiments = []

    # Iterate through the rows of the DataFrame
    for index, row in data_mcb.iterrows():
        # Initialize a variable to check if any condition is met
        condition_met = False
        current_experiment = None

        # Iterate through the list of dictionaries
        for conditions in con_list[0:]:  # Skip the first dictionary as it is for 'name'
            # Check if all conditions are met for the current row
            conditions_met = all(
                (pd.isna(value) and pd.isna(row[key])) or (row[key] == value)
                for key, value in conditions.items() if key in row.index
            )
            if conditions_met:
                current_experiment = conditions.get('experiment', None)
                break

        if current_experiment != experiment_type:
            if experiment_type:
                # If there is an ongoing experiment, record it
                experiments.append({
                    'Experiment': experiment_type,
                    'Start Time': experiment_start,
                    'End Time': row['Hours']
                })
            experiment_type = current_experiment
            experiment_start = row['Hours']

    # Add the last experiment if there is any ongoing experiment
    if experiment_type:
        experiments.append({
            'Experiment': experiment_type,
            'Start Time': experiment_start,
            'End Time': data_mcb['Hours'].iloc[-1]
        })
    experiments_renamed = experiments
    return experiments_renamed

def process_injections(excelfile, experiment_list):
    """
    links the reactor data with the GC data
    it does this by linking the start of each experiment with the corresponding injection time

    
    ----------
    DFGC = dataframe
        an dataframe wich contains the GC data, both TCD and FID
    experiment_list : list of dictionaries
        for each experiment a different dictionarie, with as keys
        End Injection Number (int)
        End Minutes (int)
        End Time
        

    Returns
    -------
    listinjections : list of dictionaries
        for each experiment a different dictionary, with as keys
        End Injection Number (int)
        End Minutes (int)
        End Time(datetime)
        experiment (str)(name of experiment)
        Start Injection Number (int)
        Start Minutes (int)
        Start Time ( datetime)

    """
    listinjections = []

    data_tcd = pd.read_excel(excelfile, header =9,skiprows=[10,11,12,13])
    data_tcd['Inject Time '] = pd.to_datetime(data_tcd['Inject Time '],errors='coerce')

    def calculate_closest_time(experiment, data_tcd):
        # Filter rows where 'Inject Time' is greater than or equal to 'Start Time'
        filtered_data = data_tcd[data_tcd['Inject Time '] >= experiment['Start Time']].copy()
        # Calculate the time difference
        filtered_data['time_difference_start'] = abs(filtered_data['Inject Time '] - experiment['Start Time'])      
        # Find the row with the minimum time difference
        closest_row_start = filtered_data.loc[filtered_data['time_difference_start'].idxmin()]
        # Extract the relevant information
        closest_datetime_start = closest_row_start['Inject Time ']
        timestamp_minutes_start = (closest_datetime_start - experiment['Start Time'])
        start_minutes_start = round(timestamp_minutes_start.seconds / 60)


        # Filter rows where 'Inject Time' is less than or equal to 'End Time'
        filtered_data_end = data_tcd[data_tcd['Inject Time '] <= experiment['End Time']].copy()
        # Calculate the time difference
        filtered_data_end['time_difference_end'] = abs(filtered_data_end['Inject Time '] - experiment['End Time'])
        # Find the row with the minimum time difference
        closest_row_end = filtered_data_end.loc[filtered_data_end['time_difference_end'].idxmin()]
        # Extract the relevant information
        closest_datetime_end = closest_row_end['Inject Time ']
        timestamp_minutes_end = (closest_datetime_end - experiment['Start Time'])
        start_minutes_end = round(timestamp_minutes_end.seconds / 60)


        return {
            'Experiment': experiment['Experiment'],
            'Start Time': experiment['Start Time'],
            'End Time': experiment['End Time'],
            'Start Minutes (Start)': start_minutes_start,
            'Start Minutes (End)': start_minutes_end,
            'Start Injection Number': closest_row_start[0],
            'End Injection Number': closest_row_end[0]
        }

    for i, experiment in enumerate(experiment_list):
        result = calculate_closest_time(experiment, data_tcd)
        listinjections.append({
            'experiment': result['Experiment'],
            'Start Time': result['Start Time'],
            'End Time': result['End Time'],
            'Start Minutes': result['Start Minutes (Start)'],
            'End Minutes': result['Start Minutes (End)'],
            'Start Injection Number': result['Start Injection Number'],
            'End Injection Number': result['End Injection Number'],
        })
 

    return listinjections

def process_excel(excelfile, result_list, second_dataframe):
    """
    is used to split the data for each experiment and link the reactor and GC data
    the splitting is based on the injection numbers
    the reactor data is added later on

    Parameters
    ----------
    DFGC : dataframe
        dataframe containing all gc data
    result_list : list of dictionaries
        return from process_injections
    second_dataframe : dataframe
        dataframe containing the reactor data

    Returns
    -------
    result_dataframes : dictionarie of dataframes
        an dictionarie, the keys are the name of an experiment
        each dataframe contains both the GC data, first and the reactor data second

    """

    df = pd.read_excel(excelfile, header = 9, skiprows=[10,11,12,13])
    df['No. '] = pd.to_numeric(df['No. '], errors='coerce')
    #df = df[~df['GC Injection'].isin(['Maximum', 'Average', 'Minimum', 'Standard Deviation', 'Relative Standard Deviation'])]
    
    result_dataframes = {}
    for value in result_list:
        experiment_name = value['experiment'].replace(' ', '_') + "_data"
        
        dhp_test = df[(df['No. '] >= value['Start Injection Number']) & 
                      (df['No. '] <= value['End Injection Number'])]
        dhp_test = dhp_test.reset_index(drop=True)

        # Adding the 'minutes reaction' column
        dhp_test['minutes reaction'] = value['Start Minutes'] + 10 * dhp_test.index
        dhp_test.insert(3, "minutes", dhp_test['minutes reaction'])
        dhp_test.drop(['minutes reaction'], axis=1, inplace=True)

        # Merge with the second dataframe based on the 'Date/Time' column
        merged_df = pd.merge_asof(dhp_test, second_dataframe,left_on='Inject Time ', right_on='Date/Time')

        # Fill missing values in specified columns with values from the last non-null row
        columns_to_fill = ['FIC-110', 'FIC-120', 'FIC-130', 'FIC-140', 'TE-310', 'RSWITCH_Val', 'PT-310', 'PT-320']
        merged_df[columns_to_fill] = merged_df[columns_to_fill].ffill()

        result_dataframes[experiment_name] = merged_df
    
    return result_dataframes

def generate_plots(dft, InputList, export_path="C:\\Users\\r0757763\\OneDrive - KU Leuven\\graphs", export_format="both"):
    """
    Generates plots of temperature, pressure and flow etc., and exports them as PNG, PDF, or both.
    """
    dft['Date/Time'] = pd.to_datetime(dft['Date/Time'], format='%d/%m/%Y %H:%M:%S')
    dft['Hours'] = (dft['Date/Time'] - dft['Date/Time'][0]).dt.total_seconds() / 3600
    dft['NORMOMILLILITER_OUT_Derivative'] = np.gradient(dft['NORMOMILLILITER_OUT'], dft['Hours'])
    dft['NORMOMILLILITER_OUT_Derivative'] = dft['NORMOMILLILITER_OUT_Derivative'].apply(lambda x: x if x > 0 else 0)

    if export_format in ["pdf", "both"]:
        pdf_path = os.path.join(export_path, "all_plots.pdf")
        pdf = PdfPages(pdf_path)

    def generate_plot(column, title, ylabel, filename=None, ylimit=None, export_path=".", export_format="both"):
        fig1, ax = plt.subplots(figsize=(14, 6))
        ax.plot(dft['Hours'], dft[column], label=column, color='blue')
        ax.set_xlabel('Hours', fontsize=14)
        ax.set_ylabel(ylabel, fontsize=14)
        ax.set_title(title)
        for exp in InputList[0::2]:
            ax.axvline(exp['Start Time'], color='gray', linestyle='--', linewidth=1)
            ax.axvline(exp['End Time'], color='red', linestyle='--', linewidth=1)
            ax.text(exp['Start Time'], 0.95 * ax.get_ylim()[1], f"{exp['Experiment']} Start", rotation=90, fontsize=8)
            ax.text(exp['End Time'], 0.95 * ax.get_ylim()[1], f"{exp['Experiment']} End", rotation=90, fontsize=8)
        ax.set_xlim(0, 25)
        if ylimit:
            ax.set_ylim(0, ylimit)
        ax.grid(True)
        plt.tight_layout()

        if export_format in ["png", "both"]:
            plt.savefig(os.path.join(export_path, filename + ".png"))
        if export_format in ["pdf", "both"]:
            pdf.savefig(fig1, bbox_inches='tight')
        plt.show()
        plt.close(fig1)

    # Plot all the defined columns
    columns_to_plot = [
        ('TIC-300_SP', 'TIC-300_SP over Time', 'Temperature'),
        ('TE-310', 'TE-310 over Time', 'Temperature'),
        ('FIC-110', 'FIC-110 over Time', 'Flow'),
        ('FI-110', 'FI-110 over Time', 'Flow'),
        ('FIC-120', 'FIC-120 over Time', 'Flow'),
        ('FI-140', 'FI-140 over Time', 'Flow'),
        ('FIC-130', 'FIC-130 over Time', 'Flow'),
        ('PT-310', 'PT-310 over Time', 'Pressure'),
        ('PT-320', 'PT-320 over Time', 'Pressure'),
        ('PT-400', 'PT-400 over Time', 'Pressure'),
        ('NORMOMILLILITER_OUT', 'NORMOMILLILITER_OUT over Time', 'NORMOMILLILITER_OUT'),
        ('NORMOMILLILITER_OUT_Derivative', 'Positive First Derivative of NORMOMILLILITER_OUT over Time', 'First Derivative'),
        
    ]
    for column, title, ylabel in columns_to_plot:
        generate_plot(column, title, ylabel, filename=column.replace('/', '_'), export_path=export_path, export_format=export_format)

    fig = plt.figure(figsize=(14, 6))
    plt.plot(dft['Hours'], dft['FIC-110'], label='FIC-110', color='blue')
    plt.plot(dft['Hours'], dft['FIC-120'], label='FIC-120', color='green')
    plt.plot(dft['Hours'], dft['FIC-130'], label='FIC-130', color='red')
    plt.plot(dft['Hours'], dft['FIC-140'], label='FIC-140', color='cyan')
    plt.ylabel('nml/min', fontsize = 14)
    plt.legend()
    ax2 = plt.gca().twinx()
    ax2.plot(dft['Hours'],dft['TIC-300_SP'], label = 'Temperature', color = 'black')
    plt.xlabel('hours', fontsize=14)
    plt.ylabel('C°', fontsize=14)
    plt.title('combined flow')
    plt.legend()
    
        
            
    for exp in InputList[0::2]:
        plt.text(exp['Start Time'], 0.9*plt.ylim()[1], exp['Experiment'])
    if export_format in ["png", "both"]:
        plt.savefig(os.path.join(export_path, "combined flow" + ".png"))
    if export_format in ["pdf", "both"]:
        pdf.savefig(fig, bbox_inches='tight')
    plt.show()
    # After all plotting:
    if export_format in ["pdf", "both"]:
        pdf.close()

def reactor_plot(dft ,name, export_path ):
    """
    generates plots of temperature, pressure and flow etc, png file

    Parameters
    ----------
    dft : dataframe
        reactor data
    
    export_path : file paht
        DESCRIPTION. gives file_path where to save

    Returns
    -------
    None.

    """
    data_mcb = dft
    data_mcb['Date/Time'] = pd.to_datetime(data_mcb['Date/Time'], format ='%d/%m/%Y %H:%M:%S')
    data_mcb['Hours'] = (data_mcb['Date/Time'] - data_mcb['Date/Time'][0]).dt.total_seconds() / 3600
    # generating new datablocks
    data_mcb['NORMOMILLILITER_OUT_Derivative'] = np.gradient(data_mcb['NORMOMILLILITER_OUT'], data_mcb['Hours'])
    data_mcb['NORMOMILLILITER_OUT_Derivative'] = data_mcb['NORMOMILLILITER_OUT_Derivative'].apply(lambda x: x if x > 0 else 0)

    def generater_plot(column, title, ylabel, filename=None, ylimit=None, export_path ="."):        
        fig, ax = plt.subplots(figsize=(14, 6))  # Create a new figure and axis
        ax.plot(data_mcb['Hours'], data_mcb[column], label=column, color='blue')
        ax.set_xlabel('Hours', fontsize=14)
        ax.set_ylabel(ylabel, fontsize=14)
        ax.set_title(f"{name} - {title}")
     
        #plt.xlim(0, 25)
        if ylimit:
            plt.ylim(0, ylimit)
        plt.grid(True)
        plt.tight_layout()
        if filename:
            full_path = os.path.join(export_path, filename)  # Combine export_path and filename
            plt.savefig(full_path)
        plt.show()
        plt.close()

    columns_to_plot = [

        # Individual values
        ('TIC-300_SP', 'TIC-300_SP over Time', 'Temperature'),
        ('TE-310', 'TE-310 over Time', 'Temperature'),
        ('FIC-110', 'FIC-110 over Time', 'Flow'),
        ('FI-110', 'FI-110 over Time', 'Flow'),
        ('FIC-120', 'FIC-120 over Time', 'Flow'),
        ('FI-140', 'FI-140 over Time', 'Flow'),
        ('FIC-130', 'FIC-130 over Time', 'Flow'),
        ('PT-310', 'PT-310 over Time', 'Pressure'),
        ('PT-320', 'PT-320 over Time', 'Pressure'),
        ('PT-400', 'PT-400 over Time', 'Pressure'),
        ('NORMOMILLILITER_OUT', 'NORMOMILLILITER_OUT over Time', 'NORMOMILLILITER_OUT'),
        ('NORMOMILLILITER_OUT_Derivative', 'Positive First Derivative of NORMOMILLILITER_OUT over Time', 'First Derivative'),
        
    ]

    for column, title, ylabel in columns_to_plot:
        generater_plot(column, title, ylabel, filename=f"{column.replace('/', '_')}{name}.png",export_path= export_path)

    
    plt.figure(figsize=(14, 6))
    plt.plot(data_mcb['Hours'], data_mcb['FIC-110'], label='FIC-110', color='blue')
    plt.plot(data_mcb['Hours'], data_mcb['FIC-120'], label='FIC-120', color='green')
    plt.plot(data_mcb['Hours'], data_mcb['FIC-130'], label='FIC-130', color='red')
    plt.plot(data_mcb['Hours'], data_mcb['FIC-140'], label='FIC-140', color='cyan')
    plt.ylabel('nml/min', fontsize = 14)
    plt.legend()
    ax2 = plt.gca().twinx()
    ax2.plot(data_mcb['Hours'],data_mcb['TIC-300_SP'], label = 'Temperature', color = 'black')
    plt.xlabel('hours', fontsize=14)
    plt.ylabel('C°', fontsize=14)
    plt.title(f'{name}combined flow')
    plt.legend()
    full_path = os.path.join(export_path,f"{name}combined_flow.png")  # Combine export_path and filename
    plt.savefig(full_path)
    plt.show()
    plt.close()



def Module1():
    """
    an integrated function that combines detect_temperature_ramp, Conditions_list, general_start_times, process_injections
    and process_excel

    Returns
    -------
    invoer : dictionary of dataframes
        an dictionarie, the keys are the name of an experiment
        each dataframe contains both the GC data, first and the reactor data second
        will be used in module 2
        

    """
    #use the detect_temperature_ramp function to detect temperature ramp
    FIC120 = float(input ("Please give the FIC-120 of the temperature ramp : "))
    FIC130 = float(input ("Please give the FIC-130 of the temperature ramp : "))
    FIC140 = float(input("Please give the FIC-140 of the temperature ramp : "))
    MinT = int(input("Please give the starting temperature of the temperature ramp : "))
    excelfile_mcr = input("Please give the excelfile with the reactor data : ")

    # Remove quotes from the file path
    excelfile_mcr = excelfile_mcr.strip('"')

    # Confirm the type of excelfile_mcr


    data = pd.read_excel(excelfile_mcr)
    detect_temperature_ramp(data, FIC120, FIC130, FIC140, MinT)
    
    #use the Conditions_list function
    conditions = input("Please give the excelfile with the conditions: ")
    conditions = conditions.strip('"')
    CON_list = Conditions_List(conditions)
    #globals()['CON_list2'] = Conditions_List(conditions)
    
    
    #use the general_start_time function
    
    starting_times = general_start_times(data, CON_list)
    #globals()['starting_times2'] = general_start_times(data, CON_list)
    
    #use the process_injections function
    excelfileGC = input ("Please give the excelfile containing the GC : ")
    excelfileGC = excelfileGC.strip('"')
    list_injections = process_injections(excelfileGC, starting_times)
    
    #use process_excel function
    invoer = process_excel(excelfileGC, list_injections, data)
    
    return invoer


def plots_Module1():
    """
    a combining function wich combines the temperature ramp detection, conditions list and the different graph functions

    Returns
    -------
    None.

    """
    excelfile_reactor = input("Please give the excelfile with the reactor data : ")
    excelfile_reactor = excelfile_reactor.strip('"')
    data = pd.read_excel(excelfile_reactor)
    export_path = input(" Pleas provide the export path: ")
    which = input("do you want to mark the end of each product, yes = A, B= No : ")
    if which == "A":
    #Tramp function
    #use the detect_temperature_ramp function to detect temperature ramp
        FIC120 = float(input ("Please give the FIC-120 of the temperature ramp : "))
        FIC130 = float(input ("Please give the FIC-130 of the temperature ramp : "))
        FIC140 = float(input("Please give the FIC-140 of the temperature ramp : "))
        MinT = int(input("Please give the MINT of the temperature ramp : "))
    
        detect_temperature_ramp(data, FIC120, FIC130, FIC140, MinT)
    
    #use Conditions_list
        conditions = input("Please give the excelfile with the conditions: ")
        conditions = conditions.strip('"')
        CON_list = Conditions_List(conditions)

    #use function general_start_hours
        hours_start_list = general_start_hours(data, CON_list)
    #decide wich function you want use
        which_function = input("Do you want to have an PDF file (A) or PNG images(B) or both (C): ")
        if which_function == "C":
            generate_plots(data, hours_start_list, export_path, export_format= "both")
        elif which_function == "A":
            generate_plots(data, hours_start_list, export_path, export_format="pdf")
        elif which_function == "B":
            generate_plots(data, hours_start_list, export_path, export_format= "png")
    else:
        name = input("what is the name: ")
        reactor_plot(data, name,export_path )

"""
functions module 2
"""
def replace_nan_with_zero(dictionary_of_dataframes):
    """
    replace n.a with 0, to avoid troubles furtheron

    Parameters
    ----------
    dictionary_of_dataframes : dictionary of dataframes
       coming from module i 

    Returns
    -------
    updated dataframe

    """
    for key, dataframe in dictionary_of_dataframes.items():
        dataframe.replace("n.a.", np.nan, inplace=True)
        dictionary_of_dataframes[key] = dataframe.fillna(0)

def calcul_som(dictf):
    for key, df in dictf.items():
        df['sumcompo'] = df.loc[:,"Argon":"C2-butene"].sum(axis=1)
    return dictf


def remove_outliers(df_dict):
    """
    A function to remove outliers based on the Interquartile Range (IQR) method.
    Outliers are removed if a data point is outside [Q1 - 3*IQR, Q3 + 3*IQR].
    
    Parameters
    ----------
    df_dict : dictionary of pandas DataFrames
        Dictionary containing dataframes with data from experiments. Each dataframe should have numerical data.

    Returns
    -------
    df_dict : dictionary of pandas DataFrames
        Updated dictionary where each dataframe has had outliers removed.
    """

    for key, df in df_dict.items():
        Q1 = df['sumcompo'].quantile(0.25)
        Q3 = df['sumcompo'].quantile(0.75)
        IQR = Q3 - Q1
        lower_bound = Q1 - 3 * IQR
        upper_bound = Q3 + 3 * IQR
        # Filter out rows with outliers in the current column
        df = df[(df['sumcompo'] >= lower_bound) & (df['sumcompo'] <= upper_bound)]

        df_dict[key] = df


    return df_dict

def normalise(dictf):
    """
    a function to normalise the obtained GC data
    where the individual composition is divided by the total composition

    Parameters
    ----------
    dictf : dictionary of dataframes
        dictionary of dataframes coming from remove outliers

    Returns
    -------
    None.

    """
    column_name = ['Carbon Dioxide', 'Hydrogen', 'Argon',  'Nitrogen', 'Carbon Monoxide', 'Methane',
                                        'Ethene', 'Ethane', 'Propylene', 'Propane', 'Isobutane',
                                        'n-Butane',  'C2-butene']

    for key, df in dictf.items():
        
        for col in column_name:
            df[col]= (df[col]/df['sumcompo'])*100
    return dictf
#normalise(RB003)


def molar_flows(DICdf, a, b):
    """
    calculates the inflow, from nml/min to percentage and mol/m³

    Parameters
    ----------
    DICdf : dictionary of dataframes
        result of module 1, contains both the reactor data and the GC data

    Returns
    -------
    updated dataframes

    """
    for key, df in DICdf.items():
        df['O2_nml_flow'] = df['FIC-110'] * 0.2095
        df['CO2_nml_flow'] = df['FIC-120']
        df['propane_nml_flow'] = df['FIC-130'] * a
        df['Argon_nml_flow'] = df['FIC-130'] * b + df['FIC-110'] * 0.00934
        df['N2_nml_flow'] = df['FIC-140'] + df['FIC-110'] * 0.7808
        df['Pcatalyst'] = (df['PT-310'] + df['PT-320']) / 2
        df['sum'] = (
            df['O2_nml_flow'] +
            df['CO2_nml_flow'] +
            df['propane_nml_flow'] +
            df['Argon_nml_flow'] +
            df['N2_nml_flow']
        )
        df['02%'] = df['O2_nml_flow'] / df['sum']
        df['CO2%'] = df['CO2_nml_flow'] / df['sum']
        df['propane%'] = df['propane_nml_flow'] / df['sum']
        df['Argon%'] = df['Argon_nml_flow'] / df['sum']
        df['N2%'] = df['N2_nml_flow'] / df['sum']
        df['O2 mol/m**3'] = df['02%'] * (df['PT-610'] + df['Pcatalyst']) * (10 ** 5) / (
                8.314 * (273.15 + df['TE-310']))
        df['CO2 mol/m**3'] = df['CO2%'] * (df['PT-610'] + df['Pcatalyst']) * (10 ** 5) / (
                8.314 * (273.15 + df['TE-310']))
        df['propane mol/m**3'] = df['propane%'] * (df['PT-610'] + df['Pcatalyst']) * (10 ** 5) / (
                8.314 * (273.15 + df['TE-310']))
        df['Argon mol/m**3'] = df['Argon%'] * (df['PT-610'] + df['Pcatalyst']) * (10 ** 5) / (
                8.314 * (273.15 + df['TE-310']))
        df['N2 mol/m**3'] = df['N2%'] * (df['PT-610'] + df['Pcatalyst']) * (10 ** 5) / (
                8.314 * (273.15 + df['TE-310']))


def out_flow(dataframes_dict, corr,IS = True):
    """
    calculates the outlfow from the gc data using an internal standard from N2, if IS = true
    has the option to correct for pressure

    Parameters
    ----------
    dataframes_dict : dictionary of dataframes
        output from inflow
    corr : boolean
        True = pressure correction
        False = no pressure correction
    IS : boolean
        True = N2 internal standard
        False = Argon  internal standard

    Returns
    -------
    dataframes_dict : dictionary of dataframes
        updated with the outflow

    """
    for key, df in dataframes_dict.items():

        column_name = ['Carbon Dioxide', 'Hydrogen', 'Argon', 'Nitrogen', 'Carbon Monoxide', 'Methane',
                                        'Ethene', 'Ethane', 'Propylene', 'Propane', 'Isobutane', 
                                        'n-Butane', 'C2-butene']
        df[column_name] = df[column_name].apply(pd.to_numeric, errors='coerce')
        if corr == True:
            
            for column in column_name:
                df[column] = df[column] / df['PT-610']
            
    for key, df in dataframes_dict.items():
        # Outlet volumetric flow
        if IS == True:
            df['IS'] = df['N2_nml_flow']/df['Nitrogen']
        else:
            
            df['IS'] = df['Argon_nml_flow'] / (df['Argon'] )
        column_Pcor = ['Carbon Dioxide', 'Hydrogen', 'Argon', 'Nitrogen',
                       'Carbon Monoxide', 'Methane', 'Ethene', 'Ethane', 'Propylene',
                       'Propane', 'Isobutane', 'n-Butane', 'C2-butene']
        for column in column_Pcor:
            new_colum_name = f'{column}_Vol_flow_nml/min'
            df[new_colum_name] = df['IS'] * (df[column] )
        df['somflow'] = df.filter(like = '_Vol_flow_nml/min').sum(axis = 1)
    return dataframes_dict

def total_atom_flow(dataframes_dict):
    """
    Calculates the total flow of C, H, and O atoms for each molecule.

    Parameters
    ----------
    dataframes_dict : dict
        Dictionary of DataFrames with molecule volume flows.

    Returns
    -------
    dict
        Updated dictionary with additional columns for C, H, and O atom flows.
    """
    for key, df in dataframes_dict.items():
        for molecule, atom_counts in Chemicals.items():
            column_nameC = f'{molecule}_C_atoms_flow'
            column_nameH = f'{molecule}_H_atoms_flow'
            column_nameO = f'{molecule}_O_atoms_flow'
            df[column_nameC] = df[f'{molecule}_Vol_flow_nml/min'] * atom_counts['C_atoms']
            df[column_nameH] = df[f'{molecule}_Vol_flow_nml/min'] * atom_counts['H_atoms']
            df[column_nameO] = df[f'{molecule}_Vol_flow_nml/min'] * atom_counts['O_atoms']

        # Sum the atom flow columns for each type of atom
        df['Total_C_atoms_flow'] = df.filter(like='_C_atoms_flow').sum(axis=1)
        df['Total_H_atoms_flow'] = df.filter(like='_H_atoms_flow').sum(axis=1)
        df['Total_O_atoms_flow'] = df.filter(like='_O_atoms_flow').sum(axis=1)

        # Optionally, remove intermediate molecule-specific atom flow columns
        for atom in ['C', 'H', 'O']:
            intermediate_columns = df.filter(like=f'_{atom}_atoms_flow').columns
            intermediate_columns = [col for col in intermediate_columns if col.startswith('Total_') is False]
            df.drop(columns=intermediate_columns, inplace=True)

            
    return dataframes_dict

def mass_balance(dataframes_dict):
    """
    Calculates the mass balance comparing general data with corresponding bypass data.
    Uses averages from 'bypass_*' DataFrames to calculate mass balances in corresponding '*_data' DataFrames.

    Parameters
    ----------
    dataframes_dict : dict
        Dictionary of DataFrames updated with total atom flows.

    Returns
    -------
    dict
        Updated dictionary with mass balance calculations.
    """
    # Dictionary to hold averages from bypass dataframes
    averages_dict = {}

    # First, compute averages from bypass dataframes
    for key, df in dataframes_dict.items():
        if key.startswith('bypass_'):
            avg_dict = {}
            for atom in ['C', 'H', 'O']:
                avg_dict[atom] = df[f'Total_{atom}_atoms_flow'].mean()
            # Strip 'bypass_' and '_data' to get the base name
            base_name = key[7:-5]
            averages_dict[base_name] = avg_dict

    # Second, apply these averages to calculate mass balances for corresponding general dataframes
    for key, df in dataframes_dict.items():
        if key.endswith('_data') and not key.startswith('bypass_'):
            base_name = key[:-5]  # Remove '_data'
            if base_name in averages_dict:
                avg_dict = averages_dict[base_name]
                for atom in ['C', 'H', 'O']:
                    avg_atom_flow = avg_dict[atom]
                    balance_column = f'{atom}-balance'
                    df[balance_column] = df[f'Total_{atom}_atoms_flow'] / avg_atom_flow if avg_atom_flow != 0 else None

    return dataframes_dict


def X_S_Y_r_bypass(dataframes_dict, M_catalyst):
    """
    a function to calculate the conversion (propane and co2), with bypass
    the selectivity (odhp, cracking and reforming)
    yield( of propylene) and space time yield
    rates of propane, co2, h2, propylene,methane, ethene, ethane and CO

    Parameters
    ----------
    dataframes_dict : dictionary of dataframes
        output from mass balance/ outflow

    Returns
    -------
    dataframes_dict : dictionary of dataframes
        updated 

    """
    for key, df in dataframes_dict.items(): #conversion bypass, propane and CO2  
        if key.endswith('_data'):
            base_nameX = key.replace('_data','')
            
            bypass_keyX = f'bypass_{base_nameX}_data'
            if bypass_keyX in dataframes_dict:
                
                bypass_avg_propane_flow = dataframes_dict[bypass_keyX]['Propane_Vol_flow_nml/min'].iloc[1:-1].mean()
                bypass_avg_CO2_flow = dataframes_dict[bypass_keyX]['Carbon Dioxide_Vol_flow_nml/min'].iloc[1:-1].mean()
                df['X_propane_bypass'] = (1-(df['Propane_Vol_flow_nml/min']/bypass_avg_propane_flow))
                df['X_CO2_bypass'] = (1-(df['Carbon Dioxide_Vol_flow_nml/min']/bypass_avg_CO2_flow))
                df['Y_propylene_bypass'] = df['Propylene_Vol_flow_nml/min']/bypass_avg_propane_flow #propylene yield
                df['STY_propylene'] = df['STY_propylene'] = df['Propylene_Vol_flow_nml/min']/M_catalyst
                df['r_propane'] = (((bypass_avg_propane_flow/10**6)*103125)/(8.314*273.15*60))*(df['X_propane_bypass']/(M_catalyst/1000))*3600
                df['r_CO2'] = (((bypass_avg_CO2_flow/10**6)*103125)/(8.314*273.15*60))*(df['X_CO2_bypass']/(M_catalyst/1000))*3600
                df['r_H2'] = (((df['Hydrogen_Vol_flow_nml/min']/10**6)*103125*60)/(8.314*273.15))/(M_catalyst/1000)
                df['r_propylene'] = (((df['Propylene_Vol_flow_nml/min']/10**6)*103125*60)/(8.314*273.15))/(M_catalyst/1000)    
                df['r_Methane'] = (((df['Methane_Vol_flow_nml/min']/10**6)*103125*60)/(8.314*273.15))/(M_catalyst/1000)
                df['r_Ethene'] = (((df['Ethene_Vol_flow_nml/min']/10**6)*103125*60)/(8.314*273.15))/(M_catalyst/1000)
                df['r_Ethane'] = (((df['Ethane_Vol_flow_nml/min']/10**6)*103125*60)/(8.314*273.15))/(M_catalyst/1000)
                df['r_CO'] = (((df['Carbon Monoxide_Vol_flow_nml/min']/10**6)*103125*60)/(8.314*273.15))/(M_catalyst/1000)
                
                #calculating selectivity
                df['CO2_consumed'] = bypass_avg_CO2_flow- df['Carbon Dioxide_Vol_flow_nml/min']
                df['FCO_from_propane'] = df['Carbon Monoxide_Vol_flow_nml/min'] - df['CO2_consumed']
                #use of bypass, big diffrence between different selectivities
                df['S_prop_odhp2'] = df['Propylene_Vol_flow_nml/min']/ (bypass_avg_propane_flow - df['Propane_Vol_flow_nml/min'])
                df['S_prop_cracking2'] = (df['Methane_Vol_flow_nml/min'] + 2* df['Ethene_Vol_flow_nml/min'])/(3*(bypass_avg_propane_flow - df['Propane_Vol_flow_nml/min']))
                df['S_prop_reforming2'] = df['FCO_from_propane']/(3*(bypass_avg_propane_flow - df['Propane_Vol_flow_nml/min']))
                df['sum2'] = df['S_prop_odhp2']+df['S_prop_cracking2']+df['S_prop_reforming2']
                
                #calculating WHSV
                df['WHSV'] = ((bypass_avg_propane_flow*10**-3)/22.4)*AM_propane*60

    return dataframes_dict

def X_S_Y_r_product(DICTF,M_catalyst):
    """
    a function to calculate the conversion (propane and co2), with product
    the selectivity (odhp, cracking and reforming)
    yield( of propylene) and space time yield
    rates of propane, co2, h2, propylene,methane, ethene, ethane and CO

    Parameters
    ----------
    DICTF : dictionary of dataframes
        output from outflow/ mass balance

    Returns
    -------
    DICTF : dictionary of dataframes
        updated

    """
    for key, df in DICTF.items():#conversion with products
        df['X_propane_productsC'] = 1-(3*df['Propane_Vol_flow_nml/min'])/df['Total_C_atoms_flow']# not possible due though addition of CO2
        df['X_H_propane_productsH'] = 1-(8*df['Propane_Vol_flow_nml/min'])/df['Total_H_atoms_flow']
        df['X_CO2_productsO'] = 1-(2*df['Carbon Dioxide_Vol_flow_nml/min'])/df['Total_O_atoms_flow']
        df['Y_propylenep'] = df['Propylene_Vol_flow_nml/min']/df['propane_nml_flow'] #propylene yield
        df['STY_propylenep'] = df['Propylene_Vol_flow_nml/min']/M_catalyst
        df['r_propanep'] = (((df['propane_nml_flow']/10**6)*103125)/(8.314*273.15*60))*(df['X_H_propane_productsH']/(M_catalyst/1000))*3600
        df['r_CO2p'] = (((df['CO2_nml_flow']/10**6)*103125)/(8.314*273.15*60))*(df['X_CO2_productsO']/(M_catalyst/1000))*3600
        df['r_H2p'] = (((df['Hydrogen_Vol_flow_nml/min']/10**6)*103125*60)/(8.314*273.15))/(M_catalyst/1000)
        df['r_propylenep'] = (((df['Propylene_Vol_flow_nml/min']/10**6)*103125*60)/(8.314*273.15))/(M_catalyst/1000)    
        df['r_Methanep'] = (((df['Methane_Vol_flow_nml/min']/10**6)*103125*60)/(8.314*273.15))/(M_catalyst/1000)
        df['r_Ethenep'] = (((df['Ethene_Vol_flow_nml/min']/10**6)*103125*60)/(8.314*273.15))/(M_catalyst/1000)
        df['r_Ethanep'] = (((df['Ethane_Vol_flow_nml/min']/10**6)*103125*60)/(8.314*273.15))/(M_catalyst/1000)
        df['r_COp'] = (((df['Carbon Monoxide_Vol_flow_nml/min']/10**6)*103125*60)/(8.314*273.15))/(M_catalyst/1000)
        #calculate selectivity
        df['CO2_consumedp'] = df['CO2_nml_flow']- df['Carbon Dioxide_Vol_flow_nml/min']
        df['FCO_from_propanep'] = df['Carbon Monoxide_Vol_flow_nml/min'] - df['CO2_consumedp']
        #selectivity based on 
        df['C_atoms_selectivity'] = df['Total_C_atoms_flow']- df['Carbon Monoxide_Vol_flow_nml/min'] -df['Carbon Dioxide_Vol_flow_nml/min'] - 3*df['Propane_Vol_flow_nml/min']
        df['denominator_selectivity'] =  df['C_atoms_selectivity'] + df['FCO_from_propanep'] 
        df['S_prop_odhp'] = 3*df['Propylene_Vol_flow_nml/min']/df['denominator_selectivity']
        df['S_prop_cracking'] = (df['Methane_Vol_flow_nml/min'] + 2* df['Ethene_Vol_flow_nml/min'])/df['denominator_selectivity']
        df['S_prop_reforming'] = df['FCO_from_propanep']/df['denominator_selectivity']
        df['sumS'] = df['S_prop_odhp']+df['S_prop_cracking']+df['S_prop_reforming']
    return DICTF

def plot_X_S_bypass(data_dict):
    """
    plots the conversion and selectivity calculated with the bypass of all the different stages in the experiment
    for the bypass calculations
    Parameters
    ----------
    data_dict : dictionary of dataframes
        from X_S_Y_r_bypass

    Returns
    -------
    graphs

    """
    for name, df in data_dict.items():
        if not name.startswith("bypass"):
            if not name.startswith("Tramp"):
                df_filtered = df.fillna(value=0)
                fig, ax1 = plt.subplots()

                ax1.scatter(df_filtered["minutes"], df_filtered["X_propane_bypass"], label="X_Propane", color = 'red' )
                ax1.scatter(df_filtered["minutes"], df_filtered["X_CO2_bypass"], label="X_CO2", color = 'black')
                ax1.set_xlabel("Minutes")
                ax1.set_ylabel("X Values", color='b')
                ax1.tick_params(axis='y', labelcolor='b')
                ax1.set_ylim(0, 0.6)

                # selectivity
                ax2 = ax1.twinx()
                selectivity_colors = ['skyblue', 'lightgreen', 'lightcoral']
                selectivity_labels = ['S_odhp', 'S_cracking', 'S_reforming']

                selectivity_data = [df['S_prop_odhp2'], df['S_prop_cracking2'], df['S_prop_reforming2']]
                ax2.stackplot(df['minutes'], *selectivity_data, labels=selectivity_labels, colors=selectivity_colors, alpha=0.5)

                ax2.set_ylabel("selectivity")
                ax2.tick_params(axis='y', labelcolor='r')
                ax2.set_ylim(0,1)

                # combined plot
                plt.title(f' Plot : {name}  with bypass ')
                lines1, labels1 = ax1.get_legend_handles_labels()
                lines2, labels2 = ax2.get_legend_handles_labels()
                plt.legend(lines1 + lines2, labels1 + labels2, loc='upper right',bbox_to_anchor=(1.5, 1), borderaxespad=0.)
                plt.show()

                # Clear the plot for the next dataframe
                plt.clf()

#plot_X_S_bypass(RB008)
def plot_X_S_products(data_dict):
    """
    plots the conversion and selectivity calculated with the bypass of all the different stages in the experiment
    for the products calculations

    Parameters
    ----------
    data_dict : dictionary of dataframes
        from X_S_Y_r_products

    Returns
    -------
    graphs

    """
    for name, df in data_dict.items():
        if not name.startswith("bypass"):
            if not name.startswith("Tramp"):
                df_filtered = df.fillna(value=0)
                fig, ax1 = plt.subplots()

                ax1.scatter(df_filtered["minutes"], df_filtered["X_H_propane_productsH"], label="X_Propane", color = 'red' )
                ax1.scatter(df_filtered["minutes"], df_filtered["X_CO2_productsO"], label="X_CO2", color = 'black')
                ax1.set_xlabel("Minutes")
                ax1.set_ylabel("X Values", color='b')
                ax1.tick_params(axis='y', labelcolor='b')
                ax1.set_ylim(0, 0.6)

                # selectivity
                ax2 = ax1.twinx()
                selectivity_colors = ['skyblue', 'lightgreen', 'lightcoral']
                selectivity_labels = ['S_odhp', 'S_cracking', 'S_reforming']

                selectivity_data = [df['S_prop_odhp'], df['S_prop_cracking'], df['S_prop_reforming']]
                ax2.stackplot(df['minutes'], *selectivity_data, labels=selectivity_labels, colors=selectivity_colors, alpha=0.5)

                ax2.set_ylabel("selectivity")
                ax2.tick_params(axis='y', labelcolor='r')
                ax2.set_ylim(0,1)

                # combined plot
                plt.title(f' Plot : {name} products')
                lines1, labels1 = ax1.get_legend_handles_labels()
                lines2, labels2 = ax2.get_legend_handles_labels()
                plt.legend(lines1 + lines2, labels1 + labels2, loc='upper right',bbox_to_anchor=(1.5, 1), borderaxespad=0.)
                plt.show()

                # Clear the plot for the next dataframe
                plt.clf()
#plot_X_S_products(RB008)
def create_gui(data,name ,output_path):
    top_columns = ['r_propane', 'r_CO2', 'r_H2', 'r_propylene', 'r_Methane', 'r_Ethane', 'r_Ethene', 'r_CO', 
                   'r_propanep', 'r_CO2p', 'r_H2p', 'r_propylenep', 'r_Methanep', 'r_Ethanep', 'r_Ethenep', 'r_COp',
                   'S_prop_odhp2', 'S_prop_cracking2', 'S_prop_reforming2', 'sum2', 
                   'S_prop_odhp', 'S_prop_cracking', 'S_prop_reforming', 'SumS']

    root = tk.Tk()
    root.title("Data Plotter")

    # Function to filter DataFrame keys
    def get_filtered_keys(data):
        return [k for k in data.keys() if not k.startswith('bypass') and not k.startswith('Tramp')]

    # Function to plot data
    def plot_data():
        col = column_combo.get()
        plot_type = plot_type_combo.get()

        fig, ax = plt.subplots()
        for df_key in get_filtered_keys(data):
            df = data[df_key]
            if col in df.columns:
                if plot_type == "Scatter Plot":
                    ax.scatter(df['minutes'], df[col], label=df_key)
                elif plot_type == "Box Plot":
                    ax.boxplot(df[col], positions=[df_keys.index(df_key)], labels=[df_key])
        
        ax.set_xlabel('Minutes')
        ax.set_ylabel(col)
        ax.set_title(f'{plot_type} of {col}_{name}')
        ax.legend()

        # Embed the plot in the tkinter window
        canvas = FigureCanvasTkAgg(fig, master=root)
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.grid(row=4, column=0, columnspan=3)
        canvas.draw()

    # Function to save the plot
    def save_plot():
        col = column_combo.get()
        plot_type = plot_type_combo.get()
        fig, ax = plt.subplots()
        for df_key in get_filtered_keys(data):
            df = data[df_key]
            if col in df.columns:
                if plot_type == "Scatter Plot":
                    ax.scatter(df['minutes'], df[col], label=df_key)
                elif plot_type == "Box Plot":
                    ax.boxplot(df[col], positions=[df_keys.index(df_key)], labels=[df_key])

        ax.set_xlabel('Minutes')
        ax.set_ylabel(col)
        ax.set_title(f'{plot_type} of {col}_{name}')
        ax.legend()

        # Save the figure
        save_path = os.path.join(output_path, f"{plot_type}_{col}_{name}.png")
        fig.savefig(save_path)
        plt.close(fig)

    df_keys = get_filtered_keys(data)
    column_combo = ttk.Combobox(root)
    column_combo.grid(row=1, column=1)

    # Update columns when DataFrames are selected
    def update_columns(*args):
        selected_dfs = [df_key for df_key in df_keys if df_key in data]
        all_cols = set(col for df_key in selected_dfs for col in data[df_key].columns if 'minutes' != col)
        ordered_cols = [col for col in top_columns if col in all_cols] + [col for col in all_cols if col not in top_columns]
        column_combo['values'] = ordered_cols

    update_columns()

    plot_types = ["Scatter Plot", "Box Plot"]
    plot_type_combo = ttk.Combobox(root, values=plot_types, state="readonly")
    plot_type_combo.grid(row=2, column=1)

    plot_button = tk.Button(root, text="Plot", command=plot_data)
    plot_button.grid(row=2, column=2)

    save_button = tk.Button(root, text="Save Plot", command=save_plot)
    save_button.grid(row=3, column=2)

    root.mainloop()


def Call_Module2(dictf):
    replace_nan_with_zero(dictf)
    calcul_som(dictf)
    remove_outliers(dictf)
    Norm = input("Do you want to normalise the GC flow , A = True, or B = False: ")
    if Norm == "A":
        normalise(dictf)
    inflow = input("what is the composition of FIC-130, A = 5% propane, B 20% propane, C other : ")
    if inflow == "A":
        molar_flows(dictf, 0.05, 0.95)
    if inflow == "B":
        molar_flows(dictf, 0.2, 0.8)
    if inflow == "C":
        prop_per = float(input("what is the percentage of propane in FIC-130? : "))
        prop = prop_per/100
        argon = 1-prop
        molar_flows(dictf, prop, argon)
    cor = input("Do you want pressure correction A = yes, B = No: ")
    i_s = input("Which internal standard do you want to use? A= N2, B = Argon: ")
    if cor == "A":
        COR = True
        if i_s == "A":
            out_flow(dictf, corr=True)
        else:
            out_flow(dictf, corr=True, IS=False)#Argon
    else:
        if i_s == "A":
            out_flow(dictf, corr= False, IS=True)
        else:
            out_flow(dictf, corr = False, IS = False)
    total_atom_flow(dictf)
    mass_balance(dictf)
    M_catalyst = float(input("please provide the mass of the catalyst used in grams : "))
    bypass = input(" bypass or products or both, A = bypass, B = product, C= both : ")
    if bypass == "A":
        X_S_Y_r_bypass(dictf, M_catalyst)

    if bypass == "B":
        X_S_Y_r_product(dictf, M_catalyst)

    else :
        X_S_Y_r_bypass(dictf, M_catalyst)

        X_S_Y_r_product(dictf, M_catalyst)

    return dictf

def AE_Module2(dictf):
    dic_new = {}
    first_name = argname('dictf')
    for key, df in dictf.items():
        if key.startswith("Tramp"):
            new_dict = {}
            df['k'] = df['r_propylenep'] / df['propane_nml_flow']
            df['ln(k)'] = np.log(df['k'])
            df['1/T'] = 1 / df['TIC-300_PV']

            X = df[['1/T']]  # Independent variable
            y = df['ln(k)']  # Dependent variable

            model = LinearRegression()

            # Fitting the model to your data
            model.fit(X, y)

            # Getting the slope (activation energy) and intercept of the linear regression line
            temperatureterm = model.coef_[0]
            intercept = model.intercept_

            R_square = model.score(X, y)
            R = 8.314

            Activation_energy = -1 * temperatureterm * R
            A = np.exp(intercept)

            new_dict['activation_energy'] = Activation_energy
            new_dict['pre_exponential_factor'] = A
            new_dict['R_squared'] = R_square

            dic_new['activation_energy'] = Activation_energy
            dic_new['pre_exponential_factor'] = A
            dic_new['R_squared'] = R_square

    return {f'{first_name}_Activation_energy': dic_new}
#RB002AE = activation_energy(RB002CAL)
def add_activation_energy(data, activation_energy):
    # Iterate through each sub-dictionary in data
    for key in data:
        # Update the sub-dictionary with the activation_energy dictionary
        data[key].update(activation_energy)
    return data
#RB002TT = add_activation_energy(RB002T, RB002AE)

def graph_Module2(dictf):
    byp = input(" bypass or products or both, A = bypass, B = product : ")
    if byp == "A":
        plot_X_S_bypass(dictf)
    if byp == "B":
        plot_X_S_products(dictf)
    export_path = input("please provide the export_path where you want to share the images : ")
    name = input("what's the name of the experiment?: ")
    #export_path = export_path.strip('"')
    create_gui(dictf,name, export_path)

def split_dataframes(dictf,fnt,M_catalyst):
    keys = ['GC', 'Reactor','Inflow','Outflow','Massbalance','X_S_Y_r_bypass', 'X_S_Y_r_products', 'Mass_catalyst']
    dic_dic = {}
    first_name = argname('dictf')
    
    for key, df in dictf.items():
        if not key.startswith("bypass"):
            new_dict = dict()
            for new_key in keys:
                if new_key == "GC":
                    new_dict[new_key] = df.loc[:, "No. ":"C2-butene"]
                    new_dict[new_key]['sumcompo'] = df['sumcompo']
                elif new_key == "Reactor":
                    new_dict[new_key] = df.loc[:, "Date/Time":"T_ramp"]
                    new_dict[new_key]['minutes'] = df['minutes']
                elif new_key == "Inflow":
                    new_dict[new_key] = df.loc[:, "O2_nml_flow":"N2 mol/m**3"]
                    new_dict[new_key]['minutes'] = df['minutes']
                elif new_key == "Outflow":
                    new_dict[new_key] = df.loc[:, "IS":"somflow"]
                    new_dict[new_key]['minutes'] = df['minutes']
                elif new_key == "Massbalance":
                    try:
                        new_dict[new_key] = df.loc[:, "Total_C_atoms_flow":"O-balance"]
                        new_dict[new_key]['minutes'] = df['minutes']
                    except:
                        continue
                elif new_key == "X_S_Y_r_bypass":
                    try:
                        new_dict[new_key] = df.loc[:, "X_propane_bypass":"WHSV"]
                        new_dict[new_key]['minutes'] = df['minutes']
                    except:
                        continue
                elif new_key == "X_S_Y_r_products":
                    try:
                        new_dict[new_key] = df.loc[:, "X_propane_productsC":"sumS"]
                        new_dict[new_key]['minutes'] = df['minutes']
                    except:
                        continue
                elif new_key == "Mass_catalyst":
                    new_dict[new_key] = M_catalyst
                dic_dic[f'{fnt}_{key}'] = new_dict

    return dic_dic
#RB002T = split_dataframes(RB002CAL, 'RB002i', M_catalyst)

def split_Module2(dictf, AE = None):
    M_catalyst = float(input("please provide the mass of the catalyst used in grams : "))
    fnt = input("what is the name of the experiment? : ")
    end_result = split_dataframes(dictf, fnt, M_catalyst)
    ae = input("Did you calculate the activation_energy, A= yes, B = No: ")
    if ae == "A":
        add_activation_energy(end_result, AE)
    return end_result
#ok = split_Module2(RB002CAL)
def save_to_excel(data_dict, output_dir='output'):
    """
    Exports nested dictionaries to multiple Excel files.
    Each sub-dictionary creates an Excel file where each DataFrame is saved in a separate sheet.
    The float values and nested dictionaries of floats are also saved in separate sheets.

    Args:
    data_dict (dict): Dictionary where each key contains another dictionary with DataFrames, floats, and possibly another dictionary with floats.
    output_dir (str): Directory to save the Excel files.
    """

    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Iterate over the main dictionary
    for key, sub_dict in data_dict.items():
        # Define file path for each sub-dictionary
        file_path = os.path.join(output_dir, f'{key}.xlsx')
        
        # Create an Excel writer for each sub-dictionary
        with pd.ExcelWriter(file_path, engine='openpyxl') as writer:
            # Iterate over items in each sub-dictionary
            for sub_key, value in sub_dict.items():
                if isinstance(value, pd.DataFrame):
                    # Write DataFrame to a sheet named after the sub_key
                    value.to_excel(writer, sheet_name=sub_key)
                elif isinstance(value, float):
                    # Handle the float by writing it to a dedicated sheet
                    df = pd.DataFrame([value], columns=['Value'])
                    df.to_excel(writer, sheet_name=f'{sub_key} Value')
                elif isinstance(value, dict):
                    # Assume this dictionary only contains floats
                    # Create a DataFrame and write it to a sheet
                    df = pd.DataFrame(list(value.items()), columns=['Metric', 'Value'])
                    df.to_excel(writer, sheet_name=f'{sub_key}')

    print(f"Excel files have been saved in '{output_dir}' directory.")

#save_to_excel(RB002TT, r'C:\Users\robbe\OneDrive - KU Leuven\2. Robbe Beirinckx Master thesis\3. Data\excelmodule2')
def save_pickle(parent_dict, output_dir="pickle_files"):
    """
    Saves each sub-dictionary in 'parent_dict' as a pickle file.
    The name of the pickle file will be the key of the sub-dictionary.
    
    :param parent_dict: Dictionary of dictionaries to be saved as pickle files.
    :param output_dir: Directory where pickle files will be saved.
    """
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for key, sub_dict in parent_dict.items():
        # Generate the file path
        file_path = os.path.join(output_dir, f"{key}.pickle")
        
        # Open the file and save the dictionary
        with open(file_path, 'wb') as file:
            pickle.dump(sub_dict, file)
    print("files have been saved as pickle file")
def save_Module2(dictf, export_path):
    print("this is for split data, please be sure that youre dataframes are split")
    save_pickle(dictf, export_path)
    save_to_excel(dictf, export_path)
#save_Module2(ok,r"C:\Users\robbe\OneDrive - KU Leuven\2. Robbe Beirinckx Master thesis\3. Data\excelmodule2")

"""
Module 3
"""
def load_Module3(file_path):
    with open(file_path, 'rb') as file:
    # Load the data from the pickle file
        out_file = pd.read_pickle(file)
    return out_file
#load_Module3(r"C:\Users\robbe\OneDrive - KU Leuven\2. Robbe Beirinckx Master thesis\3. Data\RB002_NPC_N2_NN_odhp_p6_data.pickle")

def plot_Module3(dict1, dict2, default_save_path):
    n1 = argname('dict1')
    n2 = argname('dict2')
    # Initialize the main window
    root = tk.Tk()
    root.title("Data Comparison Tool")

    # Identify common keys with DataFrame values
    common_keys = {k for k in dict1 if k in dict2 and isinstance(dict1[k], pd.DataFrame) and isinstance(dict2[k], pd.DataFrame)}

    # Drop-down for selecting the common dictionary key
    tk.Label(root, text="Select Key:").grid(row=0, column=0, padx=10, pady=10)
    key_var = tk.StringVar()
    key_dropdown = ttk.Combobox(root, textvariable=key_var)
    key_dropdown['values'] = list(common_keys)
    key_dropdown.grid(row=0, column=1)

    # Drop-down for selecting the DataFrame column (other than 'minutes' for the Y-axis)
    tk.Label(root, text="Select Y-Axis Column:").grid(row=1, column=0, padx=10, pady=10)
    column_var = tk.StringVar()
    column_dropdown = ttk.Combobox(root, textvariable=column_var)
    column_dropdown.grid(row=1, column=1)

    # Update columns when a key is selected
    def update_columns(*args):
        selected_key = key_var.get()
        if selected_key:
            columns = set(dict1[selected_key].columns).intersection(dict2[selected_key].columns)
            if 'minutes' in columns:
                columns.remove('minutes')
            column_dropdown['values'] = list(columns)

    key_var.trace('w', update_columns)

    # Radio buttons for plot type selection
    tk.Label(root, text="Select Plot Type:").grid(row=2, column=0, padx=10, pady=10)
    plot_type_var = tk.StringVar(value="scatter")
    tk.Radiobutton(root, text="Scatter Plot", variable=plot_type_var, value="scatter").grid(row=2, column=1)
    tk.Radiobutton(root, text="Box and Whisker Plot", variable=plot_type_var, value="box").grid(row=2, column=2)

    # Frame for the plot
    plot_frame = tk.Frame(root)
    plot_frame.grid(row=4, column=0, columnspan=3, pady=20)

    # Function to plot the data
    def plot_data():
        selected_key = key_var.get()
        selected_column = column_var.get()
        plot_type = plot_type_var.get()
        df1 = dict1[selected_key]
        df2 = dict2[selected_key]

        if selected_key and selected_column:
            fig, ax = plt.subplots()
            if plot_type == 'scatter' and 'minutes' in df1.columns and 'minutes' in df2.columns:
                ax.scatter(df1['minutes'], df1[selected_column], color='blue', label=f'{n1}')
                ax.scatter(df2['minutes'], df2[selected_column], color='red', label=f'{n2}')
                ax.set_xlabel('minutes')
                ax.set_ylabel(selected_column)
            elif plot_type == 'box':
                pd.concat([df1[[selected_column]], df2[[selected_column]]], axis=1, keys=[f'{n1}', f'{n2}']).boxplot(ax=ax)
            
            ax.legend()
            ax.set_title(f'{plot_type.capitalize()} Plot of {selected_column}')

            # Clear the previous plot in the frame
            for widget in plot_frame.winfo_children():
                widget.destroy()

            # Embed the plot in the tkinter frame
            canvas = FigureCanvasTkAgg(fig, master=plot_frame)
            canvas.draw()
            canvas.get_tk_widget().pack()

    # Button to plot data
    tk.Button(root, text="Plot Data", command=plot_data).grid(row=3, column=1, pady=10)

    # Function to save the plot
    def save_plot():
        initial_file_path = filedialog.asksaveasfilename(initialdir=default_save_path, defaultextension=".png", filetypes=[("PNG files", "*.png"), ("All files", "*.*")])
        if initial_file_path:
            try:
                plt.savefig(initial_file_path)
                messagebox.showinfo("Save Plot", "Plot saved successfully!")
            except Exception as e:
                messagebox.showerror("Save Plot", f"Failed to save plot: {str(e)}")

    # Button to save the plot
    tk.Button(root, text="Save Plot", command=save_plot).grid(row=3, column=2, pady=10)

    # Start the GUI event loop
    root.mainloop()

#default_save_path = r"C:/Users/robbe/OneDrive - KU Leuven/2. Robbe Beirinckx Master thesis/3. Data/excelmodule2"

#plot_Module3(RB003_odhp_p6, RB007_odhp_p6, default_save_path)



def basemodel(x, a, b,c):
    x= np.asarray(x)
    return a*np.exp(b*x) + c
columns_of_interest = ['r_propane', 'r_CO2', 'r_H2', 'r_propylene', 'r_Methane', 'r_Ethane', 'r_Ethene', 'r_CO']

def check_fit_Module3(dict1):
    columns_of_interest = ['r_propane', 'r_CO2', 'r_H2', 'r_propylene', 'r_Methane', 'r_Ethane', 'r_Ethene', 'r_CO']
    df = 'X_S_Y_r_bypass'
    result_dict = {}
    fn = argname('dict1')
    for col_name in columns_of_interest:
        X = dict1[df]['minutes'].values
        Y = dict1[df][col_name].values
        
        valid_mask = np.logical_and(Y != 0, ~np.isnan(Y))
        X_data = X[valid_mask]
        Y_data = Y[valid_mask]
                
        if len(X_data) < 2 or len(Y_data) < 2:
                print(f"Not enough data points for fitting {col_name}")
                continue
        a0 = Y_data[0]- Y_data[-1]        
                # Initial guess for parameters
        p0 = [a0, -0.1, Y_data[-1]]
                
                # Curve fitting
        try:
            popt, pcov = curve_fit(basemodel, X_data, Y_data, p0)
        except RuntimeError:
            print(f"Fit failed for {col_name} ")
            continue
                    
        a_opt, b_opt, c_opt = popt
        perr = np.sqrt(np.diag(pcov))
        fit_r = basemodel(X_data,*popt)
        mean = np.mean(Y_data)
        residuals = Y_data - fit_r
        SSE = np.sum((Y_data - fit_r)**2)
        SST = np.sum((Y_data - mean)**2)

                # Calculate R² value
        R_squared = 1 - SSE / SST
        
        stat, p_value = stats.shapiro(residuals)
        print(f'{col_name}' ,p_value>0.05)
        
        result_dict[col_name] = {
                    'optimal_parameters': popt,
                    'covariance_matrix': pcov,
                    'errors': perr,
                    'SSE': SSE,
                    'R_squared': R_squared,
                    'residuals' : residuals,
                    'residuals_normality': {
                        'normal' : p_value > 0.05 ,
                        'W-statistic': stat,
                        'p-value': p_value
                         # if p > 0.05, the data is considered normal
                         }
                }
                # Generating model data for plotting
        X_model = np.linspace(min(X_data), max(X_data), 100)
        Y_model = basemodel(X_model, a_opt, b_opt, c_opt)
                
                # Plotting
        plt.scatter(X_data, Y_data, label='Data')
        plt.plot(X_model, Y_model, color='r', label='Exponential Fit')
        plt.title(f"{col_name}-{fn} ")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.legend()
        plt.show()

                # Plot covariance matrix
        plt.imshow(pcov)
        plt.colorbar()
        plt.title(f"Covariance Matrix - {col_name}-{fn} ")
        plt.show()
                
        plt.scatter(X_data,residuals)
        plt.title(f"{col_name}-{fn}-residuals")
        plt.xlabel("index")
        plt.ylabel("residuals")
        plt.show()
    dict1['regression'] = result_dict
    return result_dict
def fit_curve(dict1, columns_of_interest, df):
    result_dict = {}
    for col_name in columns_of_interest:
        X = dict1[df]['minutes'].values
        Y = dict1[df][col_name].values
        
        valid_mask = np.logical_and(Y != 0, ~np.isnan(Y))
        X_data = X[valid_mask]
        Y_data = Y[valid_mask]
                
        if len(X_data) < 2 or len(Y_data) < 2:
                print(f"Not enough data points for fitting {col_name}")
                continue
        a0 = Y_data[0]- Y_data[-1]        
                # Initial guess for parameters
        p0 = [a0, -0.1, Y_data[-1]]
                
                # Curve fitting
        try:
            popt, pcov = curve_fit(basemodel, X_data, Y_data, p0)
        except RuntimeError:
            print(f"Fit failed for {col_name} ")
            continue
                    
        a_opt, b_opt, c_opt = popt
        perr = np.sqrt(np.diag(pcov))
        fit_r = basemodel(X_data,*popt)
        mean = np.mean(Y_data)
        residuals = Y_data - fit_r
        SSE = np.sum((Y_data - fit_r)**2)
        SST = np.sum((Y_data - mean)**2)

                # Calculate R² value
        R_squared = 1 - SSE / SST
        
        result_dict[col_name] = {
                    'optimal_parameters': popt,
                    'covariance_matrix': pcov,
                    'errors': perr,
                    'SSE': SSE,
                    'R_squared': R_squared,
                    'residuals' : residuals,
                }
    dict1['regression'] = result_dict
    return result_dict

#RB003r = fit_curve(RB003_odhp_p6, columns_of_interest)
#RB007r =fit_curve(RB007_odhp_p6, columns_of_interest)                              
#RB002r = fit_curve(RB002_odhp_p6, columns_of_interest)
def dif_sig(dic1, dic2):
  name1 = argname('dic1')
  name2 = argname('dic2')
  significant_differences = {}
  non_significant_differences = {}

  for key, value_p1 in dic1['regression'].items():
    value_p2 = dic2['regression'].get(key)
    if value_p2 is None:
      continue

    parameters_p1 = value_p1['optimal_parameters']
    errors_p1 = value_p1['errors']
    parameters_p2 = value_p2['optimal_parameters']
    errors_p2 = value_p2['errors']

    significant = False
    for i in range(len(parameters_p1)):
      min_p1 = parameters_p1[i] - 2*errors_p1[i]
      max_p1 = parameters_p1[i] + 2*errors_p1[i]
      min_p2 = parameters_p2[i] - 2*errors_p2[i]
      max_p2 = parameters_p2[i] + 2*errors_p2[i]

      if (max_p1 < min_p2 or min_p1 > max_p2):
        significant = True
        # Changed to create an inner dictionary to store all significant parameters for the key
        significant_differences.setdefault(key, {})[i] = {f'{name1} par': parameters_p1[i],
                                                          f'{name1} error': errors_p1[i],
                                                          f'{name2} par': parameters_p2[i],
                                                          f'{name2} error' : errors_p2[i]}
      else:
           non_significant_differences.setdefault(key, {})[i] = {
               f'{name1} par': parameters_p1[i],
               f'{name1} error': errors_p1[i],
               f'{name2} par': parameters_p2[i],
               f'{name2} error': errors_p2[i]
           }
        

  return significant_differences, non_significant_differences

#significant_differences, non_significant_differences = dif_sig(RB003_odhp_p6, RB007_odhp_p6)
def curve_comparison_Module3(dictf1, dictf2):
    fn1 = argname('dictf1')
    fn2 = argname('dictf2')
            
    def basemodel(x, a, b,c):
        x= np.asarray(x)
        return a*np.exp(b*x) + c
    byp = input("A = bypass, B = products : ")
    if byp == "A":
        columns_of_interest = ['r_propane', 'r_CO2', 'r_H2', 'r_propylene', 'r_Methane', 'r_Ethane', 'r_Ethene', 'r_CO']
        fit_curve(dictf1, columns_of_interest, 'X_S_Y_r_bypass')
        fit_curve(dictf2, columns_of_interest, 'X_S_Y_r_bypass')
    else:
        columns_of_interest = ['r_propanep', 'r_CO2p', 'r_H2p', 'r_propylenep', 'r_Methanep', 'r_Ethanep', 'r_Ethenep', 'r_COp']
        fit_curve(dictf1, columns_of_interest, 'X_S_Y_r_products')
        fit_curve(dictf2, columns_of_interest, 'X_S_Y_r_products')
    significant_differences, non_significant_differences = dif_sig(dictf1,dictf2)
    return significant_differences, non_significant_differences
#significant_differences, non_significant_differences = curve_comparison_Module3(RB003_odhp_p6, RB007_odhp_p6)    
