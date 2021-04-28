#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: danieleldridge

Spring 2020

RPyFR Script: A script for wrangling vibrational frequencies from Gaussian output files and computing RPFRs
over user-specified temperature ranges.

For this script to work, you must follow the following file naming format for freqchk output files (.txt):

    MoleculeID_Theory_BasisSet_IsotopicMoleculeID_note.txt

Examples:
    1. Methane_CCSD_augccpVTZ_0nosub_atom2.txt
    2. Methane_CCSD_augccpVTZ_12CH3D_atom2.txt
    3. Sulfate40H2O_B3LYP_631Gdp_33SO4_atom1.txt
    4. Thiosulfate34H2O_B3LYP_631Gdp_32S32S(16O)3_none.txt

The unsubstituted molecule must have the IsotopicMoleculeID of "0nosub" in order for the RPFRs to be
computed correctly.

Avoid special characters for basis set (and other) descriptions and use shorthand whenever possible.

Also, only place files pertaining to a single optimization/frequency calculation in the data_files folder at one time for processing.
"""

import os as os
import shutil
import pandas as pd
import numpy as np

def grab_freq(file):
    """
    This function takes a Gaussian output file from a frequency calculation (.txt, from freqchk utility) and 
    extracts the harmonic vibrational frequencies as a list item of floats.
    Note: Harmonic frequencies from calculations performed with anharmonic corrections (freq=anharmonic) 
          will be printed twice from this function because such output files list these frequencies twice.
          Thus, this function works best on frequency calculations performed without anharmonic analysis.
    """
    #read Gaussian output file
    directory = './data_files'
    with open(directory+'/' + file, "r") as f:
        file = f.read()
    #transform text in Gaussian output file into a list based on lines ('\n'):
    line_list = file.split('\n')
    #create new list comprised only of the lines containing Frequencies:
    frequency_lines_only = []
    for line in line_list:
        if line.startswith(' Freq'):
            frequency_lines_only.append(line)
    #Split each line (list item) into lists:
    frequency_lines_only_split = [line.split() for line in frequency_lines_only]
    #transform resultant list of lists into a single list:
    frequency_splice2 = []
    for x in frequency_lines_only_split:
            x = x
            frequency_splice2.extend(x)
    #Remove items from list named "Frequencies" using str.isalpha()
    #NOTE: The approach below is taken because frequency values are strings up to this point
    frequency_splice1 = []
    for y in frequency_splice2:
        if y.isalpha()==False:      
            frequency_splice1.append(y)
    #Remove '--' items to yield a list containing only frequencies
    frequency_splice = []
    for z in frequency_splice1:
        if z != '--':
            frequency_splice.append(z)
    #convert strings to floats for later computations
    frequencies =[]
    for i in frequency_splice:
        i = float(i)
        frequencies.append(i)
    #Optional printed statement that states how many frequencies were extracted.
    #print('Total number of frequencies =',len(frequencies),'; Theory =', level_theory, '; Basis Set=', basis_set)
    return(frequencies)

def process_file(file):
    """Calls upon the grab_freq function for a file and returns
    a list containing: 
    (isotopologue_name, freqs) in
    (str, list) format.
    """
    file_name = str(file.split('/')[-1][:-4])
    isotopologue_name = str(file_name.split('_')[-2]) #filename must follow correct format for this to work
    #level_theory = str(file_name.split('_')[1])
    #basis_set = str(file_name.split('_')[2])
    freqs = grab_freq(file)
    return(isotopologue_name, freqs)

#RUN FILE PROCESSING SCRIPT, GENERATING CSV FILE WITH ORDERED and LABELED FREQUENCIES
iso_output=[]
directory = './data_files'
for file in os.listdir(directory):
    if file.endswith(".txt"):
        print("Now processing: ", file)
        curr_iso_output = process_file(file)
        iso_output.append(curr_iso_output)
        file_id = file.split('_')[0]+"_"+file.split('_')[1]+"_"+file.split('_')[2] #format: MoleculeID_Theory_BasisSet

dict_output = {}
for header,freqs in iso_output:
    dict_output[header]=freqs
data_set = pd.DataFrame(dict_output)
data_set = data_set.sort_index(axis=1)
#createFolder('./'+file_id+'')
data_set.to_csv(directory+'/'+file_id+'_FREQUENCIES.csv',index='')
print('\n Harmonic vibrational frequencies: \n', data_set)

calc_RPFR_Q = input("Calculate RPFRs using these frequencies? (y/n) ")
if calc_RPFR_Q == 'y':
    #CALCULATE RPFRs AT USER-SPECIFIED TEMPERATURE INTERVAL AND INCREMENT
    # Define required constants:
    h = float(6.6260700400000000e-34)   # Planck constant (m^2*kg/s)
    k_b = float(1.3806485200000000e-23)   # Boltzmann constant (m^2*kg*s^-2*K-1)
    c = float(2.99792458000000e+10)      # Speed of light  (cm/s)   (note: wavenumbers are always given in cm^-1)
    # Set the temperature step, min and max over which to calculate:
    stepT = float(input("Enter the desired temperature step in ˚C: "))
    minT = 273.15 + float(input("Enter the desired minimum temperature in ˚C: "))
    maxT = 273.15 + float(input("Enter the desired maximum temperature in ˚C: "))+stepT
    temparray = np.arange(minT,maxT,stepT)
    l = len(temparray)
    # Create array of zeroes for the final output
    (i,j) = np.shape(data_set)
    RPFR_matrix = np.zeros(shape=(l,j+1))
    # Set up column name referencing
    names = data_set.columns
    indexes=['T (C)','T (K)'] #first two column names
    for x in range(1,j):
        indexes=np.append(indexes,names[x])
    # Nested for loops to generate output array
    for x in range(1,j):
        for T in range(0,l):
            RPFR_matrix[T,0]=np.asarray(temparray[T]-273.15)
            RPFR_matrix[T,1]=np.asarray(temparray[T])
            u = (h*c*data_set['0nosub'])/(k_b*temparray[T])
            rpf=u*np.exp(-u/2)*(1/(1-np.exp(-u)))
            u_star1=(h*c*data_set[str(names[x])])/(k_b*temparray[T])
            rpf_star1=u_star1*np.exp(-u_star1/2)*(1/(1-np.exp(-u_star1)))
            rpfr_step1 = rpf_star1/rpf #important intermediate step that avoids div. by zero errors for molecules that have lots of frequencies
            RPFR_matrix[T,x+1]=np.prod(rpfr_step1) #np.product(rpfr_step1
    # Export output to dataframe, then to csv
    RPFR_df=pd.DataFrame(RPFR_matrix,columns=indexes)
    RPFR_df.to_csv(directory+'/'+file_id+'_outputRPFRs.csv', index='', header=True, sep=',')

os.mkdir(directory+'/'+'processed_'+file_id)
processed_directory=directory+'/'+'processed_'+file_id
for file in os.listdir(directory):
    if file.startswith(file_id):
        print("Now moving: ", file)
        shutil.move(directory+'/'+file, processed_directory)