import numpy as np
import time
from math import sqrt, pi
import os, sys
from subprocess import Popen, PIPE, STDOUT

cwd = os.getcwd()

from isotope import isotope
from ranges import *
from maxwellian import *

def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def induced_data_parse(Z,A,reaction_type,directory_name,energy):
    #  The variable directory name tells us which directory (inside cwd + "/../data_master/" + str(iso) + "/") contains the data (for varying incident energies) for the observable of interest.
    #  This has the format "observable.author" so we split it accordingly here
    observable = directory_name.split('.')[0]
    citation = directory_name.split('.')[1]

    iso = isotope(Z,A,reaction_type)[0]

    #  define empty dictionary which will eventually hold the output data (note this will be structured in the same way as the dictionary holding the freya output
    output_dict = {}

    print('Begin Parsing induced Data...')
    data_begin = time.time()

    #  change the working directory to be the folder containing all of the data for this isotope

    #  set files to be a list of every file in the appropriate directory
    data_parsing_path = cwd + "/../data_master/" + str(iso) + "/" + iso+"."+directory_name
    os.chdir(data_parsing_path)
    files = os.listdir(data_parsing_path)

    #  since there will be multiple data sets for each observable, this empty dictionary will bring the name of each data set, to the corresponding general name of the observable
    #  for example, the key 'n_mult_0' will correspond to the string 'n_mult'
    #  this dictionary is used in the error calculating script so for a given data array, it can find the corresponding freya array by running the string through this dictionary
    #  it is simply easiest to define the dictionary as we go, so we know there is an entry for the name of each data array
    key_translator = {}


    #  parse the average TKE differently (since it is all in one file)
    if observable == 'TKE_bar':
        for datafile in files:
            chosenfile = datafile
        datafile = chosenfile
    else:
    #  run through every file in the list of files as defined above 
    #  set the variable 'wefoundit' to be some 0. If a file with data at the appropriate incident neutron energy is found, this variable will be set to be 1 instead.
        wefoundit = 0
        #  We will range over the different files, each of which holds data for a different incident neutron energy.
        #  Our initial threshold is $1$. I.e. if there is no file for an energy with 1 MeV of the incident energy of the FREYA events, then no data will be compared do. This threshold must change, so if there are multiple files within 1 MeV of the energy of the FREYA events, the closest one will be chosen.
        difference_threshold = 1
        for datafile in files:
            file_energy = datafile.split('_')[0]
            if isfloat(file_energy) is True:
                    if abs(float(file_energy) - float(energy)) < difference_threshold:
                        chosenfile = datafile
                        difference_threshold = abs(float(file_energy) - float(energy))
                        wefoundit=1
        if wefoundit == 0:
            chosenfile="dummy"
            print("NO DATA FILES MATCH THE INCIDENT NEUTRON ENERGY")
        datafile = chosenfile

    print("Parsing:",datafile)
    #  only process files with the appropriate beginning matching the correct isotope (in case there are irrelevant files present) 
    file_name_words = datafile.split()

    #  some of the data sets contain the variance rather than the actual values. To these sets, we attach the boolean value of 'sigma'
    #  this indication will later be used to attach this array to the appropriate data values
    #  this will eventually result in an array with the actual data, the uncertainty, and the variances 

    #  write element of the key translator for the current data set
    #  data_key = file_name_words[1]
    data_key=observable
    key_translator[data_key] = observable

    #  open the actual data file, parse the lines into the list 'lines'
    data_file_path = data_parsing_path + '/' + datafile
    infile = open(data_file_path,"r+") 
    lines = infile.readlines()

    #  define empty array
    #  length corresponds to number of "bins"
    #  width corresponds to number of things to be binned
    #  depth corresponds to the uncertainty on these values (there is sometimes "x uncertainty" as well)
    data = np.zeros((len(lines)-1,4,2))

    #  now we parse the non 3-d case (non n_A_TKE case)
    #  case 0: do nothing
    #  case 2: bins - data - 0
    #  case 3: bins - data - uncertainty

    labels = lines[0].split()
    x_label = labels[0]
    y_label = labels[1]
    z_label = None
    case = None
    for line_number in range(0, len(lines)-1):
        line = lines[line_number + 1]
        elems = line.split()
        if len(elems) == 0:
            continue
        elif len(elems) == 4:
            data[line_number, 0, 0] = elems[0]
            data[line_number, 1, 0] = elems[1]
            data[line_number, 1, 1] = elems[2] 
        elif len(elems) == 3:
            data[line_number, 0, 0] = elems[0]
            data[line_number, 1, 0] = elems[1]
            data[line_number, 1, 1] = elems[2]
        elif len(elems) == 2:
            data[line_number, 0, 0] = elems[0]
            data[line_number, 1, 0] = elems[1]
        else:
            print('Data failed to parse for: ' + str(datafile))

    data = np.array(data, dtype = np.float)

    if observable == 'TKE_bar':
        for row_number in range(0,len(data[:,0,:])):
            this_energy = data[row_number,0,0]
            if abs(float(this_energy) - float(energy)) < 1:
                chosen_data = data[row_number,1]

    #  we want to plot some things with different scales on the axes, here we define a variable for this
    if observable == 'n_spectrum':
        scale = 'log'
    else:
        scale = 'linear'

    #  if observable == 'n_spectrum':
    #      print(data)
    #      bin_differences = np.subtract(data[1:,0,0],data[:-1,0,0])
    #      areas = bin_differences * data[:-1,1,0]
    #      total_count = np.sum(areas) + bin_differences[-1] * data[-1,1,0] + 0.5
    #      print(total_count)
    #      data[:,1,0] = data[:,1,0] / total_count
    #      data[:,1,1] = data[:,1,1] / total_count
    #      print(data)


    #  the value of the output dictionary as defined above for each data key (name of data file)
    #  will be as follows:
    output_dict[str(data_key)] = [data,citation,
            x_label,
            y_label,
            str(energy),
            False,
            True,
            scale,
            z_label,
            False] # sigma indicator

    #  run through each element of the data_dictionary, and find the ones which contain only the variances for that observable
    #  these correspond to another element of the dictionary which contains the actual data, for which this original element contains the variance of
    #  we therefore combine these into one element, by writing the variance into the [:,2,0] column of the data array
    #  again, this is because the first column is the binning, the second is the data, and the third is the variance. The depth is the uncertainty. 
    #  this can be easily remembered if one recalls that the freya arrays had no additional depth (because they of course have no experimental error)
    for key in output_dict:
        key_words = key.split('.')
        if len(key_words) == 3:
            if key_words[2] == 'sigma':
                main_data_key = str( key_words[0] ) + '.' + str(key_words[1] )
                print('Sigmas transferred for:',main_data_key)
                main_data_entry = output_dict[str(main_data_key) ] 
                main_data_entry[0][:,2,0] = output_dict[key][0][:,1,0]
                main_data_entry[0][:,2,1] = output_dict[key][0][:,1,1]
                main_data_entry[9] = True
                output_dict[key] = None
        else:
            continue

    data_end = time.time()
    data_time = data_end - data_begin
    print('Time:', data_time)

    for element in output_dict:
        if output_dict[element] is None:
            print('Error in data parsing for this isotope!')

    return output_dict[str(observable)]
