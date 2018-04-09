import numpy as np
import time
from math import sqrt, pi
import os, sys
from subprocess import Popen, PIPE, STDOUT

cwd = os.getcwd()

from isotope import isotope
from ranges import *
from maxwellian import *

def data_parse(Z,A,reaction_type):

    #  use imported isotope function to get string of the isotope which corresponds to our Z,A

    iso = isotope(Z,A,reaction_type)[0]

    #  define empty dictionary which will eventually hold the output data (note this will be structured in the same way as the dictionary holding the freya output
    output_dict = {}

    print('Begin Parsing Data...')
    data_begin = time.time()

    #  assign the variable energy to be a string which can be read by the gen_par_ana script so it produces the correct events for comparison
    if reaction_type is 'induced':
        if int(A) is not 240 or 242:
            print('This isotope is not supported by Freya, or has no available data.')
        else:
            energy = 'induced'
    else:
        energy = -1
        data_parsing_path = cwd + "/../data_master/" + str(iso) + "/"

    #  change the working directory to be the folder containing all of the data for this isotope
    os.chdir(cwd+'/../data_master/'+iso)
   

    #  set files to be a list of every file in the appropriate directory
    files = os.listdir(data_parsing_path)


    #  since there will be multiple data sets for each observable, this empty dictionary will bring the name of each data set, to the corresponding general name of the observable
    #  for example, the key 'n_mult_0' will correspond to the string 'n_mult'
    #  this dictionary is used in the error calculating script so for a given data array, it can find the corresponding freya array by running the string through this dictionary
    #  it is simply easiest to define the dictionary as we go, so we know there is an entry for the name of each data array
    key_translator = {}

    #  run through every file in the list of files as defined above 
    for datafile in files:
        #  only process files with the appropriate beginning matching the correct isotope (in case there are irrelevant files present) 
        file_name_words = datafile.split('.')
        if file_name_words[0] == iso:
            #  assign the variable is3d to determine whether the data will be plotted on 3-d axes or 2d axes. 
            if file_name_words[1] == 'n_A_TKE':
                is3d = True
            else:
                is3d = False

            #  some of the data sets contain the variance rather than the actual values. To these sets, we attach the boolean value of 'sigma'
            #  this indication will later be used to attach this array to the appropriate data values
            #  this will eventually result in an array with the actual data, the uncertainty, and the variances 
            if file_name_words[3] == 'sigma':
                sigma = True
            else:
                sigma = False

            #  write element of the key translator for the current data set
            data_key = file_name_words[1]
            print('Parsing:',data_key)
            key_translator[data_key] = file_name_words[1]
            key_translator[file_name_words[1]] = file_name_words[1]

            #  open the actual data file, parse the lines into the list 'lines'
            data_file_path = data_parsing_path + datafile
            infile = open(data_file_path,"r+") 
            lines = infile.readlines()

            #  define empty array
            #  length corresponds to number of "bins"
            #  width corresponds to number of things to be binned (values, variances)
            #  depth corresponds to the uncertainty on these values (note variances sometimes have uncertainties as well)
            data = np.zeros((len(lines)-2,4,2))

            citation_line = lines[0]
            citation = citation_line[6:-2]

            if is3d == True:
                labels = lines[1].split(':')
                case = len(labels)
                #  first we deal with the 3-d case, where we are given two binnings at once (n_A_TKE for example)
                #  case corresponds to how many columns there are. This will determine what data is in what columns. 
                #  3 case: bins 1 - bins 2 - data
                #  4 case: bins 1 - bins 2 - data - uncertainty
                #  5 case: bins 1 1 - bins 1 2 - bins 2 - data - uncertainty
                #  6 case: bins 1 1 - bins 1 2 - bins 2 1 - bins 2 2 - data - uncertainty

                if case == 3 or case == 4:
                    x_label = labels[0]
                    y_label = labels[1]
                    z_label = labels[2]
                if case == 5:
                    x_label = labels[0]
                    y_label = labels[2]
                    z_label = labels[3]
                if case == 6:
                    x_label = labels[0]
                    y_label = labels[2]
                    z_label = labels[4]
                for line_number in range(0, len(lines)-2):
                    line = lines[line_number+2]
                    elems = line.split()
                    if len(elems) == 0:
                        continue
                    if len(elems) == 3:
                        data[line_number, 0, 0] = elems[0]
                        data[line_number, 1, 0] = elems[1]
                        data[line_number, 2, 0] = elems[2]
                    elif len(elems) == 4:
                        data[line_number, 0, 0] = elems[0]
                        data[line_number, 1, 0] = elems[1]
                        data[line_number, 2, 0] = elems[2]
                        data[line_number, 2, 1] = elems[3]
                    elif len(elems) == 5:
                        data[line_number, 0, 0] = float(elems[0]) + (float(elems[1]) - float(elems[0])) /2
                        data[line_number, 1, 0] = elems[2]
                        data[line_number, 2, 0] = elems[3]
                        data[line_number, 0, 1] = (float(elems[1]) - float(elems[0])) / float(elems[0])
                        data[line_number, 2, 1] = elems[4]
                    elif len(elems) == 6 and float(elems[5])>1:
                        data[line_number, 0, 0] = float(elems[0]) + (float(elems[1]) - float(elems[0]))/2
                        data[line_number, 1, 0] = float(elems[2]) + (float(elems[3]) - float(elems[2]))/2
                        data[line_number, 2, 0] = elems[4]
                        data[line_number, 0, 1] = (float(elems[1]) - float(elems[0])) / float(elems[0])
                        data[line_number, 1, 1] = (float(elems[3]) - float(elems[2])) / float(elems[2])
                        data[line_number, 2, 1] = float(elems[5])/100
                    elif len(elems) == 6:
                        data[line_number, 0, 0] = float(elems[0]) + (float(elems[1]) - float(elems[0]))/2
                        data[line_number, 1, 0] = float(elems[2]) + (float(elems[3]) - float(elems[2]))/2
                        data[line_number, 2, 0] = elems[4]
                        data[line_number, 0, 1] = (float(elems[1]) - float(elems[0])) / float(elems[0])
                        data[line_number, 1, 1] = (float(elems[3]) - float(elems[2])) / float(elems[2])
                        data[line_number, 2, 1] = elems[5]
                    else:
                        print('Data failed to parse for: ' + str(datafile))
            else:
                #  now we parse the non 3-d case (non n_A_TKE case)
                #  case 0: do nothing
                #  case 2: bins - data - 0
                #  case 3: bins - data - uncertainty

                labels = lines[1].split(':')
                x_label = labels[0]
                y_label = labels[1]
                z_label = None
                case = None
                for line_number in range(0, len(lines)-2):
                    line = lines[line_number + 2]
                    elems = line.split()
                    if len(elems) == 0:
                        continue
                    elif len(elems) == 2:
                        data[line_number, 0, 0] = elems[0]
                        data[line_number, 1, 0] = elems[1]
                    elif len(elems) == 3:
                        data[line_number, 0, 0] = elems[0]
                        data[line_number, 1, 0] = elems[1]
                        data[line_number, 1, 1] = elems[2] 
                    else:
                        print('Data failed to parse for: ' + str(datafile))
                data = np.array(data, dtype = np.float)


            if file_name_words[1] == 'mannhart':

                scale = 'log'

            elif file_name_words[1] == 'n_spectrum':
                scale = 'log'
            else:
                scale = 'linear'
            if sigma is True:
                data[:,2,0] = np.sqrt(data[:,2,0])

            #  the value of the output dictionary as defined above for each data key (name of data file)
            #  will be as follows:

            output_dict[str(data_key)] = [data,citation,
                    x_label,
                    y_label,
                    str(energy),
                    is3d,
                    sigma,
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

    return output_dict, key_translator, energy

