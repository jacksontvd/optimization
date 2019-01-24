import time
from math import sqrt, pi
import os, sys
import numpy as np
from scipy import optimize
from scipy import integrate
from scipy.integrate import simps

cwd = os.getcwd()

from ranges import * 
from gen_par_ana import gpa
from error import error
from data_parse import data_parse
from isotope import isotope

from test import *

#  define functions to block and restore printing for when freya is run many times in a row during the optimization procedure
def blockPrint():
    sys.stdout = open(os.devnull, 'w')
def enablePrint():
    sys.stdout = sys.__stdout__

def grid_run(Z,A, generate_number = None, method = None, resolution = None, **kwargs):
    print('starting grid calculation')
    #  pull the reaction type, parameters, and list of parameters
    reac_t = kwargs['reaction_type']
    parameter = kwargs['parameter']
    parameter_2 = kwargs['parameter2']
    parameters = param_list(Z,A,reac_t)

    #  pull array of ranges of parameters and index of the changing parameters
    range_array = param_ranges[parameter]
    special_index = param_ranges[parameter][2]
    range_array_2 = param_ranges[parameter_2]
    special_index_2 = param_ranges[parameter_2][2]

    #  pull array of fixed value of the two changing parameters
    fixed_value_1 = parameters[special_index]
    fixed_value_2 = parameters[special_index_2]

    resolution = np.float(resolution)

    os.chdir(cwd+'/../../freya/data_freya/')

    infile = open("inputparameters.dat","r+")
   
    #  parse appropriate data with fuction from data_parse.py
    parsed_data = data_parse(Z,A,reac_t)
    #  pull the data array out of the data_parse output list
    data_array = parsed_data[0]
    #  pull the data key dictionary from the data_parse output list
    key_translator = parsed_data[1]
    #  pull energy from the data_parse so the freya output matches the energy level of the comparison data
    Energy = parsed_data[2]
    #  pull the line number for the isotope from the isotope function
    #  this is used to rewrite the appropriate line of the parameter file in each iteration of the optimization
    #  it is also used to read off the guess values from the parameter file as it is currently
    iso = isotope(Z,A,reac_t)
    i = iso[1]
    content = infile.readlines() #reads line by line and outputs a list of each line
    line_split = content[i].split()

    opt_begin = time.time()

    print('Begin Optimizing FREYA (This may take a while...)')
    #  define function which takes in a set of parameters and returns the raw chi-squared error (total sum)
    #  see error.py for the details of this error calculation
    def err_opt(params):
        parameter = params[0]
        parameter_2 = params[1]
        parameters[special_index] = parameter
        parameters[special_index_2] = parameter_2
        return error(Z, A, parameters[0],parameters[1], parameters[2], parameters[3], parameters[4], generate_number, parsed_data, reaction_type = reac_t)[0]
        #  return test_error(Z, A, parameters[0],parameters[1], parameters[2], parameters[3], parameters[4], generate_number, parsed_data, reaction_type = reac_t)[0]

    print("calculating errors on nodes...")

    #  define the ranges to be from the minimum to the maximum values, with the difference/resolution many cells

    brute_ranges = (
                    slice(range_array[0],range_array[1],(range_array[1] - range_array[0])/resolution),
                    slice(range_array_2[0] , range_array_2[1] , (range_array_2[1] - range_array_2[0]) / resolution)
                    )
    #  brute_ranges = ((range_array[0],range_array[1]),(range_array_2[0],range_array_2[1]))

    #  block the printing for each iteration to avoid slowing the process down by taking the time to print the useless output
    #  comment this line to let it print if there is an issue with the routine

    blockPrint()
    #  x0, fval, grid, Jout = optimize.brute(err_opt , brute_ranges, full_output = True,Ns = resolution, finish = optimize.fmin)
    x0, fval, grid, Jout = optimize.brute(err_opt , brute_ranges, full_output = True, finish = optimize.fmin)
    enablePrint()

    #  set finalparams to be the 0th output of the brute routine
    finalparams = x0
    error_array = error(Z, A, parameters[0],parameters[1], parameters[2], parameters[3], parameters[4], generate_number, parsed_data, reaction_type = reac_t)
    #  error_array = test_error(Z, A, parameters[0],parameters[1], parameters[2], parameters[3], parameters[4], generate_number, parsed_data, reaction_type = reac_t)
    dof = error_array[6]
    #  set grid_values to be the final element of the output list of the brute routine
    grid_values = Jout

    if grid_values is not None:
        print("Done.")
    print("optimized parameters:",finalparams)

    ### stop optimizing

    #  define path for output to be placed
    var_path = cwd + '/../../output/grid_values/' + str(Z) + str(A) + str(Energy)

    #  create appropriate directory if it does not already exist 
    if not os.path.exists(var_path):
        os.makedirs(var_path)

    opt_end = time.time()
    print("time:",opt_end - opt_begin)

    print('Begin printing output...')
    os.chdir(var_path)
    np.savetxt(parameter + parameter_2 , grid_values)
    os.chdir(cwd)

    print("Finished.")
    print("Statistics generated in: \n"+var_path)
    return
