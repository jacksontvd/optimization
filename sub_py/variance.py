import time
from math import sqrt, pi
import os, sys, copy
import numpy as np
from scipy import optimize
from scipy import integrate

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

def probability(chi_sq_array,number):
    dof = len(chi_sq_array) - 5
    return (chi_sq_array)**(dof/2 - 1)*np.exp(-chi_sq_array/(2*number))

def variance(Z,A, generate_number = None, method = None, resolution = None, **kwargs):
    print('starting variance calculation')
    reac_t = kwargs['reaction_type']
    parameter = kwargs['parameter']
    parameters = param_ranges[str(Z)+str(A)]

    range_array = param_ranges[parameter]
    special_index = param_ranges[parameter][2]

    fixed_value = copy.copy(parameters[special_index])

    os.chdir(cwd+'/../fission_v2.0.3/data_freya/')
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

    #  define function which takes in a set of parameters and returns the weighted chi-squared error (total sum)
    #  see error.py for the details of this error calculation
    def err_opt(parameter):
        parameters[special_index] = parameter
        return error(Z, A, parameters[0],parameters[1], parameters[2], parameters[3], parameters[4], generate_number, parsed_data, reaction_type = reac_t)[0]
        #  return test_error(Z, A, parameters[0],parameters[1], parameters[2], parameters[3], parameters[4], generate_number, parsed_data, reaction_type = reac_t)[0]

    #  for grid search method initialize the brute source routine
    print("calculating errors on nodes...")

    resolution = np.float(resolution)
    #  define the ranges to be from the minimum to the maximum values, with the difference/resolution many cells
    #  brute_ranges = (slice(range_array[0],range_array[1],range_array[2]/resolution))
    brute_ranges = ((range_array[0],range_array[1]),)

    #  block the printing for each iteration to avoid slowing the process down by taking the time to print the useless output
    #  comment this line to let it print if there is an issue with the routine
    blockPrint()
    x0, fval, grid, Jout = optimize.brute(err_opt , brute_ranges, full_output = True,Ns = resolution, finish = optimize.fmin)
    enablePrint()

    #  set finalparams to be the 0th output of the brute routine
    finalparams = x0
    #  set grid_values to be the final element of the output list of the brute routine
    grid_values = Jout

    if grid_values is not None:
        print("Done.")
    print("optimized parameters:",finalparams)

    ### stop optimizing

    #  define path for output to be placed
    var_path = cwd + '/../output/variance/' + 'Z=' + str(Z) + ' A=' + str(A) + '_E=' + str(Energy) + '_' + str(method)

    #  create appropriate directory if it does not already exist 
    if not os.path.exists(var_path):
        os.makedirs(var_path)

    opt_end = time.time()
    print("time:",opt_end - opt_begin)

    print('Begin printing output...')
    os.chdir(var_path)
    np.savetxt('variance',grid_values)
    os.chdir(cwd)

    print("Finished.")
    print("Statistics generated in: \n"+var_path)
    print("Calculating variance...")

    average_error = sum(grid_values) / len(grid_values)
    big_number = average_error
    probability_array = probability(grid_values,big_number)

    def normalizer_function(x):
        index = np.searchsorted(grid,x)
        probability = probability_array[index-1]
        return probability 

    normalizer = integrate.simps(probability_array,grid)
    probability_array = np.divide(probability_array,normalizer)

    expect_array = np.multiply(probability_array,grid)
    expect2_array = np.multiply(probability_array,np.square(grid))

    average = integrate.simps(expect_array, grid)
    print("average",average)

    average2  = integrate.simps(expect2_array, grid)
    print("average2",average2)

    variance = average2 - average**2

    print("variance:",variance)
    print("events: ",generate_number)
    print(parameter)
    print("resolution:",resolution)
    print(Z,A)

    return grid_values, variance, average
