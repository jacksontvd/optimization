import time
from math import sqrt, pi
import os, sys
import numpy as np
from scipy import optimize
from scipy import integrate
from scipy.integrate import simps

cwd = os.getcwd()

sys.path.append(cwd+'/../data_master/Cf252/')
from mannhart_data import mannhart_bins, mannhart_split, mannhart_bindiff
sys.path.append(cwd)

from ranges import * 
from gen_par_ana import gpa
from error import error
from data_parse import data_parse
from isotope import isotope

#  from test import *

#  define functions to block and restore printing for when freya is run many times in a row during the optimization procedure
def blockPrint():
    sys.stdout = open(os.devnull, 'w')
def enablePrint():
    sys.stdout = sys.__stdout__

def covariance(Z,A, generate_number = None, method = None, resolution = None, **kwargs):
    print('starting covariance calculation')
    reac_t = kwargs['reaction_type']
    parameter = kwargs['parameter']
    parameter_2 = kwargs['parameter2']
    parameters = param_ranges[str(Z)+str(A)]

    range_array = param_ranges[parameter]
    special_index = param_ranges[parameter][2]
    range_array_2 = param_ranges[parameter_2]
    special_index_2 = param_ranges[parameter_2][2]

    fixed_value_1 = parameters[special_index]
    fixed_value_2 = parameters[special_index_2]

    resolution = np.float(resolution)

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
    #  define function which takes in a set of parameters and returns the raw chi-squared error (total sum)
    #  see error.py for the details of this error calculation
    def err_opt(params):
        parameter = params[0]
        parameter_2 = params[1]
        parameters[special_index] = parameter
        parameters[special_index_2] = parameter_2
        return error(Z, A, parameters[0],parameters[1], parameters[2], parameters[3], parameters[4], generate_number, parsed_data, reaction_type = reac_t)[0]
        #  return test_error(Z, A, parameters[0],parameters[1], parameters[2], parameters[3], parameters[4], generate_number, parsed_data, reaction_type = reac_t)[0]


    ### start optimizing

    print("calculating errors on nodes...")

    #  define the ranges to be from th eminimum to the maximum values, with the difference/resolution many cells

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
    #  set grid_values to be the final element of the output list of the brute routine
    grid_values = Jout

    if grid_values is not None:
        print("Done.")
    print("optimized parameters:",finalparams)

    ### stop optimizing

    #  define path for output to be placed
    var_path = cwd + '/../output/covariance/' + str(Z) + str(A) + str(Energy)

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
    print("Calculating covariance...")

    grid_rows = sum(grid_values)
    grid_total = sum(grid_rows)
    average_error = grid_total / len(grid_values)**2
    big_number = average_error

    prob_array = np.exp(-grid_values / big_number)

    rrange = np.arange(range_array[0] , range_array[1] , (range_array[1] - range_array[0])/resolution)
    rrange_2 = np.arange(range_array_2[0] , range_array_2[1] , (range_array_2[1] - range_array_2[0])/resolution)

    fixed_index_1 = np.searchsorted(rrange,fixed_value_1)
    fixed_index_2 = np.searchsorted(rrange_2,fixed_value_2)

    prob_array_1 = grid_values[:,fixed_index_2-1]
    prob_array_1 = np.exp(-prob_array_1 / np.mean(prob_array_1)).flatten()

    prob_array_2 = grid_values[fixed_index_1-1,:]
    prob_array_2 = np.exp(-prob_array_2 / np.mean(prob_array_2)).flatten()

    ###
    prob_array = np.add(prob_array,1)
    prob_array_1 = np.add(prob_array_1,1)
    prob_array_2 = np.add(prob_array_2,1)
    ###

    expect_array_1 = np.multiply(rrange,prob_array_1)
    expect_array_2 = np.multiply(rrange_2,prob_array_2)

    def normalizer_function(x,y):
        index = np.searchsorted(rrange,x)
        index2 = np.searchsorted(rrange_2,y)
        probability = prob_array[index - 1,index2 - 1]
        return probability 

    objective = normalizer_function(rrange , rrange_2[:,None])
    normalizer = simps(simps(objective,rrange),rrange_2)
    prob_array = np.divide(prob_array , normalizer)

    normalizer_1 = simps(prob_array_1 , rrange)
    expect_array_1 = np.divide(expect_array_1 , normalizer_1)

    normalizer_2 = simps(prob_array_2 , rrange_2)
    expect_array_2 = np.divide(expect_array_2 , normalizer_2)

    def pintegrand(x,y):
        index = np.searchsorted(rrange,x)
        index2 = np.searchsorted(rrange_2,y)
        probability = prob_array[index-1,index2 - 1]
        return probability * x * y

    objective = pintegrand(rrange,rrange_2[:,None])
    average = simps(simps(objective,rrange),rrange_2)
    print("average",average)

    average1 = simps(expect_array_1,rrange)
    print(kwargs['parameter'],"average: ",average1)

    average2 = simps(expect_array_2,rrange_2)
    print(kwargs['parameter2'],"average: ",average2)

    covariance = average - (average1 * average2)
    print("covariance:",covariance)

    print("events:",generate_number)
    print("resolution:",resolution)
    print(parameter , parameter_2)
    print(Z,A)

    return covariance , grid_values, rrange , rrange_2
