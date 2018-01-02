import time
from math import sqrt, pi
import os, sys
import numpy as np
from scipy import optimize
from scipy import integrate

cwd = os.getcwd()

sys.path.append(cwd+'/../data_master/Cf252/')
from mannhart_data import mannhart_bins, mannhart_split, mannhart_bindiff
sys.path.append(cwd)

from ranges import * 
from gen_par_ana import gpa
from error import error
from data_parse import data_parse
from isotope import isotope

#  define functions to block and restore printing for when freya is run many times in a row during the optimization procedure
def blockPrint():
    sys.stdout = open(os.devnull, 'w')
def enablePrint():
    sys.stdout = sys.__stdout__


def variance(Z,A, generate_number = None, method = None, resolution = None, **kwargs):
    print('starting variance calculation')
    reac_t = kwargs['reaction_type']
    parameter = kwargs['parameter']
    parameters = param_ranges[str(Z)+str(A)]
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

    #  make a zero list which will eventually contain the guess parameters from the parameter file
    #  this is a list but not an array because the optimization functions prefer lists

    opt_begin = time.time()

    print('Begin Optimizing FREYA (This may take a while...)')

    range_array = param_ranges[parameter]
    special_index = param_ranges[parameter][2]

    #  define function which takes in a set of parameters and returns the raw chi-squared error (total sum)
    #  see error.py for the details of this error calculation
    def err_opt(parameter):
        #  return error(Z, A, parameters[0],parameters[1], parameters[2], parameters[3], parameters[4], generate_number, parsed_data, reaction_type = reac_t)[0]
        parameters[special_index] = parameter
        return error(Z, A, parameters[0],parameters[1], parameters[2], parameters[3], parameters[4], generate_number, parsed_data, reaction_type = reac_t)[0]

    #  for grid search method initialize the brute source routine
    if method == 'grid': 
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

    def normalizer_function(x):
        rrange = np.arange(range_array[0] , range_array[1] , (range_array[1] - range_array[0])/resolution)
        index = np.searchsorted(rrange,x)
        error = grid_values[index-1]
        probability = np.exp(-error/10000)
        #  print("ranges:",rrange)
        #  print("value:",x)
        #  print("index:",index)
        #  print("error:",error)
        #  print("probability",probability)
        return probability 

    normalizer,int_error = integrate.quad(normalizer_function, range_array[0] , range_array[1])
    print("normalizer",normalizer)

    def integrand(x):
        rrange = np.arange(range_array[0] , range_array[1] , (range_array[1] - range_array[0])/resolution)
        index = np.searchsorted(rrange,x)
        error = grid_values[index-1]
        probability = np.exp(-error/10000)/normalizer
        #  print("ranges:",rrange)
        #  print("value:",x)
        #  print("index:",index)
        #  print("error:",error)
        #  print("probability",probability)
        return probability * x

    def integrand2(x):
        rrange = np.arange(range_array[0] , range_array[1] , (range_array[1] - range_array[0])/resolution)
        index = np.searchsorted(rrange,x)
        error = grid_values[index-1]
        probability = np.exp(-error/10000)/normalizer
        #  print("ranges:",rrange)
        #  print("value:",x)
        #  print("index:",index)
        #  print("error:",error)
        #  print("probability",probability)
        return probability * (x**2)

    average ,error_1= integrate.quad(integrand, range_array[0] , range_array[1])
    print("average",average)
    average2, error_2 = integrate.quad(integrand2, range_array[0] , range_array[1])
    print("average2",average2)
    variance = average2 - average**2
    print("variance:",variance)

    return grid_values, variance


def covariance(Z,A, generate_number = None, method = None, resolution = None, **kwargs):
    print('starting covariance calculation')
    reac_t = kwargs['reaction_type']
    parameter = kwargs['parameter']
    parameter2 = kwargs['parameter2']
    parameters = param_ranges[str(Z)+str(A)]
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

    #  make a zero list which will eventually contain the guess parameters from the parameter file
    #  this is a list but not an array because the optimization functions prefer lists

    opt_begin = time.time()

    print('Begin Optimizing FREYA (This may take a while...)')

    #       DEFINE PARAMETER RANGES

    range_array = param_ranges[parameter]
    special_index = param_ranges[parameter][2]
    range_array_2 = param_ranges[parameter2]
    special_index_2 = param_ranges[parameter2][2]

    #  define function which takes in a set of parameters and returns the raw chi-squared error (total sum)
    #  see error.py for the details of this error calculation
    def err_opt(params):
        parameter = params[0]
        parameter2 = params[1]
        #  return error(Z, A, parameters[0],parameters[1], parameters[2], parameters[3], parameters[4], generate_number, parsed_data, reaction_type = reac_t)[0]
        parameters[special_index] = parameter
        parameters[special_index_2] = parameter2
        return error(Z, A, parameters[0],parameters[1], parameters[2], parameters[3], parameters[4], generate_number, parsed_data, reaction_type = reac_t)[0]

    #  for grid search method initialize the brute source routine
    if method == 'grid': 
        resolution = np.float(resolution)
        #  define the ranges to be from th eminimum to the maximum values, with the difference/resolution many cells


        #  brute_ranges = (slice(range_array[0],range_array[1],range_array[2]/resolution))
        brute_ranges = ((range_array[0],range_array[1]),(range_array_2[0],range_array[1]))

        #  block the printing for each iteration to avoid slowing the process down by taking the time to print the useless output
        #  comment this line to let it print if there is an issue with the routine

        blockPrint()
        x0, fval, grid, Jout = optimize.brute(err_opt , brute_ranges, full_output = True,Ns = resolution, finish = optimize.fmin)
        enablePrint()

        #  set finalparams to be the 0th output of the brute routine
        finalparams = x0
        #  set grid_values to be the final element of the output list of the brute routine
        grid_values = Jout

    #  define path for output to be placed
    var_path = cwd + '/../output/variance/' + 'Z=' + str(Z) + ' A=' + str(A) + '_E=' + str(Energy) + '_' + str(kwargs['parameter']) + str(kwargs['parameter2'])

    #  create appropriate directory if it does not already exist 
    if not os.path.exists(var_path):
        os.makedirs(var_path)

    opt_end = time.time()
    print("time:",opt_end - opt_begin)

    print('Begin printing output...')
    os.chdir(var_path)
    np.savetxt('covariance',grid_values)
    os.chdir(cwd)

    print("Finished.")
    print("Statistics generated in: \n"+var_path)
    print("Calculating variance...")

    def normalizer_function(x,y):
        rrange = np.arange(range_array[0] , range_array[1] , (range_array[1] - range_array[0])/resolution)
        rrange2 = np.arange(range_array_2[0] , range_array_2[1] , (range_array_2[1] - range_array_2[0])/resolution)
        index = np.searchsorted(rrange,x)
        index2 = np.searchsorted(rrange2,y)
        error = grid_values[index - 1,index2 - 1]
        probability = np.exp(-error/10000)
        #  print("ranges:",rrange)
        #  print("value:",x)
        #  print("index:",index)
        #  print("error:",error)
        #  print("probability",probability)
        return probability 

    def y_min(x):
        return range_array_2[0]
    def y_max(x):
        return range_array_2[1]

    normalizer,int_error = integrate.dblquad(normalizer_function, range_array[0] , range_array[1],y_min , y_max)
    print("normalizer",normalizer)

    def integrand(x,y):
        rrange = np.arange(range_array[0] , range_array[1] , (range_array[1] - range_array[0])/resolution)
        rrange2 = np.arange(range_array_2[0] , range_array_2[1] , (range_array_2[1] - range_array_2[0])/resolution)
        index = np.searchsorted(rrange,x)
        index2 = np.searchsorted(rrange2,y)
        error = grid_values[index-1,index2 - 1]
        probability = np.exp(-error/10000)/normalizer
        #  print("ranges:",rrange)
        #  print("value:",x)
        #  print("index:",index)
        #  print("error:",error)
        #  print("probability",probability)
        return probability * x * y

    def integrand2(x,y):
        rrange = np.arange(range_array[0] , range_array[1] , (range_array[1] - range_array[0])/resolution)
        rrange2 = np.arange(range_array_2[0] , range_array_2[1] , (range_array_2[1] - range_array_2[0])/resolution)
        index = np.searchsorted(rrange,x)
        index2 = np.searchsorted(rrange2,y)
        error = grid_values[index-1,index2 - 1]
        probability = np.exp(-error/10000)/normalizer
        #  print("ranges:",rrange)
        #  print("value:",x)
        #  print("index:",index)
        #  print("error:",error)
        #  print("probability",probability)
        return probability * (x**2) * (y**2)

    average ,error_1= integrate.dblquad(integrand, range_array[0] , range_array[1], y_min , y_max)
    print("average",average)
    average2, error_2 = integrate.dblquad(integrand2, range_array[0] , range_array[1], y_min , y_max)
    print("average2",average2)
    covariance = average2 - average**2
    print("covariance:",covariance)
    x_array = np.arange(range_array[0] , range_array[1] , (range_array[1] - range_array[0])/resolution)
    y_array = np.arange(range_array_2[0] , range_array_2[1] , (range_array_2[1] - range_array_2[0])/resolution)

    return covariance , grid_values, x_array , y_array

