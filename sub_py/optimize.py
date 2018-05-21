import time
from math import sqrt, pi
import os, sys
import numpy as np
from scipy import optimize

cwd = os.getcwd()

from ranges import * 
from gen_par_ana import gpa
from error import *
from data_parse import data_parse
from isotope import isotope
from anneal import anneal

from test import *

#  define functions to block and restore printing for when freya is run many times in a row during the optimization procedure
def blockPrint():
    sys.stdout = open(os.devnull, 'w')
def enablePrint():
    sys.stdout = sys.__stdout__

def opt(Z,A, generate_number = None, method = None, resolution = None, **kwargs):
    print('starting')
    reac_t = kwargs['reaction_type']
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
    guesses = [0,0,0,0,0]
    guesses[0] = float(line_split[3])
    guesses[1] = float(line_split[4])
    guesses[2] = float(line_split[5])
    guesses[3] = float(line_split[6])
    guesses[4] = float(line_split[9])
    guesses = np.array(guesses)

    opt_begin = time.time()

    print('Begin Optimizing FREYA (This may take a while...)')

    #       DEFINE PARAMETER RANGES

    #  e_range = [6,12,0]
    e_range = [10,11,0]
    e_range[2] = e_range[1] - e_range[0]

    #  x_range = [1,1.5,0]
    x_range = [1.2,1.35,0]
    x_range[2] = x_range[1] - x_range[0]

    #  c_range = [1,2,0]
    c_range = [1.1,1.2,0]
    c_range[2] = c_range[1] - c_range[0]

    #  T_range = [0.5,1.5,0]
    T_range = [0.8,1,0]
    T_range[2] = T_range[1] - T_range[0]

    #  d_range = [-5,5,0]
    d_range = [-1,1,0]
    d_range[2] = d_range[1] - d_range[0]

    #  define function which takes in a set of parameters and returns the raw chi-squared error (total sum)
    #  see error.py for the details of this error calculation
    def err_opt(parameters):
        return error(Z, A, parameters[0],parameters[1], parameters[2], parameters[3], parameters[4], generate_number, parsed_data, reaction_type = reac_t)[0]
        #  return test_error(Z, A, parameters[0],parameters[1], parameters[2], parameters[3], parameters[4], generate_number, parsed_data, reaction_type = reac_t)[0]

    #  define class to give parameter bounds to optimization routine
    class MyBounds(object):
        def __init__(self, xmax=[
                e_range[1],
                x_range[1],
                c_range[1],
                T_range[1],
                d_range[1]
                ], xmin=[
                    e_range[0],
                    x_range[0],
                    c_range[0],
                    T_range[0],
                    d_range[0]
                    ] ):
            self.xmax = np.array(xmax)
            self.xmin = np.array(xmin)
        def __call__(self, **kwargs):
            x = kwargs["x_new"]
            tmax = bool(np.all(x <= self.xmax))
            tmin = bool(np.all(x >= self.xmin))
            return tmax and tmin
    mybounds = MyBounds() 

    #  for grid search method initialize the brute source routine
    if method == 'grid': 
        resolution = np.float(resolution)
        #  define the ranges to be from th eminimum to the maximum values, with the difference/resolution many cells
        brute_ranges = (slice(e_range[0],e_range[1],e_range[2]/resolution),

            slice(x_range[0] , x_range[1] , x_range[2]/resolution),
            slice(c_range[0] , c_range[1] , c_range[2]/resolution),
            slice(T_range[0] , T_range[1] , T_range[2]/resolution),
            slice(d_range[0] , d_range[1] , d_range[2]/resolution)
            )
        #  block the printing for each iteration to avoid slowing the process down by taking the time to print the useless output
        #  comment this line to let it print if there is an issue with the routine
        blockPrint()
        x0, fval, grid, Jout = optimize.brute(err_opt , brute_ranges , full_output = True, finish = optimize.fmin)
        enablePrint()

        #  set finalparams to be the 0th output of the brute routine
        finalparams = x0
        #  set grid_values to be the final element of the output list of the brute routine
        grid_values = Jout

    elif method == 'anneal':
        finalparams , this_error = anneal(err_opt , guesses)
        grid_values = None

    #  this is a hidden feature of the optimization routine, which allows the user to simply plot the data (alone and along with freya) without waiting for the optimization to happen
    elif method == 'bypass':
        print('skipping optimization')
        finalparams = np.zeros((10,0))
        res = None

    #  the default optimization routine is the stochastic method
    else: 
        #  the following lines can be commented/uncommented (leave only one uncommented) to try different optimization methods
        methods = ["BFGS","Nelder-Mead","Powell","Cobyla"]
        minimizer_kwargs = {"method": methods[int(kwargs['stochastic_type'])]}

        blockPrint()
        res = optimize.basinhopping( err_opt , guesses , minimizer_kwargs = minimizer_kwargs, 
                niter = 100, 
                accept_test = mybounds, 
                stepsize = 0.5
                )
        enablePrint()
        finalparams = res.x
        grid_values = None
    
    opt_end = time.time()

    opt_time = opt_end - opt_begin
    print('Time:',opt_time)



    #  define path for optimization output to be placed
    opt_path = cwd + '/../output/optimization/' + 'Z=' + str(Z) + ' A=' + str(A) + '_E=' + str(Energy) + '_' + str(method)

    #  create appropriate directory if it does not already exist 
    if not os.path.exists(opt_path):
        os.makedirs(opt_path)

    if finalparams is not None:
        print('Successful initial optimization.')
    else:
        print('Error in initial optimization routine!')
    print('Calling for final error estimation...')

    #  print(reac_t)
    chisq_array =  error(Z,A,None, None, None, None, None , generate_number, parsed_data , reaction_type = reac_t)
    #  chisq_array =  test_error(Z,A,finalparams[0],finalparams[1],finalparams[2],finalparams[3],finalparams[4], generate_number, parsed_data , reaction_type = reac_t)
    
    if finalparams is not None:
        print('Final Set of Parameters: ',finalparams)
    print('Final Chi-Squared Error Estimation: ',str(chisq_array[0]))
    print('These values, and further statistics/visuals generated in ',str(opt_path))

    print('Begin Generating, Analyzing, Plotting Events with Final Set of Parameters')

    freya = gpa(Z,A,Energy,'cf.plot',generate_number = generate_number)
    freya_dict = freya[0]    

    opt_end = time.time()

    print('Begin printing statistics...')
    opt_stats =open(opt_path+'/Statistics.tex', 'w+')
    opt_stats.write(
        '\\documentclass[12pt]{letter} \n \\usepackage{amsmath} ' + 
        '\n \\begin{document} '
        '\n $$\\begin{aligned} ' + 
        '\\text{Final Set of Parameters: } & ' + 
        '\n\\\\ & ' +
        'e & = & ' + str(      finalparams[0]   ) +
        '\n\\\\ & ' +
        'x & = & ' + str(        finalparams[1]     ) + 
        '\n\\\\ & ' +
        'c &= & ' + str(         finalparams[2]     ) +
        '\n\\\\ & ' +
        'c_s & = & ' + str(      finalparams[3]         ) + 
        '\n\\\\ & ' +
        'dTKE & = & ' + str(        finalparams[4]      ) + 
        '\n\\end{aligned} $$ \n' +
        '$$ \\begin{aligned} \n' +
        '\\text{Final Chi-Squared Estimation: } & ' +
        '\n\\\\ & ' +
        '\\tilde\\chi^2 & = & ' + str(    chisq_array[0]      ) +
        '\n\\\\ & ' +
        '\\chi^2 & = & ' + str(    chisq_array[2]      ) +
        '\n\\end{aligned} $$'
        '$$ \\begin{aligned} \n' +
        '\\text{Times: } & ' +
        '\n\\\\ & ' +
        '\\text{Events} & = & ' + str(        generate_number        ) +
        '\n\\\\ & ' +
        #  '\\text{Optimization Time} & = & ' + str( opt_time ) + '\\text{ sec}' +
        '\n\\end{aligned} $$'
        '\n\\end{document}')
    print("Successful.")
    print("Statistics generated (as LaTeX) in: \n"+opt_path)

    return chisq_array[0], grid_values
