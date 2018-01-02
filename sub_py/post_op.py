import time
from math import sqrt, pi
import os, sys
import matplotlib.pyplot as plt
#  from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import optimize

cwd = os.getcwd()

sys.path.append(cwd+'/../data_master/Cf252/')
from mannhart_data import mannhart_bins, mannhart_split, mannhart_bindiff
sys.path.append(cwd)

from ranges import * 
from gen_par_ana import gpa
from error import error
from data_parse import data_parse
from isotope import isotope
from plot import plot

#  define functions to block and restore printing for when freya is run many times in a row during the optimization procedure
def blockPrint():
    sys.stdout = open(os.devnull, 'w')
def enablePrint():
    sys.stdout = sys.__stdout__

def post_opt(Z,A, generate_number = None, method = None, resolution = None, **kwargs):
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

    # create directory for the plots to be placed
    data_path = cwd + '/../output/data_plots/' + 'Z=' + str(Z) + 'A=' + str(A) + '/'
    if not os.path.exists(data_path):
        os.makedirs(data_path)

    print('Begin plotting data...')

    for key in parsed_data[0]:
        element = parsed_data[0][str(key)]
        dim_status = element[5]
        data_array = element[0]

        #  recall at the beginning we defined a value to tell us later whether there is one or two sets of bins corresponding to each data point, 
        #  and therefore whether we plot it on 2 or 3 axes. 
        #  that is what determines the first separation 

        os.chdir(data_path)
        if dim_status is False:
            #  element 9 is the sigma indicator. This determines whether or not the error bars on the plot will be the variance or the uncertainty 
            #  this is reflected in the name of the pdf output of the plots
            if element[9] is True:
                print('plotting: ',key)
                plt.plot( data_array[:, 0, 0] , data_array[:, 1, 0] , color = 'r' , label = None)
                plt.plot( data_array[:, 0, 0] , data_array[:, 1, 0] , '^' , color = 'r')
                plt.errorbar( data_array[:, 0, 0] , data_array[: ,1, 0] , yerr = data_array[:, 2, 0], 
                        color = 'r', fmt = ' ' , capsize = 3, elinewidth = 1)
                plt.title(str(key)  + ' ' +  str(element[1])  + ' ' +  str(iso)  + ' ' +  '(Variance)')
                plt.xlabel( element[2] )
                plt.ylabel( element[3] )
                plt.xlim( ranges_x[key_translator[key]][0] , ranges_x[key_translator[key]][1] )
                plt.ylim( ranges_y[key_translator[key]][0] , ranges_y[key_translator[key]][1] )
                plt.xscale(element[7])
                plt.savefig(str(key) + '(variance)' + '.pdf')
                plt.close()

            else:
                print('plotting: ',key)
                plt.plot( data_array[:, 0, 0] , data_array[:, 1, 0] , color = 'r' , label = None)
                plt.plot( data_array[:, 0, 0] , data_array[:, 1, 0] , '^' , color = 'r')
                plt.errorbar( data_array[:, 0, 0] , data_array[: ,1, 0] , yerr = data_array[:, 1, 1], 
                        color = 'r', fmt = ' ' , capsize = 3, elinewidth = 1)
                plt.title(str(key)  + ' ' +  str(element[1])  + ' ' +  str(iso)  )
                plt.xlabel( element[2] )
                plt.ylabel( element[3] )
                plt.xlim( ranges_x[key_translator[key]][0] , ranges_x[key_translator[key]][1] )
                plt.ylim( ranges_y[key_translator[key]][0] , ranges_y[key_translator[key]][1] )
                plt.xscale(element[7])
                plt.savefig(str(key) + '.pdf')
                plt.close()

        elif dim_status is True:
            print('plotting: ',key)

            #  nu = np.concatenate((freya_output["light_neutrons"],freya_output["heavy_neutrons"]))
            #  A = np.concatenate((freya_output['Al'],freya_output['Ah']))
            #  TKE = np.tile(freya_output['TKE'],2)
            #
            #  plt.hist2d(A, TKE,bins=(max(A) - min(A),150),weights=nu, normed = True)
            #
            #  plt.colorbar(format='%7.2e')
            #
            #  plt.xlabel("Mass (amu)")
            #  plt.ylabel("Total Kinetic Energy (MeV)")
            #  plt.title('Neutron Yield as function of Mass, TKE')
            #
            #  plt.savefig('n_A_TKE' + '.pdf')
            #  plt.close()

            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(data_array[:,0,0] , data_array[:,1,0] , data_array[:,2,0])
            ax.set_xlabel(element[2])
            ax.set_ylabel(element[3])
            ax.set_zlabel(element[8])
            ax.set_xlim(ranges_x[key_translator[key]][0],ranges_x[key_translator[key]][1])
            ax.set_ylim(ranges_y[key_translator[key]][0],ranges_y[key_translator[key]][1])
            ax.set_zlim(ranges_z[key_translator[key]][0],ranges_z[key_translator[key]][1])
            ax.text2D(0.05, 0.95, str(key) + ' ' + str(element[1]) + ' ' + str(iso), transform=ax.transAxes)
            plt.savefig(str(key) + '.pdf')
            plt.close()

        os.chdir(cwd)

    print('Successful.')
    print('Plots generated in:',data_path)

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

    finalparams = np.zeros((10,0))
    res = None

    #  define path for optimization output to be placed
    opt_path = cwd + '/../output/optimization/' + 'Z=' + str(Z) + ' A=' + str(A) + '_E=' + str(Energy) + '_' + 'post_op'

    #  create appropriate directory if it does not already exist 
    if not os.path.exists(opt_path):
        os.makedirs(opt_path)

    print('Calling for final error estimation...')

    print(reac_t)
    chisq_array =  error(Z,A,None, None, None, None, None,generate_number, parsed_data , reaction_type = reac_t)
    
    print('Final Chi-Squared Error Estimation: ',str(chisq_array[0]))
    print('These values, and further statistics/visuals generated in ',str(opt_path))

    print('Begin Generating, Analyzing, Plotting Events with Final Set of Parameters')

    final = plot(Z, A , Energy, str(iso[0])  + '.' + 'opt', generate_number = generate_number) 
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
        '\\chi^2 & = & ' + str(    chisq_array[0]      ) +
        '\n\\end{aligned} $$'
        '$$ \\begin{aligned} \n' +
        '\\text{Times: } & ' +
        '\n\\\\ & ' +
        '\\text{Events} & = & ' + str(        generate_number        ) +
        '\n\\end{aligned} $$'
        '\n\\end{document}')
    print("Successful.")
    print("Statistics generated (as LaTeX) in: \n"+opt_path)

    print('Begin generating plots of Freya/Data comparison...')

    # plot data and freya on top of each other

    os.chdir(opt_path)

    for key in parsed_data[0]:
        element = parsed_data[0][str(key)]
        if element is None:
            continue
        else:
            print('plotting: ',key)

            dim_status = element[5]
            data_array = element[0]

            freya_data = freya_dict[key_translator[key]][0]

            freya = plt.plot(freya_data[:,0] , freya_data[:,1] , '^-' , color = 'r' , label = 'FREYA Output') 
            data = plt.plot(data_array[:,0] , data_array[:,1] , '^-' , color = 'b' , label = str(element[1]) + ' data')

            plt.xlabel(element[2])
            plt.ylabel(element[3])
            plt.xlim( ranges_x[key][0] , ranges_x[key][1])
            plt.ylim( ranges_y[key][0] , ranges_y[key][1])
            plt.xscale(element[7])

            plt.legend(("FREYA",str(element[1]) + ' data'))
            plt.title(str(key) + ' ' + str(iso[0]))
            plt.savefig(str(key) + '.pdf')
            plt.close()

    print('Begin generating plots of Chi-Squared distributions...')

    #  plot chi squared distributions
    os.chdir(opt_path)
    for key in chisq_array[1]:
        chi_dict = chisq_array[1]
        chi_x = chi_dict[key][:,0]
        chi_y = chi_dict[key][:,1]
        plt.scatter(chi_x,chi_y,color = 'b')
        plt.xlim(ranges_x[key_translator[key]][0], ranges_x[key_translator[key]][1])
        plt.xlabel('x value')
        plt.ylabel('error')
        plt.title('Chi-Squared Error for: ' + str(key) + ' (' + str(iso[0]) + ')')
        if key_translator[key] is 'mannhart': 
            plt.xscale('log')
        plt.savefig(str(key) + '_chisq.pdf')
        plt.close()
    print('Finished.')

    os.chdir(cwd)

    print("Plots generated in: \n"+opt_path)

    return 
