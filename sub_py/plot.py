import numpy as np
import matplotlib.pyplot as plt
#  from mpl_toolkits.mplot3d import Axes3D
import time
import os, sys

cwd = os.getcwd()

from gen_par_ana import gpa
from ranges import *


def plot(Z, A, Energy, output_file, **kwargs):
    #  generate events and parse into neatly organized arrays using gpa function
    freya = gpa(Z, A, Energy, output_file, e = kwargs.get('e'), 
        x = kwargs.get('x'), c = kwargs.get('c'), T = kwargs.get('T'), d = kwargs.get('d'), 
        generate_number = kwargs.get('generate_number') )
    freya_output = freya[0]
    times = freya[1]

    #  define, make path for plots to be placed
    plot_path = cwd+'/../output/'+str(output_file)
    if not os.path.exists(plot_path):
        os.makedirs(plot_path)
    print("Begin Printing Statistics...")
    plot_begin = time.time()
    if int(Energy) == -1:
        string_energy = 0

    stats =open(plot_path+'/Statistics.tex', 'w+')
    stats.write('\\documentclass[12pt]{letter} \n \\usepackage{amsmath} ' +
        '\n \\begin{document} '
        '\n $$\\begin{aligned} ' +
        '\\text{Reaction: } & ' +
        '\n\\\\ & ' +
        'Z & = & ' + str(Z) +
        '\n\\\\ & ' +
        'A & = & ' + str(A) +
        '\n\\\\ & ' +
        '\\text{Energy} & = & ' + str(string_energy) +
        '\n\\end{aligned} $$ \n' +
        '\n $$\\begin{aligned} ' +
        '\\text{Total Neutron Multiplicity Moments: } & ' +
        '\n\\\\ & ' +
        '\\bar\\nu & = & ' + str(freya_output[ 'nubar' ][0][0,1,0]) +
        '\pm' + str(freya_output[ 'nubar' ][0][0,2,0]) +
        '\n\\\\ & ' +
        '\\nu_1 & = & ' + str(freya_output[ 'nu1' ][0]) +
        '\pm' + str(freya_output[ 'nu1' ][1]) +
        '\n\\\\ & ' +
        '\\nu_2 & = & ' + str(freya_output[ 'nu2' ][0]) +
        '\pm' + str(freya_output[ 'nu2' ][1]) +
        '\n\\\\ & ' +
        '\\nu_3 & = & ' + str(freya_output[ 'nu3' ][0]) +
        '\pm' + str(freya_output[ 'nu3' ][1]) +
        '\n\\\\ & ' +
        '\\nu_4 & = & ' + str(freya_output[ 'nu4' ][0]) +
        '\pm' + str(freya_output[ 'nu4' ][1]) +
        '\n\\\\ & ' +
        '\\bar\\gamma & = & ' + str(freya_output[ 'gammabar' ][0][0,1,0]) +
        '\pm' + str(freya_output[ 'gammabar' ][0][0,2,0]) +
        '\n\\end{aligned} $$ \n' +
        '\n $$\\begin{aligned} ' +
        '\\text{Gamma statistics: } & ' +
        '\n\\\\ & ' +
        '\\text{Average photon energy} & : & ' + str(freya_output[ 'average_photon_energy' ][0][0,1,0]) +
        '\pm' + str(freya_output[ 'average_photon_energy' ][0][0,2,0]) +
        '\n\\end{aligned} $$ \n'
        '$$ \\begin{aligned} \n' +
        '\\text{Times: } & ' +
        '\n\\\\ & ' +
        '\\text{Events} & = &' + str(freya_output['number_of_events']) +
        '\n\\\\ & ' +
        '\\text{Generation Time} & = & ' + str(round(float( times[ 'gen_time' ] ),5)) + '\\text{ sec}' +
        '\n\\\\ & '
        '\\text{Parse Time} & = &' + str(round(float( times[ 'parse_time'  ] )/60,5)) + '\\text{ min}' +
        '\n\\end{aligned} $$'
        '\n\\end{document}')
    print("Successful.")
    print("Statistics generated (as LaTeX) in: \n"+plot_path)

    print("Begin Plotting...")

    for key in freya_output:
        #  give mannhart a logarithmic scaling
        if key == 'mannhart':
            scale_var = 'log'
        else:
            scale_var = 'linear'
        
        #  give n_A_TKE 3d axes
        if key == 'n_A_TKE':
            is_3d = True
        else:
            is_3d = False

        obs_master = freya_output[str(key)]

        os.chdir(plot_path)

        if isinstance(obs_master, list):
            if is_3d == False:
                observable_array = obs_master[0]

                x_lim1 = ranges_x[str(key)][0]
                x_lim2 = ranges_x[str(key)][1]
                y_lim1 = ranges_y[str(key)][0]
                y_lim2 = ranges_y[str(key)][1]

                plt.plot( observable_array[:,0] , observable_array[:,1] , color = 'r' , label = obs_master[1])
                plt.plot( observable_array[:,0] , observable_array[:,1] , '^' , color = 'r' )
                plt.errorbar( observable_array[:,0] , observable_array[:,1] , yerr = observable_array[:,2] , color = 'r' , fmt = ' ' , capsize = 3 , elinewidth = 1)
                plt.xlabel( obs_master[3] )
                plt.ylabel( obs_master[4] )
                plt.xlim( x_lim1 , x_lim2 )
                plt.ylim( y_lim1 , y_lim2 )
                plt.tick_params(direction = 'in' , labelright = True, right = True)
                plt.xscale(scale_var)
                #  plt.title( obs_master[2] )
                plt.legend()

                plt.savefig( str(key) + '.pdf' )
                plt.close()

    nu = np.concatenate((freya_output["light_neutrons"],freya_output["heavy_neutrons"]))
    A = np.concatenate((freya_output['Al'],freya_output['Ah']))
    TKE = np.tile(freya_output['TKE'],2)

    plt.hist2d(A, TKE,bins=(max(A) - min(A),150),weights=nu, normed = True)

    plt.colorbar(format='%7.2e')

    plt.xlabel("Mass (amu)")
    plt.ylabel("Total Kinetic Energy (MeV)")
    #  plt.title('Neutron Yield as function of Mass, TKE')

    plt.savefig('n_A_TKE' + '.pdf')
    plt.close()

    # Plot Fragment/Product Mass Distribution
    fragment_array = np.array( freya_output['Fragment_A'] )

    product_list = freya_output['Product_A']
    product_array = product_list[0]
    
    plt.plot(fragment_array[:,0] , fragment_array[:,1] , color = 'r' , 
            label = 'Fragment Mass') 
    plt.plot(fragment_array[:,0] , fragment_array[:,1] , '^' , color = 'r')
    plt.plot(fragment_array[:,0] , product_array[:,1] , color = 'b' , label = 'Product Mass')
    plt.plot(fragment_array[:,0] , product_array[:,1] , '^', color = 'b')
    plt.xlabel('Mass Number (A)')
    plt.ylabel('Probability')
    plt.xlim( ranges_x['Fragment_A'][0] , ranges_x['Fragment_A'][1])
    plt.ylim( ranges_y['Fragment_A'][0] , ranges_y['Fragment_A'][1])
    plt.legend()
    #  plt.title('Fragment, Product Mass Yields')
    plt.savefig('Fragment_Product_A.pdf')
    plt.close()

    os.chdir(cwd)

    if not os.path.exists(plot_path + "/raw_data/"):
        os.makedirs(plot_path+"/raw_data/")

    os.chdir(plot_path + "/raw_data/")

    for key in freya_output:
        entry = freya_output[key]
        if isinstance(entry, list):
            obs_array = entry[0]
            np.savetxt(key+".txt",obs_array)

    os.chdir(cwd)

    print("Successful.")
    plot_end = time.time()
    print('Time: ' + str(plot_end - plot_begin) + " sec")
    print("Plots Generated in: \n"+plot_path)

    return True
