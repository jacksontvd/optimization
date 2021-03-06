import time
from math import sqrt, pi
import os, sys
import matplotlib
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'
rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import numpy as np
from scipy import optimize

cwd = os.getcwd()

from ranges import * 
from gen_par_ana import gpa
from error import *
from data_parse import data_parse
from isotope import isotope
from plot import plot
from matplotlib.ticker import MaxNLocator
from maxwellian import *

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def post_opt(Z,A, generate_number = None, method = None, resolution = None, **kwargs):
    print('starting')
    reac_t = kwargs['reaction_type']
    os.chdir(cwd+'/../../freya/data_freya/')

    Energy = kwargs["energy"]
    if int(Energy) == -1:
        string_energy = 0
    else:
        string_energy = str(Energy)


    infile = open("inputparameters.dat","r+")
   
    #  parse appropriate data with fuction from data_parse.py
    parsed_data = data_parse(Z,A,reac_t,Energy)

    #  print("UHOHUOHUHOH",error(92,235,10.37 ,  1.15   ,1.3    ,0.87 ,  0.39, 10000,parsed_data,reaction_type='induced',Energy = 2.9)[0])

    #  pull the data array out of the data_parse output list
    data_array = parsed_data[0]

    #  pull the data key dictionary from the data_parse output list
    key_translator = parsed_data[1]


    #  pull the line number for the isotope from the isotope function
    #  this is used to rewrite the appropriate line of the parameter file in each iteration of the optimization
    #  it is also used to read off the guess values from the parameter file as it is currently
    iso = isotope(Z,A,reac_t)
    i = iso[1]

    content = infile.readlines() #reads line by line and outputs a list of each line
    line_split = content[i].split()

    # create directory for the plots to be placed
    data_path = cwd + '/../../output/data_plots/' + 'Z=' + str(Z) + 'A=' + str(A) + '/'
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
            print('plotting: ',key)
            plt.plot( data_array[:, 0, 0] , data_array[:, 1, 0] , color = 'r' , label = None)
            plt.plot( data_array[:, 0, 0] , data_array[:, 1, 0] , '^' , color = 'r')
            plt.errorbar( data_array[:, 0, 0],data_array[: ,1, 0],yerr = data_array[:, 1, 1], 
                    color = 'r', fmt = ' ' , capsize = 3, elinewidth = 1)
            #  plt.title(str(key)  + ' ' +  str(element[1])  + ' ' +  str(iso)  )
            plt.xlabel( element[2] )
            plt.ylabel( element[3] )
            plt.xlim( ranges_x[key_translator[key]][0] , ranges_x[key_translator[key]][1] )
            plt.ylim( ranges_y[key_translator[key]][0] , ranges_y[key_translator[key]][1] )
            plt.xscale(element[7])

            plt.savefig(str(key) + '.pdf')
            plt.close()

        elif dim_status is True:
            print('plotting: ',key)
#  #  #  #  ADD CONTOUR
        if key == "alt_n_spectrum":
            print("plotting mannhart over alternative spectrum...")

            maxwell_temp = 1.32

            #  data_array[:,1,0] = data_array[:,1,0] * MaxwellianSpectrum(maxwell_temp)
            #  data_array[:,1,1] = data_array[:,1,1] * MaxwellianSpectrum(maxwell_temp)
###
            plt.plot( data_array[:, 0, 0] , data_array[:, 1, 0] , color = 'r' , label = None)
            plt.plot( data_array[:, 0, 0] , data_array[:, 1, 0] , '^' , color = 'r')
            plt.errorbar( data_array[:, 0, 0] , data_array[: ,1, 0] , yerr = data_array[:, 1, 1], 
                    color = 'r', fmt = ' ' , capsize = 3, elinewidth = 1)
            #  plt.title(str(key)  + ' ' +  str(element[1])  + ' ' +  str(iso)  )
            plt.xlabel( element[2] )
            plt.ylabel( element[3] )
            #  plt.xlim( ranges_x[key_translator[key]][0] , ranges_x[key_translator[key]][1] )
            #  plt.ylim( ranges_y[key_translator[key]][0] , ranges_y[key_translator[key]][1] )
            plt.xscale(element[7])

            plt.savefig(str(key) + '.pdf')
            plt.close()
###

            plt.plot( data_array[:, 0, 0] , data_array[:, 1, 0] , 'r^-' , label = "Mannhart")
            plt.errorbar(data_array[:,0,0],data_array[:,1,0],
                    yerr = data_array[:, 1, 1], color = 'r', fmt = ' ' , capsize = 3, elinewidth = 1)
            spectrum_element = parsed_data[0]["n_spectrum"]
            spectrum_array = spectrum_element[0]
            data_array = spectrum_array
            total_count = np.sum(np.multiply(data_array[:,0,0],data_array[:,1,0]))*0.05
            data_array[:,1,0] = data_array[:,1,0] / total_count 
            data_array[:,1,1] = data_array[:,1,1] / total_count
            plt.plot( data_array[:, 0, 0] , data_array[:, 1, 0] , 'b^-',label = u'Göök')
            plt.errorbar(data_array[:, 0, 0],data_array[: ,1, 0] , yerr = data_array[:, 1, 1], 
                    color = 'b', fmt = ' ' , capsize = 3, elinewidth = 1)
            plt.xlabel( element[2] )
            plt.ylabel( element[3] )
            plt.xlim( ranges_x[key_translator[key]][0] , ranges_x[key_translator[key]][1] )
            plt.xscale("log")
            #  plt.yscale("log")
            lg = plt.legend(fontsize=14,numpoints=1)
            lg.draw_frame(False)
            plt.savefig("mannhart_comparison" + '.pdf')
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
    opt_path = cwd + '/../../output/optimization/' + 'Z=' + str(Z) + ' A=' + str(A) + '_E=' + str(Energy) + '_' + 'post_op'

    #  create appropriate directory if it does not already exist 
    if not os.path.exists(opt_path):
        os.makedirs(opt_path)

    print('Calling for final error estimation...')

    chisq_array =  error(Z,A,None, None, None, None, None,generate_number, parsed_data , reaction_type = reac_t, Energy = Energy)
    
    print('Final Chi-Squared Error Estimation: ',str(chisq_array[0]))
    print('These values, and further statistics/visuals generated in ',str(opt_path))

    print('Begin Generating, Analyzing, Plotting Events with Final Set of Parameters')

    freya = chisq_array[4]
    #  freya = gpa(Z,A,Energy,'cf.plot',generate_number = generate_number)
    freya_dict = freya[0]    
    times = freya[1]

    opt_end = time.time()

    print('Begin printing statistics...')
    opt_stats =open(opt_path+'/Statistics.tex', 'w+')
    opt_stats.write(
        '\\documentclass[12pt]{letter} \n \\usepackage{amsmath} ' + 
        '\n \\begin{document} '
        '\n\\[\\begin{tabular}{|c|c|c|c|c|}' + 
        '\n\hline\n' + 
        '$e_0 \\text{(1/MeV)}$ & $x$ & $c$ & $c_S$ & $d\\text{TKE} \\text{(MeV)}$ \\\\ \n' + 
        '\hline\hline\n' + 
        str(      finalparams[0]   ) +
        ' & ' +
        str(        finalparams[1]     ) + 
        ' & ' +
        str(         finalparams[2]     ) +
        ' & ' +
        str(      finalparams[3]         ) + 
        ' & ' +
        str(        finalparams[4]      ) + 
        '\n\\\\\hline' + 
        '\n \\end{tabular} \] \n' +
        '\[ \\begin{aligned} \n' +
        '\\text{Final Chi-Squared Estimation: } & ' +
        '\n\\\\ & ' +
        '\\tilde\\chi^2 & = & ' + str(    chisq_array[0]      ) +
        '\n\\\\ & ' +
        '\\chi^2 & = & ' + str(    chisq_array[2]      ) +
        '\n\\end{aligned} \]'
        '\[ \\begin{aligned} \n' +
        '\\text{Times: } & ' +
        '\n\\\\ & ' +
        '\\text{Events} & = & ' + str(        generate_number        ) +
        '\n\\end{aligned}\\]' 
        '\n\\end{document}')

    opt_stats_2 =open(opt_path+'/statistics_2.tex', 'w+')
    opt_stats_2.write('\\documentclass[12pt]{letter} \n' +
        '\\usepackage{amsmath}'+
        '\\usepackage{amsfonts}'+
        '\\usepackage{amssymb}'+
        '\\usepackage[english]{babel}'+
        '\\usepackage{graphicx}'+
        '\\usepackage{color} '+
        '\\usepackage{todonotes}'+
        '\\usepackage{cleveref}'+
        '\\usepackage{array}'+
        '\\usepackage{lipsum}\n'+
        '\\newcommand{\code}{$\mathtt{FREYA}$}'+
        '\\newcolumntype{L}{>{$}l<{$}}'+
        '\\newcolumntype{C}{>{$}c<{$}}'+
        '\\newcolumntype{R}{>{$}r<{$}}'+
        '\n \\begin{document} '+
        '\n $$\\begin{aligned} ' +
        '\\text{Reaction: } & ' +
        '\n\\\\ & ' +
        'Z & = & ' + str(Z) +
        '\n\\\\ & ' +
        'A & = & ' + str(A) +
        '\n\\\\ & ' +
        '\\text{Energy} & = & ' + str(string_energy) +
        '\n\\end{aligned} $$ \n' +
        '\n\[\\begin{tabular}{CCCC} ' +
        '\n\hline'
        '\n & \\text{ Evaluation } &\\text{\code\ Result} & \\text{C/E} \\\\'
        '\n\hline'
        '\n\\bar\\nu &' + 
#  print out data in the first column
        str(parsed_data[0]['nubar'][0][0,1,0])+ 
        '\pm '+ str(parsed_data[0]['nubar'][0][0,1,1])
        + '&' +
#  print out freya in the second column
        str(round(freya_dict[ 'nubar' ][0][0,1,0],2)) +
        '\pm ' + str(round(freya_dict[ 'nubar' ][0][0,2,0],2)) 
        + '&' +
#  print out C/E in the third column
        str(round(freya_dict[ 'nubar' ][0][0,1,0]/parsed_data[0]['nubar'][0][0,1,0],2))+ 
        '\pm '+ str(round(ratio_bar_scheme(
            freya_dict[ 'nubar' ][0][0,1,0],
            freya_dict[ 'nubar' ][0][0,2,0],
            parsed_data[0]['nubar_moments'][0][0,1,0],
            parsed_data[0]['nubar_moments'][0][0,1,1]
            ),2))+
# next line
        '\\\\' + '\n \hline'
        '\n\\nu_2 &' + 
#  print out data in the first column
        str(parsed_data[0]['nubar_moments'][0][1,1,0])+
        '\pm '+ str(parsed_data[0]['nubar_moments'][0][1,1,1])+ '&' +
#  print out freya in the second column
        str(round(freya_dict[ 'nu2' ][0],2)) +
        '\pm ' + str(round(freya_dict[ 'nu2' ][1],2)) + '&' +
#  print out C/E in the third column
        str(round(freya_dict[ 'nu2' ][0]
            /parsed_data[0]['nubar_moments'][0][1,1,0],2))+
        '\pm '+ str(round(ratio_bar_scheme(
            freya_dict[ 'nu2' ][0],
            freya_dict[ 'nu2' ][1],
            parsed_data[0]['nubar_moments'][0][1,1,0],
            parsed_data[0]['nubar_moments'][0][1,1,1]
            ),2))+
# next line
        '\\\\' + '\n \hline'
        '\n\\nu_3 &'+
#  print out data in the first column
        str(parsed_data[0]['nubar_moments'][0][2,1,0])+
        '\pm '+ str(parsed_data[0]['nubar_moments'][0][2,1,1])+ '&' +
#  print out freya in the second column
        str(round(freya_dict[ 'nu3' ][0],2)) +
        '\pm ' + str(round(freya_dict[ 'nu3' ][1],2))+ '&' +
#  print out C/E in the third column
        str(round(freya_dict[ 'nu3' ][0]
            /parsed_data[0]['nubar_moments'][0][2,1,0],2))+
        '\pm '+ str(round(ratio_bar_scheme(
            freya_dict[ 'nu3' ][0],
            freya_dict[ 'nu3' ][1],
            parsed_data[0]['nubar_moments'][0][2,1,0],
            parsed_data[0]['nubar_moments'][0][2,1,1]
            ),2))+
        '\\\\' + '\n \hline' +
        '\n\\end{tabular} \]' +
# next table
        '\n\[\\begin{tabular}{CCCC} ' +
        '\n\hline'
        '\n& \\text{ Data } &\\text{\code\ Result} & \\text{C/E} \\\\'
        '\n\hline'
        '\n\overline{N}_{\gamma} &' + 
#  print out data in the first column
        str(parsed_data[0]['gammabar'][0][0,1,0])+
        '\pm '+ str(parsed_data[0]['gammabar'][0][0,1,1])+ '&' +
#  print out freya in the second column
        str(round(freya_dict[ 'gammabar' ][0][0,1,0],2)) +
        '\pm ' + str(round(freya_dict[ 'gammabar' ][0][0,2,0],2))+ '&' +
#  print out C/E in the third column
        str(round(freya_dict[ 'gammabar' ][0][0,1,0]
            /parsed_data[0]['gammabar'][0][0,1,0],2))+
        '\pm '+ str(round(ratio_bar_scheme(
            freya_dict[ 'gammabar' ][0][0,1,0],
            freya_dict[ 'gammabar' ][0][0,2,0],
            parsed_data[0]['gammabar'][0][0,1,0],
            parsed_data[0]['gammabar'][0][0,2,0]
            ),2))+
# next line
        '\\\\' + '\n \hline'
        '\n\overline{\epsilon_\gamma} \\text{(MeV)}&' + 
#  print out data in the first column
        str(parsed_data[0]['average_photon_energy'][0][0,1,0])+
        '\pm '+ str(parsed_data[0]['average_photon_energy'][0][0,1,1])+ '&' +
#  print out freya in the second column
        str(round(freya_dict[ 'average_photon_energy' ][0][0,1,0],2)) +
        '\pm ' + str(round(freya_dict[ 'average_photon_energy' ][0][0,2,0],2))+ '&' +
#  print out C/E in the third column
        str(round(freya_dict[ 'average_photon_energy' ][0][0,1,0]
            /parsed_data[0]['average_photon_energy'][0][0,1,0],2))+
        '\pm '+ str(round(ratio_bar_scheme(
            freya_dict[ 'average_photon_energy' ][0][0,1,0],
            freya_dict[ 'average_photon_energy' ][0][0,2,0],
            parsed_data[0]['average_photon_energy'][0][0,1,0],
            parsed_data[0]['average_photon_energy'][0][0,1,1]
            ),2))+
        '\\\\' + '\n \hline'
        '\n\overline{TKE} \\text{(MeV)}&' + 
#  print out data in the first column
        str(parsed_data[0]['TKE_bar'][0][0,1,0])+
        '\pm '+ str(parsed_data[0]['TKE_bar'][0][0,1,1])+ '&' +
#  print out freya in the second column
        str(round(freya_dict[ 'TKE_bar' ][0][0,1,0],2)) +
        '\pm ' + str(round(freya_dict[ 'TKE_bar' ][0][0,2,0],2))+ '&' +
#  print out C/E in the third column
        str(round(freya_dict[ 'TKE_bar' ][0][0,1,0]
            /parsed_data[0]['TKE_bar'][0][0,1,0],2))+
        '\pm '+ str(round(ratio_bar_scheme(
            freya_dict[ 'TKE_bar' ][0][0,1,0],
            freya_dict[ 'TKE_bar' ][0][0,2,0],
            parsed_data[0]['TKE_bar'][0][0,1,0],
            parsed_data[0]['TKE_bar'][0][0,1,1]
            ),2))+
        '\\\\' + '\n \hline'
        '\n\\end{tabular} \]' +
        '$$ \\begin{aligned} \n' +
        '\\text{Times: } & ' +
        '\n\\\\ & ' +
        '\\text{Events} & = &' + str(freya_dict['number_of_events']) +
        '\n\\\\ & ' +
        '\\text{Generation Time} & = & ' + str(round(float( times[ 'gen_time' ] ),5)) + '\\text{ sec}' +
        '\n\\\\ & '
        '\\text{Parse Time} & = &' + str(round(float( times[ 'parse_time'  ] )/60,5)) + '\\text{ min}' +
        '\n\\end{aligned} $$'
        '\n\\end{document}')

    print("Successful.")
    print("Statistics generated (as LaTeX) in: \n"+opt_path)

    print('Begin generating plots of Freya/Data comparison...')

    # plot data and freya on top of each other

    os.chdir(opt_path)

    for key in parsed_data[0]:
        element = parsed_data[0][str(key)]
        if key in plot_location_dictionary.keys():
            label_locations = plot_location_dictionary[key]
        else:
            label_locations = [1,4,1]
        if element is None:
            continue
        if key_translator[key] in ['nubar','average_photon_energy','gammabar','nubar_moments','TKE_bar']:
            continue
        elif len(element[0])<2:
            print("Skipping plot with ",key," on account of a data file which seems to be too short.")
            continue
        elif key_translator[key] not in ["mannhart","n_spectrum","rest_n_spectrum","n_A_TKE"]:
            print('plotting: ',key)
            dim_status = element[5]
            data_array = element[0]
            translated_key = key_translator[key]
            freya_data = freya_dict[translated_key][0]
            ratio_array = chisq_array[3][translated_key]

            matplotlib.rcParams.update({'font.size': 14})
            plt.figure(figsize=(6,6))

            ax = plt.subplot(2,1,1,ylabel=element[3])

            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax.xaxis.set_visible(False)
            plt.plot(freya_data[:,0] , freya_data[:,1] , '^-' , color = 'r' , label = 'FREYA Output') 
            plt.errorbar( freya_data[:,0] , freya_data[:,1] , yerr = freya_data[:,2], color = 'r', fmt = ' ' , capsize = 3, elinewidth = 1)
            plt.plot(data_array[:,0,0] , data_array[:,1,0] , '^-' , color = 'b' , label = str(element[1]))
            plt.errorbar( data_array[:, 0, 0] , data_array[: ,1, 0] , yerr = data_array[:, 1, 1], color = 'b', fmt = ' ' , capsize = 3, elinewidth = 1)
            plt.ylim( ranges_y[key][0] , ranges_y[key][1])
            lg = plt.legend(("FREYA",
                str(element[1])
                ,),fontsize=14,numpoints=1,loc=label_locations[0])
            lg.draw_frame(False)
            if key in ['m_mult_smudge']:
                label = AnchoredText('$(c)$',loc=2,frameon=False)
                ax.add_artist(label)
            else:
                label = AnchoredText('$(a)$',loc=2,frameon=False)
                ax.add_artist(label)
            
            iso_label = AnchoredText('$'+str(iso[4]),frameon=False,
                    loc=label_locations[1])
            ax.add_artist(iso_label)

            ax = plt.subplot(2,1,2,xlabel=element[2],ylabel="C/E",sharex=ax)
            plt.plot(ratio_array[:,0] , ratio_array[:,1] , '^' , color = 'k') 
            plt.errorbar( ratio_array[1:None,0] , ratio_array[1:None,1] , yerr = ratio_array[1:None,2], color = 'k', fmt = ' ' , capsize = 3, elinewidth = 1)
            #  plt.ticklabel_format(style='plain',axis='x',useOffset=False)
            if key in ['n_Af','n_Af_alt','n_TKE','n_TKE_alt']:
            #  if key in ['n_mult','m_mult','n_Af','m_mult_smudge','n_TKE']:
                nonzero_ones = freya_data[np.where(np.nan_to_num(freya_data[:,1]) > 0.0001)]
                plt.xlim( min(nonzero_ones[:,0])-5 , max(nonzero_ones[:,0]))
            elif key in ['n_mult','m_mult','m_mult_smudge']: 
                plt.xlim( ranges_x[key][2] , ranges_x[key][3])
            else:
                plt.xlim( ranges_x[key][0] , ranges_x[key][1])
            if len(ranges_y[key]) > 2:
                plt.ylim( ranges_y[key][2] , ranges_y[key][3] + 0.2)
                print("fixing y limits (" + str(ranges_y[key][2])+"," +str(ranges_y[key][3]) + ") for ratio plot...")
            # Print reduced chi squared on C/E plot
            reduced_chi_sq = chisq_array[5][key]
            if key in ['m_mult_smudge']:
                label = AnchoredText('$(d)$',loc=2,frameon=False)
                ax.add_artist(label)
            else:
                label = AnchoredText('$(b)$',loc=2,frameon=False)
                ax.add_artist(label)
            cs_label = AnchoredText(r'$\chi^2_O = '+str(reduced_chi_sq)+'$',frameon=False,
                    loc=label_locations[2])
            cs_label.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
            cs_label.patch.set_alpha(0.8)
            cs_label.patch.set_edgecolor('0.8')
            ax.add_artist(cs_label)
            plt.axhline(1, color='black')

            plt.subplots_adjust(hspace=0)
            plt.savefig(str(key) + '.pdf',bbox_inches='tight')
            plt.close()

        else:
            print("Plotting log scale plots for: " + str(key) + "...")
            maxwell_temp = 1.32
            dim_status = element[5]
            data_array = element[0]
            translated_key = key_translator[key]
            freya_data = freya_dict[translated_key][0]
            ratio_array = chisq_array[3][translated_key]
            if key == "mannhart":
                ratio_array[:-1,1] = ratio_array[:-1,1] / MaxwellianSpectrum(maxwell_temp)
                ratio_array[:-1,2] = ratio_array[:-1,2] / MaxwellianSpectrum(maxwell_temp)
                data_array[:,1,0] = data_array[:,1,0] * MaxwellianSpectrum(maxwell_temp)
                #  data_array[:,1,1] = data_array[:,1,1]*MaxwellianSpectrum(maxwell_temp)
            if key == "n_spectrum":
                largest_bin = data_array[-1,0,0]
                lb_index_freya= find_nearest(freya_data[:,0],largest_bin)
                freya_bin_differences = np.subtract(
                        freya_data[1:lb_index_freya+2,0],
                        freya_data[:lb_index_freya+1,0])
                freya_areas = freya_bin_differences * freya_data[:lb_index_freya+1,1]
                smaller_area = np.sum(freya_areas) + freya_bin_differences[lb_index_freya] * freya_data[lb_index_freya,1]
                bin_differences = np.subtract(data_array[1:,0,0],data_array[:-1,0,0])
                areas = bin_differences * data_array[:-1,1,0]
                total_count = (np.sum(areas) + bin_differences[-1] * data_array[-1,1,0]) / smaller_area
                data_array[:,1,0] = data_array[:,1,0] / total_count
                data_array[:,1,1] = data_array[:,1,1] / total_count
                ratio_array[:,1] = ratio_array[:,1] * total_count
            if key == "rest_n_spectrum":
                total_count = np.sum(np.multiply(data_array[:,0,0],data_array[:,1,0]))*0.05
                data_array[:,1,0] = data_array[:,1,0] / total_count 
                data_array[:,1,1] = data_array[:,1,1] / total_count
                ratio_array[:,1] = ratio_array[:,1] * total_count
            matplotlib.rcParams.update({'font.size': 14})
### MANNHART REGULAR LOG XSCALE
            plt.figure(figsize=(6,6))
            ax = plt.subplot(2,1,1,ylabel=element[3])
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax.xaxis.set_visible(False)
            plt.plot(freya_data[:,0] , freya_data[:,1] , '^-' , color = 'r' , label = 'FREYA Output') 
            plt.errorbar( freya_data[:,0] , freya_data[:,1] , yerr = freya_data[:,2], color = 'r', fmt = ' ' , capsize = 3, elinewidth = 1)
            plt.plot(data_array[:,0,0] , data_array[:,1,0], '^-' , color = 'b' , label = str(element[1]))
            plt.errorbar( data_array[:, 0, 0] , data_array[: ,1, 0], 
                    yerr = data_array[:, 1, 1], color = 'b', fmt = ' ' , capsize = 3, elinewidth = 1)
            plt.axhline(0, color='black')
            plt.xlim( ranges_x[key][0] , ranges_x[key][1])
            plt.ylim( ranges_y[key][0] , ranges_y[key][1])
            plt.xscale('log')
            iso_label = AnchoredText('$'+str(iso[4]),frameon=False,loc=label_locations[1])
            iso_label.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
            iso_label.patch.set_alpha(0.8)
            iso_label.patch.set_edgecolor('0.8')
            ax.add_artist(iso_label)
            lg = plt.legend(("FREYA",str(element[1])),fontsize=14,numpoints=1,loc=label_locations[0])
            lg.draw_frame(False)
            label = AnchoredText('$(a)$',loc=2,frameon=False)
            ax.add_artist(label)

            ax = plt.subplot(2,1,2,xlabel=element[2],ylabel="C/E",sharex=ax)
            plt.plot(ratio_array[:,0] , ratio_array[:,1] , '^' , color = 'k') 
            plt.errorbar( ratio_array[1:,0] , ratio_array[1:,1] , 
                    yerr = ratio_array[1:,2]
                    , color = 'k', fmt = ' ' , capsize = 3, elinewidth = 1)
            plt.xlim( ranges_x[key][0] , ranges_x[key][1])
            plt.autoscale(axis = 'y')
            plt.ylim( ranges_y[key][2] , ranges_y[key][3])
            plt.xscale('log')
            plt.axhline(1, color='black')
            reduced_chi_sq = chisq_array[5][key]
            cs_label = AnchoredText(r'$\chi^2_O = '+str(reduced_chi_sq)+'$',frameon=False,loc=label_locations[2])
            cs_label.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
            cs_label.patch.set_alpha(0.8)
            cs_label.patch.set_edgecolor('0.8')
            ax.add_artist(cs_label)
            label = AnchoredText('$(b)$',loc=2,frameon=False)
            ax.add_artist(label)

            plt.subplots_adjust(hspace=0)
            plt.savefig(str(key) + "_x" + '.pdf',bbox_inches='tight')
            plt.close()

### MANNHART REGULAR LOG YSCALE

            plt.figure(figsize=(6,6))

            ax = plt.subplot(2,1,1,ylabel=element[3])
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax.xaxis.set_visible(False)
            plt.plot(freya_data[:,0] , freya_data[:,1] , '^-' , color = 'r' , label = 'FREYA Output') 
            plt.errorbar( freya_data[:,0] , freya_data[:,1] , yerr = freya_data[:,2], color = 'r', fmt = ' ' , capsize = 3, elinewidth = 1)
            plt.plot(data_array[:,0,0] , data_array[:,1,0], '^-' , color = 'b' , label = str(element[1]))
            plt.errorbar( data_array[:,0,0] , data_array[:,1,0], 
                    yerr = data_array[:,1,1], color = 'b', fmt = ' ' , capsize = 3, elinewidth = 1)
            plt.axhline(0, color='black')
            plt.xlim( ranges_x[key][0] , ranges_x[key][1])
            #  plt.ylim( ranges_y[key][0] , ranges_y[key][1])
            #  plt.ylim( 1E-6 , ranges_y[key][1])
            plt.ylim(1E-4,3)
            plt.yscale('log')
            iso_label = AnchoredText('$'+str(iso[4]),frameon=False,loc=label_locations[1])
            iso_label.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
            iso_label.patch.set_alpha(0.8)
            iso_label.patch.set_edgecolor('0.8')
            ax.add_artist(iso_label)
            lg = plt.legend(("FREYA",str(element[1])),fontsize=14,numpoints=1,loc=label_locations[0])
            lg.draw_frame(False)
            label = AnchoredText('$(c)$',loc=2,frameon=False)
            ax.add_artist(label)

            ax = plt.subplot(2,1,2,xlabel=element[2],ylabel="C/E",sharex=ax)
            plt.plot(ratio_array[:-2,0] , ratio_array[:-2,1] , '^' , color = 'k') 
            plt.errorbar( ratio_array[1:-2,0] , ratio_array[1:-2,1] , 
                    yerr = ratio_array[1:-2,2], color = 'k', fmt = ' ' , capsize = 3, elinewidth = 1)
            plt.xlim( ranges_x[key][0] , ranges_x[key][1])
            plt.autoscale(axis = 'y')
            plt.ylim( ranges_y[key][2] , ranges_y[key][3])
            plt.axhline(1, color='black')
            reduced_chi_sq = chisq_array[5][key]
            cs_label = AnchoredText(r'$\chi^2_O = '+str(reduced_chi_sq)+'$',frameon=False,loc=label_locations[2])
            ax.add_artist(cs_label)
            label = AnchoredText('$(d)$',loc=2,frameon=False)
            ax.add_artist(label)

            plt.subplots_adjust(hspace=0)
            plt.savefig(str(key) + "_y" + '.pdf',bbox_inches='tight')
            plt.close()

            continue

### plot observables with alternates (double)

        if key in ["n_TKE","n_Af"]:
            dim_status = element[5]
            data_array = element[0]
            translated_key = key_translator[key]
            freya_data = freya_dict[translated_key][0]
            ratio_array = chisq_array[3][translated_key]

            if str(key)+"_alt" in parsed_data[0]:
                print('plotting double',key,"data")
                alt_element = parsed_data[0][key+"_alt"]
                alt_data_array = alt_element[0]
                alt_ratio_array = chisq_array[3][translated_key+"_alt"]
            else:
                alt_element = element
                alt_data_array = data_array
                alt_ratio_array = ratio_array
                continue

            matplotlib.rcParams.update({'font.size': 14})
            plt.figure(figsize=(6,6))

            ax = plt.subplot(2,1,1,ylabel=element[3])

            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax.xaxis.set_visible(False)
            plt.plot(freya_data[:,0] , freya_data[:,1] , '^-' , color = 'r' , label = 'FREYA Output') 
            plt.errorbar( freya_data[:,0] , freya_data[:,1] , yerr = freya_data[:,2], color = 'r', fmt = ' ' , capsize = 3, elinewidth = 1)
            plt.plot(data_array[:,0,0] , data_array[:,1,0] , '^-' , color = 'b' , label = str(element[1]))
            plt.errorbar( data_array[:, 0, 0] , data_array[: ,1, 0] , yerr = data_array[:, 1, 1], color = 'b', fmt = ' ' , capsize = 3, elinewidth = 1)
            plt.plot(alt_data_array[:,0,0] , alt_data_array[:,1,0] , '^-' , color = 'g' , label = str(element[1]))
            plt.errorbar( alt_data_array[:, 0, 0] , alt_data_array[: ,1, 0] , yerr = alt_data_array[:, 1, 1], color = 'g', fmt = ' ' , capsize = 3, elinewidth = 1)
            plt.ylim( ranges_y[key][0] , ranges_y[key][1])
            if key == 'n_TKE':
                iso_label = AnchoredText("$"+iso[4],frameon=False,loc=label_locations[1])
            else:
                iso_label = AnchoredText("$"+iso[4],frameon=False,loc=label_locations[1])
            ax.add_artist(iso_label)
            lg = plt.legend(("FREYA",str(element[1]),str(alt_element[1])),fontsize=14,numpoints=1,loc=label_locations[0])
            lg.draw_frame(False)
            label = AnchoredText('$(a)$',loc=2,frameon=False)
            ax.add_artist(label)

            ax = plt.subplot(2,1,2,xlabel=element[2],ylabel="C/E",sharex=ax)
            plt.plot(ratio_array[:,0] , ratio_array[:,1] , '^' , color = 'b') 
            plt.errorbar( ratio_array[1:None,0] , ratio_array[1:None,1] , yerr = ratio_array[1:None,2], color = 'b', fmt = ' ' , capsize = 3, elinewidth = 1)
            plt.plot(alt_ratio_array[:,0] , alt_ratio_array[:,1] , '^' , color = 'g') 
            plt.errorbar( alt_ratio_array[1:None,0] , alt_ratio_array[1:None,1] , yerr = alt_ratio_array[1:None,2], color = 'g', fmt = ' ' , capsize = 3, elinewidth = 1)
            if len(ranges_y[key]) > 2:
                plt.ylim( ranges_y[key][2] , ranges_y[key][3] + 0.2)
                print("fixing y limits (" + str(ranges_y[key][2])+"," +str(ranges_y[key][3]) + ") for ratio plot...")
            plt.axhline(1, color='black')
            reduced_chi_sq = chisq_array[5][key]
            cs_label = AnchoredText(r'$\chi^2_O = '+str(reduced_chi_sq)+'$',frameon=False,loc=label_locations[2])
            ax.add_artist(cs_label)
            nonzero_ones = freya_data[np.where(np.nan_to_num(freya_data[:,1]) != 0)]
            plt.xlim( min(nonzero_ones[:,0])-15, max(nonzero_ones[:,0])+5)
            label = AnchoredText('$(b)$',loc=2,frameon=False)
            ax.add_artist(label)

            plt.subplots_adjust(hspace=0)
            plt.savefig(str(key)+'_double' + '.pdf',bbox_inches='tight')
            plt.close()

    print('Begin generating plots of Chi-Squared distributions...')

    #  plot chi squared distributions
    if not os.path.exists(opt_path+"/chi_sq"):
        os.makedirs(opt_path+"/chi_sq")
    os.chdir(opt_path+"/chi_sq")
    for key in chisq_array[1]:
        chi_dict = chisq_array[1]
        chi_x = chi_dict[key][:,0]
        chi_y = chi_dict[key][:,1]
        plt.scatter(chi_x,chi_y,color = 'b')
        plt.xlim(ranges_x[key_translator[key]][0], ranges_x[key_translator[key]][1])
        #  plt.xlabel()
        plt.ylabel('uncertainty',fontsize=15)
        if key_translator[key] is 'mannhart': 
            plt.xscale('log')
        plt.savefig('chisq_'+str(key)+'.pdf')
        plt.close()
    print('Finished.')

    os.chdir(cwd)

    print("Plots generated in: \n"+opt_path)

    return 
