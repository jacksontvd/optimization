import numpy as np
import time
from math import sqrt as sqrt
from math import pi as pi
import os
from scipy import stats

cwd = os.getcwd()

from gen_par_ana import gpa
from isotope import isotope
from ranges import *

def error(Z, A, e, x , c , T, d, generate_number, data, **kwargs):
    
    #  the data is brought into this function as a list by the optimization.py script
    #  assign the first element of this list (the output dictionary of data as structured in data_parse.py) to be the data_dictionary
    data_dictionary = data[0]

    #  the key translator (stripping names of data files of their numbering so they correspond directly with names of freya observables)
    #  is also a dictionary which is assigned to the name key_trans
    key_trans = data[1]

    # the string for the isotope is assigned to be iso as given by the isotope function in the isotope.py script
    reaction_t = kwargs['reaction_type']
    iso = isotope(Z,A, reaction_t)
    
    if int(iso[1]) is 9 or 10 or 11 or 7 or 5 or 3:
            Energy = '-1'
   
    #  initialize the total error to be 0. This will be added to as we calculate the error for each observable
    total_chisq = 0 
    total_clean_chisq = 0 

    #  assign the output of the generating function to be anal
    #  this is a list as given by the gpa function in gen_par_ana.py
    anal = gpa(Z, A, Energy, 'error_data.dat', e = e, x = x, c = c, T = T, d = d, generate_number = generate_number)
    anal_dict = anal[0]

    print(' Begin generating weighted chi-squared routine...')
    begin_chisq = time.time()

    #  define empty dictionary to give the arrays of chi-squared values for each bin for each data set
    #  this will be used by the optimization script to plot the errors
    chisq_dict = {}
    ratio_dict = {}
    #  define dummy error to scale data with no error for comparison with the other error contribution
    dummy_error = 0.05
    shift_error = 0.01

    keys_for_error = data_dictionary.keys()

    weights = error_weights[str(Z)+str(A)]

    for key in keys_for_error:
        if data_dictionary[key] is None:
            continue
        else:
            #  separate error calculation for n_A_TKE out from the rest of the observables
            if key_trans[key] == 'n_A_TKE':
                continue
            elif key_trans[key] == 'nubar' or key_trans[key] == 'average_photon_energy' or key_trans[key] == 'gammabar':
                print("Calculating error for:" + key)
                data_array = data_dictionary[key][0]
                data_single = data_array[0,1,0]
                single_error = data_array[0,1,1]

                anal_array = anal_dict[key_trans[key]][0]
                freya_single = anal_array[0,1,0]

                single_chisq = (freya_single - data_single)**2 / ((single_error)**2)
                chi_sq_added = single_chisq
                chi_sq_added = chi_sq_added * weights[key]
                print(key,"error:",chi_sq_added)
                total_chisq += np.nan_to_num( chi_sq_added )

                ratio = [freya_single/data_single]
            else:
                print('calculating error for: ' + key)

                #  assign the data_array to be the data array from the data_parse.py script
                #  this is of course the first element of the list for the element of the dictionary given by the current key (see data_parse for the other elements of this list)
                data_array = data_dictionary[key][0]
                data_array_depcolumn = data_array[:,1,0]

                #  the freya output array is given in a similar way to the data array, and we similarly assign 
                #  anal_array to be the corresponding freya array for the current key
                anal_array = anal_dict[key_trans[key]][0]
                anal_array_depcolumn = anal_array[:,1]

                #  initialize the chi_sq_array with zeros
                #  once this is full it will eventually be written into the chi_sq dictionary under the name of the current data file
                chi_sq_array = np.zeros( (len(anal_array)  + 1, 2) )
                dirty_chi_sq_array = np.zeros( (len(anal_array)  + 1, 2) )
                ratio = np.zeros( (len(anal_array)  + 1, 3) )
                ratio[:] = np.nan

                for element in data_array:

                    #  sort the current element of the data array into the binning of the freya array
                    row = np.searchsorted(anal_array[:,0], element[0, 0] )
                    if row >= len(anal_array_depcolumn):
                        row = row - 1

                    #  for the sorted row number of the data in freya, assign that row of the chi-squared array to have first element the same as the data
                    chi_sq_array[row,0] = element[0,0]
                    dirty_chi_sq_array[row,0] = element[0,0]
                    ratio[row,0] = element[0,0]

                    #  if we are dealing with an empty case, assure we have no contributing error
                    if element[1,0] ==0 or element[1,0] == 'NaN' or element[1,0] is None:
                        chi_sq_array[row,1] = 0
                        dirty_chi_sq_array[row,1] = None
                        ratio[row,1] = None

                    #  if we are dealing with a case with no uncertainty, assure we do not get infinite error
                    if element[1,1] == 0 or element[1,1] == 'NaN' or element[1,1] is None:
                    #  clean row
                        #  chi_sq_array[row,1] = (anal_array_depcolumn[row - 1] - element[1,0] )**2
                    #  scaled row
                        chi_sq_array[row,1] = (anal_array_depcolumn[row] - element[1,0] )**2 / ((dummy_error)**2)
                    #  dirty row
                        dirty_chi_sq_array[row,1] = (anal_array_depcolumn[row] - element[1,0] )**2 / ((dummy_error)**2)
                    #   ratio
                        this_ratio = anal_array_depcolumn[row] / element[1,0]
                        ratio[row,1] = this_ratio

                        if anal_array[row,2] is None:
                            anal_array[row,2] = 0
                        if element[1,1] is None:
                            element[1,1] = 0

                        #  potential_ratio_bar = (anal_array[row,1] * element[1,1] + anal_array[row,2] * element[1,0])/(element[1,0]**2 - element[1,1]**2)
                        a = anal_array[row,1]
                        x = anal_array[row,2]
                        b = element[1,0]
                        y = element[1,1]
                        potential_ratio_bar = np.sqrt((x**2 + a**2)/(y**2 + b**2) - a**2 / b**2)
                        ratio[row,2] = potential_ratio_bar

                    #  otherwise, take the square of the difference between the data and freya, and divide by the square of the uncertainty
                    #  this will make uncertain data matter less in the grand scheme of things than the very certain data
                    else:
                    #  clean row
                        chi_sq_array[row,1] = (anal_array_depcolumn[row] - element[1,0] )**2 / (element[1,1]**2)
                    #  shifted row
                        dirty_chi_sq_array[row,1] = (anal_array_depcolumn[row] - element[1,0] )**2 / ((element[1,1] + shift_error)**2)
                    #  ratio
                        this_ratio = anal_array_depcolumn[row] / element[1,0]
                        ratio[row,1] = this_ratio

                        if anal_array[row,2] is None:
                            anal_array[row,2] = 0
                        if element[1,1] is None:
                            element[1,1] = 0

                        #  potential_ratio_bar = (anal_array[row,1] * element[1,1] + anal_array[row,2] * element[1,0])/(element[1,0]**2 - element[1,1]**2)
                        a = anal_array[row,1]
                        x = anal_array[row,2]
                        b = element[1,0]
                        y = element[1,1]
                        potential_ratio_bar = np.sqrt((x**2 + a**2)/(y**2 + b**2) - a**2 / b**2)
                        ratio[row,2] = potential_ratio_bar

                    #  chi_sq_added = chi_sq_array[row,1]

                    #  add every individual chi-squared value to the total. This will be weighted by uncertainties, and give us a goodness of fit for freya as a whole. 
                    #  this ultimately determines whether the current set of parameters is appropriate

                    #  total_chisq +=np.nan_to_num( chi_sq_added )

                if key is 'n_Af':
                    chi_sq_added = np.mean(dirty_chi_sq_array[45:75,1])
                    #  105 - 165
                else:
                    chi_sq_added = np.mean(dirty_chi_sq_array[:,1])

                chi_sq_added = chi_sq_added * weights[key]

                print("average",key,"error: ",chi_sq_added)
                total_chisq += np.nan_to_num( chi_sq_added )

                ratio_dict[key] = ratio
                chisq_dict[key] = chi_sq_array
                clean_chisq_added = np.sum(chi_sq_array)
                total_clean_chisq += np.nan_to_num(clean_chisq_added)
    
    end_chisq = time.time()
    chisq_time = end_chisq - begin_chisq
    print('Time:', chisq_time )
    print('Total Weighted Error: ',total_chisq)
    print('Total Error: ',total_clean_chisq)

    if total_chisq is not 0:
        print('Successful.')

    return total_chisq, chisq_dict , total_clean_chisq , ratio_dict
