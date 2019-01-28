import numpy as np
import os
from copy import copy
from isotope import *
cwd = os.getcwd()

#  this script contains all of the ranges for the plots of given observables. This allows us to standardize the bounds for the plots of any observable across the board for freya, data, and error. 

parameter_labels = {}
parameter_labels['c'] = 'c'
parameter_labels['e'] = 'e_0'
parameter_labels['x'] = 'x'
parameter_labels['T'] = 'c_S'
parameter_labels['d'] = 'dTKE'

param_ranges = {}

param_ranges['e0'] = [7,12,0]
param_ranges['e'] = param_ranges['e0']

param_ranges['x'] = [1,1.5,1]

param_ranges['c'] = [1,3,2]

param_ranges['cS'] = [0.5, 1.5,3]
param_ranges['T'] = param_ranges['cS'] 

param_ranges['dTKE'] = [-5, 5,4]
param_ranges['d'] = param_ranges['dTKE']

def param_list(Z,A,reac_type):
    iso = isotope(int(Z),int(A),reac_type = reac_type)
    i = iso[1]
    infile = open(cwd+'/../inputparameters.dat','r')
    content = infile.readlines() #reads line by line and outputs a list of each line
    line = content[i]
    elems = line.split()
    plist = [float(elems[3]),float(elems[4]),float(elems[5]),float(elems[6]),float(elems[9])]
    #  print("Parameter list for Z=",Z," A=",A," is: ",plist)
    return plist

def param_var(Z,A,reac_type):
    iso = isotope(int(Z),int(A),reac_type = reac_type)
    if iso[0] == "Cf252sf":
        varlist = [1.090,0.187,0.362,0.020,0.078]
    if iso[0] == "Pu240sf":
        varlist = [0.138,0.071,0.355,0.023,0.112]
    return varlist

#  Here we hardcode the average photon energy according to FREYA (Using parameters from arXiv:1809.05587)
#  We do this, because the (incident neutron energy dependent) average photon energy data for U235(nf)+n is given as the ratio to the average photon energy for Cf, and in the process of comparing the FREYA output for U235(nf), we don't want to keep running Cf each time we compare.
cf_ape = 0.9114
#  \text{Average photon energy} & : & 0.911401600054\pm0.858108023313

#  When we add the chi squared for each observable together, we might require certain weights. I.e. one might be too many orders of magnitude larger (resp. smaller) than the others, and therefore we need to scale it down in order to properly optimize
#  Therefore we put all of these weights in this dictionary

error_weights = {}
cf252_weights = {}
cf252_weights['m_mult'] = 4
cf252_weights['m_mult_smudge'] = 0
cf252_weights['n_Af'] = 1
cf252_weights['totn_Af'] = 0
cf252_weights['n_mult'] = 100
cf252_weights['n_TKE'] = 1
cf252_weights['m_TKE'] = 1
cf252_weights['n_TKE_alt'] = 0
cf252_weights['n_Af_alt'] = 0
cf252_weights['n_A_TKE'] = 0
cf252_weights['mannhart'] = 15
cf252_weights['n_spectrum'] = 0
cf252_weights['rest_n_spectrum'] = 0
cf252_weights['product_A'] = 1
cf252_weights['nubar'] = 2
cf252_weights['nubar_moments'] = 2
cf252_weights['gammabar'] = 1
cf252_weights['average_photon_energy'] = 1
cf252_weights['TKE_bar'] = 1
error_weights['98252sf'] = cf252_weights

pu240sf_weights = {}
pu240sf_weights['m_mult'] = 1
pu240sf_weights['m_mult_smudge'] = 0
pu240sf_weights['n_Af'] = 1
pu240sf_weights['n_mult'] = 100
pu240sf_weights['n_TKE'] = 1
pu240sf_weights['m_TKE'] = 1
pu240sf_weights['n_TKE_alt'] = 1
pu240sf_weights['n_Af_alt'] = 1
pu240sf_weights['n_spectrum'] = 0
pu240sf_weights['product_A'] = 1
pu240sf_weights['nubar'] = 1
pu240sf_weights['TKE_bar'] = 1
pu240sf_weights['nubar_moments'] = 1
pu240sf_weights['gammabar'] = 2
pu240sf_weights['average_photon_energy'] = 1
error_weights['94240sf'] = pu240sf_weights

cm244_weights = {}
cm244_weights['gammabar'] = 0
cm244_weights['average_photon_energy'] = 0
cm244_weights['n_spectrum'] = 1
cm244_weights['n_mult'] = 20
cm244_weights['nubar_moments'] = 4
cm244_weights['nubar'] = 4
cm244_weights['TKE_bar'] = 4
cm244_weights['totn_Af'] = 1
cm244_weights['n_Af'] = 1
error_weights['96244sf'] = cm244_weights

u238sf_weights = {}
u238sf_weights['gammabar'] = 0
u238sf_weights['average_photon_energy'] = 0
u238sf_weights['nubar_moments'] = 10
u238sf_weights['nubar'] = 10
u238sf_weights['TKE_bar'] = 10
u238sf_weights['n_mult'] = 20
error_weights['92238sf'] = u238sf_weights

pu238sf_weights = {}
pu238sf_weights['gammabar'] = 0
pu238sf_weights['average_photon_energy'] = 0
pu238sf_weights['n_mult'] = 7
pu238sf_weights['nubar_moments'] = 1
pu238sf_weights['nubar'] = 1
error_weights['94238sf'] = pu238sf_weights

pu242sf_weights = {}
pu242sf_weights['gammabar'] = 1
pu242sf_weights['average_photon_energy'] = 1
pu242sf_weights['energy_per_photon'] = 1
pu242sf_weights['total_photon_energy'] = 1
pu242sf_weights['n_mult'] = 5
pu242sf_weights['m_Af'] = 5
pu242sf_weights['n_Af'] = 1
pu242sf_weights['nubar_moments'] = 1
pu242sf_weights['nubar'] = 1
pu242sf_weights['TKE_bar'] = 1
pu242sf_weights['n_spectrum'] = 1
pu242sf_weights['m_TKE'] = 1
error_weights['94242sf'] = pu242sf_weights

pu241nf_weights = {}
pu241nf_weights['average_photon_energy'] = 1
pu241nf_weights['gammabar'] = 1
pu241nf_weights['n_mult'] = 1
pu241nf_weights['nubar_moments'] = 1
error_weights['94241(n,f)'] = pu241nf_weights

pu242nf_weights = {}
pu242nf_weights['average_photon_energy'] = 1
pu242nf_weights['gammabar'] = 1
pu242nf_weights['n_mult'] = 1
pu242nf_weights['nubar_moments'] = 1
error_weights['94242(n,f)'] = pu242nf_weights

u239nf_weights = {}
u239nf_weights['average_photon_energy'] = 1
u239nf_weights['gammabar'] = 1
u239nf_weights['n_mult'] = 1
u239nf_weights['nubar_moments'] = 1
error_weights['92238(n,f)'] = u238nf_weights

u235nf_weights = {}
u235nf_weights['average_photon_energy'] = 1/100
u235nf_weights['gammabar'] = 1
u235nf_weights['m_Af'] = 1
u235nf_weights['m_TKE'] = 100
u235nf_weights['n_Af'] = 1
u235nf_weights['n_mult'] = 100
u235nf_weights['n_spectrum'] = 100
u235nf_weights['nubar_moments'] = 1/100
u235nf_weights['nubar'] = 10
u235nf_weights['TKE_bar'] = 1
u235nf_weights['total_photon_energy'] = 1
error_weights['92235(n,f)'] = u235nf_weights


#  The dictionary ranges_x, will take the string associated to an observable, and return the range for the $x$-axis.
ranges_x = {} 

ranges_x['A'] = [60,180]

ranges_x['Product_A'] = ranges_x['A']
ranges_x['Fragment_A'] = np.copy(ranges_x['Product_A'])

ranges_x['n_Af'] = np.copy(ranges_x['Product_A'])
ranges_x['totn_Af'] = [ranges_x['A'][0] + np.floor(ranges_x['A'][1]/3),
        ranges_x['A'][1],120,170]
ranges_x['n_Af_alt'] = np.copy(ranges_x['Product_A'])
ranges_x['m_Af'] = np.copy(ranges_x['Product_A'])
ranges_x['TKE_A'] = np.copy(ranges_x['Product_A'])

ranges_x['total_photon_energy'] = np.copy(ranges_x['Product_A'])
ranges_x['energy_per_photon'] = np.copy(ranges_x['Product_A'])

ranges_x['n_angular'] = [-1,1]

ranges_x['n_mult'] = [0,14,0,10]
ranges_x['m_mult'] = [0,70,0,22]
ranges_x['m_mult_smudge'] = [0,70,0,22]

ranges_x['m_TKE'] = [100,240]
ranges_x['n_TKE'] = [100,240]
ranges_x['n_TKE_alt'] = [100,240]
#  ranges_x['n_spectrum'] = [0,15]
ranges_x['n_spectrum'] = [0.2500E-01, 19.975]
ranges_x['rest_n_spectrum'] = [0.2500E-01, 19.975]
ranges_x['m_spectrum'] = [0,10]
#  ranges_x['mannhart'] = [0.01, 100]
#  ranges_x['mannhart'] = [0.2500E-01, 0.1280E+02]
ranges_x['mannhart'] = [0.2500E-01, 0.130E+02]
ranges_x['n_A_TKE'] = np.copy(ranges_x['A'])
ranges_x['n_TKE_A'] = np.copy(ranges_x['A'])

ranges_x['TKE_bar'] = [0,2]
ranges_x['nubar'] = [0,2]
ranges_x['nubar_moments'] = [0,4]
ranges_x['gammabar'] = [0,2]
ranges_x['average_photon_energy'] = [0,2]

ranges_y = {}

ranges_y['Product_A'] = [0,12]
ranges_y['Fragment_A'] = np.copy(ranges_y['Product_A'])

#  ranges_y['n_Af'] = [0,7,0,2.5]
ranges_y['n_Af'] = [None,None,0,2.5]
ranges_y['totn_Af'] = [0,5,0,2.5]
ranges_y['n_Af_alt'] = [0,7,0,2.5]
ranges_y['m_Af'] = [0,8]
ranges_y['TKE_A'] = [100,240]

ranges_y['total_photon_energy'] = [ 0,12 ]
ranges_y['energy_per_photon'] = [ 0,3.7 ]

ranges_y['n_angular'] = [0.2,1.2]

ranges_y['n_mult'] = [0,0.6,0,2]
# third and fourth elements are ranges for ratio of freya to data
ranges_y['m_mult'] = [0,0.25,0,2]
ranges_y['m_mult_smudge'] = [0,0.25,0,2]

ranges_y['m_TKE'] = [0,12]
ranges_y['n_TKE'] = [0,12]
ranges_y['n_TKE_alt'] = [0,12]
#  ranges_y['n_spectrum'] = [0,0.55]
ranges_y['n_spectrum'] = [None,None,None,None]
ranges_y['rest_n_spectrum'] = [None,None,None,None]
ranges_y['m_spectrum'] = [0,1.3]
ranges_y['mannhart'] = [-0.2,0.5,0,2.9]
ranges_y['n_A_TKE'] = [100, 300]
ranges_y['n_TKE_A'] = [100, 300]

ranges_y['TKE_bar'] = [100,200]
ranges_y['nubar'] = [0,10]
ranges_y['nubar_moments'] = [0,40]
ranges_y['gammabar'] = [0,50]
ranges_y['average_photon_energy'] = [0,3.7]


ranges_z = {}

ranges_z['n_A_TKE'] = [None, None]


output_keys = ranges_x.keys()


bin_number = {}

bin_number['n_A_TKE'] = 400
bin_number['n_TKE_A'] = 400
bin_number['n_TKE'] = 70
bin_number['n_TKE_alt'] = 70
bin_number['m_TKE'] = 70
bin_number['n_angular'] = 20
#  bin_number['n_spectrum'] = 400
bin_number['n_spectrum'] = 1500
bin_number['rest_n_spectrum'] = 400
bin_number['m_spectrum'] = 20
bin_number['mannhart'] = 62


bin_width = {}

for key in bin_number.keys():
    bin_width[str(key)] = (ranges_x[str(key)][1] - ranges_x[str(key)][0]) / bin_number[str(key)]


mannhart_bins = [0.2500E-01,0.4500E-01,0.6500E-01,0.8500E-01,0.1050E+00,0.1250E+00,0.1500E+00,0.1800E+00,0.2100E+00,0.2400E+00,0.2800E+00,0.3300E+00,0.3800E+00,0.4300E+00,0.4800E+00,0.5300E+00,0.5800E+00,0.6300E+00,0.6800E+00,0.7300E+00,0.7800E+00,0.8300E+00,0.8800E+00,0.9300E+00,0.1003E+01,0.1100E+01,0.1200E+01,0.1300E+01,0.1400E+01,0.1500E+01,0.1600E+01,0.1700E+01,0.1800E+01,0.1900E+01,0.2050E+01,0.2250E+01,0.2450E+01,0.2650E+01,0.2850E+01,0.3100E+01,0.3400E+01,0.3700E+01,0.4000E+01,0.4300E+01,0.4600E+01,0.4900E+01,0.5300E+01,0.5800E+01,0.6300E+01,0.6800E+01,0.7300E+01,0.7800E+01,0.8300E+01,0.8800E+01,0.9300E+01,0.9800E+01,0.1030E+02,0.1080E+02,0.1130E+02,0.1180E+02,0.1230E+02,0.1280E+02]

# empty dictionary which will hold the locations of the labels on plots
plot_location_dictionary = {}

#  Each key will correspond to an observable.
#  Each entry will be an array with three integers
#  The first (resp. second, third) integer specifies the location of the legend (resp. isotope, chi-squared error)
#  The locations are as follows:
      #  upper right 	1
      #  upper left 	2
      #  lower left 	3
      #  lower right 	4
      #  right          5
      #  center left 	6
      #  center right   7
      #  lower center   8
      #  upper center   9
      #  center 	10

#  cf:
#  plot_location_dictionary['n_mult'] = [1,6,9]
#  pu238
#  plot_location_dictionary['n_mult'] = [1,4,1]
#  pu240
#  plot_location_dictionary['n_mult'] = [1,7,1]
#  pu242
#  plot_location_dictionary['n_mult'] = [9,10,1]
#  cm
#  plot_location_dictionary['n_mult'] = [1,10,3]
#  u
plot_location_dictionary['n_mult'] = [1,10,1]

plot_location_dictionary['mannhart'] = [1,3,3]
plot_location_dictionary['n_TKE'] = [1,3,9]
plot_location_dictionary['m_mult'] = [1,7,1]
plot_location_dictionary['m_mult'] = [1,7,1]
plot_location_dictionary['m_mult_smudge'] = [1,7,1]
plot_location_dictionary['n_Af'] = [1,9,1]
plot_location_dictionary['n_spectrum'] = [1,8,1]
