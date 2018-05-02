import numpy as np

#  this script contains all of the ranges for the plots of given observables. This allows us to standardize the bounds for the plots of any observable across the board for freya, data, and error. 

parameter_labels = {}
parameter_labels['c'] = 'c'
parameter_labels['e'] = 'e_0'
parameter_labels['x'] = 'x'
parameter_labels['T'] = 'c_S'
parameter_labels['d'] = 'dTKE'

param_ranges = {}

param_ranges['c'] = [1,2,2]

param_ranges['e0'] = [7,12,0]
param_ranges['e'] = param_ranges['e0']

param_ranges['x'] = [1,1.5,1]

param_ranges['dTKE'] = [-5, 5,4]
param_ranges['d'] = param_ranges['dTKE']

param_ranges['cS'] = [0.5, 1.5,3]
param_ranges['T'] = param_ranges['cS'] 

param_ranges['98252'] = [10.37,1.27,1.18,0.87,0.52]
param_ranges['94240'] = [10.0724,1.3,3.0,0.87,-3.07119]


error_weights = {}
cf252_weights = {}
cf252_weights['m_mult'] = 10/2.8
cf252_weights['m_mult_smudge'] = 0
cf252_weights['n_Af'] = 1/100
cf252_weights['n_mult'] = 1/1.7
cf252_weights['n_TKE'] = 1/1.8
cf252_weights['n_TKE_alt'] = 0
cf252_weights['n_A_TKE'] = 0
cf252_weights['mannhart'] = 15
cf252_weights['n_spectrum'] = 0
cf252_weights['rest_n_spectrum'] = 0
cf252_weights['product_A'] = 1
cf252_weights['nubar'] = 2
cf252_weights['nubar_moments'] = 2
cf252_weights['gammabar'] = 1
cf252_weights['average_photon_energy'] = 1
error_weights['98252'] = cf252_weights
pu240sf_weights = {}
pu240sf_weights['m_mult'] = 1
pu240sf_weights['m_mult_smudge'] = 0
pu240sf_weights['n_Af'] = 1
pu240sf_weights['n_mult'] = 1
#  pu240sf_weights['n_mult'] = 10
pu240sf_weights['n_TKE'] = 1
pu240sf_weights['n_TKE_alt'] = 1
#  pu240sf_weights['n_spectrum'] = 2.5
pu240sf_weights['n_spectrum'] = 0
pu240sf_weights['product_A'] = 1
#  pu240sf_weights['nubar'] = 1/170
pu240sf_weights['nubar'] = 1
pu240sf_weights['gammabar'] = 1/18
#  pu240sf_weights['average_photon_energy'] = 1/8
pu240sf_weights['average_photon_energy'] = 1/4
error_weights['94240'] = pu240sf_weights


ranges_x = {} 

ranges_x['A'] = [60,180]

ranges_x['Product_A'] = ranges_x['A']
ranges_x['Fragment_A'] = np.copy(ranges_x['Product_A'])

ranges_x['n_Af'] = np.copy(ranges_x['Product_A'])
ranges_x['m_Af'] = np.copy(ranges_x['Product_A'])
ranges_x['TKE_A'] = np.copy(ranges_x['Product_A'])

ranges_x['total_photon_energy'] = np.copy(ranges_x['Product_A'])
ranges_x['energy_per_photon'] = np.copy(ranges_x['Product_A'])

ranges_x['n_angular'] = [-1,1]

ranges_x['n_mult'] = [0,14,0,10]
ranges_x['m_mult'] = [0,40,0,22]
ranges_x['m_mult_smudge'] = [0,40,0,22]

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

ranges_x['nubar'] = [0,2]
ranges_x['nubar_moments'] = [0,4]
ranges_x['gammabar'] = [0,2]
ranges_x['average_photon_energy'] = [0,2]

ranges_y = {}

ranges_y['Product_A'] = [0,12]
ranges_y['Fragment_A'] = np.copy(ranges_y['Product_A'])

ranges_y['n_Af'] = [0,5,0,2.5]
ranges_y['m_Af'] = [0,8]
ranges_y['TKE_A'] = [100,240]

ranges_y['total_photon_energy'] = [ 0,12 ]
ranges_y['energy_per_photon'] = [ 0,3.7 ]

ranges_y['n_angular'] = [0.2,1.2]

ranges_y['n_mult'] = [0,0.6,0,2]
# third and fourth elements are ranges for ratio of freya to data
ranges_y['m_mult'] = [0,0.25,0,2]
ranges_y['m_mult_smudge'] = [0,0.25,0,2]

ranges_y['n_TKE'] = [0,12]
ranges_y['n_TKE_alt'] = [0,12]
#  ranges_y['n_spectrum'] = [0,0.55]
ranges_y['n_spectrum'] = [None,None,None,None]
ranges_y['rest_n_spectrum'] = [None,None,None,None]
ranges_y['m_spectrum'] = [0,1.3]
#  ranges_y['mannhart'] = [0.01,2]
#  ranges_y['mannhart'] = [0.01,None]
ranges_y['mannhart'] = [-0.2,0.5,0,2.2]
ranges_y['n_A_TKE'] = [100, 300]

ranges_y['nubar'] = [0,10]
ranges_y['nubar_moments'] = [0,40]
ranges_y['gammabar'] = [0,50]
ranges_y['average_photon_energy'] = [0,3.7]


ranges_z = {}

ranges_z['n_A_TKE'] = [None, None]


output_keys = ranges_x.keys()


bin_number = {}

bin_number['n_A_TKE'] = 400
bin_number['n_TKE'] = 70
bin_number['n_TKE_alt'] = 70
bin_number['n_angular'] = 20
bin_number['n_spectrum'] = 400
bin_number['rest_n_spectrum'] = 400
bin_number['m_spectrum'] = 20
bin_number['mannhart'] = 62


bin_width = {}

for key in bin_number.keys():
    bin_width[str(key)] = (ranges_x[str(key)][1] - ranges_x[str(key)][0]) / bin_number[str(key)]


mannhart_bins = [0.2500E-01,0.4500E-01,0.6500E-01,0.8500E-01,0.1050E+00,0.1250E+00,0.1500E+00,0.1800E+00,0.2100E+00,0.2400E+00,0.2800E+00,0.3300E+00,0.3800E+00,0.4300E+00,0.4800E+00,0.5300E+00,0.5800E+00,0.6300E+00,0.6800E+00,0.7300E+00,0.7800E+00,0.8300E+00,0.8800E+00,0.9300E+00,0.1003E+01,0.1100E+01,0.1200E+01,0.1300E+01,0.1400E+01,0.1500E+01,0.1600E+01,0.1700E+01,0.1800E+01,0.1900E+01,0.2050E+01,0.2250E+01,0.2450E+01,0.2650E+01,0.2850E+01,0.3100E+01,0.3400E+01,0.3700E+01,0.4000E+01,0.4300E+01,0.4600E+01,0.4900E+01,0.5300E+01,0.5800E+01,0.6300E+01,0.6800E+01,0.7300E+01,0.7800E+01,0.8300E+01,0.8800E+01,0.9300E+01,0.9800E+01,0.1030E+02,0.1080E+02,0.1130E+02,0.1180E+02,0.1230E+02,0.1280E+02]
