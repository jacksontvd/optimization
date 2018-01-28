import os, sys, time
import numpy as np

cwd = os.getcwd()
os.chdir(cwd+'/sub_py')
sys.path.append(cwd+'/sub_py')

from isotope import isotope
from optimize import opt
from variance import variance
from covariance import covariance
from grid_plot import covariance_plot,prerun_plot

Zinput = 98
Ainput = 252
generate_number = 10000
reaction_type = 'spontaneous'
optimization_type = 'grid'
resolution = 100

def var(parameter):
    result = variance(Zinput,
            Ainput,
            generate_number,
            "grid",
            resolution,
            reaction_type = reaction_type,
            parameter = parameter)
    return result

def covar(parameter1 , parameter2):
    result = covariance(Zinput,
            Ainput,
            generate_number,
            "grid",
            resolution,
            reaction_type = reaction_type,
            parameter = parameter1,
            parameter2 = parameter2)
    return result

grid_values , e_var , e_avg = var('e')
grid_values , x_var , x_avg = var('x')
grid_values , c_var , c_avg = var('c')
grid_values , T_var , T_avg = var('T')
grid_values , d_var , d_avg = var('d')

e_x_covariance , grid_values, x_array , y_array = covar('e','x')
e_c_covariance , grid_values, x_array , y_array = covar('e','c')
e_T_covariance , grid_values, x_array , y_array = covar('e','T')
e_d_covariance , grid_values, x_array , y_array = covar('e','d')
x_c_covariance , grid_values, x_array , y_array = covar('x','c')
x_T_covariance , grid_values, x_array , y_array = covar('x','T')
x_d_covariance , grid_values, x_array , y_array = covar('x','d')
c_T_covariance , grid_values, x_array , y_array = covar('c','T')
c_d_covariance , grid_values, x_array , y_array = covar('c','d')
T_d_covariance , grid_values, x_array , y_array = covar('T','d')

print(
        e_x_covariance / (np.sqrt(e_var * x_var))
        )

print(
        e_c_covariance / (np.sqrt(e_var * c_var)),
        x_c_covariance / (np.sqrt(x_var * c_var))
        )

print(
        e_T_covariance / (np.sqrt(e_var * T_var)),
        x_T_covariance / (np.sqrt(x_var * T_var)),
        c_T_covariance / (np.sqrt(c_var * T_var))
        )

print(
        e_d_covariance / (np.sqrt(e_var * d_var)),
        x_d_covariance / (np.sqrt(x_var * d_var)),
        c_d_covariance / (np.sqrt(c_var * d_var)),
        T_d_covariance / (np.sqrt(T_var * d_var))
        )
