import os, sys, time

cwd = os.getcwd()
os.chdir(cwd+'/sub_py')
sys.path.append(cwd+'/sub_py')

from isotope import isotope
from optimize import opt
from variance import variance
from covariance import covariance
# from grid_plot import covariance_plot,prerun_plot

Zinput = 94
Ainput = 240
generate_number = 10000
reaction_type = 'spontaneous'
optimization_type = 'grid'
resolution = 10

def opti(opt_method):
    result = opt(Zinput,
        Ainput, 
        generate_number, 
        opt_method,
        resolution, 
        reaction_type = reaction_type)
    return result

def var(parameter):
    result = variance(Zinput,
            Ainput,
            generate_number,
            "grid",
            resolution,
            reaction_type = reaction_type,
            parameter = parameter)

def covar(parameter1 , parameter2):
    result = covariance(Zinput,
            Ainput,
            generate_number,
            "grid",
            resolution,
            reaction_type = reaction_type,
            parameter = parameter1,
            parameter2 = parameter2)

def plot(parameter1 , parameter2):
    result = covariance_plot(Zinput,
            Ainput, 
            generate_number, 
            "grid",
            resolution,
            reaction_type,
            parameter1,
            parameter2)
    return result

os.chdir(cwd)
os.system('./restore.sh')

if sys.argv[1] is opt:
    print(opti(optimization_type))
elif sys.argv[1] is var:
    arg = sys.argv[2]
    print(var(str(arg)))
elif sys.argv[1] is covar:
    arg1 = sys.argv[2]
    arg2 = sys.argv[3]
    print(covar(str(arg1),str(arg2)))

os.chdir(cwd)
os.system('./restore.sh')

