
import numpy as np
import os, sys

cwd = os.getcwd()
os.chdir(cwd+'/sub_py')
sys.path.append(cwd+'/sub_py')

from plot import plot
from isotope import isotope
from optimize import opt
from post_op import post_opt

if len(sys.argv) > 1:
    job_type = str(sys.argv[1])
else:
    job_type = None

if job_type == 'process':
    print("Processing...")
    Zinput=98
    Ainput= 252
    #  Zinput=94
    #  Ainput=240
    generate_number = 1000000
    method='process'
    resolution=None
    reaction_type="spontaneous"
    post_opt(Zinput , Ainput , generate_number , method , resolution , reaction_type = reaction_type , stochastic_type = 1 , iteration_number = 100)

    taskl = None

else:
    taskraw = input("What would you like to do with FREYA? (run/optimize): ")
    task = str(taskraw)
    taskl = task.lower()

if taskl == 'run' or taskl == 'run ' or taskl == ' run':
    Zinputraw = input('Z = ')
    Zinput = int(Zinputraw)
    Ainputraw = input('A = ')
    Ainput = int(Ainputraw)
    energyraw = input('Value of energy (in units of MeV, -1 for spontaneous fission): ')
    energyinput = float(energyraw)
    path_raw = input('Name of output file and location of plots: ')
    path_input = str(path_raw)
    generate_number_raw = input('How many events (Suggested 1,000,000 for high quality results): ')
    generate_number = int(generate_number_raw)
    paramraw = input("Would you like to tweak any of Freya's parameters for this run? (Yes/No): ")
    param = str(paramraw)
    paraml = param.lower()
    if paraml == 'yes':
        einput = input('e = ')
        e = float(einput)
        xinput = input('x = ')
        x = float(xinput)
        cinput = input('c = ')
        c = float(cinput)
        Tinput = input('cS = ')
        T = float(Tinput)
        dinput = input('dTKE = ')
        d = float(dinput)
        plot(Zinput, Ainput, energyinput, path_input, e = e, x = x, c = c, T = T, d = d, generate_number = generate_number)
    elif paraml == 'no':
        plot(Zinput, Ainput, energyinput, path_input, generate_number = generate_number)
    else:
        print('ERROR: Input not supported...')

elif taskl == 'optimize':
    print('Which isotope would you like to optimize with respect to?')
    Zinputraw = input('Z = ')
    Zinput = int(Zinputraw)

    Ainputraw = input('A = ')
    Ainput = int(Ainputraw)
    if Ainput == 240 or 242:
        reaction_type = input('Spontaneous or induced: ')
        reaction_type = reaction_type.lower()
    else:
        reaction_type = isotope(Z,A)[3]
    print('Which method of optimization do you prefer?')
    method_raw = input('Stochastic or Grid: ')
    method_input = str(method_raw)
    method = method_input.lower()
    generate_number_raw = input('How many events per iteration (Suggested 1,000,000 for high quality results): ')
    generate_number = int(generate_number_raw)
    optimize_bool = True
    if method == 'grid':
        resolution_raw = input('How many values should be tested within the range of each parameter: ')
        resolution = resolution_raw.lower()
    else:
        resolution = None
    if method == 'process':
        post_opt(Zinput , Ainput , generate_number , method , resolution , reaction_type = reaction_type , stochastic_type = 1 , iteration_number = 100)
        optimize_bool = False
    else:
        opt(Zinput,Ainput, generate_number, method, resolution, reaction_type = reaction_type,stochastic_type = 1, iteration_number = 100)

else:
    print('ERROR: Input is not supported')
