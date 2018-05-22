import os, sys, time

cwd = os.getcwd()
os.chdir(cwd+'/sub_py')
sys.path.append(cwd+'/sub_py')

from isotope import isotope
from optimize import opt
from variance import variance
from covariance import covariance
from grid_plot import *
from plot import plot
from covmatrix import freya_hessian

Zinput = 98
#  Zinput = 94
Ainput = 252
#  Ainput = 240
#  generate_number = 1000000
#  generate_number = 100000
generate_number = 10000
reaction_type = 'spontaneous'
optimization_type = 'grid'
#  optimization_type = 'anneal'
#  optimization_type = 'process'
stochastic_type = 0
resolution = 4
hessian_h = [0.5,0.1,0.1,0.1,0.5]
#  hessian_h = 1.e-1

def opti(opt_method):
    result = opt(Zinput,
        Ainput, 
        generate_number, 
        opt_method,
        resolution, 
        stochastic_type = stochastic_type,
        reaction_type = reaction_type)
    return result

def cluster_plot():
    result = plot(Zinput, Ainput, "-1", "plot_output",generate_number = generate_number)
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

def hess():
    result = freya_hessian(Zinput,
            Ainput,
            generate_number,
            hessian_h,
            reaction_type)
    return result

def cov_plot(parameter1 , parameter2):
    result = covariance_plot(Zinput,
            Ainput, 
            generate_number, 
            "grid",
            resolution,
            reaction_type,
            parameter1,
            parameter2)
    return result

def pre_plot():
    prerun_plot(str(Zinput),str(Ainput),cwd+'/../output/grids/'+str(Zinput)+str(Ainput)+'-1/Td',"T","d")
    prerun_plot(str(Zinput),str(Ainput),cwd+'/../output/grids/'+str(Zinput)+str(Ainput)+'-1/cT',"c","T")
    prerun_plot(str(Zinput),str(Ainput),cwd+'/../output/grids/'+str(Zinput)+str(Ainput)+'-1/cd',"c","d")
    prerun_plot(str(Zinput),str(Ainput),cwd+'/../output/grids/'+str(Zinput)+str(Ainput)+'-1/eT',"e","T")
    prerun_plot(str(Zinput),str(Ainput),cwd+'/../output/grids/'+str(Zinput)+str(Ainput)+'-1/ec',"e","c")
    prerun_plot(str(Zinput),str(Ainput),cwd+'/../output/grids/'+str(Zinput)+str(Ainput)+'-1/ed',"e","d")
    prerun_plot(str(Zinput),str(Ainput),cwd+'/../output/grids/'+str(Zinput)+str(Ainput)+'-1/ex',"e","x")
    prerun_plot(str(Zinput),str(Ainput),cwd+'/../output/grids/'+str(Zinput)+str(Ainput)+'-1/xT',"x","T")
    prerun_plot(str(Zinput),str(Ainput),cwd+'/../output/grids/'+str(Zinput)+str(Ainput)+'-1/xc',"x","c")
    prerun_plot(str(Zinput),str(Ainput),cwd+'/../output/grids/'+str(Zinput)+str(Ainput)+'-1/xd',"x","d")
    return 

def well():
    #  well_plot(str(Zinput),str(Ainput),cwd+'/output/grids/'+str(Zinput)+str(Ainput)+'-1/Td',cwd+'/output/grids/'+str(Zinput)+str(Ainput)+'-1/Te',"T","d")
    well_plot(str(Zinput),str(Ainput),
            cwd+'/../output/grids/'+str(Zinput)+str(Ainput)+'-1/cd',
            cwd+'/../output/grids/'+str(Zinput)+str(Ainput)+'-1/Td',
            "c","d","T","d")

def well_2():
    #  well_plot(str(Zinput),str(Ainput),cwd+'/output/grids/'+str(Zinput)+str(Ainput)+'-1/Td',cwd+'/output/grids/'+str(Zinput)+str(Ainput)+'-1/Te',"T","d")
    well_plot_2(str(Zinput),str(Ainput),
            cwd+'/../output/grids/'+str(Zinput)+str(Ainput)+'-1/eT',
            cwd+'/../output/grids/'+str(Zinput)+str(Ainput)+'-1/Td',
            cwd+'/../output/grids/'+str(Zinput)+str(Ainput)+'-1/cT',
            "c","d","T","d")

os.chdir(cwd+"/..")
os.system('./restore.sh')

if len(sys.argv) > 1:
    job_type = str(sys.argv[1])
else: 
    job_type = None

if job_type == 'opt':
    opti(optimization_type)
elif job_type == 'plot':
    print(cluster_plot())
elif job_type == 'var':
    print("we are here")
    arg = sys.argv[2]
    print(var(str(arg)))
elif job_type == 'covar':
    arg1 = sys.argv[2]
    arg2 = sys.argv[3]
    print(covar(str(arg1),str(arg2)))
elif job_type == 'hessian':
    print(hess())
elif job_type == 'pre_plot':
    pre_plot()
elif job_type == 'well':
    well()
elif job_type == 'well2':
    well_2()
else:
    print("ERROR: Job type not recognized.")

os.chdir(cwd+"/..")
os.system('./restore.sh')
