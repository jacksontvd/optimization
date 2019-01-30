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
from grid_run import grid_run

#  cf
#  Zinput = 98
#  cm
#  Zinput=96
#  pu
#  Zinput = 94
#  u
Zinput=92

#  Ainput= 233
#  Ainput = 238
#  Ainput = 239
#  Ainput = 240
#  Ainput = 241
#  Ainput = 242
#  Ainput = 244
#  Ainput = 252
Ainput = 235
#  Ainput = 238

#  generate_number = 1000000
generate_number = 500000
#  generate_number = 100000
#  generate_number = 10000
#  generate_number = 1

#  reaction_type = 'spontaneous'
reaction_type = 'induced'
#  thermal
#  energies = [1E-10]
#  Adams
#  energies = [0.52]
#  Boykov
#  energies = [2.9,14.7]
#  chinu
#  energies=[0.5]
#  energies = [0.5,1,1.5,2,2.5,3,3,4,5,5.5,6,7,8,9,10,11,11.5,12,13,13.5,14,15,17.5]
#  Conde
#  energies = [1.5]
#  Knitter
#  energies = [1.5]
#  Staples
#  energies = [0.5,2.5]

#  optimization_type = 'grid'
optimization_type = 'anneal'
#  optimization_type = 'bypass'
#  optimization_type = 'basinhopping'
#  optimization_type = 'process'
stochastic_type = 0

#  resolution = 4
resolution = 10

machine_epsilon = 2.22044604925E-16
#  machine_epsilon = 7./3 - 4./3 -1
hessian_h = np.array(param_list(Zinput,Ainput,reaction_type))*np.sqrt(machine_epsilon)

def opti(opt_method,neutron_energy):
    result = opt(Zinput,
        Ainput, 
        generate_number, 
        opt_method,
        resolution, 
        stochastic_type = stochastic_type,
        reaction_type = reaction_type,
        energy=neutron_energy
        )
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
    prerun_plot(str(Zinput),str(Ainput),cwd+'/../TRANSFER_FOLDER/grid_values/'+str(Zinput)+str(Ainput)+'-1/Td',"T","d")
    prerun_plot(str(Zinput),str(Ainput),cwd+'/../TRANSFER_FOLDER/grid_values/'+str(Zinput)+str(Ainput)+'-1/cT',"c","T")
    prerun_plot(str(Zinput),str(Ainput),cwd+'/../TRANSFER_FOLDER/grid_values/'+str(Zinput)+str(Ainput)+'-1/cd',"c","d")
    prerun_plot(str(Zinput),str(Ainput),cwd+'/../TRANSFER_FOLDER/grid_values/'+str(Zinput)+str(Ainput)+'-1/eT',"e","T")
    prerun_plot(str(Zinput),str(Ainput),cwd+'/../TRANSFER_FOLDER/grid_values/'+str(Zinput)+str(Ainput)+'-1/ec',"e","c")
    prerun_plot(str(Zinput),str(Ainput),cwd+'/../TRANSFER_FOLDER/grid_values/'+str(Zinput)+str(Ainput)+'-1/ed',"e","d")
    prerun_plot(str(Zinput),str(Ainput),cwd+'/../TRANSFER_FOLDER/grid_values/'+str(Zinput)+str(Ainput)+'-1/ex',"e","x")
    prerun_plot(str(Zinput),str(Ainput),cwd+'/../TRANSFER_FOLDER/grid_values/'+str(Zinput)+str(Ainput)+'-1/xT',"x","T")
    prerun_plot(str(Zinput),str(Ainput),cwd+'/../TRANSFER_FOLDER/grid_values/'+str(Zinput)+str(Ainput)+'-1/xc',"x","c")
    prerun_plot(str(Zinput),str(Ainput),cwd+'/../TRANSFER_FOLDER/grid_values/'+str(Zinput)+str(Ainput)+'-1/xd',"x","d")
    return 

#  [string,weight,parameter_index_1,parameter_index_2,parameter_name_1,parameter_name_2]
#  well1_1 = ['cd',1,   2,4,'c','dTKE']
#  well1_2 = ['Td',1,   3,4,'c_S','dTKE']
well1_1 = ['xT',1,  1,3,'x','c_S']
well1_2 = ['eT',1,  0,3,'e_0','c_S']
zoom = 50

def well():
    well_plot(str(Zinput),str(Ainput),
            cwd+'/../TRANSFER_FOLDER/grid_values/'+str(Zinput)+str(Ainput)+'-1/'+well1_1[0],
            cwd+'/../TRANSFER_FOLDER/grid_values/'+str(Zinput)+str(Ainput)+'-1/'+well1_2[0],
            well1_1,well1_2,zoom)

#  [string,weight,parameter_index_1,parameter_index_2,parameter_name_1,parameter_name_2]
if Zinput == 98:
    well2_1 = ['cT',2,      1,3,'c','c_S']
    well2_2 = ['ex',10,     0,1,'e_0','x']
    well2_3 = ['xT',0.5,    1,3,'x','c_S']
    zoom = 40
elif Zinput == 94:
    well2_1 = ['xT',3,    1,3,'x','c']
    well2_2 = ['xd',0.5,    1,4,'e_0','c_S']
    well2_3 = ['eT',4,    0,3,'e_0','c']
    zoom = 15

def well_2():
    well_plot_2(str(Zinput),str(Ainput),
            cwd+'/../TRANSFER_FOLDER/grid_values/'+str(Zinput)+str(Ainput)+'-1/'+well2_1[0],
            cwd+'/../TRANSFER_FOLDER/grid_values/'+str(Zinput)+str(Ainput)+'-1/'+well2_2[0],
            cwd+'/../TRANSFER_FOLDER/grid_values/'+str(Zinput)+str(Ainput)+'-1/'+well2_3[0],
            well2_1,well2_2,well2_3,zoom)

def grid_runn(parameter,parameter2):
    grid_run(str(Zinput),str(Ainput),generate_number,None,resolution,
            reaction_type = reaction_type,
            parameter = parameter,
            parameter2 = parameter2)


os.chdir(cwd+"/..")
os.system('./restore.sh')

if len(sys.argv) > 1:
    job_type = str(sys.argv[1])
else: 
    job_type = None
#  if len(sys.argv) > 2:
    #  n_energy = float(sys.argv[2])

if job_type == 'opt' and reaction_type == 'spontaneous':
    opti(optimization_type,-1)
elif job_type == 'opt' and reaction_type == 'induced':
    for n_energy in energies:
        opti(optimization_type,n_energy)
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
elif job_type == 'grid':
    arg1 = sys.argv[2]
    arg2 = sys.argv[3]
    grid_runn(str(arg1),str(arg2))
else:
    print("ERROR: Job type not recognized.")

os.chdir(cwd+"/..")
os.system('./restore.sh')
