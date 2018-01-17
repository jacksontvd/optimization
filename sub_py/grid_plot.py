import os, sys, time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

cwd = os.getcwd()

from covariance import covariance
from ranges import param_ranges

def covariance_plot(Zinput, Ainput, generate_number , some_type, resolution, reaction_type, parameter1, parameter2):
    covar, grid, x , y = covariance(Zinput,
            Ainput,
            generate_number,
            "grid",
            resolution,
            reaction_type = reaction_type,
            parameter = parameter1,
            parameter2 = parameter2)
    
    #  create appropriate directory if it does not already exist 

    var_path = cwd + '/../output/contours/' + 'Z=' + str(Zinput) + ' A=' + str(Ainput)

    if not os.path.exists(var_path):
        os.makedirs(var_path)

    os.chdir(var_path)
    plt.figure()
    CS = plt.contour(x,y,grid)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.title('Chi-squared error:'+parameter1+" vs. "+parameter2)

    #  labels = ['line1', 'line2','line3','line4',
               #  'line5', 'line6']
    #  for i in range(len(labels)):
        #  CS.collections[i].set_label(labels[i])
    #  plt.legend(loc='upper left')

    plt.savefig( parameter1 + parameter2 + '.pdf' )
    plt.close()
    os.chdir(cwd)
    
    return grid

def prerun_plot(Zinput, Ainput , filename, parameter1 , parameter2):
    file = open(filename,"r")
    lines = file.readlines()
    file.close()

    error_array = np.zeros((len(lines),len(lines)))

    for line_number in range(0,len(lines)):
        line = lines[line_number].strip()
        elements = line.split(' ')
        element_array = np.array(elements)
        error_array[line_number] = element_array

    #  new_error = error_array
    #  prob_array = np.exp(new_error)
    #  avg = np.sum(prob_array)
    #  prob_array = np.multiply(prob_array, 1/avg)
    #
    #  print(error_array)
    #  print(prob_array)

    range_array = param_ranges[parameter1]
    range_array_2 = param_ranges[parameter2]
    resolution  = len(lines)

    x_array = np.arange(range_array[0] , range_array[1] , (range_array[1] - range_array[0])/resolution)
    y_array = np.arange(range_array_2[0] , range_array_2[1] , (range_array_2[1] - range_array_2[0])/resolution)

    #  create appropriate directory if it does not already exist 

    var_path = cwd + '/../output/contours/' + 'Z=' + str(Zinput) + ' A=' + str(Ainput)

    if not os.path.exists(var_path):
        os.makedirs(var_path)

    os.chdir(var_path)

    CS = plt.contourf(x_array,y_array,error_array)
    #  CS = plt.contour(x_array,y_array,prob_array)
    plt.colorbar()
    #  plt.clabel(CS, inline=1, fontsize=10)

    #  fig = plt.figure()
    #  ax = fig.add_subplot(111, projection='3d')
    #  ax.contour(x_array, y_array,error_array)
    #  ax.contour(x_array, y_array,prob_array)


    plt.xlabel(parameter1)
    plt.ylabel(parameter2)

    plt.title('Chi-squared error:'+parameter1+" vs. "+parameter2)
    #  plt.title('Probability:'+parameter1+" vs. "+parameter2)

    plt.savefig(parameter1+'_vs_'+parameter2+'.pdf')
    #  plt.savefig('prob_'+parameter1+'_vs_'+parameter2+'.pdf')

    plt.close()

    os.chdir(cwd)
