import os, sys, time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import floor

cwd = os.getcwd()

from covariance import covariance
from ranges import *
from isotope import *

def zoom(array,percent,shift):
    percent = percent / 100
    number = floor(len(array)/2-percent*len(array)/2)
    return array[number+shift:-number+shift]
def zoom2d(array,percent,shift1,shift2):
    percent = percent / 100
    number = floor(len(array)/2-percent*len(array)/2)
    return array[number+shift1:-number+shift1,number+shift2:-number+shift2]

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

    range_array = param_ranges[parameter1]
    range_array_2 = param_ranges[parameter2]
    resolution  = len(lines)

    x_array = np.arange(range_array[0] , range_array[1] , (range_array[1] - range_array[0])/resolution)
    y_array = np.arange(range_array_2[0] , range_array_2[1] , (range_array_2[1] - range_array_2[0])/resolution)

    zoompercent = 90
    zoomshift1 = 0
    zoomshift2 = 0
    x_array = zoom(x_array,zoompercent,zoomshift1)
    y_array = zoom(y_array,zoompercent,zoomshift2)
    error_array = zoom2d(error_array,zoompercent,zoomshift1,zoomshift2)

    #  create appropriate directory if it does not already exist 

    var_path = cwd + '/../output/contours/' + 'Z=' + str(Zinput) + ' A=' + str(Ainput)
    if not os.path.exists(var_path):
        os.makedirs(var_path)
    os.chdir(var_path)

    plt.contourf(x_array,y_array,error_array)
    plt.colorbar()
    plt.xlabel('$'+parameter_labels[parameter1]+'$',fontsize=15)
    plt.ylabel('$'+parameter_labels[parameter2]+'$',fontsize=15)
    plt.savefig(parameter1+'_vs_'+parameter2+'.pdf')
    plt.close()
    os.chdir(cwd)

def well_plot(Zinput, Ainput , filename1, filename2, parameter1 , parameter2 , parameter3 , parameter4):
    file = open(filename1,"r")
    lines1 = file.readlines()
    file.close()
    file = open(filename2,"r")
    lines2 = file.readlines()
    file.close()

    reac_t = 'spontaneous'
    iso = isotope(Zinput,Ainput,reac_t)

    error_array1 = np.zeros((len(lines1),len(lines1)))
    error_array2 = np.zeros((len(lines2),len(lines2)))

    for line_number in range(0,len(lines1)):
        line = lines1[line_number].strip()
        elements = line.split(' ')
        element_array = np.array(elements)
        error_array1[line_number] = element_array
    for line_number in range(0,len(lines2)):
        line = lines2[line_number].strip()
        elements = line.split(' ')
        element_array = np.array(elements)
        error_array2[line_number] = element_array

    range_array1 = param_ranges[parameter1]
    range_array1_2 = param_ranges[parameter2]

    range_array2 = param_ranges[parameter1]
    range_array2_2 = param_ranges[parameter2]

    resolution  = len(lines1)

    x_array1 = np.arange(range_array1[0] , range_array1[1] , (range_array1[1] - range_array1[0])/resolution)
    y_array1 = np.arange(range_array1_2[0] , range_array1_2[1] , (range_array1_2[1] - range_array1_2[0])/resolution)

    y_array2 = np.arange(range_array2[0] , range_array2[1] , (range_array2[1] - range_array2[0])/resolution)
    x_array2 = np.arange(range_array2_2[0] , range_array2_2[1] , (range_array2_2[1] - range_array2_2[0])/resolution)

    #  create appropriate directory if it does not already exist 

    var_path = cwd + '/../output/contours/' + 'Z=' + str(Zinput) + ' A=' + str(Ainput)

    if not os.path.exists(var_path):
        os.makedirs(var_path)

    os.chdir(var_path)

    #  objective_array = error_array1 + 5 * np.fliplr(error_array2)
    def concoction(first, second):
        return  1.2 * first + 1 * second

    xs = concoction(x_array1,x_array2) - 1
    ys = concoction(y_array1,y_array2) + 1
    zs = concoction(error_array1,np.transpose(error_array2))
    zs = np.log(zs)

    plt.contourf(xs,ys,zs,cmap=plt.cm.Reds)
    plt.colorbar()
    plt.plot([1.7594],[1.43834], 'bo')
    plt.errorbar(1.7594 , 1.43834 , xerr = 1.8 , yerr = 2.36,color = 'b' , fmt = ' ' , capsize = 3 , elinewidth = 1)
    plt.xlabel('$1.2 \\times c + d$TKE(MeV)',fontsize=15)
    plt.ylabel('$1.2 \\times d$TKE (MeV) $ + c_S$',fontsize=15)
    plt.annotate("$^{252}$Cf(sf)",(2.5,7),fontsize=18)
    plt.savefig('well.pdf')
    plt.close()

    zoompercent = 50
    zoomshift1 = 2
    zoomshift2 = -1
    close_xs = zoom(xs,zoompercent,zoomshift1)
    close_ys = zoom(ys,zoompercent,zoomshift2)
    close_zs = zoom2d(zs,zoompercent,zoomshift1,zoomshift2)

    plt.contourf(close_xs,close_ys,close_zs,cmap=plt.cm.Reds)
    plt.colorbar()
    plt.plot([1.7594],[1.43834], 'bo')
    plt.xlabel('$1.2 \\times c + d$TKE (MeV)',fontsize=15)
    plt.ylabel('$1.2 \\times d$TKE (MeV) $ + c_S$',fontsize=15)
    plt.annotate("$^{252}$Cf(sf)",(2.5,7),fontsize=18)
    plt.savefig('well_closeup.pdf')
    plt.close()

    os.chdir(cwd)

def well_plot_2(Zinput, Ainput , filename1, filename2, filename3, parameter1 , parameter2 , parameter3 , parameter4):
    file = open(filename1,"r")
    lines1 = file.readlines()
    file.close()
    file = open(filename2,"r")
    lines2 = file.readlines()
    file.close()
    file = open(filename3,"r")
    lines3 = file.readlines()
    file.close()

    error_array1 = np.zeros((len(lines1),len(lines1)))
    error_array2 = np.zeros((len(lines2),len(lines2)))
    error_array3 = np.zeros((len(lines3),len(lines3)))

    for line_number in range(0,len(lines1)):
        line = lines1[line_number].strip()
        elements = line.split(' ')
        element_array = np.array(elements)
        error_array1[line_number] = element_array
    for line_number in range(0,len(lines2)):
        line = lines2[line_number].strip()
        elements = line.split(' ')
        element_array = np.array(elements)
        error_array2[line_number] = element_array
    for line_number in range(0,len(lines2)):
        line = lines3[line_number].strip()
        elements = line.split(' ')
        element_array = np.array(elements)
        error_array3[line_number] = element_array

    range_array1 = param_ranges[parameter1]
    range_array1_2 = param_ranges[parameter2]

    range_array2 = param_ranges[parameter1]
    range_array2_2 = param_ranges[parameter2]

    range_array3 = param_ranges[parameter1]
    range_array3_2 = param_ranges[parameter2]

    resolution  = len(lines1)

    x_array1 = np.arange(range_array1[0] , range_array1[1] , (range_array1[1] - range_array1[0])/resolution)
    y_array1 = np.arange(range_array1_2[0] , range_array1_2[1] , (range_array1_2[1] - range_array1_2[0])/resolution)

    y_array2 = np.arange(range_array2[0] , range_array2[1] , (range_array2[1] - range_array2[0])/resolution)
    x_array2 = np.arange(range_array2_2[0] , range_array2_2[1] , (range_array2_2[1] - range_array2_2[0])/resolution)

    y_array3 = np.arange(range_array3[0] , range_array2[1] , (range_array2[1] - range_array2[0])/resolution)
    x_array3 = np.arange(range_array3_2[0] , range_array3_2[1] , (range_array3_2[1] - range_array3_2[0])/resolution)

    #  create appropriate directory if it does not already exist 
    var_path = cwd + '/../output/contours/' + 'Z=' + str(Zinput) + ' A=' + str(Ainput)

    if not os.path.exists(var_path):
        os.makedirs(var_path)

    os.chdir(var_path)

    def concoction(first, second,third):
        return  8 * first + 1 * second + 15 * third

    xs = concoction(x_array1,x_array2,x_array3) - 1
    ys = concoction(y_array1,y_array2 ,y_array3) + 1
    #  zs = concoction(error_array1,error_array2,error_array3)
    zs = concoction(error_array1,np.transpose(error_array2),np.transpose(error_array3))
    #  zs = concoction(np.transpose(error_array1),error_array2,error_array3)

    plt.contourf(xs,ys,zs,cmap=plt.cm.Reds)
    plt.colorbar()
    #  plt.plot([50.61],[16.94], 'bo')
    #  plt.errorbar( 50.61 , 16.94 , xerr = 11.4 , yerr = 21 ,color = 'b' , fmt = ' ' , capsize = 3 , elinewidth = 1)
    plt.xlabel('$8 \\times e_0 $(MeV)$ + d$TKE (MeV)$ + 15 \\times c_s$',fontsize=15)
    plt.ylabel('$d$TKE (MeV) $ + 23 \\times c_S$',fontsize=15)
    plt.savefig('well.pdf')
    plt.close()

    zoompercent = 20
    zoomshift1 = 1
    zoomshift2 = -1
    close_xs = zoom(xs,zoompercent,zoomshift1)
    close_ys = zoom(ys,zoompercent,zoomshift2)
    close_zs = zoom2d(zs,zoompercent,zoomshift1,zoomshift2)

    plt.contourf(close_xs,close_ys,close_zs,cmap=plt.cm.Reds)
    plt.colorbar()
    #  plt.plot([50.61],[16.94], 'bo')
    plt.xlabel('$8 \\times e_0 $(MeV)$ + d$TKE (MeV)$ + 15 \\times c_s$',fontsize=15)
    plt.ylabel('$d$TKE (MeV) $ + 23 \\times c_S$',fontsize=15)
    plt.savefig('well_closeup.pdf')
    plt.close()

    os.chdir(cwd)
