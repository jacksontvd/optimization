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

    zoompercent = 80
    zoomshift1 = 0
    zoomshift2 = 0
    x_array = zoom(x_array,zoompercent,zoomshift1)
    y_array = zoom(y_array,zoompercent,zoomshift2)
    error_array = zoom2d(error_array,zoompercent,zoomshift1,zoomshift2)
    #  error_array = np.log(error_array)

    #  create appropriate directory if it does not already exist 

    var_path = cwd + '/../../output/contours/' + 'Z=' + str(Zinput) + ' A=' + str(Ainput)
    if not os.path.exists(var_path):
        os.makedirs(var_path)
    os.chdir(var_path)

    num_x_param = param_ranges[parameter1][2]
    num_y_param = param_ranges[parameter2][2]
    parameters = param_list(Zinput,Ainput,'spontaneous')
    sol_x = parameters[num_x_param]
    sol_y = parameters[num_y_param]

    plt.contourf(x_array,y_array,error_array,cmap=plt.cm.Reds)
    plt.plot(sol_x,sol_y, 'bo')
    plt.colorbar()
    plt.xlabel('$'+parameter_labels[parameter1]+'$',fontsize=15)
    plt.ylabel('$'+parameter_labels[parameter2]+'$',fontsize=15)
    plt.savefig(parameter1+'_vs_'+parameter2+'.pdf')
    plt.close()
    os.chdir(cwd)

def well_plot(Zinput, Ainput , filename1, filename2, p1_list , p2_list,zoompercent):
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

    range_array1 = param_ranges[p1_list[0][0]]
    range_array1_2 = param_ranges[p1_list[0][1]]

    range_array2 = param_ranges[p2_list[0][0]]
    range_array2_2 = param_ranges[p2_list[0][1]]

    resolution  = len(lines1)

    x_array1 = np.arange(range_array1[0] , range_array1[1] , (range_array1[1] - range_array1[0])/resolution)
    y_array1 = np.arange(range_array1_2[0] , range_array1_2[1] , (range_array1_2[1] - range_array1_2[0])/resolution)

    x_array2 = np.arange(range_array2[0] , range_array2[1] , (range_array2[1] - range_array2[0])/resolution)
    y_array2 = np.arange(range_array2_2[0] , range_array2_2[1] , (range_array2_2[1] - range_array2_2[0])/resolution)

    #  create appropriate directory if it does not already exist 

    var_path = cwd + '/../../output/contours/' + 'Z=' + str(Zinput) + ' A=' + str(Ainput)

    if not os.path.exists(var_path):
        os.makedirs(var_path)

    os.chdir(var_path)

    def concoction(first, second):
        return  p1_list[1] * first + p2_list[1] * second
    
    xs = concoction(x_array1,x_array2)
    ys = concoction(y_array1,y_array2)
    zs = concoction(error_array1,error_array2)
    zs = np.log(zs)

    parameters = param_list(Zinput,Ainput,reac_t)
    num_first_x_param = p1_list[2]
    num_second_x_param = p2_list[2]
    num_first_y_param = p1_list[3]
    num_second_y_param = p2_list[3]
    sol_x = concoction(parameters[num_first_x_param],parameters[num_second_x_param])
    sol_y = concoction(parameters[num_first_y_param],parameters[num_second_y_param])

    plt.contourf(xs,ys,zs,cmap=plt.cm.Reds)
    plt.colorbar()
    plt.plot(sol_x,sol_y, 'bo')
    plt.errorbar(sol_x,sol_y ,
            xerr = 0.01 , yerr = 0.01,color = 'b' , fmt = ' ' , capsize = 3 , elinewidth = 1)
    plt.xlabel('$'+str(p1_list[1])+'\\times '+str(p1_list[4])+' + '
            +str(p2_list[1])+'\\times '+str(p2_list[4])+'$',fontsize=15)
    plt.ylabel('$'+str(p1_list[1])+'\\times '+str(p1_list[5])+' + '
            +str(p2_list[1])+'\\times '+str(p2_list[5])+'$',fontsize=15)
    plt.annotate("$^{252}$Cf(sf)",(2.5,7),fontsize=18)
    plt.savefig('well.pdf')
    plt.close()

    zoomshift1 = 2
    zoomshift2 = -1
    close_xs = zoom(xs,zoompercent,zoomshift1)
    close_ys = zoom(ys,zoompercent,zoomshift2)
    close_zs = zoom2d(zs,zoompercent,zoomshift1,zoomshift2)

    plt.contourf(close_xs,close_ys,close_zs,cmap=plt.cm.Reds)
    plt.colorbar()
    plt.plot(sol_x,sol_y, 'bo')
    plt.xlabel('$'+str(p1_list[1])+'\\times '+str(p1_list[4])+' + '
            +str(p2_list[1])+'\\times '+str(p2_list[4])+'$',fontsize=15)
    plt.ylabel('$'+str(p1_list[1])+'\\times '+str(p1_list[5])+' + '
            +str(p2_list[1])+'\\times '+str(p2_list[5])+'$',fontsize=15)
    plt.annotate("$^{252}$Cf(sf)",(2.5,7),fontsize=18)
    plt.savefig('well_closeup.pdf')
    plt.close()

    os.chdir(cwd)

def well_plot_2(Zinput, Ainput , filename1, filename2, filename3, p1_list,p2_list,p3_list,zoompercent):
    file = open(filename1,"r")
    lines1 = file.readlines()
    file.close()
    file = open(filename2,"r")
    lines2 = file.readlines()
    file.close()
    file = open(filename3,"r")
    lines3 = file.readlines()
    file.close()

    reac_t = 'spontaneous'
    iso = isotope(Zinput,Ainput,reac_t)

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

    range_array1 = param_ranges[p1_list[0][0]]
    range_array1_2 = param_ranges[p1_list[0][1]]

    range_array2 = param_ranges[p2_list[0][0]]
    range_array2_2 = param_ranges[p2_list[0][1]]

    range_array3 = param_ranges[p3_list[0][0]]
    range_array3_2 = param_ranges[p3_list[0][1]]

    resolution  = len(lines1)

    x_array1 = np.arange(range_array1[0] , range_array1[1] , (range_array1[1] - range_array1[0])/resolution)
    y_array1 = np.arange(range_array1_2[0] , range_array1_2[1] , (range_array1_2[1] - range_array1_2[0])/resolution)

    x_array2 = np.arange(range_array2[0] , range_array2[1] , (range_array2[1] - range_array2[0])/resolution)
    y_array2 = np.arange(range_array2_2[0] , range_array2_2[1] , (range_array2_2[1] - range_array2_2[0])/resolution)

    x_array3 = np.arange(range_array3[0] , range_array3[1] , (range_array3[1] - range_array3[0])/resolution)
    y_array3 = np.arange(range_array3_2[0] , range_array3_2[1] , (range_array3_2[1] - range_array3_2[0])/resolution)

    #  create appropriate directory if it does not already exist 
    var_path = cwd + '/../../output/contours/' + 'Z=' + str(Zinput) + ' A=' + str(Ainput)

    if not os.path.exists(var_path):
        os.makedirs(var_path)

    os.chdir(var_path)

    def concoction(first,second,third):
        return  p1_list[1] *first + p2_list[1] * second + p3_list[1] * third

    xs = concoction(x_array1,x_array2,x_array3) 
    ys = concoction(y_array1,y_array2 ,y_array3) 
    zs = concoction(error_array1,error_array2,error_array3)
    #  zs = np.log(zs)

    parameters = param_list(Zinput,Ainput,reac_t)
    num_first_x_param = p1_list[2]
    num_second_x_param = p2_list[2]
    num_third_x_param = p3_list[2]
    num_first_y_param = p1_list[3]
    num_second_y_param = p2_list[3]
    num_third_y_param = p3_list[3]
    sol_x = concoction(parameters[num_first_x_param],
            parameters[num_second_x_param],parameters[num_third_x_param])
    sol_y = concoction(parameters[num_first_y_param],
            parameters[num_second_y_param],parameters[num_third_y_param])
    print("Solution:",sol_x,sol_y)

    plt.contourf(xs,ys,zs,cmap=plt.cm.Reds)
    plt.colorbar()
    plt.plot(sol_x,sol_y, 'bo')
    plt.errorbar( sol_x , sol_y ,
            xerr = 10.4725 , yerr = 1.6025 ,color = 'b' , fmt = ' ' , capsize = 3 , elinewidth = 1)
            #  xerr = 4.065 , yerr = 1.8315 ,color = 'b' , fmt = ' ' , capsize = 3 , elinewidth = 1)
            #  xerr = 0.01 , yerr = 0.01 ,color = 'b' , fmt = ' ' , capsize = 3 , elinewidth = 1)
    plt.xlabel('$'+str(p1_list[1])+'\\times '+str(p1_list[4])+
            ' + '+str(p2_list[1])+'\\times '+str(p2_list[4])+
            ' + '+str(p3_list[1])+'\\times '+str(p3_list[4])+
            '$',fontsize=15)
    plt.ylabel('$'+str(p1_list[1])+'\\times '+str(p1_list[5])+
            ' + '+str(p2_list[1])+'\\times '+str(p2_list[5])+
            ' + '+str(p3_list[1])+'\\times '+str(p3_list[5])+
            '$',fontsize=15)
    plt.savefig('well.pdf')
    plt.close()

    zoomshift1 = 5
    zoomshift2 = -2 
    close_xs = zoom(xs,zoompercent,zoomshift1)
    close_ys = zoom(ys,zoompercent,zoomshift2)
    close_zs = zoom2d(zs,zoompercent,zoomshift1,zoomshift2)

    #  plt.contourf(close_xs,close_ys,close_zs,cmap=plt.cm.Reds)
    for number in range(0,len(zs)):
        for numberr in range(0,len(zs)):
            if zs[number,numberr] > np.max(close_zs):
                zs[number,numberr] = np.max(close_zs)
    print(np.max(zs))
    plt.contourf(xs,ys,zs,cmap=plt.cm.Reds)
    plt.colorbar()
    print(np.min(zs),np.max(close_zs))
    plt.plot(sol_x,sol_y,'bo')
    plt.xlabel('$'+str(p1_list[1])+'\\times '+str(p1_list[4])+
            ' + '+str(p2_list[1])+'\\times '+str(p2_list[4])+
            ' + '+str(p3_list[1])+'\\times '+str(p3_list[4])+
            '$',fontsize=15)
    plt.ylabel('$'+str(p1_list[1])+'\\times '+str(p1_list[5])+
            ' + '+str(p2_list[1])+'\\times '+str(p2_list[5])+
            ' + '+str(p3_list[1])+'\\times '+str(p3_list[5])+
            '$',fontsize=15)
    plt.xlim(0.85 * sol_x ,1.1 * sol_x)
    plt.ylim(0.85 * sol_y ,1.1 * sol_y)
    plt.savefig('well_closeup.pdf')
    plt.close()

    os.chdir(cwd)
