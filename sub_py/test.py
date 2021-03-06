import time
from math import sqrt, pi
import os, sys
import numpy as np
from scipy import optimize
from scipy import integrate

from ranges import *

cwd = os.getcwd()

#  parameters = param_list(98,252,"spontaneous")

def test_error(Z, A, e ,x ,c ,T,d, generate_number, parsed_data, **kwargs):
    error = 1
    error += abs(e - 10.37)
    error += abs(x - 1.27)
    error += abs(c - 1.18)
    error += abs(T - 0.87)
    error += abs(d - 0.52)
    #  error += np.exp(abs(e -10.37))
    #  error += np.exp(abs(x - 1.27))
    #  error += np.exp(abs(c - 1.18))
    #  error += np.exp(abs(T - 0.87))
    #  error += np.exp(abs(d - 0.52))
    return error**2, error, None, None,None,None,10
