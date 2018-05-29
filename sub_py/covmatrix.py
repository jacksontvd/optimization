from __future__ import print_function

import numpy as np
from numpy import exp, array, linalg, zeros, isscalar
from math import gamma
from ranges import *
from error import *
from test import *
from data_parse import *

def blockPrint():
    sys.stdout = open(os.devnull, 'w')
def enablePrint():
    sys.stdout = sys.__stdout__

def probability(chi_sq_array,number,dof):
    factor = 1/(2**(dof/2)*gamma(dof/2))
    #  print(factor)
    exponential = np.exp(-chi_sq_array/(2*number))
    #  print(exponential)
    final = factor * exponential * (chi_sq_array)**(dof/2 - 1)
    #  print(final)
    return final

def log_probability(chi_sq,dof):
    factor = -np.log(gamma(dof/2)*(2**(dof/2)))
    print(factor)
    power = (1/2)*(dof-2)*np.log(chi_sq)
    print(power)
    exponential = -chi_sq/2
    print(exponential)
    return factor + power + exponential

def freya_hessian(Z,A,generate_number,h,reac_t):
    parameters = param_list(Z,A,reac_t)
    print("Using Hessian to calculate correlation matrix at the point: ",parameters)
    print("Interval used to estimate derivatives: ",h)
    parsed_data = data_parse(Z,A,reac_t)
    def objective(parameters):
        error_array = error(Z, A, parameters[0],parameters[1], parameters[2], parameters[3], parameters[4], generate_number, parsed_data, reaction_type = reac_t)
        #  error_array = test_error(Z, A, parameters[0],parameters[1], parameters[2], parameters[3], parameters[4], generate_number,None, reaction_type = reac_t)
        error_value = error_array[0]
        print("ERROR:",error_value)
        dof = error_array[6]
        print("DOF:",dof)
        log_prob = log_probability(error_value/100,dof/10)
        print("LOG_PROB:",log_prob)
        print("PROB:",np.exp(log_prob))
        return log_prob
        #  return prob

    #  blockPrint()
    result = calc_cov(objective,parameters,h)
    #  enablePrint()

    correlation_array = np.zeros((5,5))
    for i in range(0,5):
        for j in range(0,5):
            correlation_array[i,j] = result[i,j]/np.sqrt(abs(result[i,i] * result[j,j]))

    print_path = cwd + '/../../output/optimization/' + 'Z=' + str(Z) + 'A=' + str(A) + "_" + str(reac_t) + '/'
    if not os.path.exists(print_path):
        os.makedirs(print_path)

    document =open(print_path+'/matrix.tex', 'w+')
    document.write(
    '$e_0$ &'+ 
    str(round(correlation_array[0,0],3)) +'&' +
    str(round(correlation_array[0,1],3)) +'&' +
    str(round(correlation_array[0,2],3)) +'&' +
    str(round(correlation_array[0,3],3)) +'&' +
    str(round(correlation_array[0,4],3)) +
    '\\\\ \n\hline\n' +'$x$ &' +
    str(round(correlation_array[1,0],3)) +'&' +
    str(round(correlation_array[1,1],3)) +'&' +
    str(round(correlation_array[1,2],3)) +'&' +
    str(round(correlation_array[1,3],3)) +'&' +
    str(round(correlation_array[1,4],3)) +
    '\\\\ \n\hline\n' +'$c$ &' +
    str(round(correlation_array[2,0],3)) +'&' +
    str(round(correlation_array[2,1],3)) +'&' +
    str(round(correlation_array[2,2],3)) +'&' +
    str(round(correlation_array[2,3],3)) +'&' +
    str(round(correlation_array[2,4],3)) +
    '\\\\ \n\hline\n' + '$c_S$ &'+
    str(round(correlation_array[3,0],3)) +'&' +
    str(round(correlation_array[3,1],3)) +'&' +
    str(round(correlation_array[3,2],3)) +'&' +
    str(round(correlation_array[3,3],3)) + '&' +
    str(round(correlation_array[3,4],3)) +
    '\\\\ \n\hline\n' + '$d$TKE &' +
    str(round(correlation_array[4,0],3)) +'&' +
    str(round(correlation_array[4,1],3)) + '&' +
    str(round(correlation_array[4,2],3)) +'&' +
    str(round(correlation_array[4,3],3)) +'&' +
    str(round(correlation_array[4,4],3)))

    determinant = np.linalg.det(correlation_array)

    return result,correlation_array, determinant

def calc_cov(func, x, h):
    """
    Calculate the covariance matrix of the input log probability function
    using the Hessian as a gaussian approximation.

    cov = -H^{-1}

    parameters
    ---------
    func: function
        A function representing the log(probability(x))
    x: array
        the position about which to calculate the deriviatives
    h: scalar or array
        Step size, or sizes, to use for finite differencing

        If a scalar is input, it is used for all dimensions

    returns
    -------
    The covariace matrix derived from the Hessian. See calc_hess

    restrictions
    ------------
    The Hessian must be invertible
    """
    hess = calc_hess(func, x, h)
    cov = -linalg.inv(hess)
    return cov

def calc_hess(func, x, h):
    """
    calculate the Hessian matrix

       d^2 f
    ----------
     dx_i dx_j

    parameters
    ---------
    func: function
        A function than takes an array x as the argument
    x: array
        the position about which to calculate the deriviatives
    h: scalar or array
        Step size, or sizes, to use for finite differencing

        If a scalar is input, it is used for all dimensions

    returns
    -------
    The Hessian matrix calculated via finite differences
    """

    x, h = _get_x_and_h(x, h)
    ndim=x.size

    hess=zeros( (ndim,ndim) )

    for i in range(ndim):
        for j in range(i,ndim):
            if i==j:
                hess[i,i] = calc_partial_ii(func, x, i, h[i])
            else:
                hess[i,j] = calc_partial_ij(func, x, i, j, h[i], h[j])
                hess[j,i] = hess[i,j]

    return hess

def calc_partial_ij( func, xin, i, j, h, k):
    """
    calculate the second partial with respect to a single variable.
    Four function evaluations

    f(x+h,y+k) - f(x+h,y-k) - f(x-h,y+k) + f(x-h,y-k)
    -------------------------------------------------
                         4*h*k
    """
    x=xin.copy()

    # f(x+h, y+h)
    x[i] = xin[i] + h
    x[j] = xin[j] + k
    f1=func(x)

    # f(x+h, y-h)
    #x[i] = xin[i] + h
    x[j] = xin[j] - k
    f2=func(x)

    # f(x-h, y+h)
    x[i] = xin[i] - h
    x[j] = xin[j] + k
    f3=func(x)

    # f(x-h, y-h)
    #x[i] = xin[i] - h
    x[j] = xin[j] - k
    f4=func(x)

    num = f1 - f2 - f3 + f4
    denom=4.*h*k

    deriv = num/denom
    return deriv

def calc_partial_ii( func, xin, i, h ):
    """
    calculate the second partial with respect to a single variable.

    Three function evaluations instead of 4 used in calc_partial_ij

    f(x+h) - 2*f(x) + f(x-h)
    ------------------------
              h**2
    """
    x=xin.copy()

    # f(xi+h)
    x[i] = xin[i] + h
    f1=func(x)

    # f(x)
    f2=func(xin)

    # f(xi+h)
    x[i] = xin[i] - h
    f3=func(x)

    num = f1 - 2*f2 + f3
    denom=h**2

    deriv = num/denom
    return deriv

def _get_x_and_h(x, h):
    """
    get x and h as arrays.  If h is a scalar, repeat
    it to the size of x
    """

    x=array(x, ndmin=1, copy=False)
    ndim=x.size

    if isscalar(h):
        h = zeros(ndim) + h
    else:
        h=array(h, ndmin=1, copy=False)
        if h.size != x.size:
            raise ValueError("h and x have different "
                             "dimensions: %d %d" % (h.size,x.size))

    return x, h


def test():

    covtrue=array( [ [ 400.0, 0.2, 0.1  ],
                     [ 0.2,   2.0, 0.2  ],
                     [ 0.1,   0.2, 1.0] ] )

    cinvtrue=linalg.inv(covtrue)
    xmean=array( [1.0, 2.0, 3.0] )

    def calc_lnprob(x):
        xdiff = x-xmean

        chi2 = cinvtrue.dot( xdiff ).dot( xdiff )

        return exp( -0.5*chi2 )

    h=1.0e-3

    covmeas = calc_cov(calc_lnprob, xmean, h)

    print("true cov:")
    print(covtrue)
    print("meas cov:")
    print(covmeas)

    print("frac diff:")
    print( (covmeas-covtrue)/covtrue )

if __name__=="__main__":
    test()
