import numpy as np

from ranges import *

def Maxwellian(e,T):
    return np.sqrt(e)*np.exp(-e/T)*(2/np.sqrt(np.pi)/np.sqrt(T)/T)

vMaxwellian = np.vectorize(Maxwellian)
egrid = np.array(mannhart_bins)
degrid = egrid[1:]-egrid[:-1]

def MaxwellianSpectrum(T):
    egrid = mannhart_bins
    #  np.logspace(-3,2,bin_number['mannhart'])
    return vMaxwellian(egrid,T)
