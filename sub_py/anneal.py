import time
from math import sqrt, pi
import os, sys
import numpy as np
from scipy import optimize
from random import random
from error import error

cwd = os.getcwd()

def acceptance_probability(old_error , new_error , T):
    return np.exp(-(new_error-old_error)/(T))
    #  return T/(((old_error-new_error)**2 + T**2)**3)

def neighbor(old):
    sigma = old/5E2
    mu = old
    new = np.zeros((5))
    new[0] = mu[0] + sigma[0] * np.random.uniform(low=-1,high=1)
    new[1] = mu[1] + sigma[1] * np.random.uniform(low=-1,high=1)
    new[2] = mu[2] + sigma[2] * np.random.uniform(low=-1,high=1)
    new[3] = mu[3] + sigma[3] * np.random.uniform(low=-1,high=1)
    new[4] = mu[4] + sigma[4] * np.random.uniform(low=-1,high=1)
    return new

    #  define function which takes in a set of parameters and returns the raw chi-squared error (total sum)
    #  see error.py for the details of this error calculation

def anneal(objective_function,guess):
    old_error = objective_function(guess)
    T = 1.0
    T_min = 5E-2
    alpha = 0.7
    best_error = objective_function(guess)
    best_guess = guess
    #  best_guess = 0
    print("Guess:",guess,"Error:",best_error)
    while T > T_min:
        i = 1
        while i <= 50:
            new_guess = neighbor(guess)
            new_error = objective_function(new_guess)
            print("NEW GUESS-------",new_guess)
            print("NEW ERROR-------",new_error)
            ap = acceptance_probability(old_error, new_error, T)
            if ap > random():
                guess = new_guess
                old_error = new_error
                if old_error < best_error:
                    best_error = old_error
                    best_guess = guess
                    print("NEW BEST ERROR-------",best_error)
            i += 1
        T = T*alpha
    return best_guess, best_error
