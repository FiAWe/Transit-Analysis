# -*- coding: utf-8 -*-
"""
Created on Tue May 24 23:45:05 2016

Model fitting functions

@author: findlaywebb
"""

import numpy as np
import pylab as plt
import pandas as pd
from scipy.optimize import minimize
from timeit import default_timer as timer


def load(route,infile,printColumns = False):
    """Returns the data from the specified route+infile.\
    If printColumns is true the column names will be printed"""
    
    data = pd.io.parsers.read_table(route+infile)
    
    if printColumns:
        print data.columns.values
        
    return data

    
def ChiSquared(y,e,m):
    '''
    Function to calculate the chi-squared value of a
    set of data and a model.
    
    inputs:
        y (array) - measured data points
        e (array) - errors on measured data points
        m (array) - model data points
    outputs:
        X (float) - chi squared value
    '''
    diff   = y - m          # difference between data and model values
    weight = diff / e       # weighted difference, using errors
    X = sum(weight**2.0)    # sum of squares of each value
    
    return X
    
    
def ModelChiSquared(vals,model,data):
    '''
    Calculate the chi squared value for a given set of model parameters.
    
    inputs:
        vals (array)     - array containing the parameters needed to 
                            calculate model values for each data point.
        model (function) - the model to calculate the model values.
        data (array)     - data to which the model will be fitted
    outputs:
        X (float)        - chi squared value 
    '''

    # calculate model values.
    mod = model(data[0],*vals) # The '*' means fill the rest of the 
                               # function's arguments with the 
                               # values in array 'vals'
            
    # calculate chi squared for these model values and the data
    X = ChiSquared(data[1],data[2],mod)  #
    
    return X
    
#
#m = minimize(ModelChiSquared,[1,0], method = 'Nelder-Mead',args=(Model,[x,y,errs],))
#
#print 'Success = ',m['success']
#print 'Chi-squared = ',m['fun']
#print 'Best-fitting g = ',m['x'][0]
#print 'Best-fitting c = ',m['x'][1]
#
#DoF = (len(x) - 2)
#print 'Reduced chi-squared = ',m['fun']/DoF
#
## create an array of model values from our best fit params
#model = Model(x,m['x'][0],m['x'][1])  
#
## plot the data
#plt.errorbar(x,y,yerr=errs, marker='o',linestyle='none') 
## plot the model
#plt.plot(x,model,color='red')                                     
#
#plt.show()

def Bootstrap(data,params,model,n):
    '''
    Bootrap function - resample the data then refit multiple times.
    
    Inputs:
        data (array)     - data to bootstrap - list containing 
                           three arrays of [x,y,yerr]
        vals (array)     - array of initial values for each param
        model (function) - The model function to be bootstrapped
        n (int)          - The number of random samples to fit
    Outputs:
        Params (array)   - list of arrays containing the fitted values 
                           of the parameters for each random sample
    '''
    
    data = np.array(data) # make sure the data is a numpy array
    
    # create array to contain params from random samples
    Params = [[] for i in range(len(params))]
    
    # create list of n arrays of random indices between 0 and 
    # the length of the data set (-1)
    indices = np.random.randint(len(data[0]),size=(n,len(data[0])))

    # create n sets of randomly samples
    X = data[0][indices] 
    Y = data[1][indices] 
    Err = data[2][indices] 

    start = timer()
    # refit each random sample with your model
    for i in range(n):
        
        ##Show percentage complete
        if ((float(i)/n)*100)%5.0 == 0.0:
#            print (float(i)/n)*100
            print int((float(i)/n)*100),'%'
            time = timer()
            print 't = ',int(time-start), 's, i = ', i ,'/', n
        
        # chi squared fit
        mi1  = minimize(ModelChiSquared,params, method = 'Nelder-Mead',\
                            args=(model,[X[i], Y[i], Err[i]]))
        mi2  = minimize(ModelChiSquared,mi1['x'], method = 'Nelder-Mead',\
                            args=(model,[X[i], Y[i], Err[i]]))
        
        par = mi2['x']   
        
        # append best-fit params to 'Params' array
        for p in range(len(Params)):
            Params[p].append(par[p])

    return np.array(Params)

# declare a Gaussian function, which can be fit to the param distribution
def SD(x):
    '''
    Find the Standard deviation of a set of data.
    
    inputs:
        x (float) - input values 
      
    outputs:
        std (float) - standard deviation
    '''
    mean = np.mean(x)            # calculate mean
    diff = x - mean              # subtract from data
    dev = (diff **2.0) / len(x)  # square, divide by N
    std = np.sqrt(sum(dev))               # sum, square root
    
    return std
