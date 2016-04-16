#!/usr/bin/env python
"""
   for plotting the Iteration curves for the paper.
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy.random import normal
from scipy.optimize import curve_fit
from scipy.stats.distributions import  t


def func(x, a, b):
    'nonlinear function in a and b to fit to data'
    #return a * x / (b + x
    return a * np.exp(-b*x)


def fit_model(x, y): 
    initial_guess = [1.2, 0.03]
    pars, pcov = curve_fit(func, x, y, p0=initial_guess)
    alpha = 0.05 # 95% confidence interval = 100*(1-alpha)

    n = len(y)    # number of data points
    p = len(pars) # number of parameters
    dof = max(0, n - p) # number of degrees of freedom
    # student-t value for the dof and confidence level
    tval = t.ppf(1.0-alpha/2., dof) 
    for i, p,var in zip(range(n), pars, np.diag(pcov)):
        sigma = var**0.5
        print 'p{0}: {1} [{2}  {3}]'.format(i, p,p - sigma*tval,p + sigma*tval)

    return pars



# -----------------------------------------------
if __name__ == '__main__':


    # total data1
    #y = np.array([ 200., 1564., 2466., 2744., 2974., 3070.] )
    #x = np.arange(y.size) 

    # <2.0
    #y = np.array([218., 176., 160., 121., 120., ] )
    #y = np.array([218., 176., 160., 121., 120., ] )  <2.0
    y = np.array([101., 73., 60., 48., 46., 45.])   ###<1.75


    x = np.arange(y.size) 

    pars = fit_model(x, y)
    fig, ax = plt.subplots(1)

    ax.plot(x,y,'bo ')
    xfit = np.linspace(0,5)
    yfit = func(xfit, pars[0], pars[1])
    ax.plot(xfit,yfit,'b-')
    
    # total data2
    #ybar = np.array( [200., 1310., 1546., 2053., 2132., 2290.] )
    #xbar = np.arange(y.size) 


    
    # <2.0
    #ybar = np.array([294., 246., 222., 176., 165. ] )  #<20
    ybar = np.array([ 103., 93., 88., 68., 62., 58.])
    xbar = np.arange(y.size) 
    


    pars = fit_model(xbar, ybar)
    ax.plot(xbar,ybar,'ro ')
    xfitL = np.linspace(0,5)
    yfitL = func(xfit, pars[0], pars[1])
    ax.plot(xfitL,yfitL,'r-')

    plt.legend(['data','fit', 'data2', 'fit2'],loc='best')


    ax.fill_between(xfitL, yfitL, yfit, facecolor='gray', alpha=0.15)

    plt.show()
