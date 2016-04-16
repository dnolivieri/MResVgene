import matplotlib.pyplot as plt
import numpy as np
from numpy.random import normal
from scipy.optimize import curve_fit
from scipy.stats.distributions import  t

y = np.array([ 200., 1564., 2466., 2744., 2974., 3070.] )
ybar = np.array( 200., 1310., 1546., 2053., 2132., 2290.] )
x = np.arange(y.size) 


def func(x, a, b):
    'nonlinear function in a and b to fit to data'
    return a * x / (b + x)

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
    print 'p{0}: {1} [{2}  {3}]'.format(i, p,
                                  p - sigma*tval,
                                  p + sigma*tval)

#import matplotlib.pyplot as plt
plt.plot(x,y,'bo ')
xfit = np.linspace(0,5)
yfit = func(xfit, pars[0], pars[1])
plt.plot(xfit,yfit,'b-')

plt.show()
