#!/usr/bin/env python
"""
   dnolivieri:  updated ...4 dec 2015
     * used to generate the "positive" multiresolution signals.
     - based upon loci.
     http://sebastianraschka.com/Articles/2014_pca_step_by_step.html
     http://spartanideas.msu.edu/2014/04/13/implementing-a-principal-component-analysis-pca-in-python-step-by-step/
     http://sebastianraschka.com/Articles/2015_pca_in_3_steps.html
"""
import numpy as np
import sys
import itertools
from scipy import *
import struct
import re
import cPickle as pickle
import errno
import numpy as np
from sklearn.decomposition import PCA, KernelPCA

class DRmethods:
    def __init__(self, ncompnts):
        self.ncompnts = ncompnts 

    def apply_DR(self, X):
        pca = PCA(n_components = 4)
        pca.fit(X)
        Xbar = pca.transform(X)
        return pca, Xbar

## ---------------MAIN ----------------------------------
if __name__ == '__main__':

    D = DRmethods(3)

    X = np.array([[1.4160856855378861, 1.6696845258936002, 2.1381241643924174, 1.7354957705302798, 1.2572209774400114, 1.5864064108609213, 1.6741209303198727, 2.0222727516874666, 1.8764356041523746, 1.5145692005036298], 
         [1.3756005786988796, 1.7606188686001814, 2.2639998117836337, 1.5728648321219127, 0.83941751654339647, 1.3751015595632454, 1.5932825538356663, 1.4654539595008327, 1.8751080426347424, 1.5941722801239668], 
         [2.5981912683150488, 2.047003796186889, 2.2946978554977817, 2.2314350260952831, 2.0901184649431976, 2.0866771598152845, 1.8071503147470744, 2.1431257049755614, 1.8510874764250882, 2.2655455423525166], 
         [1.5480996968841654, 1.9173650757138179, 1.8578599177369277, 2.0894381520197833, 2.5536802584456124, 2.0236123682259843, 1.9094481782028025, 1.9826918109015768, 1.9310235283032458, 1.6315245218120906], 
         [1.232414281716806, 1.8672199170124479, 2.3022915715604313, 1.5677405827634416, 1.2554098582621755, 1.4622990483891356, 1.5930478996223594, 1.8530786170841327, 1.8631911036991242, 1.7093079146226813] ])

    
    print X
    """
    kpca=KernelPCA(n_components=5, kernel="rbf", fit_inverse_transform=True, gamma=10)
    kpca.fit(X)

    X_kpca = kpca.transform(X)

    
    print X_kpca

    Y = np.random.rand(10)
    print Y
    #kpca.get_covariance()
    print kpca.get_params(deep=True)
    """



    pca = PCA(n_components = 4)
    pca.fit(X)
    #print pca.get_covariance()
    Y = np.random.rand(10)
    print Y

    print "----------------------"

    Xnew= pca.transform(X)
    print "Xnew=", Xnew

    ## here is the projection.
    """
       This would be needed in the Prediction phase...
       so I would need to save the covariance matrix 
       (which is part of the 

    """
    Ynew =  pca.transform(Y)
    print "Ynew=", Ynew
