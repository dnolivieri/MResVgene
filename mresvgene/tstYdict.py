#!/usr/bin/env python
"""
   dnolivieri:  updated ...8 december 2015

      This does the prediction of Vs based 
      upon a multi-resolution training.

      - also should implement (eventually) the 
       bootstrap learning.
"""

import collections
import numpy as np
import itertools


def sliceVec( Y ):
    ybar = [ np.array(list(itertools.chain(*x[0]))) for x in Y]
    print ybar


def ListsliceVec( Y ):
    ybar=[]
    for k in range(np.power(2,3)-1):
        ybar.append(np.array([list(itertools.chain(*x[k])) for x in Y ]))


                        

    print ybar
    print 

    print ybar[1]
                            

    """
    for  i in range(3):
        print i, "....", ybar[i][0:2]
    """


#----------------------------------
if __name__ == '__main__':


    Y= [ {0: [[1.3704568655693683, 2.2154018098257131, 1.996857680129728]], 
          1: [[0.61850356498582593, 2.2575327275890844, 1.9386739528038182]], 
          2: [[2.1367939699897907, 2.0457396506803054, 1.8630343120927573]], 
          3: [[0.7005411906193626, 1.2698953365950327, 1.0038919201377712]], 
          4: [[0.54691750422357766, 2.9370043676106716, 2.75464867160112]], 
          5: [[2.9335967700369388, 2.6308882381217793, 2.3403540325813368]], 
          6: [[1.3715889242046788, 1.387549507943389, 1.3336080908015233]]},
         {0: [[1.7892697267321445, 2.0260385486123997, 1.616120881474951]], 
          1: [[2.4196413928825233, 2.2114296491889842, 1.7502767954281666]], 
          2: [[1.179929177523886, 1.7377841085031265, 1.3897847754561583]], 
          3: [[2.482604587234774, 1.9862159268860193, 2.0287189308283109]], 
          4: [[2.4572240747827041, 2.4002488500537664, 1.1042827182346882]], 
          5: [[0.7065544197233915, 1.2753565255578008, 0.96683446792500438]], 
          6: [[1.7069370797572843, 2.2741347956558218, 1.8573707926836689]]}]
    
    #sliceVec( Y )
    #ListsliceVec( Y )




    ybar=[[np.array([[ 0.69,  0.31], [ 0.75,  0.25], [0.1, 0.9]]), 
           np.array([[ 0.75,  0.25], [ 0.86,  0.14], [0.1, 0.9]]), 
           np.array([[ 0.5 ,  0.5 ], [ 0.89,  0.11], [0.1, 0.9]]),
           np.array([[ 0.94,  0.06], [ 0.95,  0.05], [0.1, 0.9]]), 
           np.array([[ 0.81,  0.19], [ 0.86,  0.14], [0.1, 0.9]]), 
           np.array([[ 0.58,  0.42], [ 0.76,  0.24], [0.1, 0.9]]), 
           np.array([[ 0.66,  0.34], [ 0.72,  0.28], [0.1, 0.9]])]]
    
    """
    p = ybar[0]
    print p
    print len(p)

    for j in range(p[0].shape[0]):
        print np.array([ list(p[k][j]) for k in range(len(p)) ])


    print 
    print
    s=[np.array([ 0.69,  0.31]), np.array([ 0.75,  0.25]), np.array([ 0.5,  0.5]), np.array([ 0.94,  0.06]), np.array([ 0.81,  0.19]), np.array([ 0.58,  0.42]), np.array([ 0.66,  0.34])]


    Xbar = np.array([ list(x) for x in s])
    print Xbar
    """



    D={'ighv': [np.array([ 0.31,  0.25,  0.5 ,  0.06,  0.19,  0.42,  0.34]),
                np.array([ 0.25,  0.14,  0.11,  0.05,  0.14,  0.24,  0.28])], 
       'igkv': [np.array([ 0.41,  0.15,  0.65 ,  0.36,  0.69,  0.42,  0.34]),
                np.array([ 0.75,  0.34,  0.81,  0.05,  0.24,  0.64,  0.78])]}



    ### This puts all on the same plane...
    print D['ighv'][1]
    print len(D['ighv'])
    for i in range( 2 ): 
        rho = [ D[x][i] for x in D.iterkeys() ]
        tau = np.array([ D[x][i].sum() for x in D.iterkeys() ])
        kmax = np.argmax(tau)
        ymax = np.max(tau)
        print rho, 
        print tau, ".....", kmax, ymax
        
    
              



