#!/usr/bin/env python
"""
   dnolivieri:  updated ...29-jan-2016
     - test the occurences...
     - 
s="RKEERKRKKERKRERERKKERKKERKKERKKERKKERKRKRDRERETERKKERQKEREER
Two implmentations:
>>> res=collections.Counter(s)
Counter({'K': 20, 'R': 20, 'E': 17, 'D': 1, 'Q': 1, 'T': 1})
>>> res.most_common(2)
[('K', 20), ('R', 20)]
>>> res.most_common(2)[0]
('K', 20)
>>> res.most_common(2)[0][1]
20

>>> import collections
>>> d=collections.defaultdict(int)
>>> for c in s:
...     d[c]+=1
>>> d
defaultdict(<type 'int'>, {'E': 17, 'D': 1, 'K': 20, 'Q': 1, 'R': 20, 'T': 1})
>>> 
#dbar = [ x[0]+x[1] for x in list(itertools.product('EDKQRT', repeat=2))]        
"""
import collections
import itertools

class Test:
    def __init__(self):
        pass
    
    def count_singles(self, s):
        res=collections.Counter(s)
        print res
        print res.most_common(2), res.most_common(2)[0][1]

    def count_sequential_doubles(self, s):

        dbar = [ s[i-1]+s[i] for i in range(len(s)) if i>0 ]
        res=collections.Counter(dbar)
        #res.most_common(2)[0][1]
        print res


#-------------------
if __name__ == '__main__':
    s="RKEERKRKKERKRERERKKERKKERKKERKKERKKERKRKRDRERETERKKERQKEREER"
    T = Test()
    T.count_singles(s)
    T.count_sequential_doubles(s)
