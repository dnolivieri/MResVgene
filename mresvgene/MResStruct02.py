#!/usr/bin/env python
"""
   dnolivieri.net:  (5 dec 2015) 
    - The Multilevel  structure of the sequence
"""
class MResStruct:
    def __init__(self, nlevels):
        self.nlevels = nlevels

    def splitN(self, str, nlevel):
        if nlevel == 0:
            return [str]
        else:
            left = str[:len(str)/2]
            right = str[len(str)/2:]
            return [] + self.splitN(left, nlevel-1) + self.splitN(right, nlevel-1)

    def get_pyramid(self, data):
        Seqs = []
        for nk in range(self.nlevels): 
            Seqs.append( self.splitN( data, nk ) )
        return Seqs


## ---------------MAIN ----------------------------------
if __name__ == '__main__':

    data = "GDTKMRISQDQLSFTRKPNRTVHISCQLSGVPLETTIVHWYQEKEGGS"
    
    M =  MResStruct( 4 ) 
    S  = M.get_pyramid(data)
    

    print "-------------"
    print "S=",S

    print
    for i in range(len(S)):
        print i, "...", S[i]
