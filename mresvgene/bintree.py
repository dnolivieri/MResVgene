#!/usr/bin/env python
"""
   dnolivieri:  updated ...4 dec 2015
     - the multiresolution tree structure for model fits.
"""

import sys
import os

class TreeNode:
    def __init__(self, data):
        self.left=None
        self.right=None
        self.data = data


    def insert(self, data):
        if data[1] < self.data[1]: 
            if self.left is None:
                self.left = TreeNode(data)
                print "L",
            else:
                self.left.insert(data)
        else:
            if self.right is None:
                self.right = TreeNode(data)
                print "R",
            else:
                self.right.insert(data)

    

    def print_tree(self):
        if self.left:
            self.left.print_tree()
        print self.data,
        if self.right:
            self.right.print_tree()






## ---------------MAIN ----------------------------------
if __name__ == '__main__':

    data = [(0,200), (240,175), (460,185)]
    root = TreeNode(data[0])
    for i in data[1:]:
        root.insert(i)

    print
    root.print_tree()
