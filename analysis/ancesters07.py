#!/usr/bin/env python

##  dnolivieri.net: (updated 30-jul-2013)
#      obtain the ancestor
#      Uses the dentropy.
#   version with names that come from newick


import numpy as np
import dendropy
from dendropy import Tree, TaxonSet

L=[]

def visit_subtree(node):
    if node is not None:
        print "actual node=",node
        for child in node.leaf_iter():
            print node, child.taxon
        print "------"
        visit_subtree(node.parent_node)
    else:
        print node


def coalesce_subtree(node):
    if node is not None:
        for child in node.leaf_iter():
            print node, node.parent_node, child.taxon
        #leaftest= lambda y: True if (y.is_leaf()) else False
        d = [x.taxon for x in node.child_nodes()]
        print "------", d, node

        if (d[0].label[0] == d[1].label[0]):
            k= t.find_node_with_taxon_label(d[0].label)
            print "k=",k
            node.remove_child(k)
            node.parent_node.collapse_clade()
            #d[1].label="a6"
            #print "after change=",node, node.parent_node
            s = [x.taxon for x in node.child_nodes()]
            #print "s=",s
            t.print_plot()

            print "node=",node
            kbar=coalesce_subtree(node.parent_node)

        else:
            print "not equal", d[0].label[0], d[1].label[0], node, type(node)
            L.append(node)
    else:
        print "beyond parent", node
        return node




def make_ancestorTree(tr, lnodes ):
    for l in lnodes:
        print l
        for child in l.leaf_iter():
            print child.taxon

        for x in l.child_nodes():
            l.remove_child(x)

    return tr

## ---------------MAIN ----------------------------------
if __name__ == '__main__':

    ### ------ tests -----------
    #s="((((a5,b6),b7),(a8,b9)),c) ;"
    #s="(((a1,(a2,a3)),b4),c);"
    #s="(((a1,b2),(a3,(a4,b5))),(a6,((b7,b8),b9))) ;"
    #t = Tree.get_from_string(s, "newick")
    ## -------------------------


    ## --- quick and dirty way of getting same names...

    t = dendropy.Tree.get_from_path('pg.nwk', 'newick')
    t.print_plot()

    tcnt=1
    for node in t.level_order_node_iter():
        if (node.taxon != None):
            node.taxon.label=node.taxon.label.split("|")[1].split(" ")[0] + str(tcnt)
            print node.level(),node.taxon.label
            tcnt+=1
            
    t.print_plot()


    cnt=0
    intNodes = []
    multifurcating = lambda x: True if (len(x.child_nodes()) == 2) else False
    for nd in t.postorder_node_iter(multifurcating):
        print "------", cnt, "-------------"
        d = [x.taxon for x in nd.child_nodes()]
        if None not in d:
            print d
            intNodes.append( nd )

        #kbar=coalesce_subtree(nd)
        #print "kbar =", kbar
        #print L
        #sbar= L[0]
        #s = [x.taxon for x in sbar.child_nodes()]
        #visit_subtree(nd)
        cnt+=1



    for nd in intNodes:
        print nd
        kbar=coalesce_subtree(nd)
        print "kbar =", kbar
        print L

    make_ancestorTree(t, L )

    t.print_plot()

