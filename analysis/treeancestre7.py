#!/usr/bin/env python
import os
import getopt, sys

import dendropy
from dendropy import treecalc
import matplotlib.pyplot as plt


"""
infile = "mam.nwk"
#order = 'mammals-igls.nwk'
order = 'lipotyphla.nwk'
"""
infile = sys.argv[1]
order  = sys.argv[2]

treem = dendropy.Tree.get_from_path(infile, "newick")

Monotreme = ['Ornithorhynchus anatinus']
Marsupial = ['Monodelphis domestica','Sarcophilus harrisii','Macropus eugenii']
Xenarthra = ['Dasypus novemcinctus','Choloepus didactylus']
Afrotheria = ['Loxodonta africana','Orycteropus afer','Procavia capensis','Echinops telfairi','Trichechus manatus','Elephantulus edwardii']
Chiroptere = ['Myotis lucifugus','Eptesicusfuscus','Pteropus vampyrus','Pteropusalecto']
Pterisodactile = ['Equus caballus','Ceratotheriumsimum']
Artiodactile = ['Camelus ferus','Bos taurus','Capra hircus','Vicugna pacos','Ovis aries','Sus scrofa']
Cetacea = ['Orcinus orca','Tursiops truncatus']
Lagomorphes = ['Oryctolagus cuniculus','Ochotona princeps']
Scandentia = ['Tupaia belangeri','Tupaia chinensis']
Soricomorphe = ['Sorex araneus']
Erinaceomorphe = ['Erinaceus europaeus']
Rodent = ['Dipodomys ordii','Octodon degus','Heterocephalus glaber','Chrysochloris asiatica','Ictidomys tridecemlineatus','Mus m2','Cavia porcellus','Microtus ochrogaster']
Carnivore = ['Odobenus rosmarus','Mustela putorius furo','Felis catus abyssinian','Canis lupus familiaris']
Primates = ['a b','Nomascus leucogenys','Otolemur agyisymbanus','Callithrix jacchus','Gorilla gorilla','Pan troglodytes','Microcebus murinus']



mix1 = ['Ornithorhynchus anatinus','Monodelphis domestica','Sarcophilus harrisii','Macropus eugenii','Dasypus novemcinctus','Loxodonta africana','Myotis lucifugus','Equus caballus','Camelus ferus','Capra hircus','Orcinus orca','Oryctolagus cuniculus','Tupaia belangeri','Sorex araneus','Erinaceus europaeus','Ictidomys tridecemlineatus','Canis lupus familiaris','Pan troglodytes']
Laurasatheria = Chiroptere + Pterisodactile + Artiodactile + Cetacea + Soricomorphe + Erinaceomorphe + Carnivore
Eutheria = Xenarthra + Afrotheria + Chiroptere + Pterisodactile + Artiodactile + Cetacea + Lagomorphes + Scandentia + Soricomorphe + Erinaceomorphe + Rodent + Carnivore + Primates
Mammals = Eutheria + Monotreme + Marsupial
Euarchontonglires =  Primates + Scandentia + Lagomorphes + Rodent
Lipotyphla = Soricomorphe + Erinaceomorphe



tt = treem.taxon_set

igls = []
tcrs = []
orden = [] # lista de labels seleccionadas
for i, t1 in enumerate(tt):
    locus =  t1.label.split('|')[3]
    especie = t1.label.split('|')[1]
    number =  t1.label.split('|')[0]

    if  especie in Lipotyphla:# <-------------------------lugar para hacer las selecciones
        orden.append(t1.label)



print 'taxones inicio ','taxones seleccionados ',len(tt),len(orden)
treem.retain_taxa_with_labels(orden)# genera el arbol seleccionado
print 'he retenido el arboles'
print treem.as_ascii_plot()
treem.write_to_path(order, 'newick')# se escribe el arbol seleccionado


print 'he ecrito el arboles'

tree = dendropy.Tree.get_from_path(order, "newick")# se lee el arbol seleccionado grabado anteriormente.
print 'he leido el arbol escrito'



pdm = treecalc.PatristicDistanceMatrix(tree)# matriz con las distancias entre nodos
print 'he hecho la matriz'

ttt = tree.taxon_set
print 'taxones en estudio ',len(ttt),len(orden)

especiesA = []
especiesB = []
L = []
N = []
Minor = []
M = []
MM = []
for i, t1 in enumerate(ttt):
    for t2 in tree.taxon_set[i+1:]:
        especieA = t1.label.split('|')[1]
        especieB = t2.label.split('|')[1]
        especiesA.append(especieA)
        especiesB.append(especieB)
        locusA = t1.label.split('|')[3]
        locusB = t2.label.split('|')[3]
        mrca = pdm.mrca(t1,t2)
        k = mrca, pdm(t1,t2),t1.label,t2.label,locusA

        if locusA == locusB and pdm(t1,t2) <= 2:

            L.append(pdm(t1,t2)) # lista de distancias de divergencia entre especies en trbv para buscar el punto de corte para el ancestro
            ss = especieA,pdm(t1,t2)
            ss2 = especieB,pdm(t1,t2)
            M.append(ss)
            MM.append(ss2)
            N.append(k)
            especiesA.append(especieA)
            especiesB.append(especieB)



print  'he seleccionado'
especies01 = []
for i in especiesA:
    if i not in especies01:
        especies01.append(i)

especies02 = []
for i in especiesB:
    if i not in especies02:
        especies02.append(i)

print especies01, '---------------------'
CC = []
AA = []
for i in especies01:
    A = []
    for k,l in M:
        if i == k:
            A.append(l)
    t = min(A)
    CC.append(t)
    AA.append(A)

BB = []
for i in especies02:
    B = []
    for k,l in MM:
        if i == k:
            B.append(l)
    t = min(A)
    BB.append(B)


#original 1.5
media = 0.7# <========================Dato basico. Se obtiene a partir del grafico. marca la zona anterior a la divergencia a estudio
print 'empiezan calculos con media en ',media
B = []
RR = []
R = []# lista de nodos dentro del area de diversificacion de especies a estudio
for i,j,l,m,n in N:
    if j <= media:
       k = i,j,l,m,n
       b =  i.distance_from_tip()
       B.append(b)
       if b <= 0.7:    # original 1.5
           R.append(k)
       else:
           RR.append(k)


q = []
for i,j,l,m,n in R:
    if i not in q:
        q.append(i)

for i in q:
    i.collapse_clade()
print 'he colapsado el arbol'
nod = tree.nodes()
qq = []
for i in q:
    if i in nod:
        qq.append(i)
listesp = []
label = []
name = []
qqq = []
count = 0
for i in qq:
    count += 1
    p = i.child_nodes()
    esp = {}
    labelb = []
    for x in p:
        aa = x.taxon
        bb = aa.label
        d = bb.split('|')[1]
        labelb.append(bb) # lista de labels por prole
        label.append(bb)# lista total de labels de proles
        if d not in esp:
           esp[d] = 1
        else:
            esp[d] +=1
    print len(labelb)
    a = p[0]
    b = a.taxon
    c= b.label
    qqq.append(c)
    h = len(p)
    e = len(esp)
    hh = c, h, e,count,esp
    name.append(hh)
    hhh = count,esp,labelb
    listesp.append(hhh)

print 'he obtenido listas de datos'
ttt = tree.taxon_set

na =[] #Lista de taxones no colsapsados
lo = [] #Lista de especies no colapsadas
for i in ttt:
    x = i.label
    if x not in label:# lista de taxones ya usados en collapse
        a =1
        d = 'one'
        b = i,a,a,d,d
        lo.append(x)
        c = x.split('|')[1]
        na.append(b)
Q = lo + qqq
P = name + na

tree = dendropy.Tree.get_from_path(order, "newick")
print 'releo el arbol original'
tree.retain_taxa_with_labels(Q) #si hay pocas especies (3 o menos) se utiliza lista Q.
tree.reroot_at_midpoint(update_splits=True)
print len(R),len(q),len(Q)
print len(lo)
tt = tree.taxon_set
yy = len(tt)


Branch = {}
for i,j,k in listesp:
    print 'branch ',i, '---->', j,
    a = str(i) + ' ' + k[0].split('|')[3]
    Branch[a] = k 

    d={}
    e={}
    for z,y in j.iteritems():
        if z in Monotreme:
            z = 'Monotreme'
            if z not in d:
                d[z] = 1
            else:
                d[z] +=1
            if z not in e:
                e[z] = y
            else:
                e[z] +=y

        if z in Marsupial:
            z='Marsupial'

            if z not in d:
                d[z] = 1

            else:
                d[z] +=1

            if z not in e:
                e[z] = y
            else:
                e[z] += y

        if z in Afrotheria:
            z='Afrotheria'

            if z not in d:
                d[z] = 1

            else:
                d[z] +=1

            if z not in e:
                e[z] = y
            else:
                e[z] +=y

        if z in Xenarthra:
            z='Xenarthra'
            if z not in d:
                d[z] = 1

            else:
                d[z] +=1

            if z not in e:
                e[z] = y
            else:
                e[z] +=y
        if z in Laurasatheria:
            z='Laurasatheria'
            if z not in d:
                d[z] = 1

            else:
                d[z] +=1
            if z not in e:
                e[z] = y
            else:
                e[z] +=y


        if z in Euarchontonglires:

            z='Euarchontonglires'

            if z not in d:
                d[z] = 1

            else:
                d[z] +=1
            if z not in e:
                e[z] = y
            else:
                e[z] +=y




    for key,value in d.iteritems():
           print key,value

    print '------------------------------------'



    for key,value in e.iteritems():
           print key,value


    print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'


ramasFile = sys.argv[3]
#ff = open('ramas.txt','r')
ff = open(ramasFile,'r')
ramas = eval(ff.read())


limp = []
for i, t1 in enumerate(tree.taxon_set):
    locus =  t1.label.split('|')[3]
    especie = t1.label.split('|')[1]

    for h,z in ramas.iteritems():
         if t1.label in z:
             b = str(h)
             print b

    for w,j,k,l,m in P:
         #if j>=6 and k >=4:
            if str(w) in t1.label:
                x = j*100/yy

                t1.label = b + '--->' + str(j) + ' (' + str(k) +')'
                limp.append(t1.label)
                if x >= 3:
                    t1.label = t1.label + '****'
                    limp.append(t1.label)

tree.retain_taxa_with_labels(limp)

print tree.as_ascii_plot()
fil = order.split('.')[0]

tree.write_to_path(fil + '01.nwk', 'newick')


branchFile=sys.argv[4]
f = open(branchFile,'w')
f.write(str(Branch))
f.close()

"""
q = especies01[0]
plt.title(q)
plt.hist(AA,bins=100)
plt.xlabel('value')
plt.ylabel('Frecuency')
plt.show()


q = especies02[1]
plt.title(q)
plt.hist(BB[1],bins=100)
plt.xlabel('value')
plt.ylabel('Frecuency')
plt.show()
"""
