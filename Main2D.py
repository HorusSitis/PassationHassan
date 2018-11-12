### Répertoire courant ###

cd /home/amorea12/Documents/T_LaSIE/PassationHassan

### Paquets à importer ###

import numpy as np
import random as rd

import pylab as pl
from pylab import *
#import matplotlib.pylab as pl

import os

#threading

import threading
import time

#multiprocessing


#stockage d'objets python

import marshal as ma

import shelve as sh

### Algorithme RSAA arrêté ###

## Codes à importer

from importlib import reload

###############################################################################################################
############################### Génération de structures périodiques aléatoires ###############################
###############################################################################################################

import RSAA_2d_ray as R2d

R2d = reload(R2d)

## Génération d'inclusions dans un réseau carré, stockage ##

cote=70#500 ; 1000
Tps=10000#1000000

# géométrie euclidienne

geo=[R2d.eucl2D,R2d.vol2D]

# lois normales pour la taille des inclusions ; paramètres #

Lois=[R2d.g_norm,R2d.g_norm]
Par=[(3,1),(1.5,5)]#[(5,0.5),(8,0.1)] pour 1000

# fractions volumiques pour les différentes phases ; marges

frac=[0.15,0.12]
delta_inc=2

L=R2d.RSAA_ph_dist2D(geo,delta_inc,Lois,Par,frac,[cote,cote],Tps,'per')

l_name='L_'+str(cote)+'_'+str(Tps)

rep_inc='IncAlea'

# stockage : listes de centres-rayons-phases et fractions volumiques #

#ma.dump(LLL, open("L_20", 'wb')) ## Sauvegarde de la liste 
#L_load = ma.load(open("L_20", "rb")) ## Rechargement de la liste

# sauvegarde de la variable maliste sous le nom "maliste" dans le fichier 'L_20_10000'
with sh.open(rep_inc+'/'+l_name) as l_sto:
    l_sto["maliste"] = L
 
# chargement de la variable maliste, on recharge les variables de dimension spatiale et du nombre de tours de boucle RSAA

#cote=
#Tps=
l_name='L_'+str(cote)+'_'+str(Tps)

with sh.open(rep_inc+'/'+l_name) as l_loa:
    L_loa = l_loa["maliste"]

## création d'une matrice contenant les inclustions : un coefficient pour un pixel ##

A_remp=R2d.Vremp2D(L_loa,R2d.eucl2D,[cote,cote],'per')

# stockage de la matrice A #

a_name='VA_'+str(cote)+'_'+str(Tps)
#répertoire ?

rep_inc='IncAlea'

with sh.open(rep_inc+'/'+a_name) as a_sto:
    a_sto["maliste"] = A_remp

#cote=
#Tps=
a_name='VA_'+str(cote)+'_'+str(Tps)

with sh.open(rep_inc+'/'+a_name) as a_loa:
    A_loa = a_loa["maliste"]

## représentation graphique ##

# Couleurs : bleu pour le fluide, gris pour une première inclusion, jaune pour une deuxième
##avec pylab
cmap=matplotlib.colors.ListedColormap(['blue','grey','yellow'])
bounds=[0,1,2,3]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
#

pl.imshow(A_remp,interpolation='none',cmap=cmap,norm=norm)
pl.axis('off')

figname=str(cote)+'f'+str(cote)+'Tps'+str(Tps)+'.png'
rep='Figures'
save_name=rep+'/'+figname

savefig(save_name)

pl.show()


##############################################################################################################################
############################### Calculs avec la POD, modèles réduits ; merci à Hassan GHRAIEB. ###############################
##############################################################################################################################

#Paquets spécifiques à POD-MOR

from fenics import *
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
#import numpy as np
from math import sqrt
import sys

#Import de codes de Hassan, pas de modification a priori. Reloads ?

cd 2DHassan/
#import calc_POD as PoD
#import Lecture_ecriture_homog as LEct
from calc_POD import *
from Lecture_ecriture_homog import *
#POD= reload(POD)
cd ../

#Imports divers : fonctions, objets etc

from DD_fun_obj import *#attention pas de chiffre au début du nom d'un paquet

import DD_fun_obj as fun_obj

fun_obj=reload(fun_obj)

### Codes éxécutés : cas d'une inclusion circulaire, le rayon du disque central est le paramètre pour la POD ###

## Etape I : réalisation des clichés, avec la méthode des éléments finis. Stockage dans snap2D/ ##
# Utilise fenics, éventuellement augmenté par dolfin #

#Caractérisation de la cellule élémentaire, tolérance ??

tol=1e-10

xinf=0.0
yinf=0.0
xsup=1.0
ysup=1.0
c_x=0.5
c_y=0.5

#determiner le domaine fixe pour interpoler la solution

domaine_fixe=Rectangle(Point(xinf,yinf),Point(xsup,ysup))
mesh_fixe=generate_mesh(domaine_fixe,80)
V_fixe=VectorFunctionSpace(mesh_fixe, "P", 2, constrained_domain=PeriodicBoundary())

plot(mesh_fixe)
plt.show()



for i in range(1,4):#attention le rayon d'un cercle doit être non nul
 r=i*0.1
 rect= Rectangle(Point(xinf,yinf),Point(xsup,ysup))
 circle = Circle(Point(c_x,c_y),r)
 domain=rect-circle
 res = 40  # Resolution of mesh
 mesh = generate_mesh(domain, res)
 #### Premier raffinement dans l'axe horizontale du cylindre 
 meshB= fun_obj.raffinemment_maillage(Point(c_x,c_y),r,mesh)
 mesh=meshB
 plot(mesh)
 plt.show()








#Stockage



## Etape II : extrapolation des clichés, domaine_fixe ##

#Attention aux codes de Hassan : erreurs ... visibles sur les figures du rapport, s'il s'agit bien des snapshots extrapolés


## Etape III : en utilisant la méthode des snapshots, calcul de la POD et des coefficients aléatoires, toujours dans domaine_fixe ##


## Etape IV : Prédictions ##

# SE1 : ... #

## etc ##










