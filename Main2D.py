##Répertoire courant

cd /home/amorea12/Documents/PassationHassan




### Paquets à importer ###

import numpy as np
import random as rd


import pylab as pl
from pylab import *
#import matplotlib.pylab as pl

import os


##Codes pour le remplissage aléatoire : algorithme RSAA

from importlib import reload

import RSAA_2d_ray as R2d
#from RSAA_2d_ray import *



R2d = reload(R2d)

#threading

import threading
import time

#multiprocessing


#stockage d'objets python

import marshal as ma

import shelve as sh



################## Exemples de remplissage ##################

L=R2d.RSAA_ph_dist2D([R2d.eucl2D,R2d.vol2D],#geom : géométrie
2,#delta : marge d'erreur
[R2d.g_norm,R2d.g_norm],#l_ray : distribution des rayons par phase
[(5,0.5),(8,0.1)],#par : paramètres respectifs des distributions
[0.15,0.12],#frac_vol : fractions volumiques demandées pour les difféérentes phases
[1000,1000],#dim : taille du système
1000000,#temps : plus grand temps d'éxécution accepté
'per')#périodicité

LL=R2d.RSAA_ph_dist2D([R2d.eucl2D,R2d.vol2D],2,[R2d.g_norm,R2d.g_norm],[(3,1),(1.5,5)],[0.15,0.12],[150,150],20000,'per')

LLL=R2d.RSAA_ph_dist2D([R2d.eucl2D,R2d.vol2D],#
2,#
[R2d.g_norm,R2d.g_norm],#
[(3,1),(1.5,5)],#
[0.15,0.12],#
[20,20],#
10000,#
'per')

#stockage

###sauvegardes : listes de centres-rayons-phases
##chargement sous un nom commun

ma.dump(LLL, open("L_20", 'wb')) ## Sauvegarde de la liste 
L_load = ma.load(open("L_20", "rb")) ## Rechargement de la liste

# sauvegarde de la variable maliste sous le nom "maliste" dans le fichier 'L_20_10000'
with sh.open('L_150_20000') as d_sto:
    d_sto["maliste"] = LL
 
# chargement de la variable maliste stockée dans le fichier 'L_20_10000'
with sh.open('L_150_10000') as d_loa:
    L_loa = d_loa["maliste"]

###--------------- représentation graphique ---------------###



#A_remp=R2d.Vremp2D(#
#L_loa,#
#R2d.eucl2D,#
#[1000,1000],#
#'per')

### plus rapide ###


A_remp=R2d.Vremp2D(#
L_loa,#
R2d.eucl2D,[150,150],'per')

# Couleurs : bleu pour le fluide, gris pour une première inclusion, jaune pour une deuxième
##avec pylab
cmap=matplotlib.colors.ListedColormap(['blue','grey','yellow'])
bounds=[0,1,2,3]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
#

#AA_loa=np.reshape(AA_loa,(50,50))#il faudra apporter une correction à R2d.Vremp2D

pl.imshow(A_remp,interpolation='none',cmap=cmap,norm=norm)
pl.axis('off')

pl.show()

#
###stockase éventuel

with sh.open('VA_150_10000') as d_sto:
    d_sto["maliste"] = A_remp

#with sh.open('VA_1000_1000000') as d_loa:
#    A_loa = d_loa["maliste"]

pl.imshow(A,interpolation='none',cmap=cmap,norm=norm)
pl.axis('off')

pl.show()













