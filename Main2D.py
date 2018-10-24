### Paquets à importer ###

import numpy as np
import random as rd

#import numpy.random as rd
#from math import *
#import pylab as pl
#from scipy import *                 
#from pylab import *

import pylab as pl
from pylab import *
#import matplotlib.pylab as pl

import os

##Répertoire courant

cd /home/amorea12/Documents/T_LaSIE/Codes

##Codes pour le remplissage aléatoire : algorithme RSAA

import RSAA_2d_ray as R2d
#from RSAA_2d_ray import *

R2d = reload(R2d)

##Exemples de remplissage

L=R2d.RSAA_ph_dist2D([eucl2D,vol2D],#geom : géométrie
2,#delta : marge d'erreur
[g_norm,g_norm],#l_ray : distribution des rayons par phase
[(5,0.5),(8,0.1)],#par : paramètres respectifs des distributions
[0.15,0.12],#frac_vol : fractions volumiques demandées pour les difféérentes phases
[500,500],#dim : taille du système
10000,#temps : plus grand temps d'éxécution accepté
'per')#périodicité

LL=R2d.RSAA_ph_dist2D([eucl2D,vol2D],2,[g_norm,g_norm],[(3,1),(1.5,5)],[0.15,0.12],[100,100],10000,'per')

##Représentation graphique

AA=remp2D(#
LL,#RSAA_ph_dist2D([eucl2D,vol2D],2,[g_norm,g_norm],[(5,0.5),(8,0.1)],[0.15,0.12],[500,500],100000,'per'),#
eucl2D,[100,100],'per')

A=remp2D(#
L,#RSAA_ph_dist2D([eucl2D,vol2D],2,[g_norm,g_norm],[(5,0.5),(8,0.1)],[0.15,0.12],[500,500],100000,'per'),#
eucl2D,[100,100],'per')

# Couleurs : bleu pour le fluide, gris pour une première inclusion, jaune pour une deuxième
##avec pylab
cmap=matplotlib.colors.ListedColormap(['blue','grey','yellow'])
bounds=[0,1,2,3]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
#

pl.imshow(AA,interpolation='none',cmap=cmap,norm=norm)
pl.axis('off')

pl.show()

#

pl.imshow(A,interpolation='none',cmap=cmap,norm=norm)
pl.axis('off')

pl.show()


