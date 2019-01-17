### ------------ Paquets à importer ------------ ###

import numpy as np
import random as rd
import pylab as pl
from pylab import *
import os

#stockage d'objets python

import marshal as ma
import shelve as sh

## Paquets spécifiques à POD-MOR ##

from fenics import *
from dolfin import *
#from mshr import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
#from math import sqrt
#from math import exp
import sys

##############################################################################################################################################################
## Etape 1demi ou lL : réalisation des fugures de microstructures périodiques correspondant aux différentes géométries existantes. Stockage dans Figures2D/ ##
##############################################################################################################################################################

#nb_lcells=
if config=='cercle unique':
 l_discs=[]
 for i in range(0,nb_lcells):
   for j in range(0,nb_lcells):
    l_discs.append(plt.Circle((0.5+i,0.5+j),0.25,color=solid_color))
 # initialisation du graohique
 fig,ax=plt.subplots()
 # taille du graphique
 ax.set_xlim((0,nb_lcells))
 ax.set_ylim((0,nb_lcells))
 # couleur de fond
 ##fig.patch.set_facecolor(fluid_color)
 ##fig.patch.set_alpha(0.7)
 ax.set_axis_bgcolor(fluid_color)
 # inclusions, placées périodiquement
 for i in range(0,nb_lcells*nb_lcells):
  ax.add_artist(l_discs[i])
 # grille
 ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))
 ax.yaxis.set_major_locator(plt.MultipleLocator(1.0))
 ax.grid(which='major', axis='x', linewidth=0.85, linestyle='-', color='0.55')
 ax.grid(which='major', axis='y', linewidth=0.85, linestyle='-', color='0.45')
 ax.set_xticklabels([])
 ax.set_yticklabels([])
 # instruction pour la sortie
 if fig_todo=='aff':
  plt.show()
 elif fig_todo=='save':
  plt.savefig("Figures2D/macro_micro"+"cellules"+str(nb_lcells)+".png")#str(int(round(100*cen[0],2)))+str(int(round(100*cen[1],2)))+str(int(round(100*r,2)))+".png")
 plt.close()














