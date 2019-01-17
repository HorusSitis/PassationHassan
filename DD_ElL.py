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

# ------------ les configurations géométriques sont :

## cercle unique : deux possibilités
### cen_snap_ray = [0.,0.]
### cen_snap_ray = [0.5,0.5]

# ------------ ou :

## deux inclusions, soit :

#config='compl' : on sous-entend que le paramètre est le rayon du disque central

#geo_p=='deuxième disque aux sommets'
#geo_p=='deuxième disque latéral'



#nb_lcells=
r=ray_snap_cen
cen=cen_snap_ray

r_c=0.25
r_v=0.15

if config=='cercle unique' and cen_snap_ray==[0.5,0.5]:
 l_discs=[]
 for i in range(0,nb_lcells):
   for j in range(0,nb_lcells):
    l_discs.append(plt.Circle((0.5+i,0.5+j),r,color=cem_color))
 # initialisation du graphique
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
 ax.grid(which='major', axis='x', linewidth=0.85, linestyle='-', color='0.45')
 ax.grid(which='major', axis='y', linewidth=0.85, linestyle='-', color='0.45')
 ax.set_xticklabels([])
 ax.set_yticklabels([])
 # instruction pour la sortie
 if fig_todo=='aff':
  plt.show()
 elif fig_todo=='save':
  plt.savefig("Figures2D/macro_micro"+"_Lsurl"+str(nb_lcells)+"unique_par"+str(int(round(100*cen[0],2)))+str(int(round(100*cen[1],2)))+str(int(round(100*r,2)))+".png")
 plt.close()
elif config=='cercle unique' and cen_snap_ray==[0.0,0.0]:
 l_discs=[]
 for i in range(0,1+nb_lcells):
   for j in range(0,1+nb_lcells):
    l_discs.append(plt.Circle((0.0+i,0.0+j),r,color=cem_color))
 # initialisation du graphique
 fig,ax=plt.subplots()
 # taille du graphique
 ax.set_xlim((0,nb_lcells))
 ax.set_ylim((0,nb_lcells))
 # couleur de fond
 ##fig.patch.set_facecolor(fluid_color)
 ##fig.patch.set_alpha(0.7)
 ax.set_axis_bgcolor(fluid_color)
 # inclusions, placées périodiquement
 for i in range(0,(1+nb_lcells)**2):
  ax.add_artist(l_discs[i])
 # grille
 ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))
 ax.yaxis.set_major_locator(plt.MultipleLocator(1.0))
 ax.grid(which='major', axis='x', linewidth=0.85, linestyle='-', color='0.45')
 ax.grid(which='major', axis='y', linewidth=0.85, linestyle='-', color='0.45')
 ax.set_xticklabels([])
 ax.set_yticklabels([])
 # instruction pour la sortie
 if fig_todo=='aff':
  plt.show()
 elif fig_todo=='save':
  plt.savefig("Figures2D/macro_micro"+"_Lsurl"+str(nb_lcells)+"unique_par"+str(int(round(100*cen[0],2)))+str(int(round(100*cen[1],2)))+str(int(round(100*r,2)))+".png")
 plt.close()
elif config=='compl' and geo_p=='deuxième disque aux sommets':
 l_cim=[]
 l_sand=[]
 for i in range(0,nb_lcells):
   for j in range(0,nb_lcells):
    l_sand.append(plt.Circle((0.5+i,0.5+j),r_c,color=sand_color))
 for i in range(0,1+nb_lcells):
   for j in range(0,1+nb_lcells):
    l_cim.append(plt.Circle((0.0+i,0.0+j),r_v,color=cem_color))
 # initialisation du graphique
 fig,ax=plt.subplots()
 # taille du graphique
 ax.set_xlim((0,nb_lcells))
 ax.set_ylim((0,nb_lcells))
 # couleur de fond
 ax.set_axis_bgcolor(fluid_color)
 # inclusions, placées périodiquement
 for i in range(0,nb_lcells*nb_lcells):
  ax.add_artist(l_sand[i])
 for i in range(0,(1+nb_lcells)**2):
  ax.add_artist(l_cim[i])
 # grille
 ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))
 ax.yaxis.set_major_locator(plt.MultipleLocator(1.0))
 ax.grid(which='major', axis='x', linewidth=0.85, linestyle='-', color='0.45')
 ax.grid(which='major', axis='y', linewidth=0.85, linestyle='-', color='0.45')
 ax.set_xticklabels([])
 ax.set_yticklabels([])
 if fig_todo=='aff':
  plt.show()
 elif fig_todo=='save':
  plt.savefig("Figures2D/macro_micro"+"_Lsurl"+str(nb_lcells)+"sommets_par"+str(int(round(100*cen[0],2)))+str(int(round(100*cen[1],2)))+str(int(round(100*r_c,2)))+str(int(round(100*r_v,2)))+".png")
 plt.close()
elif config=='compl' and geo_p=='deuxième disque latéral':
 l_cim=[]
 l_sand=[]
 for i in range(0,nb_lcells):
   for j in range(0,nb_lcells):
    l_sand.append(plt.Circle((0.5+i,0.5+j),r_c,color=sand_color))
 for i in range(0,1+nb_lcells):
   for j in range(0,nb_lcells):
    l_cim.append(plt.Circle((0.0+i,0.5+j),r_v,color=cem_color))
 # initialisation du graphique
 fig,ax=plt.subplots()
 # taille du graphique
 ax.set_xlim((0,nb_lcells))
 ax.set_ylim((0,nb_lcells))
 # couleur de fond
 ax.set_axis_bgcolor(fluid_color)
 # inclusions, placées périodiquement
 for i in range(0,nb_lcells*nb_lcells):
  ax.add_artist(l_sand[i])
 for i in range(0,(1+nb_lcells)*nb_lcells):
  ax.add_artist(l_cim[i])
 # grille
 ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))
 ax.yaxis.set_major_locator(plt.MultipleLocator(1.0))
 ax.grid(which='major', axis='x', linewidth=0.85, linestyle='-', color='0.45')
 ax.grid(which='major', axis='y', linewidth=0.85, linestyle='-', color='0.45')
 ax.set_xticklabels([])
 ax.set_yticklabels([])
 if fig_todo=='aff':
  plt.show()
 elif fig_todo=='save':
  plt.savefig("Figures2D/macro_micro"+"_Lsurl"+str(nb_lcells)+"lat_par"+str(int(round(100*cen[0],2)))+str(int(round(100*cen[1],2)))+str(int(round(100*r_c,2)))+str(int(round(100*r_v,2)))+".png")
 plt.close()






