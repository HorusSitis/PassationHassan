###############################################################################################################
############################### Génération de structures périodiques aléatoires ###############################
###############################################################################################################


### ------------ Paquets a importer ------------ ###

# paquets mathematiques
import numpy as np
import random as rd
from math import sqrt
from math import exp

# Affichage avec matplotlib

import matplotlib

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

from mpl_toolkits.mplot3d.axes3d import get_test_data
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Circle, PathPatch

from matplotlib.text import TextPath
from matplotlib.transforms import Affine2D

# Affichage avec pylab

import pylab as pl

# systeme dont calcul parallelle

import sys, os
import shelve as sh

import multiprocessing

# import des fonctions utiles a RSAA

import RSAA_3d_ray as rsa

# parametres graphiques

# cmap_conv = 'night'
cmap_conv = 'cement'

if cmap_conv == 'night':
    cmap = matplotlib.colors.ListedColormap(['blue','grey','yellow'])
elif cmap_conv == 'cement':
    cmap = matplotlib.colors.ListedColormap(['cyan','grey','orange'])
    # pour un fluide au plus deux phases solides
    fluid_color = 'cyan'
    colors_ph = ['grey','orange']

# import time


# caracteristiques de la cellule

xinf = -1.
yinf = -1.
zinf = -1.

size = 2.

xsup = xinf + size
ysup = yinf + size
zsup = zinf + size


# nombre maximal d'iterations

tps_remp = 150

# repertoire pour les listes d'inclustions

rep_inc = 'IncAlea3D'

## Fabrication de la liste des inclusions etiquetees	numfic,time = py.loadtxt("%s/rayon_ecriture.txt" %(repertoire_init), unpack=True) ##

## Pour generer les inclusions dans une cellule elementaire

# lois normales pour la taille des inclusions ; parametres #

lois_rayons = [rsa.g_norm, rsa.g_norm]
par_rayons = [(0.05,0.05), (0.1,0.02)]

# fractions volumiques pour les differentes phases ; marges

frac_vol = [0.15, 0.12]
delta_inc = 0.05

list_vol_cell_inc = rsa.RSAA_ph_eucl_cell(delta_inc, lois_rayons, par_rayons, frac_vol, [xinf, yinf, zinf], size, tps_remp, 'per')

# print(list_vol_cell_inc)

list_cell_inc = list_vol_cell_inc[0]

# Affichage de la cellule elementaire'grey','orange'

nb_lcells = 1

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# taille du graphique'grey','orange'
ax.set_xlim((xinf, xsup))
ax.set_ylim((yinf, ysup))
ax.set_zlim((zinf, zsup))
# A modifier selon la geometrie de la cellule


# couleur de fond
ax.set_axis_bgcolor(fluid_color)

# pour un parametrage des surfaces
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

# inclusions, une par periode et par cellule
for i in range(len(list_cell_inc)):

    # restitution des informations sur l'inclusion courante
    cen = list_cell_inc[i][0]
    r_sph = list_cell_inc[i][1]
    phase_index = list_cell_inc[i][2] - 1

    # coordonnees parametrees de la sphere
    x = 0.+cen[0]+r_sph*np.outer(np.cos(u), np.sin(v))## outer : produit terme à terme
    y = 0.+cen[1]+r_sph*np.outer(np.sin(u), np.sin(v))
    z = 0.+cen[2]+r_sph*np.outer(np.ones(np.size(u)), np.cos(v))

    # trace de la surface parametree x, y, z : ...
    ax.plot_surface(x, y, z, linewidth=0, antialiased=False, color=colors_ph[phase_index])

# grille
ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))
ax.yaxis.set_major_locator(plt.MultipleLocator(1.0))
ax.zaxis.set_major_locator(plt.MultipleLocator(1.0))
ax.grid(which='major', axis='x', linewidth=0.85, linestyle='-', color='0.45')
ax.grid(which='major', axis='y', linewidth=0.85, linestyle='-', color='0.45')
ax.grid(which='major', axis='z', linewidth=0.85, linestyle='-', color='0.45')
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])

# instruction pour la sortie
# if fig_todo=='aff':
plt.title('Inclusions pour une cellule unique, '+str(len(lois_rayons))+' phases solides, '+str(tps_remp)+' iterations RSAA')
plt.show()
plt.close()




# Affichage d'une microstructure periodique avec la cellule generee aleatoirement










# # inclusions, placees periodiquement
# for i in range(0,nb_lcells):
# for j in range(0,nb_lcells):
# for k in range(0,nb_lcells):
# x = 0.5+i+r_sph*np.outer(np.cos(u), np.sin(v))## outer : produit terme à terme
# y = 0.5+j+r_sph*np.outer(np.sin(u), np.sin(v))
# z = 0.5+k+r_sph*np.outer(np.ones(np.size(u)), np.cos(v))
# ax.plot_surface(x, y, z, linewidth=0, antialiased=False, color=cem_color)
# # grille
# ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))
# ax.yaxis.set_major_locator(plt.MultipleLocator(1.0))
# ax.zaxis.set_major_locator(plt.MultipleLocator(1.0))
# ax.grid(which='major', axis='x', linewidth=0.85, linestyle='-', color='0.45')
# ax.grid(which='major', axis='y', linewidth=0.85, linestyle='-', color='0.45')
# ax.grid(which='major', axis='z', linewidth=0.85, linestyle='-', color='0.45')
# ax.set_xticklabels([])
# ax.set_yticklabels([])
# ax.set_zticklabels([])
