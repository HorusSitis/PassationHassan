###############################################################################################################
############################### Génération de structures périodiques aléatoires ###############################
###############################################################################################################

# from RSAA_2d_ray import *


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

# import time


## Tache a accomplir : creer lise RSAA, matrice d'appartenance, representation graphique ##

# task = 'nothing'
task = 'listS'
# task = 'Ainc'
# task = 'graph'

# caracteristiques de la cellule

xinf = -1.
yinf = -1.
zinf = -1.

size = 2.

xsup = xinf + size
ysup = yinf + size
zsup = zinf + size

# # geometrie euclidienne
#
# geo = [rsa.eucl3D, rsa.vol3D]
#
# # lois normales pour la taille des inclusions ; parametres #
#
# lois_rayons = [rsa.g_norm, rsa.g_norm]
# par_rayons = [(0.5,3), (1.5,2)]
#
# # fractions volumiques pour les differentes phases ; marges
#
# frac_vol = [0.15, 0.12]
# delta_inc = 0.05

# nombre maximal d'iterations

tps_remp = 20

# repertoire pour les listes d'inclustions

rep_inc = 'IncAlea3D'

## Fabrication de la liste des inclusions etiquetees ##

# list_vol_inc = rsa.RSAA_ph_dist3D(geo, 2, lois_rayons, par_rayons, frac_vol, [10, 10, 10], tps_remp, 'per')
#
# list_inc = list_vol_inc[0]
# length = len(list_inc)
# array_inc = np.zeros((length, 3))
#
#
# print(list_inc)

# for i in range(length):
#     for k in range(3):
#         array_inc[i, k] = list_inc[i][k]


## Pour generer les inclusions dans une cellule elementaire

# lois normales pour la taille des inclusions ; parametres #

lois_rayons = [rsa.g_norm, rsa.g_norm]
par_rayons = [(0.05,0.05), (0.1,0.02)]

# fractions volumiques pour les differentes phases ; marges

frac_vol = [0.15, 0.12]
delta_inc = 0.05

list_vol_cell_inc = rsa.RSAA_ph_eucl_cell(delta_inc, lois_rayons, par_rayons, frac_vol, [xinf, yinf, zinf], size, tps_remp, 'per')

print(list_vol_cell_inc)

list_cell_inc = list_vol_cell_inc[0]

# Affichage de la cellule elementaire









# Affichage d'une microstructure periodique avec la cellule generee aleatoirement
