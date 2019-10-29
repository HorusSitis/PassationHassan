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

# affichage etc

import matplotlib

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

import sys, os
import shelve as sh

import multiprocessing

import pylab as pl

# import des fonctions utiles a RSAA

import RSAA_3d_ray as rsaa

# parametres graphiques

# cmap_conv = 'night'
cmap_conv = 'cement'

if cmap_conv == 'night':
    cmap = matplotlib.colors.ListedColormap(['blue','grey','yellow'])
elif cmap_conv == 'cement':
    cmap = matplotlib.colors.ListedColormap(['cyan','grey','orange'])

import time


## Tache a accomplir : creer lise RSAA, matrice d'appartenance, representation graphique ##

# task = 'nothing'
task = 'listS'
# task = 'Ainc'
# task = 'graph'

# geometrie euclidienne

geo=[eucl2D,vol2D]

# lois normales pour la taille des inclusions ; parametres #

Lois=[g_norm,g_norm]
Par=[(3,1),(1.5,5)]#[(5,0.5),(8,0.1)] pour 1000

# fractions volumiAincques pour les differentes phases ; marges

frac=[0.15,0.12]
delta_inc=2

# repertoire pour les listes d'inclustions

rep_inc='IncAlea3D'

## Fabrication de la liste des inclusions etiquetees ##
