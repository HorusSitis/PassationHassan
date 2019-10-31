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

# localisations

fig_dir = 'Figures3D/'
# fig_dir = 'IncAlea3D/'

# parametres graphiques

# fig_todo = 'aff'
fig_todo = 'save'

# cmap_conv = 'night'
cmap_conv = 'cement'

if cmap_conv == 'night':
    cmap = matplotlib.colors.ListedColormap(['blue','grey','yellow'])
elif cmap_conv == 'cement':
    cmap = matplotlib.colors.ListedColormap(['cyan','grey','orange'])
    # pour un fluide au plus deux phases solides
    fluid_color = 'cyan'
    colors_ph = ['grey','orange']

# caracteristiques de la cellule

xinf = -1.
yinf = -1.
zinf = -1.

size = 2.

xsup = xinf + size
ysup = yinf + size
zsup = zinf + size

# nombre maximal d'iterations

tps_remp = 20

# lois normales pour la taille des inclusions ; parametres #

# lois_rayons = [rsa.g_norm, rsa.g_norm]

lois_rayons = [rsa.unif, rsa.g_norm]

par_rayons = [(0.01, 0.14), (0.18, 0.02)]

# fractions volumiques pour les differentes phases ; marges

frac_vol = [0.15, 0.16]
delta_inc = 0.01

# Affichage de la cellule elementaire

nb_lcells = 3

u_res = 80
v_res = 50

### ------------ Pour generer les inclusions dans une cellule elementaire ------------ ###
# Fabrication de la liste des inclusions etiquetees

list_vol_cell_inc = rsa.RSAA_ph_cylsph_cell(delta_inc, lois_rayons, par_rayons, frac_vol, [xinf, yinf, zinf], size, tps_remp)

# print(list_vol_cell_inc)

list_cell_inc = list_vol_cell_inc[0]
array_vol = list_vol_cell_inc[1]

nb_inc = len(list_cell_inc)

print('='*60)
print('Premieres inclusions :', list_cell_inc[0:10])
print('='*60)
print('Fractions volumiques des phases solides :', array_vol)
print('Nombre d inclusions solides :', nb_inc)
print('='*60)

# ### ------------ Affichage de la la cellule generee aleatoirement ------------ ###
#
# # fig = plt.figure()
# fig_cell = plt.figure()
# # ax_cell = fig.add_subplot(111, projection='3d')
# ax_cell = fig_cell.add_subplot(111, projection='3d')
# # taille du graphique
# ax_cell.set_xlim((xinf, xsup))
# ax_cell.set_ylim((yinf, ysup))
# ax_cell.set_zlim((zinf, zsup))
# # A modifier selon la geometrie de la cellule
#
#
# # couleur de fond
# ax_cell.set_axis_bgcolor(fluid_color)
#
# # pour un parametrage des surfaces
# u = np.linspace(0, 2 * np.pi, u_res)
# v = np.linspace(0, np.pi, v_res)
#
#
# # inclusions dans la cellule elementaire
# for i in range(len(list_cell_inc)):
#
#     # restitution des informations sur l'inclusion courante
#     cen = list_cell_inc[i][0]
#     r_sph = list_cell_inc[i][1]
#     phase_index = list_cell_inc[i][2] - 1
#
#     # coordonnees parametrees de la sphere
#     x = cen[0]+r_sph*np.outer(np.cos(u), np.sin(v))## outer : produit terme à terme
#     y = cen[1]+r_sph*np.outer(np.sin(u), np.sin(v))
#     z = cen[2]+r_sph*np.outer(np.ones(np.size(u)), np.cos(v))
#
#     # trace de la surface parametree x, y, z : ...
#     ax_cell.plot_surface(x, y, z, linewidth=0, antialiased=False, color=colors_ph[phase_index])
#
#
# # grille
# ax_cell.xaxis.set_major_locator(plt.MultipleLocator(1.0))
# ax_cell.yaxis.set_major_locator(plt.MultipleLocator(1.0))
# ax_cell.zaxis.set_major_locator(plt.MultipleLocator(1.0))
# ax_cell.grid(which='major', axis='x', linewidth=0.85, linestyle='-', color='0.45')
# ax_cell.grid(which='major', axis='y', linewidth=0.85, linestyle='-', color='0.45')
# ax_cell.grid(which='major', axis='z', linewidth=0.85, linestyle='-', color='0.45')
# ax_cell.set_xticklabels([])
# ax_cell.set_yticklabels([])
# ax_cell.set_zticklabels([])
#
# # instruction pour la sortie
# if fig_todo == 'aff':
#     plt.title('Inclusions pour une cellule, '+ str(len(lois_rayons)) + ' phases solides, ' + str(tps_remp) + ' iterations RSAA')
#     plt.show()
# elif fig_todo == 'save':
#     plt.savefig(fig_dir + 'inc_cell' + '_solph'+ str(len(lois_rayons)) + '_tps' + str(tps_remp) + '.png')
#
# plt.close()
#
# ### ------------ Affichage d'une microstructure periodique avec la cellule generee aleatoirement ------------ ###
#
# # remplacer la deuxieme figure par un subplot ?
# fig_per = plt.figure()
# ax_per = fig_per.add_subplot(111, projection='3d')
# # ax_per = fig.add_subplot(112, projection='3d')
#
# # taille du graphique
# ax_per.set_xlim((xinf, xinf + nb_lcells*size))
# ax_per.set_ylim((yinf, yinf + nb_lcells*size))
# ax_per.set_zlim((zinf, zinf + nb_lcells*size))
# # A modifier selon la geometrie de la cellule
#
#
# # couleur de fond
# ax_per.set_axis_bgcolor(fluid_color)
#
# # pour un parametrage des surfaces
# u = np.linspace(0, 2 * np.pi, u_res)
# v = np.linspace(0, np.pi, v_res)
#
#
# # boucle sur les cellules
# for a in range(0,nb_lcells):
#     for b in range(0,nb_lcells):
#         for c in range(0,nb_lcells):
#
#             # inclusions, une par periode et par cellule
#             for i in range(len(list_cell_inc)):
#
#                 # restitution des informations sur l'inclusion courante
#                 cen = list_cell_inc[i][0]
#                 r_sph = list_cell_inc[i][1]
#                 phase_index = list_cell_inc[i][2] - 1
#
#                 # coordonnees parametrees de la sphere
#                 x = size*a+cen[0]+r_sph*np.outer(np.cos(u), np.sin(v))## outer : produit terme à terme
#                 y = size*b+cen[1]+r_sph*np.outer(np.sin(u), np.sin(v))
#                 z = size*c+cen[2]+r_sph*np.outer(np.ones(np.size(u)), np.cos(v))
#
#                 # trace de la surface parametree x, y, z : ...
#                 ax_per.plot_surface(x, y, z, linewidth=0, antialiased=False, color=colors_ph[phase_index])
#
#
# # grille
# ax_per.xaxis.set_major_locator(plt.MultipleLocator(1.0))
# ax_per.yaxis.set_major_locator(plt.MultipleLocator(1.0))
# ax_per.zaxis.set_major_locator(plt.MultipleLocator(1.0))
# ax_per.grid(which='major', axis='x', linewidth=0.85, linestyle='-', color='0.45')
# ax_per.grid(which='major', axis='y', linewidth=0.85, linestyle='-', color='0.45')
# ax_per.grid(which='major', axis='z', linewidth=0.85, linestyle='-', color='0.45')
# ax_per.set_xticklabels([])
# ax_per.set_yticklabels([])
# ax_per.set_zticklabels([])
#
# # instruction pour la sortie
# if fig_todo == 'aff':
#     plt.title('Inclusions pour ' + str(nb_lcells**3) + ' periodes, ' + str(len(lois_rayons)) + ' phases solides, ' + str(tps_remp) + ' iterations RSAA')
#     plt.show()
# elif fig_todo == 'save':
#     plt.savefig(fig_dir + 'inc_per3d' + str(nb_lcells) + '_solph'+ str(len(lois_rayons)) + '_tps' + str(tps_remp) + '.png')
#
# plt.close()
