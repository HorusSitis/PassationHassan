###############################################################################################################
############################### Generation de structures periodiques aleatoires ###############################
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

from matplotlib.patches import Ellipse

# from matplotlib.text import TextPath
# from matplotlib.transforms import Affine2D

# Affichage avec pylab

import pylab as pl

# systeme dont calcul parallelle

import sys, os
import shelve as sh

import multiprocessing

# import des fonctions utiles a RSAA

import RSAA_2d_ray as rsa

# localisations

fig_dir = 'Figures2D/'
# fig_dir = 'IncAlea2D/'

# parametres graphiques

fig_todo = 'aff'
# fig_todo = 'save'

# cmap_conv = 'night'
cmap_conv = 'cement'

if cmap_conv == 'night':
    cmap = matplotlib.colors.ListedColormap(['blue','grey','yellow'])
elif cmap_conv == 'cement':
    cmap = matplotlib.colors.ListedColormap(['cyan','grey','orange'])
    # pour un fluide au plus deux phases solides
    fluid_color = 'cyan'
    colors_ph = ['grey','orange']

# type d'inclusions

# type_inc = 'cer'
type_inc = 'ell'

# caracteristiques de la cellule

xinf = -1.
yinf = -1.

size = 2.

xsup = xinf + size
ysup = yinf + size

# nombre maximal d'iterations

tps_remp = 100

# lois normales pour la taille des inclusions ; parametres #

# lois_rayons = [rsa.g_norm, rsa.g_norm]

lois_rayons = [rsa.unif, rsa.g_norm]

par_rayons = [(0.005, 0.07), (0.1, 0.015)]

# fractions volumiques pour les differentes phases ; marges

frac_vol = [0.2, 0.28]
delta_inc = 0.01

# Affichage de la cellule elementaire

nb_lcells = 8

u_res = 80
v_res = 50

theta_res = 50
h_res = 50
rc_res = 25

### ------------ Pour generer les inclusions dans une cellule elementaire ------------ ###
# Fabrication de la liste des inclusions etiquetees

if type_inc == 'cer':
    list_vol_cell_inc = rsa.RSAA_ph_eucl_cell(delta_inc, lois_rayons, par_rayons, frac_vol, [xinf, yinf], size, tps_remp)
elif type_inc == 'ell':
    list_vol_cell_inc = rsa.RSAA_ph_ell_cell(delta_inc, lois_rayons, par_rayons, frac_vol, [xinf, yinf], size, tps_remp)


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

### ------------ Affichage de la la cellule generee aleatoirement ------------ ###

# l_discs=[]

# l_discs.append(plt.Circle((0.5+i,0.5+j),r,color=cem_color))

# initialisation du graphique
fig, ax=plt.subplots()
# taille du graphique
ax.set_xlim((xinf, xsup))
ax.set_ylim((yinf, ysup))
# couleur de fond

ax.set_axis_bgcolor(fluid_color)

# inclusions, placees periodiquement
# for i in range(0,nb_lcells*nb_lcells):

for i in range(nb_inc):

    if type_inc == 'cer':

        phase_index = list_cell_inc[i][2] - 1

        cen = list_cell_inc[i][0]
        ray = list_cell_inc[i][1]

        ax.add_artist(plt.Circle((cen[0], cen[1]), ray, color = colors_ph[phase_index]))

    elif type_inc == 'ell':

        phase_index = list_cell_inc[i][1] - 1

        ell = list_cell_inc[i][0]

        cen = ell[0]
        gax = ell[1][0]
        pax = gax*ell[1][1]
        theta = ell[2]

        gr_ell = Ellipse(xy = cen, width = gax, height = pax, angle = theta*360/(2*np.pi))
        gr_ell.set_facecolor(colors_ph[phase_index])
        ax.add_artist(gr_ell)



ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))
ax.yaxis.set_major_locator(plt.MultipleLocator(1.0))
ax.grid(which='major', axis='x', linewidth=0.85, linestyle='-', color='0.45')
ax.grid(which='major', axis='y', linewidth=0.85, linestyle='-', color='0.45')
ax.set_xticklabels([])
ax.set_yticklabels([])

# instruction pour la sortie
if fig_todo=='aff':
    plt.title('Inclusions pour une cellule, '+ str(len(lois_rayons)) + ' phases solides, ' + str(tps_remp) + ' iterations RSAA')
    plt.show()

plt.close()






















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
# # pour un parametrage des surfaces liees au cylindre
# theta = np.linspace(0, 2 * np.pi, theta_res)
# # h = np.linspace(0, nb_lcells, h_res)
# # rc = np.linspace(0, r_cyl, rc_res)
#
# # inclusions dans la cellule elementaire
# for i in range(len(list_cell_inc)):
#
#     # restitution des informations sur l'inclusion courante
#     cen_or_axdir = list_cell_inc[i][0]
#     ray = list_cell_inc[i][1]
#     phase_index = list_cell_inc[i][2] - 1
#
#     # inclusion courante spherique
#     if phase_index == 1:
#
#         cen = cen_or_axdir
#         r_sph = ray
#
#         # coordonnees parametrees de la sphere
#         x = cen[0]+r_sph*np.outer(np.cos(u), np.sin(v))## outer : produit terme à terme
#         y = cen[1]+r_sph*np.outer(np.sin(u), np.sin(v))
#         z = cen[2]+r_sph*np.outer(np.ones(np.size(u)), np.cos(v))
#
#         # trace de la surface parametree x, y, z : ...
#         ax_cell.plot_surface(x, y, z, linewidth=0, antialiased=False, color=colors_ph[phase_index])
#
#     # inclusion courante cylindrique
#     elif phase_index == 0:
#
#         ax_cyl = cen_or_axdir[0]
#         dir = cen_or_axdir[1]
#         r_cyl = ray
#
#         # cylindre
#         h = np.linspace(-size/2, size/2, h_res)
#         rc = np.linspace(0, r_cyl, rc_res)
#
#         if dir == 0:
#             y = ax_cyl[1] + r_cyl*np.outer(np.cos(theta),np.ones(np.size(h)))
#             z = ax_cyl[2] + r_cyl*np.outer(np.sin(theta),np.ones(np.size(h)))
#             #
#             x = np.outer(np.ones(np.size(theta)),h)
#             # disques aux extremites
#             p_ij1 = Circle((ax_cyl[1], ax_cyl[2]), r_cyl, color=colors_ph[phase_index])
#             p_ij2 = Circle((ax_cyl[1], ax_cyl[2]), r_cyl, color=colors_ph[phase_index])
#             ax_cell.add_patch(p_ij1)
#             art3d.pathpatch_2d_to_3d(p_ij1, z=xinf, zdir="x")
#             ax_cell.add_patch(p_ij2)
#             art3d.pathpatch_2d_to_3d(p_ij2, z=xinf + size, zdir="x")
#         if dir == 1:
#             x = ax_cyl[0] + r_cyl*np.outer(np.cos(theta),np.ones(np.size(h)))
#             z = ax_cyl[2] + r_cyl*np.outer(np.sin(theta),np.ones(np.size(h)))
#             #
#             y = np.outer(np.ones(np.size(theta)),h)
#             # disques aux extremites
#             p_ij1 = Circle((ax_cyl[0], ax_cyl[2]), r_cyl, color=colors_ph[phase_index])
#             p_ij2 = Circle((ax_cyl[0], ax_cyl[2]), r_cyl, color=colors_ph[phase_index])
#             ax_cell.add_patch(p_ij1)
#             art3d.pathpatch_2d_to_3d(p_ij1, z=yinf, zdir="y")
#             ax_cell.add_patch(p_ij2)
#             art3d.pathpatch_2d_to_3d(p_ij2, z=yinf + size, zdir="y")
#         if dir == 2:
#             x = ax_cyl[0] + r_cyl*np.outer(np.cos(theta),np.ones(np.size(h)))
#             y = ax_cyl[1] + r_cyl*np.outer(np.sin(theta),np.ones(np.size(h)))
#             #
#             z = np.outer(np.ones(np.size(theta)),h)
#             # disques aux extremites
#             p_ij1 = Circle((ax_cyl[0], ax_cyl[1]), r_cyl, color=colors_ph[phase_index])
#             p_ij2 = Circle((ax_cyl[0], ax_cyl[1]), r_cyl, color=colors_ph[phase_index])
#             ax_cell.add_patch(p_ij1)
#             art3d.pathpatch_2d_to_3d(p_ij1, z=zinf, zdir="z")
#             ax_cell.add_patch(p_ij2)
#             art3d.pathpatch_2d_to_3d(p_ij2, z=zinf + size, zdir="z")
#
#         # affichage des cylindres de revolution
#         ax_cell.plot_surface(x, y, z, linewidth=0, antialiased=False, color=colors_ph[phase_index])
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
#     plt.title('Inclusions pour une cellule, ' + 'cylindres et spheres, ' + str(tps_remp) + ' iterations RSAA')
#     plt.show()
# elif fig_todo == 'save':
#     plt.savefig(fig_dir + 'inc_cell' + str(nb_lcells) + '_cylsph'+ str(len(lois_rayons)) + '_tps' + str(tps_remp) + '.png')
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
#             # inclusions dans la cellule elementaire
#             for i in range(len(list_cell_inc)):
#
#                 # restitution des informations sur l'inclusion courante
#                 cen_or_axdir = list_cell_inc[i][0]
#                 ray = list_cell_inc[i][1]
#                 phase_index = list_cell_inc[i][2] - 1
#
#                 # inclusion courante spherique
#                 if phase_index == 1:
#
#                     cen = cen_or_axdir
#                     r_sph = ray
#
#                     # coordonnees parametrees de la sphere
#                     x = size*a +  cen[0]+r_sph*np.outer(np.cos(u), np.sin(v))## outer : produit terme à terme
#                     y = size*b +  cen[1]+r_sph*np.outer(np.sin(u), np.sin(v))
#                     z = size*c +  cen[2]+r_sph*np.outer(np.ones(np.size(u)), np.cos(v))
#
#                     # trace de la surface parametree x, y, z : ...
#                     ax_per.plot_surface(x, y, z, linewidth=0, antialiased=False, color=colors_ph[phase_index])
#
#                 # inclusion courante cylindrique
#                 elif phase_index == 0:
#
#                     ax_cyl = cen_or_axdir[0]
#                     dir = cen_or_axdir[1]
#                     r_cyl = ray
#
#                     # cylindre
#                     h = np.linspace(-size/2, size/2, h_res)
#                     rc = np.linspace(0, r_cyl, rc_res)
#
#                     if dir == 0:
#                         y = size*b +  ax_cyl[1] + r_cyl*np.outer(np.cos(theta),np.ones(np.size(h)))
#                         z = size*c +  ax_cyl[2] + r_cyl*np.outer(np.sin(theta),np.ones(np.size(h)))
#                         #
#                         x = size*a +  np.outer(np.ones(np.size(theta)),h)
#                         # disques aux extremites
#                         p_ij1 = Circle((size*b +  ax_cyl[1], size*c +  ax_cyl[2]), r_cyl, color=colors_ph[phase_index])
#                         p_ij2 = Circle((size*b +  ax_cyl[1], size*c +  ax_cyl[2]), r_cyl, color=colors_ph[phase_index])
#                         ax_per.add_patch(p_ij1)
#                         art3d.pathpatch_2d_to_3d(p_ij1, z=xinf, zdir="x")
#                         ax_per.add_patch(p_ij2)
#                         art3d.pathpatch_2d_to_3d(p_ij2, z=xinf + size, zdir="x")
#                     if dir == 1:
#                         x = size*a +  ax_cyl[0] + r_cyl*np.outer(np.cos(theta),np.ones(np.size(h)))
#                         z = size*c +  ax_cyl[2] + r_cyl*np.outer(np.sin(theta),np.ones(np.size(h)))
#                         #
#                         y = size*b +  np.outer(np.ones(np.size(theta)),h)
#                         # disques aux extremites
#                         p_ij1 = Circle((size*a +  ax_cyl[0], size*c +  ax_cyl[2]), r_cyl, color=colors_ph[phase_index])
#                         p_ij2 = Circle((size*a +  ax_cyl[0], size*c +  ax_cyl[2]), r_cyl, color=colors_ph[phase_index])
#                         ax_per.add_patch(p_ij1)
#                         art3d.pathpatch_2d_to_3d(p_ij1, z=yinf, zdir="y")
#                         ax_per.add_patch(p_ij2)
#                         art3d.pathpatch_2d_to_3d(p_ij2, z=yinf + size, zdir="y")
#                     if dir == 2:
#                         x = size*a +  ax_cyl[0] + r_cyl*np.outer(np.cos(theta),np.ones(np.size(h)))
#                         y = size*b +  ax_cyl[1] + r_cyl*np.outer(np.sin(theta),np.ones(np.size(h)))
#                         #
#                         z = size*c +  np.outer(np.ones(np.size(theta)),h)
#                         # disques aux extremites
#                         p_ij1 = Circle((size*a +  ax_cyl[0], size*b +  ax_cyl[1]), r_cyl, color=colors_ph[phase_index])
#                         p_ij2 = Circle((size*a +  ax_cyl[0], size*b +  ax_cyl[1]), r_cyl, color=colors_ph[phase_index])
#                         ax_per.add_patch(p_ij1)
#                         art3d.pathpatch_2d_to_3d(p_ij1, z=zinf, zdir="z")
#                         ax_per.add_patch(p_ij2)
#                         art3d.pathpatch_2d_to_3d(p_ij2, z=zinf + size, zdir="z")
#
#                     # affichage des cylindres de revolution
#                     ax_per.plot_surface(x, y, z, linewidth=0, antialiased=False, color=colors_ph[phase_index])
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
#     plt.title('Inclusions pour ' + str(nb_lcells**3) + ' periodes, cylindres et spheres, ' + str(tps_remp) + ' iterations RSAA')
#     plt.show()
# elif fig_todo == 'save':
#     plt.savefig(fig_dir + 'inc_per3d' + str(nb_lcells) + '_cylsph'+ str(len(lois_rayons)) + '_tps' + str(tps_remp) + '.png')
#
# plt.close()
