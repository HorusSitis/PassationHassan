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
    cmap = matplotlib.colors.ListedColormap(['blue', 'grey', 'yellow'])
elif cmap_conv == 'cement':
    cmap = matplotlib.colors.ListedColormap(['cyan', 'grey', 'orange'])
    # pour un fluide au plus deux phases solides
    fluid_color = 'cyan'
    colors_ph = ['grey', 'orange']

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

tps_remp = 100000

# lois normales pour la taille des inclusions ; parametres #

# lois_rayons = [rsa.g_norm, rsa.g_norm]

lois_rayons = [rsa.unif, rsa.g_norm]

par_rayons = [(0.005, 0.25), (0.1, 0.015)]

# fractions volumiques pour les differentes phases ; marges

frac_vol = [0.2, 0.28]
delta_inc = 0.005

# Affichage de la cellule elementaire

nb_lcells = 4

# Affichage de la grille pour le motif periodique

# plot_ax_per = True
plot_ax_per = False

# parametres pour le trace

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
print('Premieres inclusions :', list_cell_inc[0:5])
print('='*60)
print('Fractions volumiques des phases solides :', array_vol)
print('Nombre d inclusions solides :', nb_inc)
print('='*60)

# sys.exit('%'*20+' Fin de la generation RSAA '+'%'*20)
### ------------ Affichage de la la cellule generee aleatoirement ------------ ###

# l_discs=[]

# l_discs.append(plt.Circle((0.5+i,0.5+j),r,color=cem_color))

# initialisation du graphique
fig, ax_cell=plt.subplots()
# taille du graphique
ax_cell.set_xlim((xinf, xsup))
ax_cell.set_ylim((yinf, ysup))

# couleur de fond
ax_cell.set_axis_bgcolor(fluid_color)

# inclusions, cellule unique
for i in range(nb_inc):

    if type_inc == 'cer':

        phase_index = list_cell_inc[i][2] - 1

        cen = list_cell_inc[i][0]
        ray = list_cell_inc[i][1]

        ax_cell.add_artist(plt.Circle((cen[0], cen[1]), ray, color = colors_ph[phase_index]))

    elif type_inc == 'ell':

        phase_index = list_cell_inc[i][1] - 1

        ell = list_cell_inc[i][0]

        cen = ell[0]
        gax = ell[1][0]
        pax = gax*ell[1][1]
        theta = ell[2]

        gr_ell = Ellipse(xy = cen, width = gax, height = pax, angle = theta*360/(2*np.pi))
        gr_ell.set_facecolor(colors_ph[phase_index])
        gr_ell.set_edgecolor(colors_ph[phase_index])
        ax_cell.add_artist(gr_ell)



ax_cell.xaxis.set_major_locator(plt.MultipleLocator(1.0))
ax_cell.yaxis.set_major_locator(plt.MultipleLocator(1.0))
ax_cell.grid(which='major', axis='x', linewidth=0.85, linestyle='-', color='0.45')
ax_cell.grid(which='major', axis='y', linewidth=0.85, linestyle='-', color='0.45')
ax_cell.set_xticklabels([])
ax_cell.set_yticklabels([])

# instruction pour la sortie
if fig_todo=='aff':
    if type_inc == 'cer':
        plt.title('Inclusions pour une cellule, Gamma_sf cercles, '+ str(len(lois_rayons)) + ' phases solides, ' + str(tps_remp) + ' iterations RSAA')
    elif type_inc == 'ell':
        plt.title('Inclusions pour une cellule, Gamma_sf ellipses, '+ str(len(lois_rayons)) + ' phases solides, ' + str(tps_remp) + ' iterations RSAA')
    plt.show()
elif fig_todo == 'save':
    if type_inc == 'cer':
        plt.savefig(fig_dir + 'inc_cell2d' + str(nb_lcells) + '_cer'+ str(len(lois_rayons)) + '_tps' + str(tps_remp) + '.png')
    elif type_inc == 'ell':
        plt.savefig(fig_dir + 'inc_cell2d' + str(nb_lcells) + '_ell'+ str(len(lois_rayons)) + '_tps' + str(tps_remp) + '.png')



plt.close()

# print('%'*20+' Fin de l affichage; cellule unique '+'%'*20)

### ------------ Affichage d'une microstructure periodique avec la cellule generee aleatoirement ------------ ###

# remplacer la deuxieme figure par un subplot ?
fig_per, ax_per = plt.subplots()

# taille du graphique
ax_per.set_xlim((xinf, xinf + nb_lcells*size))
ax_per.set_ylim((yinf, yinf + nb_lcells*size))

# couleur de fond
ax_per.set_axis_bgcolor(fluid_color)

# boucle sur les periodes
for a in range(0,nb_lcells):
    for b in range(0,nb_lcells):

        # inclusions, cellule translatee
        for i in range(nb_inc):

            if type_inc == 'cer':

                phase_index = list_cell_inc[i][2] - 1

                cen = size*np.array([a, b]) + list_cell_inc[i][0]
                ray = list_cell_inc[i][1]

                ax_per.add_artist(plt.Circle((cen[0], cen[1]), ray, color = colors_ph[phase_index]))

            elif type_inc == 'ell':

                phase_index = list_cell_inc[i][1] - 1

                ell = list_cell_inc[i][0]

                cen = size*np.array([a, b]) + ell[0]
                gax = ell[1][0]
                pax = gax*ell[1][1]
                theta = ell[2]

                gr_ell = Ellipse(xy = cen, width = gax, height = pax, angle = theta*360/(2*np.pi))
                gr_ell.set_edgecolor(colors_ph[phase_index])
                gr_ell.set_facecolor(colors_ph[phase_index])
                ax_per.add_artist(gr_ell)

if plot_ax_per == True:

    ax_per.xaxis.set_major_locator(plt.MultipleLocator(1.0))
    ax_per.yaxis.set_major_locator(plt.MultipleLocator(1.0))
    ax_per.grid(which='major', axis='x', linewidth=0.85, linestyle='-', color='0.45')
    ax_per.grid(which='major', axis='y', linewidth=0.85, linestyle='-', color='0.45')
    ax_per.set_xticklabels([])
    ax_per.set_yticklabels([])

# instruction pour la sortie
if fig_todo=='aff':
    if type_inc == 'cer':
        plt.title('Inclusions pour ' + str(nb_lcells**2) + ' periodes, cercles, ' + str(tps_remp) + ' iterations RSAA')
    elif type_inc == 'ell':
        plt.title('Inclusions pour ' + str(nb_lcells**2) + ' periodes, ellipses, ' + str(tps_remp) + ' iterations RSAA')
    plt.show()
elif fig_todo == 'save':
    if type_inc == 'cer':
        plt.savefig(fig_dir + 'inc_per2d' + str(nb_lcells) + '_cer'+ str(len(lois_rayons)) + '_tps' + str(tps_remp) + '.png')
    elif type_inc == 'ell':
        plt.savefig(fig_dir + 'inc_per2d' + str(nb_lcells) + '_ell'+ str(len(lois_rayons)) + '_tps' + str(tps_remp) + '.png')

plt.close()
