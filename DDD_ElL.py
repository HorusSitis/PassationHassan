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

## Paquets spécifiques à la 3d ##

from mpl_toolkits.mplot3d.axes3d import get_test_data
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Circle, PathPatch

##############################################################################################################################################################
## Etape 1demi ou lL : réalisation des fugures de microstructures périodiques correspondant aux différentes géométries existantes. Stockage dans Figures3D/ ##
##############################################################################################################################################################
#import numpy as np
#import matplotlib.pyplot as plt
#from matplotlib.patches import Circle, PathPatch
from matplotlib.text import TextPath
from matplotlib.transforms import Affine2D
# This import registers the 3D projection, but is otherwise unused.
#from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
#import mpl_toolkits.mplot3d.art3d as art3d
# ------------ les configurations géométriques sont :

## inclusion unique : deux possibilités
### 'sphère unique'
### 'cylindre unique'

# ------------ ou :

## deux inclusions, soit :

#config='compl_ss' : on sous-entend que le paramètre est le rayon de la sphère centrale
##geo_p==''

#config='compl_sc'

##geo_p=='rayon des quarts de cylindres variable'
##geo_p=='rayon de la sphère variable'



#nb_lcells=
#r=ray_snap_cen
# cen=cen_snap_ray
#cen_snap_ray=0.35

r_cyl=0.5
r_sph=0.7
r_sph_v=0.2

# ------------ code à éxécuter pour réaliser et exploiter les figures ------------ #


if config=='sph_un':
    cen=cen_snap_ray
    # initialisation du graphique
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # taille du graphique
    ax.set_xlim((xinf, xinf + nb_lcells*(xsup - xinf)))
    ax.set_ylim((xinf, xinf + nb_lcells*(xsup - xinf)))
    ax.set_zlim((xinf, xinf + nb_lcells*(xsup - xinf)))
    # couleur de fond
    ##fig.patch.set_facecolor(fluid_color)
    ##fig.patch.set_alpha(0.7)
    # ax.set_axis_bgcolor(fluid_color)
    # pour un paramétrage des surfaces
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    # inclusions, placées périodiquement
    for i in range(0,nb_lcells):
        for j in range(0,nb_lcells):
            for k in range(0,nb_lcells):
                x = 0.5*(xsup + xinf)+i*(xsup - xinf)+r_sph*np.outer(np.cos(u), np.sin(v))## outer : produit terme à terme
                y = 0.5*(xsup + xinf)+j*(xsup - xinf)+r_sph*np.outer(np.sin(u), np.sin(v))
                z = 0.5*(xsup + xinf)+k*(xsup - xinf)+r_sph*np.outer(np.ones(np.size(u)), np.cos(v))
                ax.plot_surface(x, y, z, linewidth=0, antialiased=False, color=cem_color)
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
    if fig_todo=='aff':
        plt.show()
    elif fig_todo=='save':
        # plt.savefig("Figures3D/macro_micro"+"_Lsurl"+str(nb_lcells)+"unique_par"+str(int(round(100*cen[0],2)))+str(int(round(100*cen[1],2)))+str(int(round(100*r_sph,2)))+".png")
        plt.savefig("Figures3D/macro_micro"+"_Lsurl"+str(nb_lcells)+'_'+config+'_l_article'+".png")
    plt.close()
    #############################################################################################################
    #############################################################################################################
elif config=='cyl_un':
    # initialisation du graphique
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # taille du graphique
    ax.set_xlim((xinf, xinf + nb_lcells*(xsup - xinf)))
    ax.set_ylim((xinf, xinf + nb_lcells*(xsup - xinf)))
    ax.set_zlim((xinf, xinf + nb_lcells*(xsup - xinf)))
    # couleur de fond
    # ax.set_axis_bgcolor(fluid_color)
    # pour un paramétrage des surfaces
    theta = np.linspace(0, 2 * np.pi, 100)
    h = np.linspace(0, nb_lcells, 100)
    rc = np.linspace(0,r_cyl,100)
    # inclusions, placées périodiquement
    ## cylindres
    for i in range(0,nb_lcells):
        for j in range(0,nb_lcells):
            x = 0.5*(xsup + xinf)+i*(xsup - xinf)+r_cyl*np.outer(np.cos(theta),np.ones(np.size(h)))
            y = np.outer(np.ones(np.size(theta)),h)
            z = 0.5*(xsup + xinf)+j*(xsup - xinf)+r_cyl*np.outer(np.sin(theta),np.ones(np.size(h)))
            ax.plot_surface(x, y, z, linewidth=0, antialiased=False, color=cem_color)
    ## disques
    for i in range(0,nb_lcells):
        for j in range(0,nb_lcells):
            p_ij1 = Circle((0.5*(xsup + xinf)+i*(xsup - xinf), 0.5*(xsup + xinf)+j*(xsup - xinf)), r_cyl, color=cem_color)
            p_ij2 = Circle((0.5*(xsup + xinf)+i*(xsup - xinf), 0.5*(xsup + xinf)+j*(xsup - xinf)), r_cyl, color=cem_color)
            ax.add_patch(p_ij1)
            art3d.pathpatch_2d_to_3d(p_ij1, z=0, zdir="y")
            ax.add_patch(p_ij2)
            art3d.pathpatch_2d_to_3d(p_ij2, z=nb_lcells, zdir="y")
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
    if fig_todo=='aff':
        plt.show()
    elif fig_todo=='save':
        plt.savefig("Figures3D/macro_micro"+"_Lsurl"+str(nb_lcells)+"unique_cyl"+str(int(round(100*r_cyl,2)))+".png")
    plt.close()
    #############################################################################################################
    #############################################################################################################
elif config=='2sph':
    # initialisation du graphique
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # taille du graphique
    ax.set_xlim((xinf, xinf + nb_lcells*(xsup - xinf)))
    ax.set_ylim((xinf, xinf + nb_lcells*(xsup - xinf)))
    ax.set_zlim((xinf, xinf + nb_lcells*(xsup - xinf)))
    # couleur de fond
    # ax.set_axis_bgcolor(fluid_color)
    # pour un paramétrage des surfaces
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    # inclusions, placées périodiquement
    for i in range(0,nb_lcells):
        for j in range(0,nb_lcells):
            for k in range(0,nb_lcells):
                xc = 0.5*(xsup + xinf)+i*(xsup - xinf)+r_sph*np.outer(np.cos(u), np.sin(v))## outer : produit cartésien
                yc = 0.5*(xsup + xinf)+j*(xsup - xinf)+r_sph*np.outer(np.sin(u), np.sin(v))
                zc = 0.5*(xsup + xinf)+k*(xsup - xinf)+r_sph*np.outer(np.ones(np.size(u)), np.cos(v))
                ax.plot_surface(xc, yc, zc, linewidth=0, antialiased=False, color=sand_color)
    for i in range(0,1+nb_lcells):
        for j in range(0,1+nb_lcells):
            for k in range(0,1+nb_lcells):
                xv = xinf+i*(xsup - xinf)+r_sph_v*np.outer(np.cos(u), np.sin(v))
                yv = xinf+j*(xsup - xinf)+r_sph_v*np.outer(np.sin(u), np.sin(v))
                zv = xinf+k*(xsup - xinf)+r_sph_v*np.outer(np.ones(np.size(u)), np.cos(v))
                ax.plot_surface(xv, yv, zv, linewidth=0, antialiased=False, color=cem_color)
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
    if fig_todo=='aff':
        plt.show()
    elif fig_todo=='save':
        plt.savefig("Figures3D/macro_micro"+"_Lsurl"+str(nb_lcells)+"diag_par"+str(int(round(100*cen[0],2)))+str(int(round(100*cen[1],2)))+str(int(round(100*r_sph,2)))+".png")
    plt.close()
    #############################################################################################################
    #############################################################################################################
elif config=='cylsph':
    # initialisation du graphique
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # taille du graphique
    ax.set_xlim((xinf, xinf + nb_lcells*(xsup - xinf)))
    ax.set_ylim((xinf, xinf + nb_lcells*(xsup - xinf)))
    ax.set_zlim((xinf, xinf + nb_lcells*(xsup - xinf)))
    # couleur de fond
    ##fig.patch.set_facecolor(fluid_color)
    ##fig.patch.set_alpha(0.7)
    # ax.set_axis_bgcolor(fluid_color)
    # pour un paramétrage des surfaces liées à la sphère
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    # inclusions, placées périodiquement
    for i in range(0,nb_lcells):
        for j in range(0,nb_lcells):
            for k in range(0,nb_lcells):
                x = 0.5*(xsup + xinf)+i*(xsup - xinf)+r_sph*np.outer(np.cos(u), np.sin(v))## outer : produit cartésien
                y = 0.5*(xsup + xinf)+j*(xsup - xinf)+r_sph*np.outer(np.sin(u), np.sin(v))
                z = 0.5*(xsup + xinf)+k*(xsup - xinf)+r_sph*np.outer(np.ones(np.size(u)), np.cos(v))
                ax.plot_surface(x, y, z, linewidth=0, antialiased=False, color=sand_color)
    # pour un paramétrage des surfaces liées au cylindre
    theta = np.linspace(0, 2 * np.pi, 100)
    h = np.linspace(xinf, xinf + nb_lcells*(xsup - xinf), 100)
    rc = np.linspace(0,r_cyl,100)
    # inclusions, placées périodiquement
    ## cylindres
    for i in range(0,1+nb_lcells):
        for j in range(0,1+nb_lcells):
            x = xinf+i*(xsup - xinf)+r_cyl*np.outer(np.cos(theta),np.ones(np.size(h)))
            y = np.outer(np.ones(np.size(theta)),h)
            z = xinf+j*(xsup - xinf)+r_cyl*np.outer(np.sin(theta),np.ones(np.size(h)))
            ax.plot_surface(x, y, z, linewidth=0, antialiased=False, color=cem_color)
    ## disques
    for i in range(0,nb_lcells):
        for j in range(0,nb_lcells):
            p_ij1 = Circle((xinf+i*(xsup - xinf), xinf+j*(xsup - xinf)), r_cyl, color=cem_color)
            p_ij2 = Circle((xinf+i*(xsup - xinf), xinf+j*(xsup - xinf)), r_cyl, color=cem_color)
            ax.add_patch(p_ij1)
            art3d.pathpatch_2d_to_3d(p_ij1, z=xinf, zdir="y")
            ax.add_patch(p_ij2)
            art3d.pathpatch_2d_to_3d(p_ij2, z=(xsup - xinf)*nb_lcells, zdir="y")
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
    if fig_todo=='aff':
        plt.show()
    elif fig_todo=='save':
        # plt.savefig("Figures3D/macro_micro"+"_Lsurl"+str(nb_lcells)+"sph_cyl_par"+str(int(round(100*cen[0],2)))+str(int(round(100*cen[1],2)))+str(int(round(100*r_sph,2)))+".png")
        plt.savefig("Figures3D/macro_micro"+"_Lsurl"+str(nb_lcells)+'_'+config+'_l_article'+".png")
    plt.close()
