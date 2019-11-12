# -*- coding: utf-8 -*-
### Une commande possible dans le terminal ###

#--- mpirun -np 8 python3 Main3D.py ---#
#--- affiche npfois 'pas encore' fait avec l'etape IV ---#

# Attention : on execute parallelement

### Paquets a importer ###

# from fenics import *
# from dolfin import *
# from mshr import *
# import matplotlib.pyplot as plt
import numpy as np
# from math import sqrt
# import sys

# Calcul parallele

import multiprocessing

# Performances

import time

# Stockage d'objets python

import marshal as ma
import shelve as sh

# Fonctions maison

from DDD_fun_obj import *

## Paquets specifiques a la 3d ##

from mpl_toolkits.mplot3d.axes3d import get_test_data
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Circle, PathPatch

##########################################################
### ------------ Code a lire : conditions ------------ ###
##########################################################

list_rho_appr = np.linspace(0.05, 0.40, 8)
# list_rho_appr = np.linspace(0,1, 0.45, 8)

list_rho_test = np.linspace(0.11, 0.44, 4)

# Choix de la resolution du maillage : nombre de noeuds par cote du cube

res_gmsh=10
if typ_msh=='gms':
    res=res_gmsh

# -------------------- Geometrie du probleme -------------------- #

config='sph_un'
# config='cyl_un'
# config='cylsph'
# config='2sph'

### inclusions simples
if config=='sph_un':

    dom_fixe="am"
    geo_p='ray'#'cen'#
    ##
    conf_mess='sphere unique'

    if geo_p=='ray':
        geo_mess='rayon variable'
        cen_snap_ray=[0.5,0.5,0.5]
    elif geo_p=='cen':
        geo_mess='centre variable'
        ray_snap_cen=0.35
        csr_list=[[0.5,0.5,0.3+0.05*k] for k in range(1,1+Nsnap)]

elif config=='cyl_un':

    dom_fixe="am"
    geo_p='ray'#'axe'#
    ##
    conf_mess='cylindre unique'
    if geo_p=='ray':
        geo_mess='rayon variable'
        cen_snap_ray=[0.5,0.,0.5]
        top_snap_ray=[0.5,0.5]
    elif geo_p=='axe':
        asr_list=[[0.5,0.3+0.05*k] for k in range(1,1+Nsnap)]

### inclusions composees
elif config=='2sph':
    conf_mess='deux spheres'
    dom_fixe="solid"#"am"#

    geo_p='ray'
    geo_mess='rayon de la sphere centrale variable'
    ##

elif config=='cylsph':
    conf_mess='un cylindre et une sphere'
    dom_fixe="ray_min"#"am"#"solid"#

    geo_p='ray_sph'#'ray_cyl'#'ray_linked'
    ###ray_sph pour le test diff_ion###
    ##
    if geo_p=='ray_cyl':
        geo_mess='rayon du cylindre variable'
        ## utilisation du domaine fixe avec annulation du rayon du cylindre dans le fichier general ##
        fixe_comp='cyl_sph'#'sph_un'#'ray_min'#
    elif geo_p=='ray_sph':
        geo_mess='rayon de la sphere variable'
    elif geo_p=='ray_linked':
        geo_mess='rayons lies'
