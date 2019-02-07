# -*- coding: utf-8 -*-
### Une commande possible dans le terminal ###

#--- mpirun -np 8 python3 Main3D.py ---#
#--- affiche npfois 'pas encore' fait avec l'étape IV ---#

# Attention : on éxécute parallèlement 

### Paquets à importer ###

from fenics import *
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
import sys

# Calcul parallèle

import multiprocessing

# Performances

import time

# Stockage d'objets python

import marshal as ma
import shelve as sh

# Fonctions maison

from DDD_fun_obj import *

## Paquets spécifiques à la 3d ##

from mpl_toolkits.mplot3d.axes3d import get_test_data
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Circle, PathPatch

##########################################################
### ------------ Code à lire : conditions ------------ ###
##########################################################

E_=False
E_lL=False

EI=False

EII=False

EIII=False
EIV=True

res_fixe=6
fixe_aff=False
#res=6
#slices_cyl=5
Nsnap=8
rempUsnap='par8'#'seq'
#c_x=0.5
#c_y=0.5
#c_z=0.5
#r=0.35#pour une réalisation unique
npas_err=20
fig_todo='save'#'aff'

typ_msh='gms'#''
## choix du domaine fixe
dom_fixe='0001'#'0001'#''

# nom de l'appareil utilisé pour générer les données enregistrées
computer='MECALAC_29x8'##'T1700_35x8'#



repertoire_parent="Res3D/"


## -------------------- Etape I -------------------- ##

D_k=1.0
Nsnap=8
npas_err=20
ordo='Ordr'#'Nordr'

# Parallélisation du calcul des snapshots

parallelize=True

# Choix de la résolution du maillage : nombre de noeuds par côté du cube

res_gmsh=10

if typ_msh=='gms':
 res=res_gmsh

# Choix du paramètre géométrique : variable et message de sortie _mess

## Cas d'une sphère unique
config='sph_un'
if config=='sph_un':
 conf_mess='sphère unique'

geo_p='ray'#
if geo_p=='ray':
 geo_mess='rayon variable'
cen_snap_ray=[0.5,0.5,0.5]
#geo_p='centre'
#ray_snap_cen=5
#csr_list=[[0.5,0.5,0.3+0.05*k] for k in range(1,1+Nsnap)]
### l'inclusion solide ne rencontre pas les bords du domaine, sous peine d'une incompatibilité entre les mailles de faces opposées (?)

## Cas d'un cylindre unique : axe parallèle à Oy
#config='cyl_un'

if config=='cyl_un':
 conf_mess='cylindre unique'
#geo_p='ray'
#geo_p='axe'
asr_list=[[0.5,0.3+0.05*k] for k in range(1,1+Nsnap)]
### point auquel l'axe rencontre le plan Oxz. L'inclusion solide ne rencontre pas les bords du domaine.

## Cas d'un cylindre périodique aux arètes et une sphère unique au centre
#config='compl'
#geo_p='sphh'
#geo_p='sph_cyl'

## ------------ Etape lL 1demi : Affichage de midrostructures périodiques ------------ ##

nb_lcells=5
cem_color='grey'#'r'
sand_color='orange'
fluid_color='cyan'

if E_lL :
 exec(open("DDD_ElL.py").read())

# Exécution

snap_done=True #-------------------> pour calculer les snapshots seulement si nécessaire

if EI :
 exec(open("DDD_EI.py").read())

##computer=
### Pour les étapes qui suivent, on peut choisir l'ordinateur qui a effectué le calcul des snapshots physiques

## -------------------- Etape II -------------------- ##

exsnap_done=False

if EII :
 exec(open("DDD_EII.py").read())

## -------------------- Etape III -------------------- ##

from PO23D import *
#rempUsnap='par8'#'seq'

if EIII :
 exec(open("DDD_EIII.py").read())

N_mor=2#99,9%
N_mor=4#99,99%
N_mor=6#99,999%
N_mor=8#100%

## -------------------- Etape IV -------------------- ##

N_mor=4
r_nouv=0.33#0.22#0.33#0.44

# La mesure du temps d'éxécution doit se faire avec l'option 'save' de fig_todo

ind_res=True#False###----------> on précise la résolution du maillage, qui apparaît ou non dans le fichier contenant Phi

#0.22 : 0.011% d'erreur, tps d'éxécution ~ 1''/70''
#0.33 : 0.012% d'erreur, tps d'éxécution ~ 1''/60''

if EIV :
 exec(open("DDD_EIV.py").read())




















