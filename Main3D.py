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

EI=True
snap_done=False
EII=False
exsnap_done=False

EIII=False
EIV=False



# nom de l'appareil utilisé pour générer les données enregistrées
computer='T1700_35x8'#'MECALAC_29x8'#

# paramètres pour l'éxécution des étapes : affichage, tests de périodicité etc


fig_todo='save'
typ_msh='gms'#''
D_k=1.0
Nsnap=8
npas_err=20

typ_sol="bic_cyr"#"default"#seulement si res=10##
#moy_null=True

ordo='Ordr'#'Nordr'
# apprentissage : calcul parallèle ou séquentiel, prise en compte de la résolution

gen_snap='seq'#'par8'#''seq_par'#

# répertoire pour les résultats

repertoire_parent="Res3D/"


## -------------------- Etape I -------------------- ##

# Parallélisation du calcul des snapshots

parallelize=True

# Choix de la résolution du maillage : nombre de noeuds par côté du cube

res_gmsh=50

if typ_msh=='gms':
 res=res_gmsh

config='sph_un'#'cylsph'#'2sph'#'cyl_un'#

### inclusions simples
if config=='sph_un':
 dom_fixe="am"
 geo_p='ray'#'cen'#
 ##
 conf_mess='sphère unique'
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
### inclusions composées
elif config=='2sph':
 conf_mess='deux sphères'
 dom_fixe="solid"#"am"#
 geo_p='ray'
 geo_mess='rayon de la sphère centrale variable'
 ##
elif config=='cylsph':
 conf_mess='un cylindre et une sphère'
 dom_fixe="ray_min"#"am"#"solid"#
 geo_p='ray_sph'#'ray_cyl'#'ray_linked'###ray_sph pour le test diff_ion##
 ##
 if geo_p=='ray_cyl':
  geo_mess='rayon du cylindre variable'
  fixe_comp='cyl_sph'##utilisation du domaine fixe avec annulation du rayon du cylindre dans le fichier général##'sph_un'#'ray_min'#
 elif geo_p=='ray_sph':
  geo_mess='rayon de la sphère variable'
 elif geo_p=='ray_linked':
  geo_mess='rayons liés'

## ------------ Etape lL 1demi : Affichage de microstructures périodiques ------------ ##

nb_lcells=5
cem_color='grey'#'r'
sand_color='orange'
fluid_color='cyan'

if E_lL :
 exec(open("DDD_ElL.py").read())

# Exécution

if EI :
 exec(open("DDD_EI.py").read())

##computer=
### Pour les étapes qui suivent, on peut choisir l'ordinateur qui a effectué le calcul des snapshots physiques

## -------------------- Etape II -------------------- ##

if EII :
 exec(open("DDD_EII.py").read())

## -------------------- Etape III -------------------- ##

from PO23D import *
#rempUsnap='par8'#'seq'

if EIII :
 exec(open("DDD_EIII.py").read())

## Sphère ou cylindre unique
#N_mor=4#99,9%
#N_mor=5#99,99%

## -------------------- Etape IV -------------------- ##

N_mor=3
moy_null=False
r_nouv=0.44#0.33#0.22#0.11#

# La mesure du temps d'éxécution doit se faire avec l'option 'save' de fig_todo

ind_fixe=True##-----------> dom_fixe devant le 'Phi'
ind_res=True#False###----------> on précise la résolution du maillage, qui apparaît ou non dans le fichier contenant Phi

#0.22 : 0.011% d'erreur, tps d'éxécution ~ 1''/70''
#0.33 : 0.012% d'erreur, tps d'éxécution ~ 1''/60''

if EIV :
 exec(open("DDD_EIV.py").read())




















