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

import multiprocessing

from DDD_fun_obj import *

##########################################################
### ------------ Code à lire : conditions ------------ ###
##########################################################

E_=False
E_lL=True

EI=False
EII=False

EIII=False
EIV=False


res_fixe=6
fixe_aff=False
res=6
slices_cyl=5
Nsnap=8
rempUsnap='par8'#'seq'
c_x=0.5
c_y=0.5
c_z=0.5
#r=0.35#pour une réalisation unique
npas_err=20
fig_todo='aff'




## -------------------- Etape I -------------------- ##

D_k=1.0
Nsnap=8
npas_err=20

# Parallélisation du calcul des snapshots

parallelize=True

# Choix du paramètre géométrique

## Cas d'une sphère unique
config='sphère unique'

geo_p='rayon'
cen_snap_ray=[0.5,0.5,0.5]
#geo_p='centre'
ray_snap_cen=0.35
#csr_list=[[0.5,0.5,0.3+0.05*k] for k in range(1,1+Nsnap)]
### l'inclusion solide ne rencontre pas les bords du domaine, sous peine d'une incompatibilité entre les mailles de faces opposées (?)

## Cas d'un cylindre unique : axe parallèle à Oy
#config='cylindre unique'

#geo_p='rayon'
#geo_p='axe'
#csr_list=[[0.5,0.3+0.05*k] for k in range(1,1+Nsnap)]
### point auquel l'axe rencontre le plan Oxz. L'inclusion solide ne rencontre pas les bords du domaine.

## Cas d'un cylindre périodique aux arètes et une sphère unique au centre
#config='compl'
#geo_p=='rayon de la sphère variable'
#geo_p=='rayon du cylindre variable'

## ------------ Etape lL 1demi : Affichage de midrostructures périodiques ------------ ##

nb_lcells=1
cem_color='grey'#'r'
sand_color='orange'
fluid_color='cyan'

if E_lL :
 exec(open("DDD_ElL.py").read())

# Exécution

if EI :
 exec(open("DDD_EI.py").read())

## -------------------- Etape II -------------------- ##

if EII :
 exec(open("DDD_EII.py").read())

## -------------------- Etape III -------------------- ##

from PO23D import *
#rempUsnap='par8'#'seq'

if EIII :
 exec(open("DDD_EIII.py").read())


## -------------------- Etape IV -------------------- ##























