### ------------ Paquets à importer ------------ ###

# paquets mathématiques
import numpy as np
import random as rd
from math import sqrt
from math import exp

# affichage etc

import pylab as pl
#from pylab import *

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

# système

import os
import sys

from importlib import reload

# calcul parallèle

import threading
import time
import multiprocessing

import subprocess

# stockage d'objets python

import marshal as ma
import shelve as sh

## Paquets spécifiques à POD-MOR ##

from fenics import *
from dolfin import *
from mshr import *

# Fonctions 2D

from DD_fun_obj import *
#import DD_fun_obj as F2d
#F2d=reload(F2d)

##########################################################
### ------------ Code à lire : conditions ------------ ###
##########################################################

E_=False
E_lL=False

EI=False
snap_done=True
EII=False
exsnap_done=True

EIII=False
EIV=True

#res_fixe=20
#res=100#pour les sommets#20
res_gmsh=100

fixe_aff=False
fig_todo='save'

### ------------ Etape 0 : Génération de microstructures périodiques aléatoires ------------ ###

if E_ :
 from DD_E0 import *

#########################################################
### ----------------- Etapes I à IV ----------------- ###
#########################################################

# choix du type de maillage

typ_msh='gms'#''
print(typ_msh)
if typ_msh=='gms':
 res=res_gmsh
 res_fixe=res_gmsh

# nom de l'appareil utilisé pour générer les données enregistrées
computer='MECALAC_29x8'##'T1700_35x8'##



repertoire_parent="Res2D/"

### ------------ Exécution des étapes demandées en préambule, imports spécifiques ------------ ###


from LEc import *

D_k=1.0
Nsnap=8
npas_err=50
ordo='Ordr'#'Nordr'

gen_snap='par8'#'seq'

# Choix du paramètre géométrique

## Cas d'un disque unique
config='cer_un'
if config=='cer_un':
 conf_mess='disque unique'

geo_p='ray'
if geo_p=='ray':
 geo_mess='rayon variable'

#emplacement du disque unique : ici, centré ou aux sommets
cen_snap_ray=[0.,0.]#[0.5,0.5]#

### Disque centré ou aux sommets

if cen_snap_ray==[0.5,0.5]:
 conf_mess=conf_mess+" centré"
 mention=""
elif cen_snap_ray==[0.,0.]:
 conf_mess=conf_mess+" aux sommets"
 mention="_som"
 config=config+mention

### geo_p='centre'
#ray_snap_cen=0.25
#csr_list=[[0.05*k,0.5] for k in range(1,1+Nsnap)]


## Cas de deux inclusions périodiques : une inclusion centrale et deux ou quatre latérales par cellule
#config='compl'
## le rayon du disque central est variable
#geo_p='deuxième disque aux sommets'
#geo_p='deuxième disque latéral'
#

## ------------ Etape lL 1demi : Affichage de microstructures périodiques ------------ ##

nb_lcells=1
cem_color='grey'#'r'
sand_color='orange'
fluid_color='cyan'

if E_lL :
 exec(open("DD_ElL.py").read())

## ---------- Etape I, mêmes paramètres que pour 1demi ---------- ##

# Exécution

if EI :
 exec(open("DD_EI.py").read())

## ---------- Etape II ---------- ##

if EII :
 exec(open("DD_EII.py").read())

## ---------- Etape III ---------- ##

from PO23D import *
rempUsnap='par8'#'seq'

# On soustrait éventuellement la moyenne des snapshots interpolés à chaque snapshot
moy_mod='soust_m'#'' par défaut

if EIII :
 exec(open("DD_EIII.py").read())

## ---------- Etape IV ---------- ##

N_mor=5
r_nouv=0.22#0.44#0.33#

# La mesure du temps d'éxécution doit se faire avec l'option 'save' de fig_todo

if EIV :
 exec(open("DD_EIV.py").read())













