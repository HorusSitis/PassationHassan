### ------------ Paquets à importer ------------ ###

import subprocess

import numpy as np
import random as rd
import pylab as pl
from pylab import *
#import matplotlib.pylab as pl
import os

from importlib import reload

#threading

import threading
import time
import multiprocessing

#multiprocessing

import multiprocessing
#from joblib import Parallel, delayed

#stockage d'objets python

import marshal as ma
import shelve as sh

## Paquets spécifiques à POD-MOR ##

from fenics import *
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from math import sqrt
from math import exp
import sys

# Fonctions 2D

from DD_fun_obj import *
#import DD_fun_obj as F2d
#F2d=reload(F2d)

##########################################################
### ------------ Code à lire : conditions ------------ ###
##########################################################

E_=False
E_lL=True

EI=False
EII=False

EIII=False
EIV=False

res_fixe=20
res=20

fixe_aff=False
fig_todo='save'

### ------------ Etape 0 : Génération de microstructures périodiques aléatoires ------------ ###

if E_ :
 from DD_E0 import *

#########################################################
### ----------------- Etapes I à IV ----------------- ###
#########################################################

### ------------ Exécution des étapes demandées en préambule, imports spécifiques ------------ ###

repertoire_parent="Res2D/"
from LEc import *

D_k=1.0
Nsnap=8
npas_err=20

# Choix du paramètre géométrique

## Cas d'un disque unique
#config='cercle unique'

#geo_p='rayon'
cen_snap_ray=[0.5,0.5]#[0.,0.]#

#geo_p='centre'
ray_snap_cen=0.25
#csr_list=[[0.05*k,0.5] for k in range(1,1+Nsnap)]


## Cas de deux inclusions périodiques : une inclusion centrale et deux ou quatre latérales par cellule
config='compl'
## le rayon du disque central est variable
#geo_p='deuxième disque aux sommets'
geo_p='deuxième disque latéral'
#

## ------------ Etape lL 1demi : Affichage de midrostructures périodiques ------------ ##

nb_lcells=1
cem_color='grey'#'r'
sand_color='orange'
fluid_color='cyan'

if E_lL :
 exec(open("DD_ElL.py").read())

## ---------- Etape I, mêmes paramètres que pour 1demi ---------- ##

if EI :
 exec(open("DD_EI.py").read())

## ---------- Etape II ---------- ##

if EII :
 exec(open("DD_EII.py").read())

## ---------- Etape III ---------- ##

from PO23D import *
rempUsnap='par8'#'seq'

if EIII :
 exec(open("DD_EIII.py").read())


## ---------- Etape IV ---------- ##















