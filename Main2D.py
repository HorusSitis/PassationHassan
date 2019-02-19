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
exsnap_done=False

EIII=False
EIV=True



#res_fixe=20
#res=100#pour les sommets#20
res_gmsh=100

fixe_aff=False
fig_todo=''

### ------------ Etape 0 : Génération de microstructures périodiques aléatoires ------------ ###

if E_ :
 from DD_E0 import *

#########################################################
### ----------------- Etapes I à IV ----------------- ###
#########################################################

### ------------------ Important : degré pour la résolution par éléments finis ------------------ ###
VFS_degree=2#3#
## degré 2 : comme en dimension 3, permet d'éviter les erreurs de périodicité pour des pas qui nen sont pas de la forme 2"n, où n est un diviseur de 100 ##

config='compl'#'cer_un'#

if config=='cer_un':
 test_snap='test'#'testbis'#''#'ntest'#
 dom_fixe="am"
 ##
 geo_p='ray'#'cen'#
 cen_snap_ray=[0.,0.]#[0.5,0.5]#
 ##
 conf_mess='disque unique'
 ##
 if geo_p=='ray':
  geo_mess='rayon variable'
 ### geo_p='centre'
 #ray_snap_cen=0.25
 #csr_list=[[0.05*k,0.5] for k in range(1,1+Nsnap)]
 ##
 if cen_snap_ray==[0.5,0.5]:
  conf_mess=conf_mess+" centré"
  mention=""
 elif cen_snap_ray==[0.,0.]:
  conf_mess=conf_mess+" aux sommets"
  mention="_som"
  config=config+mention
elif config=='compl':
 test_snap=''#'test'#
 ##
 dom_fixe="solid"#"am"#
 ##
 geo_p='diag'#'hor'#
 ##
 if geo_p=='diag':
  cen_snap_ray=[0.,0.]
 elif geo_p=='hor':
  cen_snap_ray=[0.,0.5]
 ##
 conf_mess='deux disques par période'
 mess_prefix=' rayon central variable'
 mention=''
 if geo_p=='diag':
  geo_mess='alignés en diagonale, '+mess_prefix
 elif geo_p=='hor':
  geo_mess='alignés horizontalement, '+mess_prefix

# choix du type de maillage

typ_msh='gms'#''
print(typ_msh)
if typ_msh=='gms':
 res=res_gmsh
 res_fixe=res_gmsh

# nom de l'appareil utilisé pour générer les données enregistrées
computer='T1700_35x8'##'MECALAC_29x8'##



repertoire_parent="Res2D/"

### ------------ Exécution des étapes demandées en préambule, imports spécifiques ------------ ###


from LEc import *

D_k=1.0
Nsnap=8
npas_err=50
ordo='Ordr'#'Nordr'

gen_snap='par8'#'seq'

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
##moy_mod='soust_m'#'' par défaut

if EIII :
 exec(open("DD_EIII.py").read())

##99,99% :
### compl diag N_mor=3 pour "solid", 4 pour "am"
### compl hor N_mor=4

## ---------- Etape IV ---------- ##

N_mor=4#cer_ _som ; 
#r_nouv=0.11#0.22#0.44#0.33#0.10#

# La mesure du temps d'éxécution doit se faire avec l'option 'save' de fig_todo

if EIV :
 ## initialisation : chargement des résultats des étapes précédentes
 phi_name=test_snap+'Phi'+dom_fixe+'_dim'+str(Nsnap)+'_'+config+'_'+geo_p+'_deg'+str(VFS_degree)+'_'+"res"+str(res)+'_'+ordo+'_'+computer
 with sh.open(repertoire_parent+phi_name) as phi_loa:
  Phi_prime_v = phi_loa["maliste"]
 ## Appel du script de l'étape IV
 ## boucle pour l'éxécution du modèle réduit sur le jeu de parapètres de test
 for r_nouv in [0.11,0.22,0.33,0.44]:
  exec(open("DD_EIV.py").read())













