### ------------ Paquets a importer ------------ ###

# paquets mathematiques
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

# systeme

import os
import sys

from importlib import reload

# calcul parallele

import threading
import time
import multiprocessing

import subprocess

# stockage d'objets python

import marshal as ma
import shelve as sh

## Paquets specifiques a POD-MOR ##

from fenics import *
from dolfin import *
from mshr import *

# Fonctions 2D

from DD_fun_obj import *
#import DD_fun_obj as F2d
#F2d=reload(F2d)

##########################################################
### ------------ Code a lire : conditions ------------ ###
##########################################################

E_=False
E_lL=False

EI=False
snap_done=True
EII=False
exsnap_done=True
test_Dhom=True


EIII=False
EIV=True
EIVfixe=False



res_gmsh=100

fixe_aff=False
fig_todo=''

### ------------ Etape 0 : Generation de microstructures periodiques aleatoires ------------ ###

if E_ :
 from DD_E0 import *

#########################################################
### ----------------- Etapes I a IV ----------------- ###
#########################################################

### ------------------ Important : degre pour la resolution par elements finis ------------------ ###
VFS_degree=2#3#
## degre 2 : comme en dimension 3, permet d'eviter les erreurs de periodicite pour des pas qui nen sont pas de la forme 2"n, ou n est un diviseur de 100 ##

config='cer_un'#'compl'#

if config=='cer_un':
 test_snap='i_per'#'solid_1'#''#
 ### on exclut 'solid_2' ###
 dom_fixe="am"#"multiray"#"ray_min"#
 ##
 geo_p='ray'#'cen'#
 cen_snap_ray=[0.5,0.5]#[0.,0.]#
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
  conf_mess=conf_mess+" centre"
  mention=""
  ### 'i_per' : Nrom=... ; 'solid_1': Nrom=... ; ### ???
 elif cen_snap_ray==[0.,0.]:
  conf_mess=conf_mess+" aux sommets"
  mention="_som"
  config=config+mention
  ### 'i_per' : Nrom=2 ; 'solid_1': Nrom=5 ; 'solid_2': Nrom=2 ---> ne fonctionne pas
elif config=='compl':
 test_snap='solid_1'#'solid_2'#''#
 ##
 dom_fixe="solid"#"am"#
 ##
 geo_p='hor'#'diag'#
 ##
 if geo_p=='diag':
  cen_snap_ray=[0.,0.]
  ### 'solid_1': Nrom=3 ;
 elif geo_p=='hor':
  cen_snap_ray=[0.,0.5]
 ##
 conf_mess='deux disques par periode'
 mess_prefix=' rayon central variable'
 mention=''
 if geo_p=='diag':
  geo_mess='alignes en diagonale, '+mess_prefix
 elif geo_p=='hor':
  geo_mess='alignes horizontalement, '+mess_prefix

# choix du type de maillage

typ_msh='gms'#''#
#print(typ_msh)
if typ_msh=='gms':
 res=res_gmsh
 res_fixe=res_gmsh

# raffinement de maillages

Nrefine=1
#crow=(1/res_gmsh)*1e-1
typ_refi='front'#'vol'#
lg_crow=-1
crow=2*10**(lg_crow)

if typ_refi=='vol':
 refi_mess='Couronne : '+str(crow)
elif typ_refi=='front':
 refi_mess='Surface'

# nom de l'appareil utilise pour generer les donnees enregistrees
computer='MECALAC_29x8'##'T1700_35x8'##

# apprentissage : calcul parallele ou sequentiel, prise en compte de la resolution

gen_snap='par8'#'seq'#'seq_par'#

# repertoire pour les resultats

repertoire_parent="Res2D/"

### ------------ Execution des etapes demandees en preambule, imports specifiques ------------ ###


from LEc import *

D_k=1.0
Nsnap=8
deb=1# par defaut##5# pour les tests cer_un et diag##6# pour hor##
npas_err=50
ordo='Ordr'#'Nordr'


## ------------ Etape lL 1demi : Affichage de microstructures periodiques ------------ ##

nb_lcells=1
cem_color='grey'#'r'
sand_color='orange'
fluid_color='cyan'

if E_lL :
 exec(open("DD_ElL.py").read())

## ---------- Etape I, memes parametres que pour 1demi ---------- ##

# Execution

if EI :
 exec(open("DD_EI.py").read())

## ---------- Etape II ---------- ##

if EII :
 exec(open("DD_EII.py").read())

## ---------- Etape III ---------- ##

from PO23D import *

# On soustrait eventuellement la moyenne des snapshots interpoles a chaque snapshot
##moy_mod='soust_m'#'' par defaut

if EIII :
 exec(open("DD_EIII.py").read())

##99,99% :
### compl diag N_mor=3 pour "solid", 4 pour "am"
### compl hor N_mor=4

## ---------- Etape IV ---------- ##

N_mor=2
r_nouv=0.33#0.44#0.22#0.11#0.10#

# La mesure du temps d'execution doit se faire avec l'option 'save' de fig_todo

if EIV :
 exec(open("DD_EIV.py").read())

if EIVfixe:
 exec(open("DD_EIV_fixe.py").read())
