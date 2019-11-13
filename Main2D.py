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

# Parametres

from DD_par import *


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

### ------------ Execution des etapes demandees en preambule, imports specifiques ------------ ###


from LEc import *




## ------------ Etape lL 1demi : Affichage de microstructures periodiques ------------ ##



if E_lL :
    exec(open("DD_ElL.py").read())

## ---------- Etape I, memes parametres que pour 1demi ---------- ##

# Execution

if EI :
    exec(open("DD_EI.py").read())

## ---------- Etape II ---------- ##

if EII :
    exec(open("DD_EII.py").read())

if exsnap_done and test_Dhom :
    exec(open("DD_EIIintgrad"))

## ---------- Etape III ---------- ##

from PO23D import *



if EIII :
    exec(open("DD_EIII.py").read())



## ---------- Etape IV ---------- ##



# La mesure du temps d'execution doit se faire avec l'option 'save' de fig_todo

if EIV :
    exec(open("DD_EIV.py").read())

# if EIVfixe:
#  exec(open("DD_EIV_fixe.py").read())
