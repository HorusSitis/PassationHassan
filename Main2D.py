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
from joblib import Parallel, delayed

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
EI=True
EII=False
EIII=False
EIV=False

res_fixe=40
res=40

#c_x=0.5
#c_y=0.5
#r=0.35

fixe_aff=True
fig_todo='aff'

### ------------ Etape 0 : Génération de microstructures périodiques aléatoires ------------ ###

if E_ :
 from DD_E0 import *

#########################################################
### ----------------- Etapes I à IV ----------------- ###
#########################################################

### ------------ Exécution des étapes demandées en préambule, imports spécifiques ------------ ###

## ---------- Etape I ---------- ##

repertoire_parent="Res2D/"
from LEc import *

D_k=1.0
Nsnap=8
rempUsnap='par8'#'seq'
npas_err=20

if EI :
 exec(open("DD_EI.py").read())

## ---------- Etape II ---------- ##

if EII :
 exec(open("DD_EII.py").read())

## ---------- Etape III ---------- ##

from PO23D import *

if EIII :
 exec(open("DD_EIII.py").read())


## ---------- Etape IV ---------- ##















