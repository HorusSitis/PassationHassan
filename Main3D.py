# -*- coding: utf-8 -*-
### Une commande possible dans le terminal ###

#--- mpirun -np 8 python3 Main3D.py ---#
#--- affiche npfois 'pas encore' fait avec l'étape IV ---#

# Attention : on éxécute parallèlement 


##########################################################
### ------------ Code à lire : conditions ------------ ###
##########################################################

E_=False
EI=True
EII=False
EIII=False
EIV=False


res_fixe=6
fixe_aff=False
res=6
Nsnap=8
rempUsnap='par8'#'seq'
c_x=0.5
c_y=0.5
c_z=0.5
#r=0.35#pour une réalisation unique
npas_err=20
fig_todo='aff'

### Répertoire courant ###

#cd /home/amorea12/Documents/T_LaSIE/PassationHassan

### Paquets à importer ###

##############################################################################################################################
############################### Calculs avec la POD, modèles réduits ; merci à Hassan GHRAIEB. ###############################
##############################################################################################################################

#Paquets spécifiques à POD-MOR

from fenics import *
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
import sys

import multiprocessing

from DDD_fun_obj import *
#import DDD_fun_obj as F3d
#from importlib import reload
#F3d=reload(F3d)

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################


## ---------- Etape I ---------- ##

repertoire_parent="Res2D/"
from LEc import *

D_k=1.0
Nsnap=8
#snapshots='par8'#'seq'
npas_err=20

# Parallélisation du calcul des snapshots

parallelize=True

# Choix du paramètre géométrique

#geo_p='rayon'
csr_list=[[0.5,0.5,0.3+0.05*k] for k in range(1,1+Nsnap)]## l'inclusion solide ne rencontre pas les bords du domaine, sous peine d'une incompatibilité entre les mailles de faces opposées (?)
geo_p='centre'

# Exécution

if EI :
 exec(open("DDD_EI.py").read())




## ---------- Etape II ---------- ##

if EII :
 exec(open("DDD_EII.py").read())

## ---------- Etape III ---------- ##

from PO23D import *
#rempUsnap='par8'#'seq'

if EIII :
 exec(open("DDD_EIII.py").read())


## ---------- Etape IV ---------- ##























