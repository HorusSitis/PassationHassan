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

from DD_pars import *

##########################################################
### ------------ Code a lire : conditions ------------ ###
##########################################################

E_=False
E_lL=False

EI=False
snap_done=True

EII=False
exsnap_done=True
test_Dhom=False

EIII=True

EIV=False
Report = True

# EIVfixe=False

# res_gmsh=100

### ------------ Etape 0 : Generation de microstructures periodiques aleatoires ------------ ###

if E_ :
    from DD_E0 import *

#########################################################
### ----------------- Etapes I a IV ----------------- ###
#########################################################

### ------------ Execution des etapes demandees en preambule, imports specifiques ------------ ###

from LEc import *

## ------------ Etape lL 1demi : Affichage de microstructures periodiques ------------ ##

if E_lL :
    exec(open("DD_ElL.py").read())

## ---------- Etape I, memes parametres que pour 1demi ---------- ##

if EI :
    exec(open("DD_EI.py").read())

## ---------- Etape II : int_grad eventuellement sur dom_fixe ---------- ##

if EII :
    exec(open("DD_EII.py").read())

if exsnap_done and test_Dhom :
    exec(open("DD_EIIintgrad.py").read())

## ---------- Etape III ---------- ##

from PO23D import *

if EIII :
    exec(open("DD_EIII.py").read())

## ---------- Etape IV ---------- ##

# La mesure du temps d'execution doit se faire avec l'option 'save' de fig_todo

if EIV :

    nom_fichier_pg='Perf3D/' + 'pg_' + computer + 'res' + str(res_gmsh) + config + geo_p + 'Nmor' + str(N_mor)
    nom_fichier_gr='Perf3D/' + 'gr_' + computer + 'res' + str(res_gmsh) + config + geo_p + 'Nmor' + str(N_mor)

    registre_pg=open(nom_fichier_pg + '.txt','w')
    registre_gr=open(nom_fichier_gr + '.txt','w')

    ## en tete du tableau de resultats et performances
    registre_pg.write('\\'+'begin{tabular}')
    registre_pg.write('{|c|c||c|c|c||c|}')
    registre_pg.write('\n')
    registre_pg.write('\\'+'hline'+'\n')

    registre_pg.write('\\'+'rowcolor{'+'lightgray'+'}')
    registre_pg.write('\\'+'('+'\\'+'tilde{'+'\\'+'rho'+'}'+'\\'+')'+'&')
    registre_pg.write('Nodes'+'&')

    registre_pg.write('\\'+'('+'\\'+'frac{'+'\\'+'int'+'\\'+'nabla'+'\\'+'chi_{rom}}{|'+'\\'+'Omega|}'+'\\'+')'+'&')
    registre_pg.write('\\'+'('+'\\'+'frac{'+'\\'+'int'+'\\'+'nabla'+'\\'+'chi_{fem}}{|'+'\\'+'Omega|}'+'\\'+')'+'&')
    registre_pg.write('\\'+'('+'Err'+'\\'+')'+'&')

    registre_pg.write('\\'+'('+'\\'+'mathcal{G}^{rom}'+'\\'+')'+'\\'+'\\'+'\n')

    registre_pg.write('\\'+'hline'+'\n')

    ## en tete du tableau de performances temporelles relatives
    registre_gr.write('\\'+'begin{tabular}')
    registre_gr.write('{|c|c||c|c|c|c|}')
    registre_gr.write('\n')
    registre_gr.write('\\'+'hline'+'\n')

    registre_gr.write('\\'+'rowcolor{'+'lightgray'+'}')
    registre_gr.write('\\'+'('+'\\'+'tilde{'+'\\'+'rho'+'}'+'\\'+')'+'&')
    registre_gr.write('Nodes'+'&')

    registre_gr.write('\\'+'('+'t_{'+'\\'+'phi^{nouv}}/t_{ROM}'+'\\'+')'+'&')
    registre_gr.write('\\'+'('+'t_{Ab}/t_{ROM}'+'\\'+')'+'&')
    registre_gr.write('\\'+'('+'t_{solve}/t_{ROM}'+'\\'+')'+'&')
    registre_gr.write('\\'+'('+'t_{D^{hom}}/t_{ROM}'+'\\'+')'+'\\'+'\\'+'\n')
    # registre_gr.write('\\'+'('+'t_{fem}'+'\\'+')'+'&')

    registre_gr.write('\\'+'hline'+'\n')

    for rho in list_rho_test:
    # for rho in [0.33]:

        r_nouv=rho
        exec(open("DD_EIV.py").read())


    ## fin des deux tableaux
    registre_pg.write('\\'+'end{tabular}')
    registre_gr.write('\\'+'end{tabular}')

    registre_pg.close()
    registre_gr.close()

    nom_tab_latex_pg = '../GitLab/rom_diffeo_dhom/latex_article/' + 'pg_' + config + '_res' + str(res_gmsh) + '_raydeb_o' + str(int(100*2*list_rho_appr[0]))
    nom_tab_latex_gr = '../GitLab/rom_diffeo_dhom/latex_article/' + 'gr_' + config + '_res' + str(res_gmsh) + '_raydeb_o' + str(int(100*2*list_rho_appr[0]))

    os.rename(nom_fichier_pg + '.txt', nom_tab_latex_pg + '.tex')
    os.rename(nom_fichier_gr + '.txt', nom_tab_latex_gr + '.tex')

# if EIVfixe:
#  exec(open("DD_EIV_fixe.py").read())
