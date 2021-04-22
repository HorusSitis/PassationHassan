# -*- coding: utf-8 -*-
### Une commande possible dans le terminal ###

#--- mpirun -np 8 python3 Main3D.py ---#
#--- affiche npfois 'pas encore' fait avec l'etape IV ---#

# Attention : on execute parallelement

### Paquets a importer ###

from fenics import *
from dolfin import *
### --- from mshr import * --- ###
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, log

import sys
import os

# Calcul parallele

import multiprocessing

# Performances

import time

# Stockage d'objets python

import marshal as ma
import shelve as sh

# affichage etc

import pylab as pl

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
#
## Paquets specifiques a la 3d ##

from mpl_toolkits.mplot3d.axes3d import get_test_data
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Circle, PathPatch

# nom de l'appareil utilise pour generer les donnees enregistre_pges
computer='MECALAC_29x8'#'T1700_35x8'#

# Parametres geometriques : topologie et longueurs

# from DDD_geoset import *
exec(open('DDD_geoset.py', encoding='utf-8').read())

# Fonctions maison

exec(open('DDD_fun_obj.py', encoding='utf-8').read())

##########################################################
### ------------ Code a lire : conditions ------------ ###
##########################################################

E_ = False
E_lL = True

EI = False
mesh_appr_done = True
snap_done = False
err_per_calc = False

mesh_ex_done = False
EII = False
exsnap_done = False
test_Dhom = True

EIII = False

EIV = False
Report = True

EV = False





## ------------ Etape lL 1demi : Affichage de microstructures periodiques ------------ ##

if E_lL :
    exec(open("DDD_ElL.py").read())

## -------------------- Etape I -------------------- ##

# Execution

if EI :
    exec(open("DDD_EI.py").read())

##computer=
### Pour les etapes qui suivent, on peut choisir l'ordinateur qui a effectue le calcul des snapshots physiques

## -------------------- Etape II -------------------- ##

if EII :
    exec(open("DDD_EII.py").read())

## -------------------- Etape III -------------------- ##

# from PO23D import *
exec(open("PO23D.py").read())

#rempUsnap='par8'#'seq'
# seuil_ener = 99.99

if EIII :
    exec(open("DDD_EIII.py").read())

## -------------------- Etape IV -------------------- ##

# registre_N_mor_name = 'Perf3D/' + 'N_mor_' + 'ener_nu10E' + expo + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)

if EIV:

    exec(open('tools_chi.py', encoding='utf-8').read())

    tab_N_mor = np.load(registre_N_mor_name + '.npy')
    N_mor = tab_N_mor[0]

    # tableaux contenant les performances du ROM
    arr_int_grad_fem = np.zeros(len(list_rho_test))
    arr_int_grad_rom = np.zeros(len(list_rho_test))

    if config == 'cylsph':
        arr_int_grad_yy_fem = np.zeros(len(list_rho_test))
        arr_int_grad_yy_rom = np.zeros(len(list_rho_test))

    arr_nodes = np.zeros(len(list_rho_test))

    arr_err_rel = np.zeros(len(list_rho_test))
    arr_var_rel = np.zeros(len(list_rho_test))
    arr_var_rel_yy = np.zeros(len(list_rho_test))
    arr_var_rel_chi = np.zeros(len(list_rho_test))

    ## pour l'instant : on en compte pas le temps de maillage, identique pour FEM et ROM

    arr_t = np.zeros((len(list_rho_test), 6))

    # fichiers pour enregistrer textuellement les performances du ROM
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
    # registre_pg.write('\\'+'('+'\\'+'mathcal{G}^{rom}_{-'+'\\'+'phi}'+'\\'+')'+'&')
    # registre_pg.write('\\'+'('+'\\'+'mathcal{G}^{rom}_{solve}'+'\\'+')'+'\\'+'\\'+'\n')

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

    for i in range(len(list_rho_test)):

        rho = list_rho_test[i]
        r_nouv=rho
        exec(open("DDD_EIV.py").read())


    ## fin des deux tableaux
    registre_pg.write('\\'+'end{tabular}')
    registre_gr.write('\\'+'end{tabular}')

    registre_pg.close()
    registre_gr.close()

    nom_tab_latex_pg = '../GitLab/rom_diffeo_dhom/latex_article/' + 'pg_' + config + '_res' + str(res_gmsh) + '_raydeb_o' + str(int(100*2*list_rho_appr[0]))
    nom_tab_latex_gr = '../GitLab/rom_diffeo_dhom/latex_article/' + 'gr_' + config + '_res' + str(res_gmsh) + '_raydeb_o' + str(int(100*2*list_rho_appr[0]))

    os.rename(nom_fichier_pg + '.txt', nom_tab_latex_pg + '.tex')
    os.rename(nom_fichier_gr + '.txt', nom_tab_latex_gr + '.tex')

    # sauvegarde des performances
    np.save(registre_perf_num['int_grad_fem'] + '.npy', arr_int_grad_fem)
    np.save(registre_perf_num['int_grad_rom'] + '.npy', arr_int_grad_rom)
    if config == 'cylsph':
        np.save(registre_perf_num['int_grad_yy_fem'] + '.npy', arr_int_grad_yy_fem)
        np.save(registre_perf_num['int_grad_yy_rom'] + '.npy', arr_int_grad_yy_rom)

    np.save(registre_perf_num['nodes'] + '.npy', arr_nodes)
    np.save(registre_perf_num['err_rel'] + '.npy', arr_err_rel)
    np.save(registre_perf_num['var_rel'] + '.npy', arr_var_rel)
    if config == 'cylsph':
        np.save(registre_perf_num['var_rel_yy'] + '.npy', arr_var_rel_yy)
    np.save(registre_perf_num['t'] + '.npy', arr_t)

    np.save(registre_perf_num['var_rel_chi'] + '.npy', arr_var_rel_chi)

## -------------------- Etape V -------------------- ##

figures_repository = 'Figures3D/'
figures_repository = '../GitLab/rom_diffeo_dhom/figures_interp/'


if EV :
    exec(open("DDD_EV.py").read())
