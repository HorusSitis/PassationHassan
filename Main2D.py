# -*- coding: utf-8 -*-
### Une commande possible dans le terminal ###

#--- mpirun -np 8 python3 Main3D.py ---#
#--- affiche npfois 'pas encore' fait avec l'etape IV ---#

# Attention : on execute parallelement

### ------------ Paquets a importer ------------ ###

# paquets mathematiques
import numpy as np
import random as rd
from math import sqrt
from math import exp

# affichage etc

import pylab as pl

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
### --- from mshr import * --- ###

# Parametres

# from DD_pars import *
exec(open('DD_pars.py', encoding = 'utf8').read())

# Fonctions 2D

exec(open('DD_fun_obj.py', encoding='utf-8').read())

##########################################################
### ------------ Code a lire : conditions ------------ ###
##########################################################

E_ = False
E_lL = False

EI = False
mesh_appr_done = True
snap_done = True

err_eval = False

mesh_ex_done = False
EII = True
exsnap_done = False

test_Dhom = False

EIII = False

EIV = False
Report = True

# EIVfixe=False

EV = False

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
    exec(open('DD_ElL.py', encoding='utf-8').read())

## ---------- Etape I, memes parametres que pour 1demi ---------- ##

if EI :
    exec(open('DD_EI.py', encoding='utf-8').read())

## ---------- Etape II : int_grad eventuellement sur dom_fixe ---------- ##

if EII :
    exec(open('DD_EII.py', encoding='utf-8').read())

if exsnap_done and test_Dhom :
    exec(open('DD_EIIintgrad.py', encoding='utf-8').read())

## ---------- Etape III ---------- ##

# from PO23D import *
exec(open('PO23D.py', encoding='utf-8').read())

if EIII :
    exec(open('DD_EIII.py', encoding='utf-8').read())

if config == 'cer_un':
    if mention == '':
        N_mor = 2
elif config == 'compl':
    if geo_p == 'diag':
        N_mor = 3
    elif geo_p == 'hor':
        N_mor = 3

## ---------- Etape IV ---------- ##

# La mesure du temps d'execution doit se faire avec l'option 'save' de fig_todo

if EIV :

    # tableaux contenant les performances du ROM
    arr_int_grad_fem = np.zeros(len(list_rho_test))
    arr_int_grad_rom = np.zeros(len(list_rho_test))

    arr_nodes = np.zeros(len(list_rho_test))

    arr_err_rel = np.zeros(len(list_rho_test))
    arr_var_rel = np.zeros(len(list_rho_test))
    arr_var_rel_yy = np.zeros(len(list_rho_test))

    ## pour l'instant : on en compte pas le temps de maillage, identique pour FEM et ROM

    ## performances temporelles
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
        exec(open('DD_EIV.py', encoding='utf-8').read())


    ## fin des deux tableaux
    registre_pg.write('\\'+'end{tabular}')
    registre_gr.write('\\'+'end{tabular}')

    registre_pg.close()
    registre_gr.close()

    # sauvegarde des performances
    np.save(registre_perf_num['int_grad_fem'] + '.npy', arr_int_grad_fem)
    np.save(registre_perf_num['int_grad_rom'] + '.npy', arr_int_grad_rom)
    if config == 'compl':
        np.save(registre_perf_num['int_grad_yy_fem'] + '.npy', arr_int_grad_yy_fem)
        np.save(registre_perf_num['int_grad_yy_rom'] + '.npy', arr_int_grad_yy_rom)

    np.save(registre_perf_num['nodes'] + '.npy', arr_nodes)
    np.save(registre_perf_num['err_rel'] + '.npy', arr_err_rel)
    np.save(registre_perf_num['var_rel'] + '.npy', arr_var_rel)
    if config == 'compl':
        np.save(registre_perf_num['var_rel_yy'] + '.npy', arr_var_rel_yy)
    np.save(registre_perf_num['t'] + '.npy', arr_t)

    # nom_tab_latex_pg = '../GitLab/rom_diffeo_dhom/latex_article/' + 'pg_' + config + '_res' + str(res_gmsh) + '_raydeb_o' + str(int(100*2*list_rho_appr[0]))
    # nom_tab_latex_gr = '../GitLab/rom_diffeo_dhom/latex_article/' + 'gr_' + config + '_res' + str(res_gmsh) + '_raydeb_o' + str(int(100*2*list_rho_appr[0]))

    nom_tab_latex_pg = '../GitLab/rom_diffeo_dhom/latex_article/' + 'meshr_pg_' + config + '_res' + str(res_gmsh)# + '_raydeb_o' + str(int(100*2*list_rho_appr[0]))
    nom_tab_latex_gr = '../GitLab/rom_diffeo_dhom/latex_article/' + 'meshr_gr_' + config + '_res' + str(res_gmsh)# + '_raydeb_o' + str(int(100*2*list_rho_appr[0]))

    if config == 'compl':
        nom_tab_latex_pg = nom_tab_latex_pg + '_rayp' + str(int(round(100*ray_p,2)))
        nom_tab_latex_gr = nom_tab_latex_gr + '_rayp' + str(int(round(100*ray_p,2)))

    os.rename(nom_fichier_pg + '.txt', nom_tab_latex_pg + '.tex')
    os.rename(nom_fichier_gr + '.txt', nom_tab_latex_gr + '.tex')


## -------------------- Etape V -------------------- ##

figures_repository = 'Figures3D/'
figures_repository = '../GitLab/rom_diffeo_dhom/figures_interp/'


if EV :
    exec(open("DD_EV.py").read())
