# -*- coding: utf-8 -*-
### Une commande possible dans le terminal ###

#--- mpirun -np 8 python3 Main3D.py ---#
#--- affiche npfois 'pas encore' fait avec l'etape IV ---#

# Attention : on execute parallelement

### Paquets a importer ###

from fenics import *
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt

import sys
import os

# Calcul parallele

import multiprocessing

# Performances

import time

# Stockage d'objets python

import marshal as ma
import shelve as sh

## Paquets specifiques a la 3d ##

from mpl_toolkits.mplot3d.axes3d import get_test_data
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Circle, PathPatch

# Parametres geometriques : topologie et longueurs

from DDD_geoset import *

# Fonctions maison

exec(open('DDD_fun_obj.py', encoding='utf-8').read())

##########################################################
### ------------ Code a lire : conditions ------------ ###
##########################################################

E_=False
E_lL=False

EI=True
mesh_appr_done = True
snap_done=True

mesh_ex_done = False
EII=False
exsnap_done=False
test_Dhom=False

EIII=False

EIV=False
Interpolation=True
Report=True


# nom de l'appareil utilise pour generer les donnees enregistre_pges
computer='MECALAC_29x8'#'T1700_35x8'#


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

from PO23D import *
#rempUsnap='par8'#'seq'

if EIII :
    exec(open("DDD_EIII.py").read())

## -------------------- Etape IV -------------------- ##

if EIV and res_gmsh!=50:

    if Interpolation:

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

        for rho in list_rho_test:

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

    else:
        exec(open("DDD_EIV_fixe.py").read())

# # -------------- res = 50 : conditions pour l'execution du modele reduit -------------- #
#
#
# l_conf_g=[('2sph','ray'),('cylsph','ray_sph'),('cylsph','ray_cyl')]
#
# dico_conf_g={}
# dico_conf_g[('2sph','ray')]=([20.40,20.88,20.40,17.73],2,[0.11,0.22,0.33,0.44],4)# dernier parametre : nombre de calculs deja faits
# dico_conf_g[('cylsph','ray_sph')]=([18.54,18.79,18.47,16.60],2,[0.11,0.22,0.33,0.44],2)
# dico_conf_g[('cylsph','ray_cyl')]=([19.78,18.08,14.16,9.86],4,[0.11,0.22,0.33,0.44],4)
#
# if EIV and res_gmsh==50:
#  for conf_g in l_conf_g:
#   config=conf_g[0]
#   geo_p=conf_g[1]
#   #
#   ### inclusions composees
#   if config=='2sph':
#    conf_mess='deux spheres'
#    dom_fixe="solid"
#    geo_p='ray'
#    geo_mess='rayon de la sphere centrale variable'
#    ##config + geo_p + 'Nmor' + str(N_mor)
#   if config=='cylsph':
#    conf_mess='un cylindre et une sphere'
#    dom_fixe="ray_min"
#    ##
#    if geo_p=='ray_cyl':
#     geo_mess='rayon du cylindre variable'
#    elif geo_p=='ray_sph':
#     geo_mess='rayon de la sphere variable'
#   ###
#   l_tm=dico_conf_g[conf_g][0]
#   N_mor=dico_conf_g[conf_g][1]
#   l_rho=dico_conf_g[conf_g][2]
#   deb=dico_conf_g[conf_g][3]
#   ##
#   nom_fichier='Perf3D/'+computer+'_res'+str(res_gmsh)+config+geo_p+'Nmor'+str(N_mor)+'.txt'
#   registre_pg=open(nom_fichier,'w')
#   #
#   for i in range(deb,len(l_rho)):
#     print('rayon',str(l_rho[i]),'maillage :',l_tm[i],'secondes')
#     print('##############################################################################')
#     #
#     r_nouv=l_rho[i]
#     t_meshing=l_tm[i]
#     with open("DDD_EIV.py",'r') as op:
#      exec(op.read())
#   registre_pg.close()
