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

# Fonctions maison

from DDD_fun_obj import *

## Paquets specifiques a la 3d ##

from mpl_toolkits.mplot3d.axes3d import get_test_data
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Circle, PathPatch

# Parametres geometriques : topologie et longueurs

from DDD_geoset import *

##########################################################
### ------------ Code a lire : conditions ------------ ###
##########################################################

E_=False
E_lL=False

EI=False
snap_done=False

EII=False
exsnap_done=False
test_Dhom=False
EIII=False

EIV=True
Interpolation=True
Report=True


# nom de l'appareil utilise pour generer les donnees enregistrees
computer='MECALAC_29x8'#'T1700_35x8'#

# Choix de la resolution du maillage : nombre de noeuds par cote du cube

res_gmsh=10
if typ_msh=='gms':
    res=res_gmsh

# parametres pour l'execution des etapes : affichage, tests de periodicite etc

fig_todo='save'
typ_msh='gms'#''
D_k=1.0

Nsnap=len(list_rho_appr)

npas_err=20
typ_sol="bic_cyr"#"default"#seulement si res=10##
ordo='Ordr'#'Nordr'

Nrefine=1
crow=(1/res_gmsh)*1e-1
typ_refi='vol'#'front'#

# apprentissage : calcul parallele ou sequentiel, prise en compte de la resolution

gen_snap='par8'
# gen_snap='par4'
# gen_snap='seq'
# gen_snap='seq_par'

# repertoire pour les resultats
repertoire_parent="Res3D/"



## -------------------- Etape I -------------------- ##

## ------------ Etape lL 1demi : Affichage de microstructures periodiques ------------ ##

nb_lcells=5
cem_color='grey'
sand_color='orange'
fluid_color='cyan'

if E_lL :
    exec(open("DDD_ElL.py").read())

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

N_mor=3
t_meshing=1.89
r_nouv=0.33#0.44#0.22#0.11#

# La mesure du temps d'execution doit se faire avec l'option 'save' de fig_todo

ind_fixe=True##-----------> dom_fixe devant le 'Phi'
ind_res=True#False###----------> on precise la resolution du maillage, qui apparaît ou non dans le fichier contenant Phi

if EIV and res_gmsh!=50:

    if Interpolation:

        nom_fichier='Perf3D/' + computer + 'res' + str(res_gmsh) + config + geo_p + 'Nmor' + str(N_mor)
        registre=open(nom_fichier + '.txt','w')

        registre.write('\\'+'begin{tabular}')
        registre.write('{|c|c||c|c|c||c|c|c||c|c|}')
        registre.write('\n')
        registre.write('\\'+'hline'+'\n')

        registre.write('\\'+'('+'\\'+'tilde{'+'\\'+'rho'+'}'+'\\'+')'+'&')
        registre.write('Nodes'+'&')

        registre.write('\\'+'('+'{'+'\\'+'frac{'+'\\'+'int'+'\\'+'nabla'+'\\'+'chi}{'+'\\'+'Omega}}_{ROM}'+'\\'+'&')
        registre.write('\\'+'('+'{'+'\\'+'frac{'+'\\'+'int'+'\\'+'nabla'+'\\'+'chi}{'+'\\'+'Omega}}_{FOM}'+'\\'+'&')
        registre.write('\\'+'('+'Err'+')'+'\\'+'&')

        registre.write('\\'+'('+'t_{build}'+'\\'+')'+'&')
        registre.write('\\'+'('+'t_{solve}'+'\\'+')'+'&')
        registre.write('\\'+'('+'t_{FEM}'+'\\'+')'+'&')

        registre.write('\\'+'('+'\\'+'mathcal{G}_{rom}'+'\\'+')'+'&')
        registre.write('\\'+'('+'\\'+'mathcal{G}_{rom-sol}'+'\\'+')'+'\\'+'\\'+'\n')

        for rho in list_rho_test:

            r_nouv=rho
            exec(open("DDD_EIV.py").read())


        registre.write('\\'+'end{tabular}')

        registre.close()
        nom_tab_latex = 'LaTeXArticle/'+ config + '_res' + str(res_gmsh) + '_raydeb_o' + str(int(100*2*list_rho_appr[0]))# geo_p + 'Npod' + str(N_mor) +
        os.rename(nom_fichier + '.txt', nom_tab_latex + '.tex')

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
#   registre=open(nom_fichier,'w')
#   #
#   for i in range(deb,len(l_rho)):
#     print('rayon',str(l_rho[i]),'maillage :',l_tm[i],'secondes')
#     print('##############################################################################')
#     #
#     r_nouv=l_rho[i]
#     t_meshing=l_tm[i]
#     with open("DDD_EIV.py",'r') as op:
#      exec(op.read())
#   registre.close()
