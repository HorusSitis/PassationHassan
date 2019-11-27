### ------------ Paquets a importer ------------ ###

import numpy as np

import os
import sys

# from fenics import *
from dolfin import *
from mshr import *
# import matplotlib.pyplot as plt
import numpy as np

##########################################################
### ------------ Code a lire : conditions ------------ ###
##########################################################

fixe_aff=False

# fig_todo=''
# fig_todo='aff'
fig_todo='save'

import time


mesh_repository = 'maillages_per/2D/'


### ------------ Implementation du domaine periodique ------------ ###

tol=1e-10

xinf=-1.0
yinf=-1.0
# xinf=0.
# yinf=0.
xsup=1.0
ysup=1.0

xyinfsup = [[xinf, yinf], [xsup, ysup]]

# determiner le domaine fixe pour interpoler la solution

dimension=2

class PeriodicBoundary(SubDomain):
    # Left boundary is 'target domain' G
    def inside(self, x, on_boundary):
        return on_boundary and not(near(x[0],xsup,tol) or near(x[1],ysup,tol))## merci a Arnold Douglas
    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        for i in range(dimension):
            if near(x[i],xsup,tol):
                y[i]=xinf
            else:
                y[i]=x[i]

### ------------------ Important : degre pour la resolution par elements finis ------------------ ###
VFS_degree=2
# VFS_degree=3
## degre 2 : comme en dimension 3, permet d'eviter les erreurs de periodicite pour des pas qui nen sont pas de la forme 2'n, ou n est un diviseur de 100 ##

# config='cer_un'
config='compl'

if config=='cer_un':
    # pour editer les fichiers de maillages
    mesh_prefix = 'maillage_trou2D_'
    # espace des snapshots
    test_snap='i_per'
    # test_snap = ''
    # test_snap = 'solid_1'
    dom_fixe= 'am'
    ##
    geo_p='ray'
    # geo_p = 'cen'
    cen_snap_ray = [(xinf + xsup)/2., (yinf + ysup)/2.]
    # cen_snap_ray = [xinf, yinf]
    ##
    conf_mess='disque unique'
    ##
    if geo_p=='ray':
        geo_mess='rayon variable'
    ### geo_p='centre'
    ##
    if cen_snap_ray == [(xinf + xsup)/2., (yinf + ysup)/2.]:
        conf_mess = conf_mess+' centre'
        mention = ''
    elif cen_snap_ray == [xinf, yinf]:
        conf_mess = conf_mess+' aux sommets'
        mention = '_som'
    config = config + mention
    # pour des procedures communes aa toutes les configurations
    ray_p = 0.
elif config=='compl':
    ray_p = 0.3
    # espace des snapshots
    test_snap = 'solid_1'#''#
    ##
    dom_fixe='solid'
    # dom_fixe = 'am'
    ##
    geo_p='hor'
    # geo_p='diag'
    # pour editer les fichiers de maillages
    mesh_prefix = 'maillage_trous2D_'+geo_p+'_'
    ##
    if geo_p=='diag':
        cen_snap_ray=[xinf, yinf]
    elif geo_p=='hor':
        cen_snap_ray=[xinf, (yinf + ysup)/2.]
    ##
    conf_mess='deux disques par periode'
    mess_prefix=' rayon central variable'
    mention=''
    if geo_p=='diag':
        geo_mess='alignes en diagonale, '+mess_prefix
    elif geo_p=='hor':
        geo_mess='alignes horizontalement, '+mess_prefix

### ------------ Important : liste des rayons pour l'apprentissage et les tests ------------ ###

N_snap = 8

if xinf == 0.:
    if config != 'compl' or (config == 'compl' and geo_p == 'diag'):
        rho_appr_min = 0.05
        # rho_appr_min = 0.1
        rho_appr_max = 0.4
        # rho_appr_max = 0.45
        list_rho_test = np.linspace(0.11, 0.44, 4)
    elif config == 'compl' and geo_p == 'hor':
        rho_appr_min = 0.01
        rho_appr_max = 0.028
        list_rho_test = np.linspace(0.04, 0.1, 0.2, 0.3)


if config != 'compl' or (config == 'compl' and geo_p == 'diag'):
    rho_appr_min = 0.1
    rho_appr_max = 0.8
    list_rho_test = np.linspace(0.15, 0.85, 8)
elif config == 'compl' and geo_p == 'hor':
    rho_appr_min = 0.07
    rho_appr_max = 0.56
    list_rho_test = np.linspace(0.105, 0.595, 8)

list_rho_appr = np.linspace(rho_appr_min, rho_appr_max, N_snap)

# choix du type de maillage

res_gmsh = 100

typ_msh='gms'
# typ_msh=''

if typ_msh=='gms':
    res=res_gmsh
    res_fixe=res_gmsh

# raffinement de maillages : ocouronnes, volumes

Nrefine=1

# crow=(1/res_gmsh)*1e-1

typ_refi='front'
# typ_refi='vol'

lg_crow=-1
crow=2*10**(lg_crow)

if typ_refi=='vol':
    refi_mess='Couronne : '+str(crow)
elif typ_refi=='front':
    refi_mess='Surface'

# nom de l'appareil utilise pour generer les donnees enregistrees
computer='MECALAC_29x8'##'T1700_35x8'##

# repertoire pour les resultats

repertoire_parent='Res2D/'

### ------------ Execution des etapes demandees en preambule, imports specifiques ------------ ###

from LEc import *

## ------------ Etape lL 1demi : Affichage de microstructures periodiques ------------ ##

nb_lcells = 1

cem_color = 'grey'
sand_color = 'orange'
fluid_color = 'cyan'

## ---------- Etape I, memes parametres que pour 1demi ---------- ##

# solveur choisi

# typ_sol = 'lu'
# typ_sol = 'bic_cyr'
typ_sol = 'kr_null_vect'

## ---------- Etape II ---------- ##

D_k=1.0

# important :
deb = 1

# par defaut
npas_err = 50
ordo = 'Ordr'
# ordo='Nordr'

# apprentissage : calcul parallele ou sequentiel, prise en compte de la resolution

gen_snap = 'par8'
# gen_snap = 'par4'
# gen_snap = 'seq'

# # non prepare dans EI-II
# gen_snap = 'seq_par'

## ---------- Etape III ---------- ##

seuil_ener_pour=99.99

## ---------- Etape IV ---------- ##

N_mor=2
