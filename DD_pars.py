### ------------ Paquets a importer ------------ ###

import numpy as np

import os
import sys

# from fenics import *
from dolfin import *
# ### --- from mshr import * --- ###
# import matplotlib.pyplot as plt
import numpy as np

##########################################################
### ------------ Code a lire : conditions ------------ ###
##########################################################

fixe_aff=False

# fig_todo=''
fig_todo='aff'
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

config='cer_un'
# config='compl'

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
    # pour des procedures communes a toutes les configurations
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

if config == 'compl':
    if geo_p == 'diag':
        ray_p = 0.3
    elif geo_p == 'hor':
        ray_p = 0.25


## ------------ Porosite ------------ ##

size = (xsup - xinf)
cell_vol = size**dimension

def epsilon_p(r, config,  r_f):

    # size = (xsup - xinf)
    # cell_vol = size**dimension

    if config!='compl':
        fluid_vol = (cell_vol - pi*r**2)
    else:
        fluid_vol = cell_vol - pi*(r**2+r_f**2)

    epsilon_p = fluid_vol/cell_vol

    return epsilon_p


### ------------ Important : liste des rayons pour l'apprentissage et les tests ------------ ###

N_snap = 8

# if xinf == 0.:
#     if config != 'compl' or (config == 'compl' and geo_p == 'diag'):
#         rho_appr_min = 0.05
#         # rho_appr_min = 0.1
#         rho_appr_max = 0.4
#         # rho_appr_max = 0.45
#         list_rho_test = np.linspace(0.11, 0.44, 4)
#     elif config == 'compl' and geo_p == 'hor':
#         rho_appr_min = 0.01
#         rho_appr_max = 0.028
#         list_rho_test = np.linspace(0.04, 0.1, 0.2, 0.3)


if config != 'compl' or (config == 'compl' and geo_p == 'diag'):
    rho_appr_min = 0.1
    rho_appr_max = 0.8
    list_rho_test = np.linspace(0.15, 0.75, 7)
elif config == 'compl' and geo_p == 'hor':
    rho_appr_min = 0.07
    rho_appr_max = 0.56
    list_rho_test = np.linspace(0.105, 0.595, 8)

list_rho_appr = np.linspace(rho_appr_min, rho_appr_max, N_snap)

# --------------------- Important : pas des maillages --------------------- #
res_gmsh = 25
# ------------------------------------------ #

typ_msh='gms'
# typ_msh=''

if typ_msh=='gms':
    res=res_gmsh
    res_fixe=res_gmsh

# ------------------------------------------ #

# raffinement de maillages : ocouronnes, volumes
Nrefine = 1

# crow=(1/res_gmsh)*1e-1

typ_refi = 'front'
# typ_refi='vol'

lg_crow=-1
crow=2*10**(lg_crow)

if typ_refi=='vol':
    refi_mess='Couronne : '+str(crow)
elif typ_refi=='front':
    refi_mess='Surface'

# ------------------------------------------ #
# ------------------------------------------ #

# nom de l'appareil utilise pour generer les donnees enregistrees
computer='MECALAC_29x8'##'T1700_35x8'##

# repertoire pour les resultats

repertoire_parent='Res2D/'

### ------------ Execution des etapes demandees en preambule, imports specifiques ------------ ###

from LEc import *

## ------------ Etape lL 1demi : Affichage de microstructures periodiques ------------ ##

nb_lcells = 8

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

## -------------------- Etape III -------------------- ##

if config == 'cer_un':
    rg_perf_fact = ''
    ray_fix = 0
else:
    ray_fix = ray_p
    rg_perf_fact = '_rayf' + str(int(round(100*ray_fix, 2)))

if config == 'sph_un' or config == 'cyl_un':
    rg_perf_fact = ''
else:
    rg_perf_fact = '_rayf' + str(int(round(100*ray_fix, 2)))

seuil_ener_pour = 99.99

from math import log

nu = 1 - seuil_ener_pour/100
nu_log = log(nu)/log(10)
expo = str(int(round(nu_log, 0)))

registre_N_mor_name = 'Perf2D/' + 'N_mor_' + 'ener_nu10E' + expo + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)

## ---------- Etape IV ---------- ##



## ---------- Etape V ---------- ##



## ------------ Registre pour les performances des tests ------------ ##





registre_perf_num = dict()

registre_perf_num['int_grad_fem'] = 'Perf2D/' + 'IG_fem_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)
registre_perf_num['int_grad_rom'] = 'Perf2D/' + 'IG_rom_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)

if config == 'compl':
    registre_perf_num['int_grad_yy_fem'] = 'Perf2D/' + 'IG_yy_fem_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)
    registre_perf_num['int_grad_yy_rom'] = 'Perf2D/' + 'IG_yy_rom_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)


registre_perf_num['nodes'] = 'Perf2D/' + 'nodes_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)
registre_perf_num['err_rel'] = 'Perf2D/' + 'err_rel_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)
registre_perf_num['var_rel'] = 'Perf2D/' + 'var_rel_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)
registre_perf_num['var_rel_yy'] = 'Perf2D/' + 'var_rel_yy_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)
registre_perf_num['t'] = 'Perf2D/' + 'tps_exec_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)

registre_perf_num['var_rel_chi'] = 'Perf2D/' + 'var_rel_chi_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)

## ------------ Registres de snapshots et vecteurs POD ------------ ##

## Solutions physiques vectorisees
l_name='Lchi_'+str(N_snap)+'_'+config+'_'+geo_p+'_'+'sur'+str(res)+'_'+ordo+'_'+computer

## Snapshots vectorises utilisables pour la POD
u_name = 'Usnap_' + dom_fixe + '_' + str(N_snap) + '_' + config + '_' + geo_p + '_' + 'res' + str(res) + '_' + ordo + '_' + computer

## Chargement de la base POD complete
phi_name='Phi'+dom_fixe+'_dim'+str(N_snap)+'_'+config+'_'+geo_p+'_'+'res'+str(res)+'_'+ordo+'_'+computer
