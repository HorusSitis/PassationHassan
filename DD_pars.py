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

### ------------ Implementation du domaine periodique ------------ ###

tol=1e-10

xinf=0.0
yinf=0.0
xsup=1.0
ysup=1.0

# determiner le domaine fixe pour interpoler la solution

dimension=2

class PeriodicBoundary(SubDomain):
    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return on_boundary and not(near(x[0],xsup,tol) or near(x[1],ysup,tol))## merci a Arnold Douglas
    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        for i in range(dimension):
            if near(x[i],1.0,tol):
                y[i]=0.0
            else:
                y[i]=x[i]

### ------------------ Important : degre pour la resolution par elements finis ------------------ ###
VFS_degree=2
# VFS_degree=3
## degre 2 : comme en dimension 3, permet d'eviter les erreurs de periodicite pour des pas qui nen sont pas de la forme 2"n, ou n est un diviseur de 100 ##

config='cer_un'
# config='compl'

if config=='cer_un':
    test_snap='i_per'
    # test_snap = ''
    # test_snap = 'solid_1'
    ### on exclut 'solid_2' ###
    dom_fixe= "am"
    # dom_fixe = "multiray"
    # dom_fixe = "ray_min"
    ##
    geo_p='ray'
    # geo_p = "cen"
    cen_snap_ray=[0.5,0.5]
    # cen_snap_ray = [0.,0.]
    ##
    conf_mess='disque unique'
    ##
    if geo_p=='ray':
        geo_mess='rayon variable'
    ### geo_p='centre'
    # ray_snap_cen=0.25
    # csr_list=[[0.05*k,0.5] for k in range(1,1+Nsnap)]
    ##
    if cen_snap_ray==[0.5,0.5]:
        conf_mess=conf_mess+" centre"
        mention=""
        ### 'i_per' : Nrom=... ; 'solid_1': Nrom=... ; ### ???
    elif cen_snap_ray==[0.,0.]:
        conf_mess=conf_mess+" aux sommets"
        mention="_som"
        config=config+mention
        ### 'i_per' : Nrom=2 ; 'solid_1': Nrom=5 ; 'solid_2': Nrom=2 ---> ne fonctionne pas
elif config=='compl':
    test_snap='solid_1'#'solid_2'#''#
    ##
    dom_fixe="solid"
    # dom_fixe = "am"
    ##
    geo_p='hor'
    # geo_p='diag'
    ##
    if geo_p=='diag':
        cen_snap_ray=[0.,0.]
        ### 'solid_1': Nrom=3 ;
    elif geo_p=='hor':
        cen_snap_ray=[0.,0.5]
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

list_rho_appr = np.linspace(rho_appr_min, rho_appr_max, N_snap)


# choix du type de maillage

res_gmsh=100

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

repertoire_parent="Res2D/"

### ------------ Execution des etapes demandees en preambule, imports specifiques ------------ ###

from LEc import *


## ------------ Etape lL 1demi : Affichage de microstructures periodiques ------------ ##

nb_lcells=1
cem_color='grey'
sand_color='orange'
fluid_color='cyan'

## ---------- Etape I, memes parametres que pour 1demi ---------- ##
## ---------- Etape II ---------- ##

D_k=1.0
# Nsnap=8

# important :
deb=1

# par defaut##5# pour les tests cer_un et diag##6# pour hor##
npas_err=50
ordo='Ordr'
# ordo='Nordr'
# apprentissage : calcul parallele ou sequentiel, prise en compte de la resolution

gen_snap = 'par8'
# gen_snap = 'par4'
#
# gen_snap = 'seq'

# # non prepare dans EI-II
# gen_snap = 'seq_par'

## ---------- Etape III ---------- ##

seuil_ener_pour=99.99

# On soustrait eventuellement la moyenne des snapshots interpoles a chaque snapshot
## moy_mod='soust_m'#'' par defaut


## 99,99% :
### compl diag N_mor=3 pour "solid", 4 pour "am"
### compl hor N_mor=4

## ---------- Etape IV ---------- ##

# cer_un 99,99% ; ...
N_mor=2

# r_nouv=0.11
# r_nouv=0.22
r_nouv=0.33
# r_nouv=0.44
