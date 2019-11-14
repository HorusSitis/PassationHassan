# -*- coding: utf-8 -*-
### Une commande possible dans le terminal ###

#--- mpirun -np 8 python3 Main3D.py ---#
#--- affiche npfois 'pas encore' fait avec l'etape IV ---#

# Attention : on execute parallelement

### Paquets a importer ###

# from fenics import *
from dolfin import *
from mshr import *
# import matplotlib.pyplot as plt
import numpy as np

##########################################################
### ------------ Code a lire : conditions ------------ ###
##########################################################

# configuration du domaine periodique

tol=1e-10

xinf=0.0
yinf=0.0
zinf=0.0
xsup=1.0
ysup=1.0
zsup=1.0

dimension=3

if res_gmsh==10:
    lw=0.27
elif res_gmsh==20:
    lw=0.15
elif res_gmsh==50:
    lw=0.01

r_s_0=0.15
r_v_0=0.15
r_c_0=0.15

r_min=0.05

class PeriodicBoundary(SubDomain):
    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return on_boundary and not(near(x[0],xsup,tol) or near(x[1],ysup,tol) or near(x[2],zsup,tol))
    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        for i in range(dimension):
            if near(x[i],1.0,tol):
                y[i]=0.0
            else:
                y[i]=x[i]


### ------------ Important : liste des rayons pour l'apprentissage et les tests ------------ ###

N_snap = 8

rho_appr_min = 0.05
# rho_appr_min = 0.1
rho_appr_max = 0.4
# rho_appr_max = 0.45

list_rho_appr = np.linspace(rho_appr_min, rho_appr_max, N_snap)
list_rho_test = np.linspace(0.11, 0.44, 4)


# Choix de la resolution du maillage : nombre de noeuds par cote du cube

res_gmsh=20

typ_msh='gms'
# typ_msh=''

if typ_msh=='gms':
    res=res_gmsh

# -------------------- Geometrie du probleme -------------------- #

config='sph_un'
# config='cyl_un'
# config='cylsph'
# config='2sph'

### inclusions simples
if config=='sph_un':

    dom_fixe="am"
    geo_p='ray'#'cen'#
    ##
    conf_mess='sphere unique'

    if geo_p=='ray':
        geo_mess='rayon variable'
        cen_snap_ray=[0.5,0.5,0.5]
    elif geo_p=='cen':
        geo_mess='centre variable'
        ray_snap_cen=0.35
        csr_list=[[0.5,0.5,0.3+0.05*k] for k in range(1,1+Nsnap)]

elif config=='cyl_un':

    dom_fixe="am"
    geo_p='ray'#'axe'#
    ##
    conf_mess='cylindre unique'
    if geo_p=='ray':
        geo_mess='rayon variable'
        cen_snap_ray=[0.5,0.,0.5]
        top_snap_ray=[0.5,0.5]
    elif geo_p=='axe':
        asr_list=[[0.5,0.3+0.05*k] for k in range(1,1+Nsnap)]

### inclusions composees
elif config=='2sph':
    conf_mess='deux spheres'
    dom_fixe="solid"#"am"#

    geo_p='ray'
    geo_mess='rayon de la sphere centrale variable'
    ##

elif config=='cylsph':
    conf_mess='un cylindre et une sphere'
    dom_fixe="ray_min"#"am"#"solid"#

    geo_p='ray_sph'#'ray_cyl'#'ray_linked'
    ###ray_sph pour le test diff_ion###
    ##
    if geo_p=='ray_cyl':
        geo_mess='rayon du cylindre variable'
        ## utilisation du domaine fixe avec annulation du rayon du cylindre dans le fichier general ##
        fixe_comp='cyl_sph'#'sph_un'#'ray_min'#
    elif geo_p=='ray_sph':
        geo_mess='rayon de la sphere variable'
    elif geo_p=='ray_linked':
        geo_mess='rayons lies'
