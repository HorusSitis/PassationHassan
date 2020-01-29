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

# Choix de la resolution du maillage : nombre de noeuds par cote du cube

res_gmsh=10
# res_gmsh=20
# res_gmsh=25

typ_msh='gms'
# typ_msh=''

if typ_msh=='gms':
    res=res_gmsh

# configuration du domaine periodique

tol=1e-10

xinf=-1.0
yinf=-1.0
zinf=-1.0
xsup=1.0
ysup=1.0
zsup=1.0

xyzinfsup = [[xinf, yinf, zinf], [xsup, ysup, zsup]]

dimension=3

# mesh_repository = 'maillages_per/2D/'
mesh_repository = 'maillages_per/' + str(dimension) + 'D/'

if res_gmsh==5:
    lw=0.27
elif res_gmsh==10:
    lw=0.15
elif res_gmsh == 20 or res_gmsh == 25:
    lw=0.01



class PeriodicBoundary(SubDomain):
    # Left boundary is 'target domain' G
    def inside(self, x, on_boundary):
        return on_boundary and not(near(x[0],xsup,tol) or near(x[1],ysup,tol) or near(x[2],zsup,tol))
    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        for i in range(dimension):
            if near(x[i],xsup,tol):
                y[i]=xinf
            else:
                y[i]=x[i]



# parametres pour l'execution des etapes : affichage, tests de periodicite etc

fig_todo = 'aff'
# fig_todo = 'save'

typ_msh='gms'#''
D_k=1.0

# N_snap=len(list_rho_appr)

N_snap = 8

npas_err=10
# typ_sol='bic_cyr'#'default'#seulement si res=10##


typ_sol = 'kr_null_vect'

ordo='Ordr'#'Nordr'

Nrefine=1
crow=(1/res_gmsh)*1e-1
typ_refi='vol'#'front'#

# apprentissage : calcul parallele ou sequentiel, prise en compte de la resolution

# gen_snap='par8'
# gen_snap='par4'
# gen_snap='par2'
gen_snap='seq'
# gen_snap='seq_par'

# repertoire pour les resultats
repertoire_parent='Res3D/'



# -------------------- Geometrie du probleme -------------------- #

config = 'sph_un'
# config = 'cyl_un'
# config = 'cylsph'
# config = '2sph'

## pour le rom
# ray_fix = 0.36

ray_fix = 0.39
## pour des exemples avec MEF
# ray_fix = 0.5

### ------------ Important : liste des rayons pour l'apprentissage et les tests ------------ ###

N_snap = 8

## apprentissage
if config == 'sph_un' or ray_fix <= 0.5:

    rho_appr_min = 0.1
    # rho_appr_min = 0.2
    rho_appr_max = 0.8
    # rho_appr_max = 0.9

    list_rho_appr = np.linspace(rho_appr_min, rho_appr_max, N_snap)

## test
rho_diff = (rho_appr_max - rho_appr_min)/(N_snap - 1)
rho_test_min = rho_appr_min + 0.5*rho_diff
rho_test_max = rho_appr_max - 0.5*rho_diff

list_rho_test = np.array([0.35])
# list_rho_test = np.array([0.35, 0.65])
# list_rho_test = np.linspace(0.11, 0.44, 4)
list_rho_test = np.linspace(rho_test_min, rho_test_max, N_snap - 1)

### ------------ Parametres geometriques : configurations ------------ ###

### inclusions simples
if config == 'sph_un':

    dom_fixe = 'am'
    geo_p = 'ray'#'cen'#
    ##
    conf_mess = 'sphere unique'

    if geo_p == 'ray':
        geo_mess='rayon variable'
        cen_snap_ray = [0.5,0.5,0.5]
    elif geo_p == 'cen':
        geo_mess = 'centre variable'
        ray_snap_cen = 0.35
        csr_list = [[0.5,0.5,0.3+0.05*k] for k in range(1,1+Nsnap)]

elif config == 'cyl_un':

    dom_fixe = 'am'
    geo_p = 'ray'#'axe'#
    ##
    conf_mess = 'cylindre unique'
    if geo_p == 'ray':
        geo_mess = 'rayon variable'
        cen_snap_ray = [0.5,0.,0.5]
        top_snap_ray = [0.5,0.5]
    elif geo_p == 'axe':
        asr_list = [[0.5,0.3+0.05*k] for k in range(1,1+Nsnap)]

### inclusions composees
elif config == '2sph':
    conf_mess = 'deux spheres'
    dom_fixe = 'solid'

    geo_p = 'ray'
    geo_mess = 'rayon de la sphere centrale variable'
    ##

elif config == 'cylsph':
    conf_mess = 'un cylindre et une sphere'
    dom_fixe = 'ray_min'

    geo_p = 'ray_sph'
    # geo_p = 'ray_cyl'

    ##
    if geo_p=='ray_cyl':
        geo_mess='rayon du cylindre variable'
        ## utilisation du domaine fixe avec annulation du rayon du cylindre dans le fichier general ##
        fixe_comp='cylsph'#'sph_un'#'ray_min'#
    elif geo_p=='ray_sph':
        geo_mess='rayon de la sphere variable'
    elif geo_p=='ray_linked':
        geo_mess='rayons lies'


## ------------ Porosite ------------ ##

size = (xsup - xinf)
cell_vol = size**dimension

def epsilon_p(r, config, geo_p, r_f):

    if config == 'sph_un':
        fluid_vol = cell_vol - 4/3*pi*r**3
    elif config=='cyl_un':
        fluid_vol = cell_vol - size*pi*r**2
    elif config=='2sph':
        fluid_vol = cell_vol - 4/3*pi*(r**3+r_f**3)
    elif config=='cylsph' :
        if geo_p == 'ray_sph':
            fluid_vol = cell_vol - 4/3*pi*r**3 - size*pi*r_f**2
        elif geo_p == 'ray_cyl':
            fluid_vol = cell_vol - 4/3*pi*r_f**3 - size*pi*r**2

    epsilon_p = fluid_vol/cell_vol

    return epsilon_p

## ------------ Pour les noms de fichiers de maillages ------------ ##

if config == 'sph_un':
    mesh_prefix = 'cubesphere_periodique_triangle_'
elif config == '2sph':
    mesh_prefix = 'cube2sph_periodique_triangle_'
elif config == 'cylsph':
    mesh_prefix = 'cubecylsph_periodique_triangle_'

## ------------ Etape lL 1demi : Affichage de microstructures periodiques ------------ ##

nb_lcells = 5
cem_color = 'grey'
sand_color = 'orange'
fluid_color = 'cyan'

## -------------------- Etape III -------------------- ##

if config == 'sph_un' or config == 'cyl_un':
    rg_perf_fact = ''
else:
    rg_perf_fact = '_rayf' + str(int(round(100*ray_fix, 2)))

seuil_ener = 99.99

nu = 1 - seuil_ener/100
nu_log = log(nu)/log(10)
expo = str(int(round(nu_log, 0)))

registre_N_mor_name = 'Perf3D/' + 'N_mor_' + 'ener_nu10E' + expo + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)


## -------------------- Etape IV -------------------- ##


# La mesure du temps d'execution doit se faire avec l'option 'save' de fig_todo

ind_fixe = True ##-----------> dom_fixe devant le 'Phi'
ind_res = True #False ###----------> on precise la resolution du maillage, qui apparait ou non dans le fichier contenant Phi

## -------------------- Etape V -------------------- ##



## ------------ Registre pour les performances des tests ------------ ##

# if config == 'sph_un' or config == 'cyl_un':
#     rg_perf_fact = ''
# else:
#     rg_perf_fact = '_rayf' + str(int(round(100*ray_fix, 2)))

registre_perf_num = dict()

registre_perf_num['int_grad_fem'] = 'Perf3D/' + 'IG_fem_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)
registre_perf_num['int_grad_rom'] = 'Perf3D/' + 'IG_rom_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)

if config == 'cylsph':
    registre_perf_num['int_grad_yy_fem'] = 'Perf3D/' + 'IG_yy_fem_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)
    registre_perf_num['int_grad_yy_rom'] = 'Perf3D/' + 'IG_yy_rom_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)


registre_perf_num['nodes'] = 'Perf3D/' + 'nodes_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)
registre_perf_num['err_rel'] = 'Perf3D/' + 'err_rel_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)
registre_perf_num['var_rel'] = 'Perf3D/' + 'var_rel_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)
registre_perf_num['var_rel_yy'] = 'Perf3D/' + 'var_rel_yy_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)
registre_perf_num['t'] = 'Perf3D/' + 'tps_exec_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)

registre_perf_num['var_rel_chi'] = 'Perf3D/' + 'var_rel_chi_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)

## ------------ Registres de snapshots et vecteurs POD ------------ ##

## Solutions physiques vectorisees
l_name='Lchi_'+str(N_snap)+'_'+config+'_'+geo_p+'_'+'sur'+str(res)+'_'+ordo+'_'+computer

## Snapshots vectorises utilisables pour la POD
u_name = 'Usnap_' + dom_fixe + '_' + str(N_snap) + '_' + config + '_' + geo_p + '_' + 'res' + str(res) + '_' + ordo + '_' + computer

## Chargement de la base POD complete
phi_name='Phi'+dom_fixe+'_dim'+str(N_snap)+'_'+config+'_'+geo_p+'_'+'res'+str(res)+'_'+ordo+'_'+computer

# if ind_res:
#     phi_name='Phi'+dom_fixe+'_dim'+str(N_snap)+'_'+config+'_'+geo_p+'_'+'res'+str(res)+'_'+ordo+'_'+computer
# elif ind_fixe:
#     phi_name='Phi'+dom_fixe+'_dim'+str(N_snap)+'_'+config+'_'+geo_p+'_'+ordo+'_'+computer
# else:
#     phi_name='Phi'+'_dim'+str(N_snap)+'_'+config+'_'+geo_p+'_'+ordo+'_'+computer
