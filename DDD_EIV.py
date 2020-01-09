
# -*- coding: Latin-1 -*-
#################################################################################################
## Etape IV : Predictions. Choisir les parametres du probleme a resoudre par le modele reduit.#

# mesh_repository = 'maillages_per/' + str(dimension) + 'D/'
# if config == 'sph_un':
#     mesh_prefix = 'cubesphere_periodique_triangle_'
# elif config == '2sph':
#     mesh_prefix = 'cube2sph_periodique_triangle_'
# elif config == 'cylsph':
#     mesh_prefix = 'cubecylsph_periodique_triangle_'

## maillage du domaine fixe

if dom_fixe == 'am':
    mesh_f_name = 'cube_periodique_triangle' + '_' + dom_fixe + '_sur' + str(res_gmsh) + 'fixe'
## inclusions multiples, unique rayon variable
elif dom_fixe == 'solid':
    # mesh_fixe_prefix ='cube' + config + '_periodique_triangle_'
    mesh_fixe_prefix = mesh_prefix
    if config == '2sph':
        mesh_f_name = mesh_fixe_prefix + 'fixe_som' + '_rayp' + str(int(round(100*ray_fix,2))) + '_sur' + str(res_gmsh)
elif dom_fixe == 'ray_min':
    # utilisation du domaine fixe avec annulation du rayon du cylindre dans le fichier general
    fixe_comp = True
    if config == 'cylsph':
        if geo_p == 'ray_sph':
            mesh_f_name = mesh_prefix + 'rayc' + str(int(round(100*rho_appr_min,2))) + '_rayp' + str(int(round(100*ray_fix,2))) + '_sur' + str(res_gmsh)
        elif geo_p == 'ray_cyl':
            mesh_f_name = mesh_prefix + 'rayc' + str(int(round(100*ray_fix,2))) + '_rayp' + str(int(round(100*rho_appr_min,2))) + '_sur' + str(res_gmsh)

mesh_fixe = Mesh(mesh_repository + mesh_f_name + '.xml')

# fonctions test du domaine fixe
V_fixe = VectorFunctionSpace(mesh_fixe,'P',2,constrained_domain=PeriodicBoundary())

# Performances

import time

### ------------ Etapes reproduites : dependances directes de Main3D ------------ ###

# nu = 1 - seuil_ener/100
# nu_log = log(nu)/log(10)
# expo = str(int(round(nu_log, 0)))
#
# registre_N_mor_name = 'Perf3D/' + 'N_mor_' + 'ener_nu10E' + expo + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)
# np.save(registre_Nmor_name + '.npy')

# tab_N_mor = np.load(registre_N_mor_name + '.npy')
# N_mor = tab_N_mor[0]

nb_modes = N_mor

# --------------------- SE0 : maillage et fonctions tests du domaine fixe --------------------- #

start_mesh = time.time()

if config == 'sph_un':
    mesh_nouv_name = mesh_prefix + 'rayc' + str(int(round(100*r_nouv,2))) + '_sur' + str(res)
elif config == '2sph':
    mesh_nouv_name = mesh_prefix + 'rayc' + str(int(round(100*r_nouv,2))) + '_rayp' + str(int(round(100*ray_fix,2))) + '_sur' + str(res)
elif config == 'cylsph' and geo_p == 'ray_sph':
    mesh_nouv_name = mesh_prefix + 'rayc' + str(int(round(100*r_nouv,2))) + '_rayp' + str(int(round(100*ray_fix,2))) + '_sur' + str(res)
elif config == 'cylsph' and geo_p == 'ray_cyl':
    mesh_nouv_name = mesh_prefix + 'rayc' + str(int(round(100*ray_fix,2))) + '_rayp' + str(int(round(100*r_nouv,2))) + '_sur' + str(res)

# generation des fichiers avec gmsh
creer_maill_per_gpar(config, geo_p, xyzinfsup, r_nouv, ray_fix, res_gmsh)

# appelle le maillage reconverti depuis maillages_per/3D
mesh_nouv = Mesh(mesh_repository + mesh_nouv_name + '.xml')

# fin de la generation du maillage courant
end_mesh = time.time()

t_meshing = end_mesh - start_mesh

# creation de l'espace des champs admissibles
V_nouv = VectorFunctionSpace(mesh_nouv, 'P', 2, constrained_domain=PeriodicBoundary())


# --------------------- SE1 : projection de la base POD sur le nouveau domaine --------------------- #

## On initialise le temps de calcul ##

start=time.time()

## Taille du maillage du domaine fixe ##

nb_noeuds_fixe = V_fixe.dim()

with sh.open(repertoire_parent+phi_name) as phi_loa:
    Phi_prime_v = phi_loa['maliste']

## Creation de la base POD tronquee, sous forme vectorielle

Phi_mor=Phi_prime_v[:,range(0,nb_modes)]

## Extrpolation des fonctions de la base POD pour former le modele reduit defini sur V_nouv

nb_noeuds_nouv=V_nouv.dim()
Phi_nouv_v=np.zeros((nb_noeuds_nouv,nb_modes))

list_pod_nouv=[]
phi_nouv=Function(V_nouv)
phi_fixe=Function(V_fixe)

for n in range(0,nb_modes):
    phi_fixe.vector().set_local(Phi_mor[:,n])
    # extrapolation du snapshot au domaine fixe
    phi_fixe.set_allow_extrapolation(True)
    phi_n_nouv=interpolate(phi_fixe,V_nouv)
    # on range le vecteur de POD interpolee dans la matrice Phi_nouv_v
    Phi_nouv_v[:,n]=phi_n_nouv.vector().get_local()

## Stockage de la matrice du modele reduit

## On enregistre et imprime le temps d'execution de SE1

end=time.time()

t_phi_nouv=end-start

print('se1 faite ',t_phi_nouv,' secondes')

# --------------------- SE2 : resolution du modele reduit --------------------- #

## On reinitialise le temps de calcul ##

start=time.time()

## On ecrit les deux tenseurs qui comportent les coefficients de l'equation du modele reduit : ceux-ci dependent des vecteurs de la base POD projetee

if config=='sph_un' or config=='cyl_un':
    Coeff=calc_Ab_3D(V_nouv, mesh_nouv, Phi_nouv_v, r_nouv, cen_snap_ray, nb_modes, config)
else:
    Coeff=calc_Ab_compl_3D(mesh_repository + mesh_nouv_name, Phi_nouv_v, nb_modes)

A=Coeff[0]
b=Coeff[1]

end=time.time()

t_int_Ab = end - start

# print('A :',A,'b :',b)
## On resoud le modele reduit

start = time.time()

a_nouv=np.linalg.solve(A.T,-b)

## On enregistre et imprime le temps d'execution de SE2

end=time.time()

t_rom_linear=end-start

print('se2 faite', t_int_Ab + t_rom_linear, 'secondes')

# --------------------- SE3 : calcul du nouveau champ de vecteurs, affichage --------------------- #

## On reinitialise le temps de calcul ##

start=time.time()

## On initialise et affiche le champ chi_nouv

chi_nouv_v=np.dot(Phi_nouv_v,a_nouv)
chi_nouv=Function(V_nouv)
chi_nouv.vector().set_local(chi_nouv_v)

plot(chi_nouv, linewidth=lw)
plt.title('Rho = 0,'+str(int(round(100*r_nouv,2))),fontsize=40)
if fig_todo=='aff':
    plt.show()
elif fig_todo=='save':
    plt.savefig('Figures3D/sol_rom'+str(int(round(100*r_nouv,2)))+'_sur'+str(N_snap)+config+'_'+geo_p+'res'+str(res)+'.png')
plt.close()

## Exploitation du champ ainsi obtenu
r=r_nouv
rho=r_nouv

## a faire avec une instruction conditionnelle
if err_per_calc:
    # Affichage des valeurs et erreurs de la solution periodique, quelle que soit la configuration
    if config=='sph_un' or config=='cyl_un':
        err_per_gr(cen_snap_ray,r_nouv,chi_nouv,npas_err,fig_todo)
    elif config=='2sph':
        err_per_gr_compl(config,r_v_0,chi_nouv,npas_err,fig_todo)
    elif config=='cylsph':
        if geo_p=='ray_sph':
            err_per_gr_compl(config,r_c_0,chi_nouv,npas_err,fig_todo)
        elif geo_p=='ray_cyl':
            err_per_gr_compl(config,r_nouv,chi_nouv,npas_err,fig_todo)


# Tenseur de diffusion homogeneise

## Integrale de chi sur le domaine fluide
T_chi_rom=np.zeros((3,3))
for k in range(0,3):
    for l in range(0,3):
        T_chi_rom[k,l]=assemble(grad(chi_nouv)[k,l]*dx)
T_chi_rom_omega = T_chi_rom/cell_vol
## Integrale de l'identite sur le domaine fluide
### Calcul de la porosite
porosity = epsilon_p(r_nouv, config, geo_p, ray_fix)
### Integration du terme constant du coefficient d diffusion, sur le domaine fluide
D = porosity*np.eye(3)
## Calcul et affichage du tenseur Dhom
Dhom_kMOR=D_k*(D + T_chi_rom_omega.T)

print('Coefficient Dhom_k11 '+conf_mess+', '+geo_mess+' valeur '+str(rho)+' MOR :',Dhom_kMOR[0,0])

## On enregistre et imprime le temps d'execution de SE3

end=time.time()

t_rom_Dhom=end-start

print('se3 faite ',t_rom_Dhom,' secondes')

# --------------------- SE4 : comparaison avec la methode des elements finis --------------------- #

## On reinitialise le temps de calcul ##

start=time.time()

## On reinitialise le champ chi_nouv pour la methode des elements finis

#res=20
cen_snap_ray=[(xinf + xsup)/2., (zinf + zsup)/2., (zinf + zsup)/2.]

if config=='sph_un':
    chi_nouv=snapshot_sph_per(cen_snap_ray,r_nouv,res,typ_sol)
elif config=='cyl_un':
    chi_nouv=snapshot_cyl_per(top_snap_ray,r_nouv,res,typ_sol)
elif config=='2sph':
    if geo_p=='ray':
        chi_nouv = snapshot_compl_per(r_nouv, ray_fix, config, res_gmsh)
elif config == 'cylsph':# or config == '2sph':
    chi_compl = snapshot_compl_per(r_nouv, ray_fix, config, res_gmsh)

## Exploitation du champ ainsi obtenu
rho=r_nouv
r=r_nouv

plot(chi_nouv, linewidth=lw)
plt.title('Rho = 0,'+str(int(round(100*r_nouv,2))),fontsize=40)
if fig_todo=='aff':
    plt.show()
elif fig_todo=='save':
    plt.savefig('Figures3D/sol_rom'+str(int(round(100*r_nouv,2)))+'_sur'+str(N_snap)+config+'_'+geo_p+'res'+str(res)+'.png')
plt.close()

## a faire avec une instruction conditionnelle
if err_per_calc:
    # Affichage des valeurs et erreurs de la solution periodique, quelle que soit la configuration
    if config=='sph_un' or config=='cyl_un':
        #err_per_ind_01(chi_n,cen,r,npas_err)
        err_per_gr(cen_snap_ray,r_nouv,chi_nouv,npas_err,fig_todo)
    elif config=='2sph':
        err_per_gr_compl(config,r_v_0,chi_nouv,npas_err,fig_todo)
    elif config=='cylsph':
        if geo_p=='ray_sph':
            err_per_gr_compl(config,r_c_0,chi_nouv,npas_err,fig_todo)
        elif geo_p=='ray_cyl':
            err_per_gr_compl(config,r_nouv,chi_nouv,npas_err,fig_todo)

# Tenseur de diffusion homogeneise

## Integrale de chi sur le domaine fluide
T_chi_fom = np.zeros((3, 3))
for k in range(0, 3):
    for l in range(0, 3):
        T_chi_fom[k,l] = assemble(grad(chi_nouv)[k, l]*dx)
T_chi_fom_omega = T_chi_fom/cell_vol

## Integrale de l'identite sur le domaine fluide : voir ce qui precede avec la porosite
print('Noeuds :', V_nouv.dim())
print('Porosite :', porosity)

## Calcul et affichage du tenseur Dhom
Dhom_kMEF=D_k*(D+T_chi_fom_omega.T)
print('Coefficient Dhom_k11 '+conf_mess+', '+geo_mess+' valeur '+str(rho)+ ' MEF :',Dhom_kMEF[0,0])

## Comparaison
err_rel_dhom = 100*(Dhom_kMOR[0,0] - Dhom_kMEF[0,0])/Dhom_kMEF[0,0]
print('Erreur relative Dhom MEF-MOR :', err_rel_dhom , ' pourcent')

err_rel_ig = 100*(T_chi_rom_omega[0,0] - T_chi_fom_omega[0,0])/T_chi_fom_omega[0,0]
print('Erreur relative int_grad MEF-MOR :', err_rel_ig , ' pourcent')

var_rel_ig = (np.sqrt(2*(T_chi_rom_omega[0,0]**2 + T_chi_fom_omega[0,0]**2)/((T_chi_rom_omega[0,0] + T_chi_fom_omega[0,0])**2) - 1))
var_rel_ig_yy = (np.sqrt(2*(T_chi_rom_omega[1,1]**2 + T_chi_fom_omega[1,1]**2)/((T_chi_rom_omega[1,1] + T_chi_fom_omega[1,1])**2) - 1))

## On enregistre et imprime le temps d'execution de SE4

end=time.time()

t_fem=end - start

print('se4 faite ', t_fem, ' secondes')

##############################################################################
########################## Evaluation de la methode ##########################
##############################################################################

if Report :
    print('#'*78)
    print('#'*26+' Evaluation de la methode '+'#'*26)
    # print('########################## Evaluation de la methode ##########################')
    print('#'*78)
    print('Resultats '+conf_mess+', '+geo_mess+' valeur '+str(r_nouv)+' :')
    print('#'*78)    #
    ## Porosite ##
    print('Porosite :',porosity)
    ## Generation du maillage : temps d'execution ... ##
    print('Maillage :',t_meshing,'secondes')
    ## Dhom, erreur, et temps d'execution ##
    print('%'*78)
    print('Coefficient Dhom_k11 MOR :',Dhom_kMOR[0,0])
    print('Coefficient Dhom_k11 MEF :',Dhom_kMEF[0,0])
    print('Erreur relative MEF-MOR :',err_rel_dhom,'pourcent')
    print('%'*78)
    print('Coefficient int_grad11 MOR :',T_chi_rom_omega[0,0])
    print('Coefficient int_grad11 MEF :',T_chi_fom_omega[0,0])
    print('Erreur relative MEF-MOR :',err_rel_ig,'pourcent')
    print('%'*78)
    ## Temps de calcul a evaluer ##
    t_rom = t_phi_nouv + t_int_Ab + t_rom_linear + t_rom_Dhom
    print('Tps_ROM :', t_rom)
    print('-'*78)
    print('Tps_phi_nouv :',t_phi_nouv,'secondes')
    print('Tps_int_Ab :', t_int_Ab, 'secondes')
    print('Tps_solve :', t_rom_linear, 'secondes')
    print('Tps_dhom :',t_rom_Dhom, 'secondes')
    print('='*78)
    print('Tps_FEM :', t_fem, 'secondes')
    ## Nombre de noeuds ##
    print('Noeuds :',V_nouv.dim())
    ## Rapports de temps de calcul, sans unite ##
    print('%'*78)
    R_rom=(t_phi_nouv+t_int_Ab+t_rom_linear+t_rom_Dhom+t_meshing)/(t_fem+t_meshing)
    R_interpolation=t_phi_nouv/(t_fem+t_meshing)
    print('Gain de temps maillage EF compris :', 1./R_rom)#,'sans unite')
    # print('Contribution de l interpolation :',R_interpolation)#,'sans unite')
    # print('Difference :',R_rom-R_interpolation)#,'sans unite')
    print('='*78)
    R_rom_Nmaill = (t_phi_nouv+t_int_Ab+t_rom_linear+t_rom_Dhom)/t_fem
    R_interpolation_Nmaill = t_phi_nouv/t_fem
    Rdiff_rom = R_rom_Nmaill - R_interpolation_Nmaill
    R_rom_solve = t_rom_linear/t_fem
    print('Gain de temps sans maillage :', 1./R_rom_Nmaill)#,'sans unite')
    # print('Interpolation sans maillage :', 1./R_interpolation_Nmaill)#,'sans unite')
    print('ROM seul sans maillage :', 1./Rdiff_rom)#,'sans unite')
    print('Resolution du ROM seule :', round(0.0001/R_rom_solve, 2), '10E4')
    #
    print('#'*78)
    print('#'*78)


##############################################################################
##############################################################################

## ecriture des resultats dans le tableau _pg

# registre_pg.write(str(2*r_nouv)+'&')
registre_pg.write(str(r_nouv)+'&')
registre_pg.write(str(V_nouv.dim())+'&')

registre_pg.write(str(round(T_chi_rom_omega[0,0], 4))+'&')
registre_pg.write(str(round(T_chi_fom_omega[0,0], 4))+'&')
registre_pg.write(str(round(err_rel_ig, 2))+'\\'+'%'+'&')

registre_pg.write(str(round(1./R_rom_Nmaill, 2))+'\\'+'\\'+'\n')
# registre_pg.write(str(round(1./Rdiff_rom, 2))+'&')
# registre_pg.write(str(round(0.0001/R_rom_solve, 2))+'\\'+'('+'\\'+'cdot 10^4'+'\\'+')'+'\\'+'\\'+'\n')

registre_pg.write('\\'+'hline'+'\n')

## ecriture des ersultats dans le tableau _gr

# registre_gr.write(str(2*r_nouv)+'&')
registre_gr.write(str(r_nouv)+'&')
registre_gr.write(str(V_nouv.dim())+'&')

registre_gr.write(str(round(100*(t_phi_nouv/t_rom), 2))+'\\'+'%'+'&')
registre_gr.write(str(round(100*(t_int_Ab/t_rom), 2))+'\\'+'%'+'&')
registre_gr.write(str(round(100*(t_rom_linear/t_rom), 4))+'\\'+'%'+'&')
registre_gr.write(str(round(100*(t_rom_Dhom/t_rom), 2))+'\\'+'%'+'\\'+'\\'+'\n')
# registre_pg.write(str(round(t_fem, 2))+'s'+'&')


registre_gr.write('\\'+'hline'+'\n')


##############################################################################
##############################################################################
##############################################################################
##############################################################################


## ecriture des resultats dans les vecteurs arr_

arr_int_grad_fem[i] = T_chi_fom_omega[0,0]
arr_int_grad_rom[i] = T_chi_rom_omega[0,0]

if config == 'cylsph':
    arr_int_grad_yy_fem[i] = T_chi_fom_omega[0,0]
    arr_int_grad_yy_rom[i] = T_chi_rom_omega[0,0]

arr_nodes[i] = V_nouv.dim()

arr_err_rel[i] = err_rel_ig
arr_var_rel[i] = var_rel_ig
arr_var_rel_yy[i] = var_rel_ig_yy


arr_t[i, 0] = t_fem

t_rom = t_phi_nouv + t_int_Ab + t_rom_linear + t_rom_Dhom
arr_t[i, 1] = t_rom

arr_t[i, 2] = t_phi_nouv
arr_t[i, 3] = t_int_Ab
arr_t[i, 4] = t_rom_linear
arr_t[i, 5] = t_rom_Dhom
