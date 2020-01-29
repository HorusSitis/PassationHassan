
# -*- coding: Latin-1 -*-
#################################################################################################
## Etape IV : Predictions. Choisir les parametres du probleme a resoudre par le modele reduit. ##
#################################################################################################

# from tools_chi import *

# Definition du domaine Omega_fixe :
if dom_fixe == 'am':
    mesh_fixe_name = 'maillage_fixe2d_am'#.xml'
elif config == 'compl':
    mesh_fixe_name = 'maillage_trous2D_'+geo_p+'_fixe'#.xml'


# print(mesh_fixe_name)
mesh_fixe = Mesh(mesh_repository + mesh_fixe_name + '.xml')

V_fixe = VectorFunctionSpace(mesh_fixe, 'P', VFS_degree, constrained_domain=PeriodicBoundary())

# --------------------- SE0 : maillage et fonctions tests du domaine fixe --------------------- #

start_mesh = time.time()

creer_maill_per_gpar(config, geo_p, mention, xyinfsup, rho, ray_p)

if config == 'compl':
    mesh_name = mesh_prefix + str(int(round(100*rho,2))) + '_rayp' + str(int(round(100*ray_p,2)))
else:
    mesh_name = mesh_prefix + mention + str(int(round(100*rho,2)))

print('%'*80)
print('Maillage rom : ' + mesh_name)
print('%'*80)

# appelle le maillage reconverti depuis maillages_per/2D
mesh_nouv_rom = Mesh(mesh_repository + mesh_name + '.xml')

end_mesh = time.time()

# pointe vers le meme maillage que la fonction snapshot_ de DD_fun_obj
V_nouv = VectorFunctionSpace(mesh_nouv_rom, 'P', VFS_degree, constrained_domain=PeriodicBoundary())

# Performances
tps_mesh = end_mesh - start_mesh

### ------------ Etapes reproduites : dependances directes de Main3D ------------ ###

nb_modes = N_mor

# --------------------- SE1 : projection de la base POD sur le nouveau domaine --------------------- #

## On initialise le temps de calcul ##

start = time.time()

## Taille du maillage du domaine fixe ##

nb_noeuds_fixe = V_fixe.dim()

## Chargement de la base POD complete

with sh.open(repertoire_parent+phi_name) as phi_loa:
    Phi_prime_v = phi_loa['maliste']

## Creation de la base POD tronquee, sous forme vectorielle

Phi_mor=Phi_prime_v[:,range(0,nb_modes)]

## Extrapolation des fonctions de la base POD pour former le modele reduit defini sur V_nouv

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
    # affichage des modes extrapoles
    plot(phi_n_nouv)
    # plt.title('Phi '+str(n+1)+' sur Omega_nouv',fontsize=30)
    if fig_todo=='aff':
        plt.title('Phi '+str(n+1)+' sur Omega_nouv',fontsize=30)
        plt.show()
    elif fig_todo=='save':
        plt.savefig('Figures2D/phi_nouv_'+str(n+1)+'_'+config+'_'+geo_p+'.png')
    plt.close()
    # on range le vecteur de POD interpolee dans la matrice Phi_nouv_v
    Phi_nouv_v[:,n]=phi_n_nouv.vector().get_local()

## On enregistre et imprime le temps d'execution de SE1

end=time.time()

t_phi_nouv=end-start

print('se1 faite ',t_phi_nouv,' secondes')

# --------------------- SE2 : resolution du modele reduit --------------------- #

## On reinitialise le temps de calcul ##

start=time.time()

## On ecrit les deux tenseurs qui comportent les coefficients de l'equation du modele reduit : ceux-ci dependent des vecteurs de la base POD projetee

if test_snap=='i_per':
    Coeff=calc_Ab_2D(V_nouv,mesh_nouv_rom,Phi_nouv_v,r_nouv,cen_snap_ray,nb_modes)
else:
    Coeff=calc_Ab_compl(V_nouv,mesh_nouv_rom,Phi_nouv_v,nb_modes,test_snap)
    print('modele reduit complexe utilise')

A=Coeff[0]
b=Coeff[1]

end=time.time()

t_int_Ab = end - start

## On resoud le modele reduit

start = time.time()

a_nouv=np.linalg.solve(A.T,-b)

## On enregistre et imprime le temps d'execution de SE2

end=time.time()

t_rom_linear=end-start

print('se2 faite', t_int_Ab + t_rom_linear, 'secondes')
# print(A,b,a_nouv)

# --------------------- SE3 : calcul du nouveau champ de vecteurs, affichage --------------------- #

## On reinitialise le temps de calcul ##

start=time.time()

## On initialise et affiche le champ chi_nouv_rom

chi_nouv_rom_v=np.dot(Phi_nouv_v,a_nouv)
chi_nouv_rom=Function(V_nouv)
chi_nouv_rom.vector().set_local(chi_nouv_rom_v)

plot(chi_nouv_rom)#, linewidth=0.55)
plt.title('Solution ROM', fontsize=30)#'Rho = 0,'+str(int(round(100*r_nouv,2))),fontsize=40)
if fig_todo=='aff':
    plt.show()
else:
    plt.savefig('Figures2D/solROM_'+config+'_'+geo_p+str(int(round(100*r_nouv,2)))+'.png')
plt.close()

## Exploitation du champ ainsi obtenu
r=r_nouv
rho=r_nouv

# Affichage des valeurs et erreurs de la solution periodique, quelle que soit la configuration

# Tenseur de diffusion homogeneise
## Integrale de chi sur le domaine fluide
T_chi_rom=np.zeros((2,2))
for k in range(0,2):
    for l in range(0,2):
        T_chi_rom[k,l]=assemble(grad(chi_nouv_rom)[k,l]*dx)
## Integrale de l'identite sur le domaine fluide
if config!='compl':
    por=(1-pi*r**2)
else:
    por=1-pi*(r**2+0.15**2)
D=por*np.eye(2)
## Calcul et affichage du tenseur Dhom
Dhom_kMOR=D_k*(D+T_chi_rom.T)
print('Coefficient Dhom_k11 '+conf_mess+', '+geo_mess+' valeur '+str(rho)+' MOR :',Dhom_kMOR[0,0])

## On enregistre et imprime le temps d'execution de SE3

end=time.time()

t_rom_Dhom=end-start

print('se3 faite ',t_rom_Dhom,' secondes')

# --------------------- SE4 : comparaison avec la methode des elements finis --------------------- #

## On reinitialise le temps de calcul ##

start=time.time()

## On reinitialise le champ chi_nouv_full pour la methode des elements finis


if test_snap=='i_per':
    chi_nouv_full=snapshot_circ_per(cen_snap_ray,r_nouv,res)
else:
    # chi_nouv_full=snapshot_compl_per(geo_p,r_nouv,cen_snap_ray,mention,test_snap)
    chi_nouv_full=snapshot_compl_per(geo_p,r_nouv, cen_snap_ray, test_snap, ray_p)# mention,

## Exploitation du champ ainsi obtenu
rho=r_nouv
r=r_nouv

# Affichage des valeurs et erreurs de la solution periodique, quelle que soit la configuration
#err_per_ind_01(chi_n,cen,r,npas_err)
for npas_test in []:#30,40]:#7,15,16,60,125,250]:
    err_per_gr(cen_snap_ray,r_nouv,chi_nouv_full,npas_test,fig_todo)

# Tenseur de diffusion homogeneise
## Integrale de chi sur le domaine fluide
T_chi_fom=np.zeros((2,2))
for k in range(0,2):
    for l in range(0,2):
        T_chi_fom[k,l]=assemble(grad(chi_nouv_full)[k,l]*dx)
## Integrale de l'identite sur le domaine fluide
if config!='compl':
    por=(1-pi*r**2)
else:
    por=1-pi*(r**2+0.15**2)
D=por*np.eye(2)
print('Noeuds :',V_nouv.dim())
print('Porosite :',por)
## Calcul et affichage du tenseur Dhom
Dhom_kMEF=D_k*(D+T_chi_fom.T)
print('Coefficient Dhom_k11 '+conf_mess+', '+geo_mess+' valeur '+str(rho)+ ' MEF :',Dhom_kMEF[0,0])

## Sortie graphique

plot(chi_nouv_full)#, linewidth=0.55)
plt.title('Solution EF',fontsize=30)
if fig_todo=='aff':
    plt.show()
else:
    plt.savefig('Figures2D/solFEM_'+config+'_'+geo_p+str(int(round(100*r_nouv,2)))+'.png')
plt.close()

## Comparaison

err_rel_dhom=100*(Dhom_kMOR[0,0]-Dhom_kMEF[0,0])/Dhom_kMEF[0,0]
print('Erreur relative Dhom MEF-MOR :', err_rel_dhom , ' pourcent')

err_rel_ig=100*(T_chi_rom[0,0]-T_chi_fom[0,0])/T_chi_fom[0,0]
print('Erreur relative int_grad MEF-MOR :', err_rel_ig , ' pourcent')

var_rel_ig = 100*(np.sqrt(2*(T_chi_rom[0,0]**2 + T_chi_fom[0,0]**2)/((T_chi_rom[0,0] + T_chi_fom[0,0])**2) - 1))
if config == 'compl':
    var_rel_ig_yy = 100*(np.sqrt(2*(T_chi_rom[1,1]**2 + T_chi_fom[1,1]**2)/((T_chi_rom[1,1] + T_chi_fom[1,1])**2) - 1))



## On enregistre et imprime le temps d'execution de SE4

end=time.time()

t_fem=end-start

print('se4 faite ',t_fem,' secondes')

chi_nouv_rom.set_allow_extrapolation(True)
# chi_nouv_rom_prime = interpolate(chi_nouv_rom, V_nouv)
chi_nouv_full_prime = interpolate(chi_nouv_full, V_nouv)

var_rel_chi = 100*var_rel_func(chi_nouv_rom, chi_nouv_full_prime)

##############################################################################
########################## Evaluation de la methode ##########################
##############################################################################

if Report :
    print('#'*78)
    print('#'*26+' Evaluation de la methode '+'#'*26)
    print('#'*78)
    print('Resultats '+conf_mess+', '+geo_mess+' valeur '+str(r_nouv)+' :')
    print('#'*78)    #
    ## Porosite ##
    print('Porosite :',por)
    ## Dhom, erreur, et temps d'execution ##
    print('%'*78)
    print('Coefficient Dhom_k11 MOR :',Dhom_kMOR[0,0])
    print('Coefficient Dhom_k11 MEF :',Dhom_kMEF[0,0])
    print('Erreur relative MEF-MOR :',err_rel_dhom,'pourcent')
    print('%'*78)
    print('Coefficient int_grad11 MOR :',T_chi_rom[0,0])
    print('Coefficient int_grad11 MEF :',T_chi_fom[0,0])
    print('Erreur relative MEF-MOR :',err_rel_ig,'pourcent')
    print('%'*78)
    print('Variance relative MEF-MOR IG :',var_rel_ig,'pourcent')
    print('Variance relative MEF-MOR chi :',var_rel_chi,'pourcent')
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
    R_rom_Nmaill = (t_phi_nouv+t_int_Ab+t_rom_linear+t_rom_Dhom)/t_fem
    R_rom_solve = t_rom_linear/t_fem
    print('Gain de temps :', 1./R_rom_Nmaill)
    # print('Resolution du ROM seule :', round(0.0001/R_rom_solve, 2), '10E4')
    #
    print('#'*78)
    print('#'*78)

##############################################################################
##############################################################################

## ecriture des resultats dans le tableau _pg

# registre_pg.write(str(2*r_nouv)+'&')
registre_pg.write(str(r_nouv)+'&')
registre_pg.write(str(V_nouv.dim())+'&')

registre_pg.write(str(round(T_chi_rom[0,0], 4))+'&')
registre_pg.write(str(round(T_chi_fom[0,0], 4))+'&')
registre_pg.write(str(round(err_rel_ig, 3))+'\\'+'%'+'&')

registre_pg.write(str(round(1./R_rom_Nmaill, 2))+'\\'+'\\'+'\n')

registre_pg.write('\\'+'hline'+'\n')

## ecriture des ersultats dans le tableau _gr

registre_gr.write(str(r_nouv)+'&')
registre_gr.write(str(V_nouv.dim())+'&')

registre_gr.write(str(round(100*(t_phi_nouv/t_rom), 2))+'\\'+'%'+'&')
registre_gr.write(str(round(100*(t_int_Ab/t_rom), 2))+'\\'+'%'+'&')
registre_gr.write(str(round(100*(t_rom_linear/t_rom), 4))+'\\'+'%'+'&')
registre_gr.write(str(round(100*(t_rom_Dhom/t_rom), 2))+'\\'+'%'+'\\'+'\\'+'\n')


registre_gr.write('\\'+'hline'+'\n')

## ecriture des resultats dans les vecteurs arr_

arr_int_grad_fem[i] = T_chi_fom[0,0]
arr_int_grad_rom[i] = T_chi_rom[0,0]

arr_nodes[i] = V_nouv.dim()

arr_err_rel[i] = err_rel_ig
arr_var_rel[i] = var_rel_ig

arr_var_rel_chi[i] = var_rel_chi

arr_t[i, 0] = t_fem

t_rom = t_phi_nouv + t_int_Ab + t_rom_linear + t_rom_Dhom
arr_t[i, 1] = t_rom

arr_t[i, 2] = t_phi_nouv
arr_t[i, 3] = t_int_Ab
arr_t[i, 4] = t_rom_linear
arr_t[i, 5] = t_rom_Dhom
