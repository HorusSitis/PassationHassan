
# -*- coding: Latin-1 -*-
#################################################################################################
## Etape IV : Predictions. Choisir les parametres du probleme a resoudre par le modele reduit.#

# maillage du domaine fixe

mesh_dir="maillages_per/3D/"

## inclusions simples ou rayons lies

if dom_fixe=="am":
    mesh_f_name=mesh_dir+"cube_periodique_triangle"+"_"+dom_fixe+"_sur"+str(res_gmsh)+"_fixe.xml"
## inclusions multiples, unique rayon variable
elif dom_fixe=="solid":
    mesh_fixe_prefix=mesh_dir+"cube"+config+"_periodique_triangle_"
    if config=='2sph':
        mesh_f_name=mesh_fixe_prefix+"fixe"+str(int(round(100*r_v_0,2)))+"sur"+str(res_gmsh)+".xml"
    elif config=='cylsph':
        ## rayon du cylindre aux aretes ou de la sphere centrale fixes a 0.15 ##
        if geo_p=='ray_sph':
            mesh_f_name=mesh_fixe_prefix+str(int(round(100*r_c_0,2)))+"fixe"+"sur"+str(res_gmsh)+".xml"
        elif geo_p=='ray_cyl':
            if fixe_comp=='cyl_sph':
                mesh_f_name=mesh_fixe_prefix+"fixe"+str(int(round(100*r_s_0,2)))+"sur"+str(res_gmsh)+".xml"
            elif fixe_comp=='sph_un':
                mesh_f_name=mesh_dir+"cubesphere_periodique_triangle_"+str(int(round(100*r_s_0,2)))+"sur"+str(res_gmsh)+".xml"
elif dom_fixe=="ray_min":
    # utilisation du domaine fixe avec annulation du rayon du cylindre dans le fichier general
    fixe_comp=True
    if config=='cylsph':
        if geo_p=='ray_sph':
            mesh_f_name=mesh_dir+"cube"+config+"_periodique_triangle_"+str(int(round(100*r_c_0,2)))+str(int(round(100*r_min,2)))+"sur"+str(res_gmsh)+".xml"
        elif geo_p=='ray_cyl':
            mesh_f_name=mesh_dir+"cube"+config+"_periodique_triangle_"+str(int(round(100*r_min,2)))+str(int(round(100*r_s_0,2)))+"sur"+str(res_gmsh)+".xml"

mesh_fixe=Mesh(mesh_f_name)

# fonctions test du domaine fixe

V_fixe=VectorFunctionSpace(mesh_fixe,'P',2,constrained_domain=PeriodicBoundary())

# Performances

import time

### ------------ Etapes reproduites : dependances directes de Main3D ------------ ###

nb_modes=N_mor

if config=='sph_un':
    mesh_n_name=mesh_dir+"cubesphere_periodique_triangle_"+str(int(round(100*r_nouv,2)))+"sur"+str(res_gmsh)
elif config=='cyl_un':
    mesh_n_name=mesh_dir+"cubecylindre_periodique_triangle_"+str(int(round(100*r_nouv,2)))+"sur"+str(res_gmsh)
if config=='2sph':
    mesh_n_name=mesh_dir+"cube"+config+"_periodique_triangle_"+str(int(round(100*r_nouv,2)))+str(int(round(100*r_v_0,2)))+"sur"+str(res_gmsh)
elif config=='cylsph':
    if geo_p=='ray_sph':
        mesh_n_name=mesh_dir+"cube"+config+"_periodique_triangle_"+str(int(round(100*r_c_0,2)))+str(int(round(100*r_nouv,2)))+"sur"+str(res_gmsh)
    elif geo_p=='ray_cyl':
        mesh_n_name=mesh_dir+"cube"+config+"_periodique_triangle_"+str(int(round(100*r_nouv,2)))+str(int(round(100*r_s_0,2)))+"sur"+str(res_gmsh)

mesh_nouv=Mesh(mesh_n_name+".xml")

V_nouv=VectorFunctionSpace(mesh_nouv, "P", 2, constrained_domain=PeriodicBoundary())


# --------------------- SE1 : projection de la base POD sur le nouveau domaine --------------------- #

## On initialise le temps de calcul ##

start=time.time()

## Taille du maillage du domaine fixe ##

nb_noeuds_fixe=V_fixe.dim()

## Chargement de la base POD complete

if ind_res:
    phi_name='Phi'+dom_fixe+'_dim'+str(N_snap)+'_'+config+'_'+geo_p+'_'+"res"+str(res)+'_'+ordo+'_'+computer
elif ind_fixe:
    phi_name='Phi'+dom_fixe+'_dim'+str(N_snap)+'_'+config+'_'+geo_p+'_'+ordo+'_'+computer
else:
    phi_name='Phi'+'_dim'+str(N_snap)+'_'+config+'_'+geo_p+'_'+ordo+'_'+computer


with sh.open(repertoire_parent+phi_name) as phi_loa:
    Phi_prime_v = phi_loa["maliste"]

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
    Coeff=calc_Ab_3D(V_nouv,mesh_nouv,Phi_nouv_v,r_nouv,cen_snap_ray,nb_modes,config)
else:
    Coeff=calc_Ab_compl_3D(mesh_n_name,Phi_nouv_v,nb_modes)

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
plt.title("Rho = 0,"+str(int(round(100*r_nouv,2))),fontsize=40)
if fig_todo=='aff':
    plt.show()
elif fig_todo=='save':
    plt.savefig("Figures3D/sol_rom"+str(int(round(100*r_nouv,2)))+"_sur"+str(N_snap)+config+'_'+geo_p+"res"+str(res)+".png")
plt.close()

## Exploitation du champ ainsi obtenu
r=r_nouv
rho=r_nouv

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
        #print("Integrale ",assemble(grad(chi_nouv)[k,l]*dx))
## Integrale de l'identite sur le domaine fluide
### Calcul de la porosite
if config=='sph_un':
    por=1-4/3*pi*r_nouv**3
elif config=='cyl_un':
    por=1-pi*r_nouv**2
elif config=='2sph':
    por=1-4/3*pi*(r_nouv**3+r_v_0**3)
elif config=='cylsph':
    if geo_p=='ray_sph':
        r_s=r_nouv
        r_c=r_c_0
    elif geo_p=='ray_cyl':
        r_s=r_s_0
        r_c=r_nouv
    por=1-4/3*pi*r_s**3-pi*r_c**2
### Integration du terme constant du coefficient d diffusion, sur le domaine fluide
D=por*np.eye(3)
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

## On reinitialise le champ chi_nouv pour la methode des elements finis

#res=20
if config=='sph_un':
    chi_nouv=snapshot_sph_per(cen_snap_ray,r_nouv,res,typ_sol)
elif config=='cyl_un':
    chi_nouv=snapshot_cyl_per(top_snap_ray,r_nouv,res,typ_sol)
elif config=='2sph':
    if geo_p=='ray':
        chi_nouv=snapshot_compl_per(r_nouv,r_v_0,config,res_gmsh,typ_sol)
elif config=='cylsph':
    if geo_p=='ray_sph':
        chi_nouv=snapshot_compl_per(r_nouv,r_c_0,config,res_gmsh,typ_sol)
    elif geo_p=='ray_cyl':
        chi_nouv=snapshot_compl_per(r_s_0,r_nouv,config,res_gmsh,typ_sol)

## Exploitation du champ ainsi obtenu
rho=r_nouv
r=r_nouv

plot(chi_nouv, linewidth=lw)
plt.title("Rho = 0,"+str(int(round(100*r_nouv,2))),fontsize=40)
if fig_todo=='aff':
    plt.show()
elif fig_todo=='save':
    plt.savefig("Figures3D/sol_rom"+str(int(round(100*r_nouv,2)))+"_sur"+str(N_snap)+config+'_'+geo_p+"res"+str(res)+".png")
plt.close()

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
T_chi_fom=np.zeros((3,3))
for k in range(0,3):
    for l in range(0,3):
        T_chi_fom[k,l]=assemble(grad(chi_nouv)[k,l]*dx)
## Integrale de l'identite sur le domaine fluide : voir ce qui precede avec la porosite
print('Noeuds :',V_nouv.dim())
print('Porosite :',por)

## Calcul et affichage du tenseur Dhom
Dhom_kMEF=D_k*(D+T_chi_fom.T)
print('Coefficient Dhom_k11 '+conf_mess+', '+geo_mess+' valeur '+str(rho)+ ' MEF :',Dhom_kMEF[0,0])

## Comparaison

err_rel_dhom=100*(Dhom_kMOR[0,0]-Dhom_kMEF[0,0])/Dhom_kMEF[0,0]
print('Erreur relative Dhom MEF-MOR :', err_rel_dhom , ' pourcent')

err_rel_ig=100*(T_chi_rom[0,0]-T_chi_fom[0,0])/T_chi_fom[0,0]
print('Erreur relative int_grad MEF-MOR :', err_rel_ig , ' pourcent')

## On enregistre et imprime le temps d'execution de SE4

end=time.time()

t_fem=end-start

print('se4 faite ',t_fem,' secondes')

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
    print('Porosite :',por)
    ## Generation du maillage : temps d'execution ... ##
    print('Maillage :',t_meshing,'secondes')
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
    ## Temps de calcul a evaluer ##
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
    print('Gain de temps sans maillage :', 1./R_rom_Nmaill)#,'sans unite')
    # print('Interpolation sans maillage :', 1./R_interpolation_Nmaill)#,'sans unite')
    print('ROM seul sans maillage :', 1./Rdiff_rom)#,'sans unite')
    #
    print('#'*78)
    print('#'*78)


##############################################################################
##############################################################################


registre.write(str(2*r_nouv)+'&')
registre.write(str(V_nouv.dim())+'&')

registre.write(str(round(T_chi_rom[0,0], 4))+'&')
registre.write(str(round(T_chi_fom[0,0], 4))+'&')
registre.write(str(round(err_rel_ig, 2))+'\\'+'%'+'&')

registre.write(str(round(t_phi_nouv, 2))+'s'+'&')
registre.write(str(round(t_int_Ab, 2))+'s'+'&')
registre.write(str(round(t_rom_linear, 2))+'s'+'&')
registre.write(str(round(t_rom_Dhom, 2))+'s'+'&')
registre.write(str(round(t_fem, 2))+'s'+'&')

registre.write(str(round(1./R_rom_Nmaill, 2))+'&')
registre.write(str(round(1./Rdiff_rom, 2))+'\\'+'\\'+'\n')

registre.write('\\'+'hline'+'\n')
