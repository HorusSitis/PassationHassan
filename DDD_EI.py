#######################################################################################################################################
## Etape I : realisation des cliches, avec la methode des elements finis. Calcul du tenseur d'homogeneisation. Stockage dans snap2D/ ##
#######################################################################################################################################


## Boucle pour la creation des snapshots, avec un parametre pouvant etre le rayon d'une inclusion circulaire, ou l'emplacement de son centre

# Pour avoir des fonctions "top-level" a paralleliser

## Sphere unique

#if geo_p=='ray':
cen_snap_ray=[0.5,0.5,0.5]
def snap_sph_ray(N_par):
    # rho : on utilise la liste d'apprentissage definie dans DDD_geoset
    rho = list_rho_appr[N_par]

    chi_r=snapshot_sph_per(cen_snap_ray,rho,res_gmsh,typ_sol)
    chi_r_v=chi_r.vector().get_local()

    return([N_par,chi_r_v])

# #if geo_p=='cen':
# #ray_snap_cen=0.25
# #csr_list=[[0.5,0.5,0.05*k] for k in range(1,1+N_snap)]
# #c_par : parametre scalaire pour la position du centre
# def snap_sph_cen(c_par):
#     cen_snap_ray=csr_list[c_par-1]
#     chi_c=snapshot_sph_per(cen_snap_ray,ray_snap_cen,res_gmsh)
#     chi_c_v=chi_c.vector().get_local()
#     return([c_par,chi_c_v])

## Cylindre unique

#if geo_p=='ray':
axe_snap_ray=[0.5,0.5]
def snap_cyl_ray(N_par):
    # rho : on utilise la liste d'apprentissage definie dans DDD_geoset
    rho = list_rho_appr[N_par]

    chi_r=snapshot_cyl_per(axe_snap_ray,rho,res_gmsh,typ_sol)
    chi_r_v=chi_r.vector().get_local()
    return([r_par,chi_r_v])

# #if geo_p=='axe':
# #ray_snap_axe=0.25
# asr_list=[[0.5,0.3+0.05*k] for k in range(1,1+N_snap)]
# #c_par : parametre scalaire pour la position du centre
# def snap_cyl_axe(c_par):
#     cen_snap_ray=csr_list[c_par-1]
#     chi_c=snapshot_sph_per(axe_snap_ray,ray_snap_axe,res)
#     chi_c_v=chi_c.vector().get_local()
#     return([c_par,chi_c_v])



r_s_0=0
r_c_0=0
r_v_0=0

if config=='2sph':
    r_v_0=0.15
elif config=='cylsph':
    if geo_p=='ray_cyl':
        r_s_0=0.15
    elif geo_p=='ray_sph':
        r_c_0=0.15


def snap_compl_ray(N_par):
    # rho : on utilise la liste d'apprentissage definie dans DDD_geoset
    rho = list_rho_appr[N_par]
    ## deux spheres ##
    if geo_p=='ray':
        chi_compl=snapshot_compl_per(rho,r_v_0,config,res_gmsh,typ_sol)
    ## un cylindre et une sphere ##
    elif geo_p=='ray_sph':
        chi_compl=snapshot_compl_per(rho,r_c_0,config,res_gmsh,typ_sol)
    elif geo_p=='ray_cyl':
        chi_compl=snapshot_compl_per(r_s_0,rho,config,res_gmsh,typ_sol)
    ## on vectorise la fonction calculee par MEF ##
    chi_compl_v=chi_compl.vector().get_local()
    ## on renvoie un vecteur etiquete, utilisable avec l'option 'par8' ##
    return([rho_par,chi_compl_v])






# ------------------------- Snapshots, conditionnellement ------------------------- #

if not snap_done:

    # -------- Calcul des snapshots, sous forme vectorielle, avec des etiquettes -------- #
    ### Generation parallele des snapshots ###
    if gen_snap=='par8' or gen_snap=='par4':
        if gen_snap=='par8':
            nproc=8
        else:
            nproc=4
        pool=multiprocessing.Pool(processes=nproc)
        if config=='sph_un':
            if geo_p=='ray':
                list_chi_n_v=pool.map(snap_sph_ray,(N for N in range(0,N_snap)))
            # elif geo_p=='cen':
            #     list_chi_n_v=pool.map(snap_sph_cen,(N for N in range(0,N_snap)))
        elif config=='cyl_un':
            if geo_p=='ray':
                list_chi_n_v=pool.map(snap_cyl_ray,(N for N in range(0,N_snap)))
            # elif geo_p=='axe':
            #     list_chi_n_v=pool.map(snap_cyl_axe,(N for N in range(0,N_snap)))
            else:
                list_chi_n_v=pool.map(snap_compl_ray,(N for N in range(0,N_snap)))
    ### Generation sequentielle des snapshots, pour des tests de la methode des elements finis ###
    elif gen_snap=='seq':
        start=time.time()
        list_chi_n_v=[]
        for n in range(1,1+N_snap):
            print(n)
            if config=='sph_un':
                if geo_p=='ray':
                    list_chi_n_v.append(snap_sph_ray(n))
                elif geo_p=='cen':
                    list_chi_n_v.append(snap_sph_cen(n))
            elif config=='cyl_un':
                if geo_p=='ray':
                    list_chi_n_v.append(snap_cyl_ray(n))
                elif geo_p=='axe':
                    list_chi_n_v.append(snap_cyl_axe(n))
            else:
                list_chi_n_v.append(snap_compl_ray(n))
        end=time.time()
        print('temps EF : ',end-start,' secondes')
    # ### Generation parallele pour chaque snapshot, pour de gros maillages ###
    # elif gen_snap=='seq_par':
    #     list_chi_n_v=[]
    # -------- enregistrement des fonctions vectorisees dans une liste -------- #
    # Construction de la liste des snapshots vectorises : cas d'un parametre geometrique definissant un ordre - lien avec la porosite ; ou non.
    list_chi_v=[]
    if geo_p=='ray' or config=='compl':
        for n in range(0,N_snap):
            for i in range(0,N_snap):
                if list_chi_n_v[i][0]==n:
                    chi_n_v=list_chi_n_v[i][1]
                    list_chi_v.append(chi_n_v)
    else:
        for i in range(0,N_snap):
            chi_n_v=list_chi_n_v[i][1]
            list_chi_v.append(chi_n_v)

    # Liste des snapshots : sauvegarde, on precise l'identite de la machine qui a effectue le calcul
    l_name='Lchi_'+str(N_snap)+'_'+config+'_'+geo_p+'_'+"sur"+str(res)+'_'+ordo+'_'+computer
    # sauvegarde de la liste des solutions indexees calculees avec la methode des elements finis
    with sh.open(repertoire_parent+l_name) as l_sto:
        l_sto["maliste"] = list_chi_v

# Matrice des snapshots : plus tard, voir l'etape II

else :
    l_name='Lchi_'+str(N_snap)+'_'+config+'_'+geo_p+'_'+"sur"+str(res)+'_'+ordo+'_'+computer
    with sh.open(repertoire_parent+l_name) as l_loa:
        list_chi_v = l_loa["maliste"]

# --------------------------------------------------------------------------------- #

# Exploitation des solution du probleme aux elements finis
res=res_gmsh

if res==10:
    lw=0.27
elif res==20:
    lw=0.15
elif res==50:
    lw=0.01

for n in range(0,N_snap):
    # Extraction du snapshot de rang n
    chi_n_v=list_chi_v[n]
    r=list_rho_appr[n]
    # On cree un maillage pour reecrire les snapshots sous la forme de fonctions

    if config=='sph_un':
        if geo_p=='ray':
            cen=cen_snap_ray
            # r=list_rho_appr[n]
            if typ_msh=='gms':
                mesh_name="cubesphere_periodique_triangle_"+str(int(round(100*r,2)))+"sur"+str(res)
            else:
                mesh=creer_maill_sph(cen,r,res)
        elif geo_p=='cen':
            r=ray_snap_cen
            mesh=creer_maill_sph(csr_list[n],r,res)

    elif config=='cyl_un':## avec gmsh
        if geo_p=='ray':
            # r=list_rho_appr[n]
            mesh_name="cubecylindre_periodique_triangle_"+str(int(round(100*r,2)))+"sur"+str(res)
        #elif geo_p=='axe':
            #r=ray_snap_ax
            #mesh=creer_maill_cyl(acr_list[n-1],r,res)

    elif config=='2sph':
        # r=n*0.05
        r_s=r
        r_v=r_v_0
        mesh_name="cube"+config+"_periodique_triangle_"+str(int(round(100*r_s,2)))+str(int(round(100*r_v_0,2)))+"sur"+str(res)

    elif config=='cylsph':
        # r=n*0.05
        if geo_p=='ray_cyl':
            r_c=r
            r_s=r_s_0
        elif geo_p=='ray_sph':
            r_c=r_c_0
            r_s=r
        mesh_name="cube"+config+"_periodique_triangle_"+str(int(round(100*r_c,2)))+str(int(round(100*r_s,2)))+"sur"+str(res)

    print("maillages_per/3D/"+mesh_name+".xml")
    mesh=Mesh("maillages_per/3D/"+mesh_name+".xml")
    V_n=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
    # On restitue la forme fonctionnelle du snapshot courant
    chi_n=Function(V_n)
    chi_n.vector().set_local(chi_n_v)
    # Representation graphique
    plot(chi_n, linewidth=lw)
    plt.tight_layout(pad=0)
    if r<0.1:
        plt.title("Rho = 0,0"+str(int(round(100*r,2))), fontsize=40)
    else:
        plt.title("Rho = 0,"+str(int(round(100*r,2))),fontsize=40)
    if fig_todo=='aff':
        plt.show()
    else:
        plt.savefig("Figures3D/sol_"+str(n)+"_sur"+str(N_snap)+config+'_'+geo_p+"res"+str(res)+".png")
    plt.close()
    # Affichage des valeurs et erreurs de la solution periodique, quelle que soit la configuration
    #err_per_ind_01(chi_n,cen,r,npas_err)
    if config=='cyl_un' and geo_p=='ray':
        cen=[0.5,0.,0.5]# on triche un peu : on prend une face prevee d'une demie-sphere au lieu d'une face privee du disque frontal du cylindre
    if config=='2sph':
        err_per_gr_compl(config,r_v,chi_n,npas_err,fig_todo)
    elif config=='cylsph':
        err_per_gr_compl(config,r_c,chi_n,npas_err,fig_todo)
    else:
        err_per_gr(cen,r,chi_n,npas_err,fig_todo)
    # Tenseur de diffusion homogeneise
    ## Integrale de khi sur le domaine fluide
    T_chi=np.zeros((3,3))
    for k in range(0,3):
        for l in range(0,3):
            T_chi[k,l]=assemble(grad(chi_n)[k,l]*dx)
    #print(T_chi)
    ## Integrale de l'identite sur le domaine fluide
    if config=='sph_un':
        por=1-4/3*pi*r**3
    elif config=='cyl_un':
        por=1-pi*r**2
    elif config=='2sph':
        por=1-4/3*pi*(r_s**3+r_v**3)
    elif config=='cylsph' :
        por=1-4/3*pi*r_s**3-pi*r_c**2
    D=por*np.eye(3)
    ## Calcul et affichage du tenseur Dhom
    Dhom_k=D_k*(D+T_chi.T)
    #print(('Tenseur Dhom_k',Dhom_k))
    print("Noeuds",V_n.dim())
    print("Porosite :",por)
    print('Coefficient Dhom_k11EF, snapshot '+str(n)+", "+conf_mess+', '+geo_mess+" :",Dhom_k[0,0])
    integ=assemble(chi_n[1]*dx)
    print('Valeur moyenne : ',integ)
    ## Anisotropie
    mod_diag=max(abs(Dhom_k[0,0]),abs(Dhom_k[1,1]),abs(Dhom_k[2,2]))
    mod_ndiag=0
    for i in range(0,3):
        for j in range(0,i):
            if abs(Dhom_k[i,j])>mod_ndiag:
                mod_ndiag=abs(Dhom_k[i,j])
        for j in range(i+1,3):
            if abs(Dhom_k[i,j])>mod_ndiag:
                mod_ndiag=abs(Dhom_k[i,j])
    #print("Anisotropie : ",mod_ndiag/mod_diag)
