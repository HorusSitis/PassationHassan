#######################################################################################################################################
## Etape I : realisation des cliches, avec la methode des elements finis. Calcul du tenseur d'homogeneisation. Stockage dans snap2D/ ##
#######################################################################################################################################


if typ_msh=='gms':
    res_fixe=res_gmsh
    if dom_fixe=="am":
        mesh_f_name="maillages_per/2D/maillage_fixe2D_am.xml"
    elif config=='compl':
        mesh_f_name="maillages_per/2D/maillage_trous2D_"+geo_p+"_fixe.xml"

## Boucle pour la creation des snapshots, avec un parametre pouvant etre le rayon d'une inclusion circulaire, ou l'emplacement de son centre ##
# Calcule aussi le tenseur de diffusion homogeneise #


## Cercle unique

#if geo_p=='ray':
#cen_snap_ray=[0.5,0.5]
def snap_circ_ray(N_par):
    rho=list_rho_appr[N_par]

    if test_snap=='i_per':
        chi_r=snapshot_circ_per(cen_snap_ray,rho,res)
    else:
        chi_r=snapshot_compl_per(geo_p,rho,cen_snap_ray,mention,test_snap)
    chi_r_v=chi_r.vector().get_local()
    return([N_par,chi_r_v])

# #if geo_p=='cen':
# #ray_snap_cen=0.25
# #csr_list=[[0.5,0.3+0.05*k] for k in range(1,1+N_snap)]
# #c_par : parametre scalaire pour la position du centre
# def snap_circ_cen(c_par):
#     cen_snap_ray=csr_list[c_par-1]
#     chi_c=snapshot_circ_per(cen_snap_ray,ray_snap_cen,res)
#     chi_c_v=chi_c.vector().get_local()
#     return([c_par,chi_c_v])

def snap_compl_ray(N_par):

    rho=list_rho_appr[N_par]

    chi_compl=snapshot_compl_per(geo_p,rho,cen_snap_ray,mention,test_snap)
    chi_compl_v=chi_compl.vector().get_local()

    return([N_par,chi_compl_v])

# ------------------------- Snapshots, conditionnellement ------------------------- #

if not snap_done:

    # Calcul des snapshots, sous forme vectorielle
    if gen_snap=='par8':
        # Generation parallele des snapshots
        pool=multiprocessing.Pool(processes=8)
        if config=='cer_un':
            if geo_p=='ray':
                list_chi_n_v=pool.map(snap_circ_ray,(n for n in range(0,N_snap)))
            elif geo_p=='cen':
                list_chi_n_v=pool.map(snap_circ_cen,(n for n in range(0,N_snap)))
            elif config=='cer_un_som':
                list_chi_n_v=pool.map(snap_circ_ray,(n for n in range(0,N_snap)))
            elif config=='compl':
                list_chi_n_v=pool.map(snap_compl_ray,(n for n in range(0,N_snap)))

    # Generation sequentielle des snapshots
    elif gen_snap=='seq':
        start=time.time()
        list_chi_n_v=[]
        for n in range(0,N_snap):
            if geo_p=='ray':
                list_chi_n_v.append(snap_circ_ray(n))
            elif geo_p=='cen':
                list_chi_n_v.append(snap_circ_cen(n))
            elif config=='cer_un_som':
                list_chi_n_v.append(snap_circ_ray(n))
            elif config=='compl':
                list_chi_n_v.append(snap_compl_ray(n))

        end=time.time()
        print('resolution EF : ',end-start,' secondes')

    ### utilise en 3D, pour le calcul parallele d'un snapshot individuel
    #
    ## enregistrement des donnees dans une liste
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
    l_name='Lchi_'+str(N_snap)+'_'+config+'_'+geo_p+'_deg'+str(VFS_degree)+'_'+ordo+'_'+computer
    # sauvegarde de la liste des solutions indexees calculees avec la methode des elements finis
    with sh.open(repertoire_parent+l_name) as l_sto:
        l_sto["maliste"] = list_chi_v
    # Matrice des snapshots : plus tard, voir l'etape II

else :
    l_name='Lchi_'+str(N_snap)+'_'+config+'_'+geo_p+'_deg'+str(VFS_degree)+'_'+ordo+'_'+computer
    with sh.open(repertoire_parent+l_name) as l_loa:
        list_chi_v = l_loa["maliste"]


# --------------------------------------------------------------------------------- #

# for rho in list_rho_appr:
for n in range(0,N_snap):
    # Extraction du snapshot de rang n
    chi_n_v=list_chi_v[n]
    # On cree un maillage pour reecrire les snapshots sous la forme de fonctions
    mesh_directory="maillages_per/2D/"

    rho=list_rho_appr[n]
    r = rho

    if config=='cer_un':
        if geo_p=='ray':
            cen=cen_snap_ray
        if typ_msh=='gms':
            mesh_name="maillage_trou2D_"+str(int(round(100*r,2)))
            mesh=Mesh(mesh_directory+mesh_name+".xml")
            plot(mesh)
            plt.tight_layout(pad=0)
            if fig_todo=='aff':
                plt.show()
            elif fig_todo=='save' and r==0.25:
                plt.savefig('Figures2D/maillage_gmsh_per_'+config+geo_p+'_ray'+str(int(round(100*r,2)))+'png')
            plt.close()
        else:
            mesh=creer_maill_circ(cen,r,res)
    elif config=='cer_un_som':
        r=n*0.05
        mention="_som"
        if typ_msh=='gms':
            mesh_name="maillage_trou2D"+mention+"_"+str(int(round(100*r,2)))
            mesh=Mesh(mesh_directory+mesh_name+".xml")
            plot(mesh)
            plt.tight_layout(pad=0)
            if fig_todo=='aff':
                plt.show()
            elif fig_todo=='save' and (r==0.25):
                plt.savefig('Figures2D/maillage_gmsh_per_'+config+geo_p+'_ray'+str(int(round(100*r,2)))+".png")
            plt.close()
        else:
            mesh=creer_maill_circ([c_x,c_y],r,res)

    else:
        r_fixe=0.15

        mesh_name="maillage_trous2D_"+geo_p+"_"+str(int(round(100*rho,2)))
        mesh=Mesh(mesh_directory+mesh_name+".xml")
        plot(mesh)
        plt.tight_layout(pad=0)
        if fig_todo=='aff':
            plt.show()
        elif fig_todo=='save' and (r==0.25):
            plt.savefig('Figures2D/maillage_gmsh_per_'+config+geo_p+'_ray'+str(int(round(100*r,2)))+".png")
        plt.close()

    V_n=VectorFunctionSpace(mesh, 'P', VFS_degree, constrained_domain=PeriodicBoundary())
    # On restitue la forme fonctionnelle du snapshot courant
    chi_n=Function(V_n)
    print(V_n.dim())
    print(len(chi_n_v))

    chi_n.vector().set_local(chi_n_v)
    # Figures et erreurs
    plot(chi_n)
    # if n==0 and rho_appr_min <1:
    #     plt.title("Rho = 0,0"+str(int(round(100*r,1))), fontsize=40)
    # else:
    #     plt.title("Rho = 0,"+str(int(round(100*r,2))),fontsize=40)
    #plot(grad(chi_n)[:,0]
    #plot(grad(chi_n)[:,1]
    if fig_todo=='aff':
        plt.show()
    elif fig_todo=='save':
        plt.savefig("Figures2D/sol_"+str(n+1)+"_sur"+str(N_snap)+config+'_'+geo_p+".png")
    plt.close()
    if config!='compl':
        err_per_gr(cen_snap_ray,r,chi_n,npas_err,fig_todo)
    else:
        err_per_gr(cen_snap_ray,r_fixe,chi_n,npas_err,fig_todo)
    ##
    # Tenseur de diffusion homogeneise
    ## Integrale de chi sur le domaine fluide
    T_chi=np.zeros((2,2))
    for k in range(0,2):
        for l in range(0,2):
            T_chi[k,l]=assemble(grad(chi_n)[k,l]*dx)
    ## Integrale de l'identite sur le domaine fluide
    if config!='compl':
        por=(1-pi*r**2)
    else:
        por=1-pi*(r**2+0.15**2)
    D=por*np.eye(2)
    ## Calcul et affichage du tenseur Dhom
    Dhom_k=D_k*(D+T_chi.T)
    print("Porosite :", por)
    print('Coefficient Dhom_k11, snapshot '+str(n+1)+", "+conf_mess+', '+geo_mess+" :",Dhom_k[0,0])
    ## Anisotropie
    mod_diag=min(abs(Dhom_k[0,0]),abs(Dhom_k[1,1]))
    mod_ndiag=0
    for i in range(0,2):
        for j in range(0,i):
            if abs(Dhom_k[i,j])>mod_ndiag:
                mod_ndiag=abs(Dhom_k[i,j])
        for j in range(i+1,2):
            if abs(Dhom_k[i,j])>mod_ndiag:
                mod_ndiag=abs(Dhom_k[i,j])
    print("Anisotropie I- : ",max(abs(Dhom_k[0,0]),abs(Dhom_k[1,1]))/mod_diag-1)
    print("Anisotropie II- : ",mod_ndiag/mod_diag)
    # Affichage des composantes scalaires : solution
    if config=='cer_un' and geo_p=='ray':
        fig_chi(cen_snap_ray,r,chi_n,fig_todo)
