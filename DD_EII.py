#####################################################################################################################################
######################################### Etape II : extrapolation des cliches, domaine_fixe ########################################
#####################################################################################################################################



# Chargement de la liste des snapshots physiques

l_name='Lchi_'+str(N_snap)+'_'+config+'_'+geo_p+'_deg'+str(VFS_degree)+'_'+ordo+'_'+computer

with sh.open(repertoire_parent+l_name) as l_loa:
    list_chi_v = l_loa["maliste"]

# Extrapolation au domaine Omega_fixe :
if dom_fixe=="am":
    mesh_fixe=Mesh("maillages_per/2D/maillage_fixe2D_am.xml")
elif dom_fixe=="multiray":
    mesh_fixe=Mesh("maillages_per/2D/maillage_fixe2d_"+dom_fixe+".xml")
elif config=='compl':
    mesh_fixe=Mesh("maillages_per/2D/maillage_trous2D_"+geo_p+"_fixe.xml")
elif dom_fixe=="ray_min":
    if config=='cer_un':
        mesh_fixe=Mesh('maillages_per/2D/maillage_trou2D_5.xml')

V_fixe = VectorFunctionSpace(mesh_fixe, 'P', VFS_degree, constrained_domain=PeriodicBoundary())

plot(mesh_fixe)
if fig_todo=='aff':
    plt.show()
plt.close()
# Extrapolation des solutions du probleme physique

def extra_snap(n):

    r = list_rho_appr[n]

    # chargement du snapshot courant
    chi_n_v=list_chi_v[n]

    # mise sous forme d'une fonction EF
    if config!='compl':
        if cen_snap_ray==[0.,0.]:
            mention="_som"
        else:
            mention=""
    if typ_msh=='gms':
        mesh_name="maillages_per/2D/maillage_trou2D"+mention+"_"+str(int(round(100*r,2)))+".xml"
        print(mesh_name)
        mesh=Mesh(mesh_name)
    else:
        mesh_name="maillages_per/2D/maillage_trous2D_"+geo_p+"_"+str(int(round(100*r,2)))+".xml"
        mesh=Mesh(mesh_name)
    V_n=VectorFunctionSpace(mesh, 'P', VFS_degree, constrained_domain=PeriodicBoundary())
    chi_n=Function(V_n)
    chi_n.vector().set_local(chi_n_v)
    # extrapolation du snapshot au domaine fixe
    chi_n.set_allow_extrapolation(True)
    chi_n_prime=interpolate(chi_n,V_fixe)
    ## on range le snapshot dans une liste
    #list_snap.append(chi_n_prime)
    chi_n_prime_v=chi_n_prime.vector().get_local()
    return([n,chi_n_prime_v])


# Constitution de la matrice des snapshots

if not exsnap_done:
    # production de la matrice des  snapshots

    # preparation du calcul parallele gros grains
    pool=multiprocessing.Pool(processes=8)
    list_chi_n_prime_v=pool.map(extra_snap,(n for n in range(0,N_snap)))

    # Remplissage de la matrice
    nb_noeuds=V_fixe.dim()
    Usnap=np.zeros((nb_noeuds,N_snap))
    for n in range(0,N_snap):
        for i in range(0,N_snap):
            if list_chi_n_prime_v[i][0]==n:
                Usnap[:,n]=list_chi_n_prime_v[i][1]

    # Stochage de la matrice des snapshots
    u_name='Usnap_'+dom_fixe+str(N_snap)+'_'+config+'_'+geo_p+'_deg'+str(VFS_degree)+'_'+ordo+'_'+computer
    #
    with sh.open(repertoire_parent+u_name) as u_sto:
        u_sto["maliste"] = Usnap

else:
    # Chargement de la matrice des snapshots
    u_name='Usnap_'+dom_fixe+str(N_snap)+'_'+config+'_'+geo_p+'_deg'+str(VFS_degree)+'_'+ordo+'_'+computer
    print(repertoire_parent+u_name)
    with sh.open(repertoire_parent+u_name) as u_loa:
        Usnap = u_loa["maliste"]



# Representations graphiques

list_snap=[]
for n in range(0,N_snap):
    chi_prime=Function(V_fixe)
    chi_prime.vector().set_local(Usnap[:,n])
    # remplissage de la liste de fonctions
    list_snap.append(chi_prime)

cen=cen_snap_ray
for n in range(0,N_snap):

    r = list_rho_appr[n]

    chi_prime_n=list_snap[n]

    if fig_todo=='aff' or fig_todo=='save':
        # Affichage des valeurs de la solution interpolee
        plot(chi_prime_n)
        plt.title("Snapshot "+str(n+1),fontsize=40)
        if fig_todo=='aff':
            plt.show()
        else:
            plt.savefig("Figures2D/snap_"+str(n+1)+"_sur"+str(N_snap)+config+'_'+geo_p+".png")
        plt.close()
    else:
        print('pffrrh !')

    # Affichage des composantes scalaires : interpolee
    if config=='cer_un' and geo_p=='ray' and fig_todo!='':
        fig_chi(cen_snap_ray,r,chi_prime_n,fig_todo)

# # test possible de Dhom calcule directement sur des interpolees : a ecrire sur un autre fichier ?
# if not test_Dhom:
#     sys.exit('fin de l etape II, sans tests d integration sur des sousdomaines')#-------------------------------------------------------------------------------------------------------------------------------------------
