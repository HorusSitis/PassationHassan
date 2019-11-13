#####################################################################################################################################
######################################### Etape II : extrapolation des cliches, domaine_fixe ########################################
#####################################################################################################################################

### ------------ Reproduire eventuellement pour des etapes ulterieures. Laisser seulement dans DD_fun_obj ? ------------ ###

import time

tol=1e-10

xinf=0.0
yinf=0.0
xsup=1.0
ysup=1.0

#determiner le domaine fixe pour interpoler la solution

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

# Chargement de la liste des snapshots physiques

l_name='Lchi_'+str(Nsnap)+'_'+config+'_'+geo_p+'_deg'+str(VFS_degree)+'_'+ordo+'_'+computer

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

    if geo_p=='hor':
        r=0.01+0.04*n
    else:
        r=0.05*n

    # chargement du snapshot courant
    chi_n_v=list_chi_v[n-1]

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
    list_chi_n_prime_v=pool.map(extra_snap,(n for n in range(1,1+Nsnap)))

    # Remplissage de la matrice
    nb_noeuds=V_fixe.dim()
    Usnap=np.zeros((nb_noeuds,Nsnap))
    for n in range(1,1+Nsnap):
        for i in range(0,Nsnap):
            if list_chi_n_prime_v[i][0]==n:
                Usnap[:,n-1]=list_chi_n_prime_v[i][1]

    # Stochage de la matrice des snapshots
    u_name='Usnap_'+dom_fixe+str(Nsnap)+'_'+config+'_'+geo_p+'_deg'+str(VFS_degree)+'_'+ordo+'_'+computer
    #
    with sh.open(repertoire_parent+u_name) as u_sto:
        u_sto["maliste"] = Usnap

else:
    # Chargement de la matrice des snapshots
    u_name='Usnap_'+dom_fixe+str(Nsnap)+'_'+config+'_'+geo_p+'_deg'+str(VFS_degree)+'_'+ordo+'_'+computer
    print(repertoire_parent+u_name)
    with sh.open(repertoire_parent+u_name) as u_loa:
        Usnap = u_loa["maliste"]



# Representations graphiques

list_snap=[]
for n in range(1,1+Nsnap):
    chi_prime=Function(V_fixe)
    chi_prime.vector().set_local(Usnap[:,n-1])
    # remplissage de la liste de fonctions
    list_snap.append(chi_prime)

cen=cen_snap_ray
for n in range(1,1+Nsnap):

    if geo_p=='hor':
        r=0.01+0.04*n
    else:
        r=0.05*n
    chi_prime_n=list_snap[n-1]

    if fig_todo=='aff' or fig_todo=='save':
        # Affichage des valeurs de la solution interpolee
        plot(chi_prime_n)
        plt.title("Snapshot "+str(n),fontsize=40)
        if fig_todo=='aff':
            plt.show()
        else:
            plt.savefig("Figures2D/snap_"+str(n)+"_sur"+str(Nsnap)+config+'_'+geo_p+".png")
        plt.close()
    else:
        print('pffrrh !')

    # Affichage des composantes scalaires : interpolee
    if config=='cer_un' and geo_p=='ray' and fig_todo!='':
        fig_chi(cen_snap_ray,r,chi_prime_n,fig_todo)

# test possible de Dhom calcule directement sur des interpolees : a ecrire sur un autre fichier ?
if not test_Dhom:
    sys.exit('fin de l etape II, sans tests d integration sur des sousdomaines')#-------------------------------------------------------------------------------------------------------------------------------------------
