#####################################################################################################################################
######################################### Etape II : extrapolation des cliches, domaine fixe ########################################
#####################################################################################################################################

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


if exsnap_done:
    # Chargement de la matrice des snapshots
    u_name='Usnap_'+dom_fixe+str(N_snap)+'_'+config+'_'+geo_p+'_deg'+str(VFS_degree)+'_'+ordo+'_'+computer
    print(repertoire_parent+u_name)
    with sh.open(repertoire_parent+u_name) as u_loa:
        Usnap = u_loa["maliste"]
else:
    sys.exit('snapshots non encore generes')


# Representations graphiques

list_snap=[]
for n in range(0,N_snap):
    chi_prime=Function(V_fixe)
    chi_prime.vector().set_local(Usnap[:,n])
    # remplissage de la liste de fonctions
    list_snap.append(chi_prime)

cen=cen_snap_ray
for n in range(0,N_snap):

    if geo_p=='hor':
        r=0.01+0.04*n
    else:
        r=0.05*n
    chi_prime_n=list_snap[n]

    if fig_todo=='aff' or fig_todo=='save':
        # Affichage des valeurs de la solution interpolee
        plot(chi_prime_n)
        plt.title("Snapshot "+str(n),fontsize=40)
        if fig_todo=='aff':
            plt.show()
        else:
            plt.savefig("Figures2D/snap_"+str(n)+"_sur"+str(N_snap)+config+'_'+geo_p+".png")
        plt.close()
    else:
        print('pffrrh !')

    # Affichage des composantes scalaires : interpolee
    if config=='cer_un' and geo_p=='ray' and fig_todo!='':
        fig_chi(cen_snap_ray,r,chi_prime_n,fig_todo)










# test possible de Dhom calcule directement sur des interpolees : a ecrire sur un autre fichier ?
if not test_Dhom:
    sys.exit('fin de l etape II, sans tests d integration sur des sousdomaines')#-------------------------------------------------------------------------------------------------------------------------------------------

nb_noeuds_dom_fixe=V_fixe.dim()

n_mp_refi=8
fig_mesh=True

def width(i):
    return crow/i**2

def f_testDhom(n):

    rho=list_rho_appr[n]
    r = rho
    chi_prime_n=list_snap[n]

    por_prime=1

    ## Creation du maillage fixe local
    if dom_fixe=="am":
        mesh_fixe=Mesh("maillages_per/2D/maillage_fixe2D_am.xml")
    elif config=='compl':
        mesh_fixe=Mesh("maillages_per/2D/maillage_trous2D_"+geo_p+"_fixe.xml")
    elif dom_fixe=="ray_min":
        if config=='cer_un':
            mesh_fixe=Mesh('maillages_per/2D/maillage_trou2D_5.xml')
    elif dom_fixe=="multiray":
        mesh_fixe=Mesh("maillages_per/2D/maillage_fixe2d_"+dom_fixe+".xml")
    V_fixe=VectorFunctionSpace(mesh_fixe, 'P', VFS_degree, constrained_domain=PeriodicBoundary())

    ## Interpolation post-snapshot sur le maillage physique
    mesh_name="maillages_per/2D/maillage_trou2D"+mention+"_"+str(int(round(100*r,2)))+".xml"
    mesh=Mesh(mesh_name)
    V_n=VectorFunctionSpace(mesh,'P',VFS_degree,constrained_domain=PeriodicBoundary())
    chi_postprime_n=Function(V_n)
    chi_prime_n.set_allow_extrapolation(True)
    chi_postprime_n=interpolate(chi_prime_n,V_n)

    ## calcul de int_grad
    T_chi_postprime=np.zeros((2,2))
    for k in range(0,2):
        for l in range(0,2):
            T_chi_postprime[k,l]=assemble(grad(chi_postprime_n)[k,l]*dx)
    ## Integration sur le domaine fluide avec le maillage fixe
    print("Avec le maillage fixe")

    # Raffinement du maillage
    start=time.time()
    # Boucle pour des rafinements successifs
    for i in range(1,1+Nrefine):
        print('raffinement',i)
        markers = MeshFunction("bool", mesh_fixe, mesh_fixe.topology().dim())
        markers.set_all(False)
        for c in cells(mesh_fixe):
            if typ_refi=='vol':
                for f in facets(c):
                    if ((f.midpoint()[0]-cen_snap_ray[0])**2+(f.midpoint()[1]-cen_snap_ray[1])**2<=(r*(1+width(i)))**2) and ((f.midpoint()[0]-cen_snap_ray[0])**2+(f.midpoint()[1]-cen_snap_ray[1])**2>=(r*(1-width(i)))**2):
                        markers[c]=True
                for v in vertices(c):
                    if (v.point().x()-cen_snap_ray[0])**2+(v.point().y()-cen_snap_ray[1])**2>=(r*(1-width(i)))**2 and (v.point().x()-cen_snap_ray[0])**2+(v.point().y()-cen_snap_ray[1])**2<=(r*(1+width(i)))**2:
                        markers[c]=True
            elif typ_refi=='front':
                list_sgn=[]
                for v in vertices(c):
                    if (v.point().x()-cen_snap_ray[0])**2+(v.point().y()-cen_snap_ray[1])**2>r**2:
                        list_sgn.append(1)
                    elif (v.point().x()-cen_snap_ray[0])**2+(v.point().y()-cen_snap_ray[1])**2==r**2:
                        list_sgn.append(0)
                    else:
                        list_sgn.append(-1)
                # on marque les cellules qui coupent la frontière du domaine fluide virtuel
                if list_sgn==[1,1,1] or list_sgn==[-1,-1,-1]:
                    markers[c]=False
                else:
                    markers[c]=True
        # raffinement à l'étape i
        mesh_fixe=refine(mesh_fixe, markers, redistribute=True)
    end=time.time()
    tps_refi=end-start

    # Interpolation sur le maillage raffine
    start=time.time()
    if Nrefine>0:
        print('interpolation rayon',str(int(round(100*r,2))))
        V_fixe=VectorFunctionSpace(mesh_fixe, 'P', VFS_degree, constrained_domain=PeriodicBoundary())
        chi_prime_n.set_allow_extrapolation(True)
        chi_prime_n=interpolate(chi_prime_n,V_fixe)
    #
    nb_noeuds_refi=V_fixe.dim()
    end=time.time()
    tps_interp=end-start
    if Nrefine>0 and fig_mesh and n==2:
        plot(mesh_fixe)
        plt.title('Maillage raffine '+str(Nrefine)+' fois rayon '+str(int(round(100*r,2)))+'x10e-2')
        plt.show()
        plt.close()

    # Creation du domaine d'integration sur le maillage fixe
    class DomPhysFluide(SubDomain):
        def inside(self, x, on_boundary):
            return True if ((x[0]-0.5)**2+(x[1]-0.5)**2>=r**2) else False
    dom_courant=DomPhysFluide()
    subdomains=MeshFunction('size_t',mesh_fixe,mesh_fixe.topology().dim())
    subdomains.set_all(1)
    dom_courant.mark(subdomains,12829)
    dxf=Measure("dx", domain=mesh_fixe, subdomain_data=subdomains)

    # Calcul de int_grad en restreignant sur plusieurs domaines
    T_chi_restr_prime=np.zeros((2,2))
    co_T_chi_restr_prime=np.zeros((2,2))
    for k in range(0,2):
        for l in range(0,2):
            T_chi_restr_prime[k,l]=assemble(grad(chi_prime_n)[k,l]*dxf(12829))
            co_T_chi_restr_prime[k,l]=assemble(grad(chi_prime_n)[k,l]*dxf(1))

    return([tps_refi,tps_interp,T_chi_restr_prime,co_T_chi_restr_prime,T_chi_postprime,nb_noeuds_refi])

print('Liste des snapshots :',u_name)

pool=multiprocessing.Pool(processes=n_mp_refi)
list_refi_interp=pool.map(f_testDhom,(n for n in range(0,N_snap)))

nom_fichier='Res2D/'+'testDhom'+config+geo_p+'_'+dom_fixe+'_crown10'+str(lg_crow)+'_Nrefi'+str(Nrefine)+'.txt'
registre=open(nom_fichier,'w')

for n in range(0,N_snap):

    # calcul de rho d'apres les parametres geometriques exploitables avec multiprocessing
    rho=list_rho_appr[n]
    r = rho

    if geo_p=='ray':
        por=1-pi*r**2
    por_prime=1

    ## Interpolation post-snapshot sur le maillage physique
    mesh_name="maillages_per/2D/maillage_trou2D"+mention+"_"+str(int(round(100*r,2)))+".xml"
    mesh=Mesh(mesh_name)
    V_n=VectorFunctionSpace(mesh,'P',VFS_degree,constrained_domain=PeriodicBoundary())
    chi_postprime_n=Function(V_n)
    chi_prime_n.set_allow_extrapolation(True)
    chi_postprime_n=interpolate(chi_prime_n,V_n)

    T_chi_postprime=np.zeros((2,2))
    for k in range(0,2):
        for l in range(0,2):
            T_chi_postprime[k,l]=assemble(grad(chi_postprime_n)[k,l]*dx)

    ## Chargement des resultats d'integration des tenseurs d'homogeneisation sur les differents maillages
    [tps_refi,tps_interp,T_chi_restr_prime,co_T_chi_restr_prime,T_chi_postprime,nb_noeuds_refi]=list_refi_interp[n]

    ## Integration sur le domaine virtuel : carre unite eventuellement prive d'un point
    T_chi_prime=np.zeros((2,2))
    for k in range(0,2):
        for l in range(0,2):
            T_chi_prime[k,l]=assemble(grad(chi_prime_n)[k,l]*dx)

    ## les trois coefficients a comparer
    ###
    D=por*np.eye(2)
    integr_k_postprime=D_k*(D+T_chi_postprime.T)
    Dhom_k_postprime=integr_k_postprime*(1/por_prime)#por en dimensionnel
    ###
    integr_k_restr_prime=D_k*(D+T_chi_restr_prime.T)
    Dhom_k_restr_prime=integr_k_restr_prime*(1/por_prime)
    ###
    D_moy=D_k*(D+0.5*(T_chi_restr_prime.T-co_T_chi_restr_prime.T))*(1/por_prime)
    ###
    D_prime=por_prime*np.eye(2)
    integr_k_prime=D_k*(D_prime+T_chi_prime.T)
    Dhom_k_prime=integr_k_prime*(1/por_prime)
    ##
    print('##################################################################')
    print('Geometrie : '+conf_mess+', '+geo_mess+', '+str(int(round(100*r,2))))
    print('Domaine fixe : '+dom_fixe)
    print('Noeuds du maillage fixe :',nb_noeuds_dom_fixe)
    print('------------------------------------------------------------------')
    print(refi_mess)
    print('Nombre de tours de raffinement :',Nrefine)
    print('Noeuds du maillage raffine :',nb_noeuds_refi)
    print('Raffinemenent du maillage :',tps_refi,'secondes')
    print('Interpolation sur le maillage raffine :',tps_interp,'secondes')
    print('------------------------------------------------------------------')
    print()
    print('DUsnap fixe restreint au domaine courant :',Dhom_k_restr_prime[0,0])#,'porosite',por)
    print('DUsnap physique :',Dhom_k_postprime[0,0])#,'porosite',por)
    print('DUsnap domaine fixe moyenne :',D_moy[0,0])
    print('------------------------------------------------------------------')
    print('Erreur relative :',100*(Dhom_k_restr_prime[0,0]-Dhom_k_postprime[0,0])/Dhom_k_postprime[0,0],'pourcent')
    print('##################################################################')
    #
    registre.write('##################################################################'+'\n')
    registre.write('Rho : '+str(int(round(100*r,2)))+'\n')
    registre.write('Domaine fixe : '+dom_fixe+'\n')
    registre.write('Noeuds du maillage fixe : '+str(nb_noeuds_dom_fixe)+'\n')
    registre.write('------------------------------------------------------------------'+'\n')
    registre.write(refi_mess+'\n')
    registre.write('Nombre de tours de raffinement : '+str(Nrefine)+'\n')
    registre.write('Noeuds du maillage raffine : '+str(nb_noeuds_refi)+'\n')
    registre.write('Raffinemenent du maillage : '+str(tps_refi)+' secondes'+'\n')
    registre.write('Interpolation sur le maillage raffine : '+str(tps_interp)+' secondes'+'\n')
    registre.write('------------------------------------------------------------------'+'\n')
    #registre.write('DUsnap fixe : '+str(Dhom_k_prime[0,0],'porosite fixe'+str(por_prime)+'\n')
    registre.write('Tenseur virtuel : '+str(T_chi_prime)+'\n')
    registre.write('DUsnap fixe restreint au domaine courant : '+str(Dhom_k_restr_prime[0,0])+'\n')
    registre.write('DUsnap physique : '+str(Dhom_k_postprime[0,0])+'\n')
    registre.write('DUsnap domaine fixe moyenne : '+str(D_moy[0,0])+'\n')
    registre.write('------------------------------------------------------------------'+'\n')
    registre.write('Erreur relative : '+str(100*(Dhom_k_restr_prime[0,0]-Dhom_k_postprime[0,0])/Dhom_k_postprime[0,0])+' pourcent'+'\n')
    registre.write('##################################################################'+'\n')


registre.close()
