#####################################################################################################################################
######################################### Etape II : extrapolation des cliches, domaine_fixe ########################################
#####################################################################################################################################


# mesh_dir='maillages_per/3D/'
# mesh_repository = 'maillages_per/' + str(dimension) + 'D/'

### ------------ Nommage du domaine fixe ------------ ###

# if config == 'sph_un':
#     mesh_prefix = 'cubesphere_periodique_triangle_'
# elif config == '2sph':
#     mesh_prefix = 'cube2sph_periodique_triangle_'
# elif config == 'cylsph':
#     mesh_prefix = 'cubecylsph_periodique_triangle_'

## inclusions simples ou rayons lies

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

### ------------ Generation du maillage fixe avec Gmsh ------------ ###





if not mesh_ex_done:

    xinf = xyzinfsup[0][0]
    yinf = xyzinfsup[0][1]
    zinf = xyzinfsup[0][2]
    xsup = xyzinfsup[1][0]
    ysup = xyzinfsup[1][1]
    zsup = xyzinfsup[1][2]

    xcen = (xinf + xsup)/2.
    ycen = (yinf + ysup)/2.
    zcen = (zinf + zsup)/2.

    ## on selectionne le nom du fichier canevas

    if config == 'sph_un':
        fichier_sansentete = open(mesh_repository + 'cube_periodique_triangle' + '_' + dom_fixe + '_sansgpar' + '.txt', 'r')
    elif config == '2sph':
        fichier_sansentete = open(mesh_repository + mesh_prefix + 'fixe_som' + '_sansgpar' + '.txt', 'r')
    elif config == 'cylsph':
        fichier_sansentete = open(mesh_repository + mesh_prefix + 'sansgpar' + '.txt', 'r')

    ## fichier a editer et remplacer par un script .geo

    gen_mesh = open(mesh_repository + mesh_f_name + '.txt', 'w')

    if config == 'sph_un' and dom_fixe == 'am':
        print('Pas de rayon')

    elif config == '2sph' and dom_fixe == 'solid':
        gen_mesh.write('S = ' + str(ray_fix) + ';' + '\n' + '\n')

    elif config == 'cylsph' and dom_fixe == 'ray_min':
        # centre de la maille
        gen_mesh.write('xc = ' + str(xcen) + ';' + '\n')
        gen_mesh.write('yc = ' + str(ycen) + ';' + '\n')
        gen_mesh.write('zc = ' + str(zcen) + ';' + '\n' + '\n')
        # rayons
        if geo_p == 'ray_sph':
            gen_mesh.write('R = ' + str(rho_appr_min) + ';' + '\n')
            gen_mesh.write('S = ' + str(ray_fix) + ';' + '\n' + '\n')
        elif geo_p == 'ray_cyl':
            gen_mesh.write('R = ' + str(ray_fix) + ';' + '\n')
            gen_mesh.write('S = ' + str(rho_appr_min) + ';' + '\n' + '\n')

    # bornes du domaine

    gen_mesh.write('xmin = ' + str(xinf) + ';' + '\n')
    gen_mesh.write('ymin = ' + str(yinf) + ';' + '\n')
    gen_mesh.write('zmin = ' + str(zinf) + ';' + '\n')
    gen_mesh.write('xmax = ' + str(xsup) + ';' + '\n')
    gen_mesh.write('ymax = ' + str(ysup) + ';' + '\n')
    gen_mesh.write('zmax = ' + str(zsup) + ';' + '\n' + '\n')

    # resolution

    gen_mesh.write('step = ' + str(round(1./res, 2)) + ';' + '\n' + '\n')

    # utilisation du canevas

    for line in fichier_sansentete :
        gen_mesh.write(line)

    gen_mesh.close()
    fichier_sansentete.close()

    ### Conversion du fichier texte obtenu et stockage du maillage en fichiers .xml
    os.rename(mesh_repository + mesh_f_name + '.txt', mesh_repository + mesh_f_name + '.geo')

    ### on genere le maillage avec gmsh
    os.system('gmsh -' + str(dimension) + ' ' + mesh_repository + mesh_f_name+ '.geo')

    if fig_todo == 'aff':
        ## affichage du maillage genere selon les parametres d'entree
        os.system('gmsh ' + mesh_repository + mesh_f_name + '.msh')

    ### on convertit le maillage avec dolfin
    os.system('dolfin-convert ' + mesh_repository + mesh_f_name + '.msh' + ' ' + mesh_repository + mesh_f_name + '.xml')




# sys.exit()







### ------------ Generation du maillage fixe avec FEniCS et definition de VFS ------------ ###

print('Maillage fixe :', mesh_repository + mesh_f_name)
mesh_fixe = Mesh(mesh_repository + mesh_f_name + '.xml')

# fonctions test du domaine fixe

V_fixe = VectorFunctionSpace(mesh_fixe,'P',2,constrained_domain=PeriodicBoundary())

### ------------ Etapes reproduites : dependances directes de Main3D ------------ ###

# Chargement de la liste des snapshots physiques

# l_name='Lchi_' + str(N_snap) + '_'+config + '_' + geo_p + '_' + 'sur' + str(res) + '_' + ordo + '_' + computer
print(l_name)
with sh.open(repertoire_parent + l_name) as l_loa:
    list_chi_v = l_loa['maliste']


# Extrapolation au domaine Omega_fixe : inclusion spherique de rayon 0.0001, chi_prime defini sur ce domaine

def extra_snap(n):
    # valeur de rho^j
    r = list_rho_appr[n]
    # chargement du snapshot courant
    chi_n_v=list_chi_v[n]
    # mise sous forme d'une fonction EF


    # generation du maillage avec FEniCS

    if config == 'sph_un':
        nom_fichier_avecgpar = mesh_prefix + 'rayc' + str(int(round(100*r,2))) + '_sur' + str(res)
    elif config == '2sph':
        nom_fichier_avecgpar = mesh_prefix + 'rayc' + str(int(round(100*r,2))) + '_rayp' + str(int(round(100*ray_fix,2))) + '_sur' + str(res)
    elif config == 'cylsph' and geo_p == 'ray_sph':
        nom_fichier_avecgpar = mesh_prefix + 'rayc' + str(int(round(100*r,2))) + '_rayp' + str(int(round(100*ray_fix,2))) + '_sur' + str(res)
    elif config == 'cylsph' and geo_p == 'ray_cyl':
        nom_fichier_avecgpar = mesh_prefix + 'rayc' + str(int(round(100*ray_fix,2))) + '_rayp' + str(int(round(100*r,2))) + '_sur' + str(res)
    mesh_name = nom_fichier_avecgpar

    mesh = Mesh(mesh_repository + mesh_name + '.xml')

    # exploitation du maillage
    V_n = VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
    chi_n = Function(V_n)
    chi_n.vector().set_local(chi_n_v)
    # extrapolation du snapshot au domaine fixe
    chi_n.set_allow_extrapolation(True)
    chi_n_prime = interpolate(chi_n, V_fixe)
    # verification sur chi_n
    ###
    ## on range le snapshot dans une liste
    chi_n_prime_v = chi_n_prime.vector().get_local()

    ##
    return([n, chi_n_prime_v])



# Constitution de la matrice des snapshots

if not exsnap_done:

    pool=multiprocessing.Pool(processes=8)
    list_chi_n_prime_v=pool.map(extra_snap, (n for n in range(0, N_snap)))

    # Remplissage de la matrice
    nb_noeuds = V_fixe.dim()
    Usnap = np.zeros((nb_noeuds, N_snap))
    for n in range(0, N_snap):
        for i in range(0, N_snap):
            if list_chi_n_prime_v[i][0] == n:
                Usnap[:, n]=list_chi_n_prime_v[i][1]
    # Stockage de la matrice des snapshots
    # u_name = 'Usnap_' + dom_fixe + '_' + str(N_snap) + '_' + config + '_' + geo_p + '_' + 'res' + str(res) + '_' + ordo + '_' + computer
    #
    with sh.open(repertoire_parent+u_name) as u_sto:
        u_sto['maliste'] = Usnap
else:
    # Chargement de la matrice des snapshots
    # u_name = 'Usnap_' + dom_fixe + '_' + str(N_snap) + '_' + config + '_' + geo_p + '_' + 'res' + str(res) + '_' + ordo + '_' + computer
    with sh.open(repertoire_parent + u_name) as u_loa:
        Usnap = u_loa['maliste']


# Representations graphiques

list_snap = []
for n in range(0, N_snap):
    chi_prime=Function(V_fixe)
    chi_prime.vector().set_local(Usnap[:, n - 1])
    # remplissage de la liste de fonctions
    list_snap.append(chi_prime)

#cen=cen_snap_ray
for n in range(0, N_snap):

    chi_prime_n = list_snap[n]
    r = list_rho_appr[n]

    porosity = epsilon_p(r, config, geo_p, ray_fix)
    if config == 'sph_un' or config == 'cube2sph':
        por_fix = epsilon_p(0., config, geo_p, ray_fix)
    elif config == 'cylsph':
        por_fix = epsilon_p(rho_appr_min, config, geo_p, ray_fix)

    if config == 'sph_un':
        nom_fichier_avecgpar = mesh_prefix + 'rayc' + str(int(round(100*r,2))) + '_sur' + str(res)
    elif config == '2sph':
        nom_fichier_avecgpar = mesh_prefix + 'rayc' + str(int(round(100*r,2))) + '_rayp' + str(int(round(100*ray_fix,2))) + '_sur' + str(res)
    elif config == 'cylsph' and geo_p == 'ray_sph':
        nom_fichier_avecgpar = mesh_prefix + 'rayc' + str(int(round(100*r,2))) + '_rayp' + str(int(round(100*ray_fix,2))) + '_sur' + str(res)
    elif config == 'cylsph' and geo_p == 'ray_cyl':
        nom_fichier_avecgpar = mesh_prefix + 'rayc' + str(int(round(100*ray_fix,2))) + '_rayp' + str(int(round(100*r,2))) + '_sur' + str(res)
    mesh_name = nom_fichier_avecgpar

    # Affichage des valeurs de la solution interpolee
    plot(chi_prime_n, linewidth=lw)
    plt.title('Snapshot ' + str(n + 1), fontsize=30)
    if fig_todo == 'aff':
        plt.show()
    else:
        plt.savefig('Figures3D/snap_interp'+dom_fixe+'_'+str(n)+'_sur'+str(N_snap)+config+'_'+geo_p+'res'+str(res)+'.png')
    plt.close()

    # --------- Affichage des valeurs et erreurs de la solution periodique, quelle que soit la configuration --------- #

    ## ------------------------ Configurations complexes ------------------------ ##
    if test_Dhom and (config=='cylsph' or config=='2sph'):
        mesh = Mesh(mesh_repository + mesh_name + '.xml')
        ##
        V_n=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
        chi_n=Function(V_n)
        chi_prime_n.set_allow_extrapolation(True)
        chi_n=interpolate(chi_prime_n,V_n)
        ##
        T_chi=np.zeros((3,3))
        for k in range(0,3):
            for l in range(0,3):
                T_chi[k,l]=assemble(grad(chi_n)[k,l]*dx)
        T_chi_omega = T_chi/cell_vol
        ##
        # if geo_p=='ray':
        #     class DomPhysFluide(SubDomain):
        #         def inside(self, x, on_boundary):
        #             return True if (x[0]**2+x[1]**2+x[2]**2>=r**2) else False
        # elif geo_p=='ray_sph':
        #     class DomPhysFluide(SubDomain):
        #         def inside(self, x, on_boundary):
        #             return True if (x[0]**2+x[1]**2+x[2]**2>=r_s**2) else False
        # elif geo_p=='ray_cyl':
        #     class DomPhysFluide(SubDomain):
        #         def inside(self, x, on_boundary):
        #             return True if (x[0]**2+x[2]**2>=r_c**2 and (1-x[0])**2+x[2]**2>=r_c**2 and x[0]**2+(1-x[2])**2>=r_c**2 and (1-x[0])**2+(1-x[2])**2>=r_c**2) else False
        # dom_courant = DomPhysFluide()
        # subdomains = MeshFunction('size_t', mesh_fixe, mesh_fixe.topology().dim())
        # subdomains.set_all(1)
        # dom_courant.mark(subdomains,12829)
        # dxf = Measure('dx', domain=mesh_fixe, subdomain_data=subdomains)
        # T_chi_restr_prime=np.zeros((3,3))
        # for k in range(0,3):
        #     for l in range(0,3):
        #         T_chi_restr_prime[k,l]=assemble(grad(chi_prime_n)[k,l]*dxf(12829))
        # T_chi_restr_prime_omega=T_chi_restr_prime/cell_vol
        # ##
        # T_chi_prime=np.zeros((3,3))
        # for k in range(0,3):
        #     for l in range(0,3):
        #         T_chi_prime[k,l]=assemble(grad(chi_prime_n)[k,l]*dxf)
        # T_chi_prime_omega=T_chi_prime/cell_vol
        ## les trois coefficients a comparer
        ###
        D = porosity*np.eye(3)
        integr_k_postprime = D_k*(D+T_chi_omega.T)
        Dhom_k_postprime = integr_k_postprime#*(1/por_fix)
        # ###
        # D = porosity*np.eye(3)
        # integr_k_restr_prime = D_k*(D+T_chi_restr_prime_omega.T)
        # Dhom_k_restr_prime = integr_k_restr_prime#*(1/por_fix)
        # ###
        # D_prime = por_fix*np.eye(3)
        # integr_k_prime = D_k*(D_prime+T_chi_prime_omega.T)
        # Dhom_k_prime = integr_k_prime#*(1/por_fix)
        ##
    ## ------------------------ Cylindre unique, test ------------------------ ##
    elif test_Dhom and config == 'cyl_un':
        mesh_name = 'cubecylindre_periodique_triangle_'+str(int(round(100*r,2)))+'sur'+str(res)
        mesh = Mesh(mesh_repository + mesh_name + '.xml')
        ##
        V_n=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
        chi_n=Function(V_n)
        chi_prime_n.set_allow_extrapolation(True)
        chi_n=interpolate(chi_prime_n,V_n)
        ##
        T_chi=np.zeros((3,3))
        for k in range(0,3):
            for l in range(0,3):
                T_chi[k,l]=assemble(grad(chi_n)[k,l]*dx)
        T_chi_omega = T_chi/cell_vol
        ## a remplacer ##
        porosity = epsilon_p(r, config, geo_p, r_f)
        D=porosity*np.eye(3)
        integr_k_postprime=D_k*(D+T_chi_omega.T)
        Dhom_k_postprime=integr_k_postprime#*(1/por_fix)#por en dimensionnel

    ## ------------------------ Sphere unique, test ------------------------ ##
    elif test_Dhom and config == 'sph_un':
        ##
        # mesh_name = 'cubesphere_periodique_triangle_'+str(int(round(100*r,2)))+'sur'+str(res_gmsh)
        print('Verification : maillage courant', mesh_name)
        mesh = Mesh(mesh_repository + mesh_name + '.xml')
        ##
        V_n = VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
        chi_n=Function(V_n)
        chi_prime_n.set_allow_extrapolation(True)
        chi_n=interpolate(chi_prime_n, V_n)
        ##
        T_chi=np.zeros((3,3))
        for k in range(0,3):
            for l in range(0,3):
                T_chi[k,l]=assemble(grad(chi_n)[k,l]*dx)
        T_chi_omega=T_chi/cell_vol
        ##
        class DomPhysFluide(SubDomain):
            def inside(self, x, on_boundary):
                return True if (x[0]**2+x[1]**2+x[2]**2>=r**2) else False
        dom_courant = DomPhysFluide()
        subdomains = MeshFunction('size_t',mesh_fixe,mesh_fixe.topology().dim())
        subdomains.set_all(1)
        dom_courant.mark(subdomains,12829)
        ##
        dxf = Measure('dx', domain = mesh_fixe, subdomain_data=subdomains)
        T_chi_restr_prime=np.zeros((3,3))
        for k in range(0,3):
            for l in range(0,3):
                T_chi_restr_prime[k,l]=assemble(grad(chi_prime_n)[k,l]*dxf(12829))
        T_chi_restr_prime_omega=T_chi_restr_prime/cell_vol
        ##
        T_chi_prime=np.zeros((3,3))
        for k in range(0,3):
            for l in range(0,3):
                T_chi_prime[k,l]=assemble(grad(chi_prime_n)[k,l]*dxf)
        T_chi_prime_omega=T_chi_prime/cell_vol
        ## les trois coefficients a comparer
        ###
        D = porosity*np.eye(3)
        integr_k_postprime = D_k*(D+T_chi_omega.T)
        Dhom_k_postprime = integr_k_postprime#*(1/por_fix)
        ###
        D = porosity*np.eye(3)
        integr_k_restr_prime = D_k*(D+T_chi_restr_prime_omega.T)
        Dhom_k_restr_prime = integr_k_restr_prime#*(1/por_fix)
        ###
        D_prime = por_fix*np.eye(3)
        integr_k_prime = D_k*(D_prime+T_chi_prime_omega.T)
        Dhom_k_prime = integr_k_prime#*(1/por_fix)
        ##

    # --------- Sortie en shell : resultats en IG et Dhom des differentes methodes d'integration --------- #
    if test_Dhom:

        print('#'*78)
        print('Geometrie : ' + conf_mess + ', ' + geo_mess + ', ' + str(int(round(100*r, 2))))
        print('%'*78)
        if config == 'sph_un':
            print('Porosite fixe', por_fix)
            print('DUsnap fixe :', Dhom_k_prime[0,0])
            print('-'*78)
            print('DUsnap fixe restreint au domaine courant :', Dhom_k_restr_prime[0,0])
            print('%'*78)
        print('Porosite', porosity)
        print('IG physique :', T_chi_omega[0,0])
        print('DUsnap physique :', Dhom_k_postprime[0,0])
        print('%'*78)
    ###

    # --------- Erreur de periodicite --------- #
    if err_per_calc:
        # Dans ce cas, les faces du cube sont entieres
        ##err_per_gr(cen,r,chi_prime_n,npas_err,fig_todo)
        if config == 'cylsph' and geo_p == 'ray_sph':
            err_per_gr_compl(config, ray_fix, chi_prime_n, npas_err, fig_todo)
