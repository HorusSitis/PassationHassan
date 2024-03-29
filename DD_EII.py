#####################################################################################################################################
######################################### Etape II : extrapolation des cliches, domaine_fixe ########################################
#####################################################################################################################################

print('='*100)
print('Etape II, dimension 2')
print('='*100)

# Chargement de la liste des snapshots physiques

# l_name='Lchi_'+str(N_snap)+'_'+config+'_'+geo_p+'_deg'+str(VFS_degree)+'_'+ordo+'_'+computer

with sh.open(repertoire_parent+l_name) as l_loa:
    list_chi_v = l_loa['maliste']
# sys.exit('!'*100)
# Extrapolation au domaine Omega_fixe :
if dom_fixe == 'am':
    mesh_fixe_name = 'maillages_per/2D/maillage_fixe2d_am'#.xml'
elif config == 'compl':
    mesh_fixe_name = 'maillages_per/2D/maillage_trous2D_'+geo_p+'_fixe'#.xml'

# generation du maillage fixe
if not mesh_ex_done:

    xinf = xyinfsup[0][0]
    yinf = xyinfsup[0][1]
    xsup = xyinfsup[1][0]
    ysup = xyinfsup[1][1]

    xcen = (xinf + xsup)/2.
    ycen = (yinf + ysup)/2.
    # elif mention == '_som':
    #     #

    ### on cree le fichier qui code le maillage pour gmsh
    fichier_sansentete = open(mesh_fixe_name + '_sansxyinfsup' + '.txt', 'r')

    gen_mesh = open(mesh_fixe_name + '.txt', 'w')

    gen_mesh.write('pas_mesh = ' + str(round(1./res_gmsh, 2)) + ';' + '\n')

    if geo_p == 'ray_min':
        gen_mesh.write('R = ' + str(rho_appr_min) + ';' + '\n')
        if mention == '':
            gen_mesh.write('cx = ' + str(xcen) + ';' + '\n')
            gen_mesh.write('cy = ' + str(ycen) + ';' + '\n')

    if config == 'compl':
        gen_mesh.write('S = ' + str(ray_p) + ';' + '\n')
        if geo_p == 'hor':
            gen_mesh.write('cy = ' + str((yinf + ysup)/2.) + ';' + '\n')

    gen_mesh.write('xmin = ' + str(xinf) + ';' + '\n')
    gen_mesh.write('ymin = ' + str(yinf) + ';' + '\n')
    gen_mesh.write('xmax = ' + str(xsup) + ';' + '\n')
    gen_mesh.write('ymax = ' + str(ysup) + ';' + '\n')

    for line in fichier_sansentete :
        gen_mesh.write(line)

    gen_mesh.close()
    fichier_sansentete.close()

    ### Conversion du fichier texte obtenu et stockage du maillage en fichiers .xml
    os.rename(mesh_fixe_name + '.txt', mesh_fixe_name + '.geo')

    ### on genere le maillage avec gmsh
    os.system('gmsh -' + str(dimension) + ' ' + mesh_fixe_name+ '.geo')

    if fig_todo == 'aff':
        ## affichage du maillage genere selon les parametres d'entree
        os.system('gmsh ' + mesh_fixe_name + '.msh')

    ### on convertit le maillage avec dolfin
    os.system('dolfin-convert ' + mesh_fixe_name + '.msh' + ' ' + mesh_fixe_name + '.xml')


mesh_fixe = Mesh(mesh_fixe_name + '.xml')

# if config == 'compl':
#     mesh_name = mesh_prefix + str(int(round(100*r,2))) + '_rayp' + str(int(round(100*ray_p,2)))
# else:
#     mesh_name = mesh_prefix + mention + str(int(round(100*r,2)))
#
#
# mesh = Mesh(mesh_repository + mesh_name + '.xml')
# sys.exit()
V_fixe = VectorFunctionSpace(mesh_fixe, 'P', VFS_degree, constrained_domain=PeriodicBoundary())

plot(mesh_fixe, linewidth=lw_bis)
if fig_todo=='aff':
    plt.show()
elif fig_todo == 'save':
    plt.savefig('Figures2D/maillage_gmsh_per_fix_'+config+geo_p+'_res'+str(res_gmsh)+'.png')
plt.close()
# Extrapolation des solutions du probleme physique

def extra_snap(n):
    r = list_rho_appr[n]

    # chargement du snapshot courant
    chi_n_v=list_chi_v[n]

    # mise sous forme d'une fonction EF
    if config!='compl':
        if cen_snap_ray == [xinf, yinf]:
            mention = '_som'
        else:
            mention = ''

    # chargement du maillage genere avec DD_EI
    if config == 'compl':
        mesh_name = mesh_prefix + str(int(round(100*r,2))) + '_rayp' + str(int(round(100*ray_p,2)))
    else:
        mesh_name = mesh_prefix + mention + str(int(round(100*r,2)))

    mesh = Mesh(mesh_repository + mesh_name + '.xml')
    V_n=VectorFunctionSpace(mesh, 'P', VFS_degree, constrained_domain=PeriodicBoundary())

    chi_n=Function(V_n)
    chi_n.vector().set_local(chi_n_v)
    # extrapolation du snapshot au domaine fixe
    chi_n.set_allow_extrapolation(True)
    chi_n_prime=interpolate(chi_n,V_fixe)
    #
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
    with sh.open(repertoire_parent+u_name) as u_sto:
        u_sto['maliste'] = Usnap

else:
    # Chargement de la matrice des snapshots
    print(repertoire_parent+u_name)
    with sh.open(repertoire_parent+u_name) as u_loa:
        Usnap = u_loa['maliste']



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

    # Affichage des valeurs de la solution interpolee
    plot(chi_prime_n, linewidth=lw_bis)
    if fig_todo=='aff':
        plt.title('Snapshot '+str(n+1),fontsize=40)
        plt.show()
    else:
        plt.savefig('Figures2D/snap_'+str(n+1)+'_sur'+str(N_snap)+config+'_'+geo_p+'_res'+str(res_gmsh)+'.png')
    plt.close()

    # Affichage des composantes scalaires : interpolee
    if config=='cer_un' and geo_p=='ray' and fig_todo!='':
        fig_chi(cen_snap_ray,r,chi_prime_n,fig_todo)
