# -*- coding: utf-8 -*-

# from fenics import *
#
# from dolfin import *
# from mshr import *
# import matplotlib.pyplot as plt
# import numpy as np
# from math import sqrt
# from math import exp
# import sys

# tol=1e-10
# xinf=0.0
# yinf=0.0
# zinf=0.0
# xsup=1.0
# ysup=1.0
# zsup=1.0
#
# typ_msh='gms'

############################# Une classe de sous domaines pour tous les calculs : structure periodique avec repetition de la cellule elementaire #############################

# dimension=3
#
# class PeriodicBoundary(SubDomain):
#  # Left boundary is "target domain" G
#  def inside(self, x, on_boundary):
#   return on_boundary and not(near(x[0],xsup,tol) or near(x[1],ysup,tol) or near(x[2],zsup,tol))
#  # Map right boundary (R) to left boundary (G), top to bottom, back to front
#  def map(self, x, y):
#   for i in range(dimension):
#    if near(x[i],1.0,tol):
#     y[i]=0.0
#    else:
#     y[i]=x[i]

############################# Pour creer des maillages, avec des familles de cellules elementaires #############################

def crown(r):#epaisseur de la couronne dans laquelle le maillage est raffine
    return r+tol#*(1+0.2*exp(-r**2))#1.2*r

def raffinement_maillage_sph_per(cen,r,mesh):
    markers = MeshFunction("bool", mesh, mesh.topology().dim())
    markers.set_all(False)
    # on cree une liste des centres des inclusions voisines de la cellule elementaire
    l_cen=[]
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                l_cen.append([cen[0]+i,cen[1]+j,cen[2]+k])
    for c in cells(mesh):
        for f in facets(c):
            # Raffinement autour de l'inclusion periodique
            for cen_per in l_cen:
                if (sqrt((f.midpoint()[0]-cen_per[0])**2+(f.midpoint()[1]-cen_per[1])**2+(f.midpoint()[2]-cen_per[2])**2)<=crown(r)):
                    markers[c] = True
            # Raffinement aux bords du domaine
            if any([f.midpoint()[k]==0 or f.midpoint()[k]==1 for k in range(0,3)]):
                markers[c]=True
    mesh=refine(mesh, markers, redistribute=True)
    return mesh

def creer_maill_sph(cen,r,res):#valable quel que soit la position de l'inclusion : centre, choisi aleatoirement. 1.2*r<0.5.
    # Creation de la cellule elementaire avec inclusion
    box=Box(Point(0,0,0),Point(1,1,1))#Box(Point(-1,-1,-1),Point(2,2,2))
    domain=box
    l_sph=[]
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                l_sph.append(Sphere(Point(cen[0]+i,cen[1]+j,cen[2]+k),r))
    for sph_per in l_sph:
        domain=domain-sph_per
    # Creation du permier maillage
    mesh=generate_mesh(domain,res)
    # On raffine le long du bord de l'inclusion
    mesh=raffinement_maillage_sph_per(cen,r,mesh)
    #
    return(mesh)

def raffinement_maillage_cyl_per(top,r,mesh):
    markers = MeshFunction("bool", mesh, mesh.topology().dim())
    markers.set_all(False)
    # on cree une liste des centres des inclusions voisines de la cellule elementaire
    l_top=[]
    for i in range(-1,2):
        for j in range(-1,2):
            l_top.append([top[0]+i,top[1]+j])
    for c in cells(mesh):
        for f in facets(c):
            # Raffinement autour de l'inclusion periodique
            for top_per in l_top:
                if (sqrt((f.midpoint()[0]-top_per[0])**2+(f.midpoint()[2]-top_per[1])**2)<=crown(r)):
                    markers[c] = True
    # Raffinement aux bords du domaine
            if any([f.midpoint()[k]==0 or f.midpoint()[k]==1 for k in range(0,3)]):
                markers[c]=True
    mesh=refine(mesh, markers, redistribute=True)
    return mesh

def creer_maill_cyl(top,r,slices_cyl,res):#valable quel que soit la position de l'inclusion : centre, choisi aleatoirement. 1.2*r<0.5.
    # Creation de la cellule elementaire avec inclusion
    box=Box(Point(0,0,0),Point(1,1,1))
    domain=box
    l_cyl=[]
    for i in range(-1,2):
        for j in range(-1,2):
            l_cyl.append(Cylinder(Point(top[0]+i,0.,top[1]+j),Point(top[0]+i,1.,top[1]+j),r,slices_cyl))
    print(len(l_cyl))
    for cyl_per in l_cyl:
        domain=domain-cyl_per
    domain=box-l_cyl[4]
    # Creation du permier maillage
    mesh=generate_mesh(domain,res)
    # On raffine le long du bord de l'inclusion
    #mesh=raffinement_maillage_cyl_per(top,r,mesh)
    #
    return(mesh)

############################# Pour creer des maillages avec gmsh, a partir de fichiers .txt en guise de canevas #############################

def creer_maill_per_gpar(config, geo_p, xyzinfsup, rho, ray_fix, res):

    xinf = xyzinfsup[0][0]
    yinf = xyzinfsup[0][1]
    zinf = xyzinfsup[0][2]
    xsup = xyzinfsup[1][0]
    ysup = xyzinfsup[1][1]
    zsup = xyzinfsup[1][2]

    xcen = (xinf + xsup)/2.
    ycen = (yinf + ysup)/2.
    zcen = (zinf + zsup)/2.

    # if config == 'sph_un':
    #     mesh_prefix = 'cubesphre_periodique_triangle_'
    # elif config == '2sph':
    #     mesh_prefix = 'cube2sph_periodique_triangle_'
    # elif config == 'cylsph':
    #     mesh_prefix = 'cubecylsph_periodique_triangle_'

    ### on cree le fichier qui code le maillage pour gmsh
    mesh_name = mesh_prefix + 'sansgpar'
    fichier_sansentete = open(mesh_repository + mesh_name + '.txt', 'r')

    if config == 'sph_un':
        nom_fichier_avecgpar = mesh_prefix + 'rayc' + str(int(round(100*rho,2))) + '_sur' + str(res)
    elif config == '2sph':
        nom_fichier_avecgpar = mesh_prefix + 'rayc' + str(int(round(100*rho,2))) + '_rayp' + str(int(round(100*ray_fix,2))) + '_sur' + str(res)
    elif config == 'cylsph' and geo_p == 'ray_sph':
        nom_fichier_avecgpar = mesh_prefix + 'rayc' + str(int(round(100*rho,2))) + '_rayp' + str(int(round(100*ray_fix,2))) + '_sur' + str(res)
    elif config == 'cylsph' and geo_p == 'ray_cyl':
        nom_fichier_avecgpar = mesh_prefix + 'rayc' + str(int(round(100*ray_fix,2))) + '_rayp' + str(int(round(100*rho,2))) + '_sur' + str(res)

    gen_mesh = open(mesh_repository + nom_fichier_avecgpar + '.txt', 'w')

    gen_mesh.write('/'*64 + '\n')
    gen_mesh.write('/'*10 + ' Fichier de maillage rempli automatiquement ' + '/'*10 + '\n')
    gen_mesh.write('/'*64 + '\n' + '\n')

    # parametres geometriques

    if config == 'sph_un':
        gen_mesh.write('R = ' + str(rho) + ';' + '\n' + '\n')
    if config == 'cube2sph':
        gen_mesh.write('R = ' + str(rho) + ';' + '\n')
        gen_mesh.write('S = ' + str(ray_fix) + ';' + '\n' + '\n')
    if config == 'cyl_sph' and geo_p == 'ray_sph':
        gen_mesh.write('R = ' + str(rho) + ';' + '\n')
        gen_mesh.write('S = ' + str(ray_fix) + ';' + '\n' + '\n')
    if config == 'cyl_sph' and geo_p == 'ray_cyl':
        gen_mesh.write('R = ' + str(ray_fix) + ';' + '\n')
        gen_mesh.write('S = ' + str(rho) + ';' + '\n' + '\n')

    # centre de la maille

    gen_mesh.write('xc = ' + str(xcen) + ';' + '\n')
    gen_mesh.write('yc = ' + str(ycen) + ';' + '\n')
    gen_mesh.write('zc = ' + str(zcen) + ';' + '\n' + '\n')

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
    os.rename(mesh_repository + nom_fichier_avecgpar + '.txt', mesh_repository + nom_fichier_avecgpar + '.geo')

    ### on genere le maillage avec gmsh
    os.system('gmsh -' + str(dimension) + ' ' + mesh_repository + nom_fichier_avecgpar + '.geo')

    if fig_todo == 'aff':
        ## affichage du maillage genere selon les parametres d'entree
        os.system('gmsh ' + mesh_repository + nom_fichier_avecgpar + '.msh')

    ### on convertit le maillage avec dolfin
    os.system('dolfin-convert ' + mesh_repository + nom_fichier_avecgpar + '.msh' + ' ' + mesh_repository + nom_fichier_avecgpar + '.xml')

    ### pas de sortie pour la procedure
    return()

############################# Pour creer des snapshots, inclusion circulaire periodique unique #############################

def snapshot_sph_per(cen,r,res,typ_sol):

    # mesh_prefix = 'cubesphere_periodique_triangle_'
    nom_fichier_avecgpar = mesh_prefix + 'rayc' + str(int(round(100*r,2))) + '_sur' + str(res)
    # creer_maill_per_gpar(config, geo_p, xyzinfsup, r, 0., res)
    mesh_s_r = Mesh(mesh_repository + nom_fichier_avecgpar + '.xml')

    # On pose et on resoud le probleme aux elements finis
    V = VectorFunctionSpace(mesh_s_r, 'P', 2, constrained_domain = PeriodicBoundary())

    ## On definit l'interface fluide-solide, periodique a geometrie spherique
    l_cen = []
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                l_cen.append([cen[0]+i,cen[1]+j,cen[2]+k])
    class inclusion_periodique(SubDomain):
        def inside(self,x,on_boundary):
            return (on_boundary and any([between((x[0]-c[0]), (-r-tol, r+tol)) for c in l_cen]) and any([between((x[1]-c[1]), (-r-tol, r+tol)) for c in l_cen]) and any([between((x[2]-c[2]), (-r-tol, r+tol)) for c in l_cen]))#points de la frontiere du systeme compris dans la boule de centre et rayons cen et r, pour la norme infinie

    ## Utilisation des classes definies precedemment : mesure de la limite du domaine fluide
    Gamma_sf = inclusion_periodique()
    boundaries = MeshFunction("size_t", mesh_s_r, mesh_s_r.topology().dim()-1)
    print('Facettes : ',mesh_s_r.num_edges())
    # On attribue une valeur par defaut aux frontieres du domaine fluide, qui concerne plus particulierement l'interface fluide-fluide
    boundaries.set_all(1)
    Gamma_sf.mark(boundaries, 7)
    ds = Measure("ds")(subdomain_data=boundaries)
    num_ff = 1
    num_sphere = 7

    ## On resoud le probleme faible, avec une condition de type Neumann au bord de l'obstacle
    normale = FacetNormal(mesh_s_r)
    nb_noeuds = V.dim()
    u = TrialFunction(V)
    v = TestFunction(V)
    a = tr(dot((grad(u)).T, grad(v)))*dx
    l = -dot(normale,v)*ds(num_sphere)

    ### Resolution

    u1 = Function(V)
    A = assemble(a)
    L = assemble(l)

    ## solveur de krylov pour un probleme bien pose, dans un hyperplan
    if typ_sol == 'kr_null_vect':

        solver = dolfin.PETScKrylovSolver("cg")

        solver.parameters["absolute_tolerance"] = 1e-6
        solver.parameters["relative_tolerance"] = 1e-6
        solver.parameters["maximum_iterations"] = 5000
        solver.parameters["error_on_nonconvergence"] = True
        solver.parameters["monitor_convergence"] = True
        solver.parameters["report"] = True

        solver.set_operator(A)

        null_vec = dolfin.Vector(u1.vector())
        V.dofmap().set(null_vec, 1.0)
        null_vec *= 1.0/(dolfin.norm(null_vec, norm_type = 'l2'))

        null_space = dolfin.VectorSpaceBasis([null_vec])
        dolfin.as_backend_type(A).set_nullspace(null_space)
        null_space.orthogonalize(L)

        # Resolution et restitution de chi
        solver.solve(u1.vector(), L)
        chi = u1

    ## autres solveurs : probleme mal pose a piori
    else:
        if typ_sol=='default':
            solve(a==L,u)
        elif typ_sol=='bic_cyr':
            u = Function(V)
            solver_correction = KrylovSolver("bicgstab", "amg")
            solver_correction.parameters["absolute_tolerance"] = 1e-6
            solver_correction.parameters["relative_tolerance"] = 1e-6
            solver_correction.parameters["maximum_iterations"] = 5000
            solver_correction.parameters["error_on_nonconvergence"]= True
            solver_correction.parameters["monitor_convergence"] = True
            solver_correction.parameters["report"] = True
            # # AA=assemble(a)
            # # LL=assemble(L)
            # # solver_correction.solve(AA,u.vector(),LL)
            # A=assemble(a)
            # L=assemble(l)
            solver_correction.solve(A,u.vector(),L)
        # Annulation de la valeur moyenne
        porosity=1-(4/3)*pi*r**3
        moy_u_x=assemble(u[0]*dx)/porosity
        moy_u_y=assemble(u[1]*dx)/porosity
        moy_u_z=assemble(u[2]*dx)/porosity
        moy=Function(V)
        moy=Constant((moy_u_x,moy_u_y,moy_u_z))
        print("Valeur moyenne de u :",[moy_u_x,moy_u_y,moy_u_z])
        moy_V=interpolate(moy,V)
        moy_Vv=moy_V.vector().get_local()
        u_v=u.vector().get_local()
        chi_v=u_v-moy_Vv
        chi=Function(V)
        chi.vector().set_local(chi_v)
        chi=u

    ### Resultat : snapshot
    return(chi)

def snapshot_cyl_per(top,r,res,typ_sol):### ------------------> resolution : avec gmsh
    t_x,t_z=top[0],top[1]
    mesh_name="maillages_per/3D/cubecylindre_periodique_triangle_"+str(int(round(100*r,2)))+"sur"+str(res)+".xml"
    mesh_c_r=Mesh(mesh_name)
    # On pose et on resoud le probleme aux elements finis
    V=VectorFunctionSpace(mesh_c_r, 'P', 2, constrained_domain=PeriodicBoundary())
    ## On definit l'interface fluide-solide, periodique a geometrie spherique
    l_top=[]
    for i in range(-1,2):
        for j in range(-1,2):
            l_top.append([top[0]+i,top[1]+j])
    class inclusion_periodique(SubDomain):
        def inside(self,x,on_boundary):
            return (on_boundary and any([between((x[0]-t[0]), (-r-tol, r+tol)) for t in l_top]) and any([between((x[2]-t[1]), (-r-tol, r+tol)) for t in l_top]))#points de la frontiere du dysteme compris dans ... pour la norme infinie
    ## Utilisation des classes definies precedemment : mesure de la limite du domaine fluide
    Gamma_sf = inclusion_periodique()
    boundaries = MeshFunction("size_t", mesh_c_r, mesh_c_r.topology().dim()-1)
    print('Facettes : ',mesh_c_r.num_edges())
    # On attribue une valeur par defaut aux frontieres du domaine fluide, qui concerne plus particulierement l'interface fluide-fluide
    boundaries.set_all(1)
    Gamma_sf.mark(boundaries, 7)
    ds = Measure("ds")(subdomain_data=boundaries)
    num_ff=1
    num_cyl=7
    ## On resoud le probleme faible, avec une condition de type Neumann au bord de l'obstacle
    normale = FacetNormal(mesh_c_r)
    nb_noeuds=V.dim()
    u = TrialFunction(V)
    v = TestFunction(V)
    a=tr(dot((grad(u)).T, grad(v)))*dx
    l=-dot(normale,v)*ds(num_cyl)

    ### Resolution
    u1 = Function(V)
    A = assemble(a)
    L = assemble(l)

    ## solveur de krylov pour un probleme bien pose, dans un hyperplan
    if typ_sol == 'kr_null_vect':

        solver = dolfin.PETScKrylovSolver("cg")

        solver.parameters["absolute_tolerance"] = 1e-6
        solver.parameters["relative_tolerance"] = 1e-6
        solver.parameters["maximum_iterations"] = 5000
        solver.parameters["error_on_nonconvergence"] = True
        solver.parameters["monitor_convergence"] = True
        solver.parameters["report"] = True

        solver.set_operator(A)

        null_vec = dolfin.Vector(u1.vector())
        V.dofmap().set(null_vec, 1.0)
        null_vec *= 1.0/(dolfin.norm(null_vec, norm_type = 'l2'))

        null_space = dolfin.VectorSpaceBasis([null_vec])
        dolfin.as_backend_type(A).set_nullspace(null_space)
        null_space.orthogonalize(L)

        # Resolution et restitution de chi
        solver.solve(u1.vector(), L)
        chi = u1

    ## probleme mal pose a priori
    else:
        if typ_sol=='default':
            solve(a==l,u)
        elif typ_sol=='bic_cyr':
            u=Function(V)
            solver_correction = KrylovSolver("bicgstab", "amg")
            solver_correction.parameters["absolute_tolerance"] = 1e-6
            solver_correction.parameters["relative_tolerance"] = 1e-6
            solver_correction.parameters["maximum_iterations"] = 5000
            solver_correction.parameters["error_on_nonconvergence"]= True
            solver_correction.parameters["monitor_convergence"] = True
            solver_correction.parameters["report"] = True
            # # AA=assemble(a)
            # # LL=assemble(L)
            # # solver_correction.solve(AA,u.vector(),LL)
            # A=assemble(a)
            # L=assemble(l)
            solver_correction.solve(A,u.vector(),L)
        ## Annulation de la valeur moyenne
        porosity=1-pi*r**2
        moy_u_x=assemble(u[0]*dx)/porosity
        moy_u_y=assemble(u[1]*dx)/porosity
        moy_u_z=assemble(u[2]*dx)/porosity
        moy=Function(V)
        moy=Constant((moy_u_x,moy_u_y,moy_u_z))
        print("Valeur moyenne de u :",[moy_u_x,moy_u_y,moy_u_z])
        moy_V=interpolate(moy,V)
        moy_Vv=moy_V.vector().get_local()
        u_v=u.vector().get_local()
        chi_v=u_v-moy_Vv
        chi=Function(V)
        chi.vector().set_local(chi_v)
        chi=u

    ### Resultat : snapshot
    return(chi)

def snapshot_compl_per(rho, ray_fix, config, res, typ_sol):### ------------------> resolution : avec gmsh

    ## convention : nommage du fichier .xml de maillage
    # mesh_prefix = 'cube2sph_periodique_triangle_'
    # mesh_prefix = 'cubecylsph_periodique_triangle_'
    if config == '2sph':
        nom_fichier_avecgpar = mesh_prefix + 'rayc' + str(int(round(100*rho,2))) + '_rayp' + str(int(round(100*ray_fix,2))) + '_sur' + str(res)
    elif config == 'cylsph' and geo_p == 'ray_sph':
        nom_fichier_avecgpar = mesh_prefix + 'rayc' + str(int(round(100*rho,2))) + '_rayp' + str(int(round(100*ray_fix,2))) + '_sur' + str(res)
    elif config == 'cylsph' and geo_p == 'ray_cyl':
        nom_fichier_avecgpar = mesh_prefix + 'rayc' + str(int(round(100*ray_fix,2))) + '_rayp' + str(int(round(100*rho,2))) + '_sur' + str(res)

    ## creation du maillage
    # creer_maill_per_gpar(config, geo_p, xyzinfsup, rho, ray_fix, res)
    mesh = Mesh(mesh_repository + nom_fichier_avecgpar + '.xml')

    ## VFS pour le probleme variationnel
    V = VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
    print('Noeuds : ',V.dim())

    ## On definit la bordure du domaine, sur laquelle integrer le second membre "L" de l'equation en dimension finie
    boundaries = MeshFunction('size_t', mesh, mesh_repository + nom_fichier_avecgpar + "_facet_region" + ".xml")
    print('Facettes : ', mesh.num_edges())
    ds = Measure("ds")(subdomain_data=boundaries)

    ## Marquage des bordures pour la condition de Neumann
    num_solid_boundary=1700
    print("Numero de la frontiere physique sf :",num_solid_boundary)

    ## On resoud le probleme faible, avec une condition de type Neumann au bord de l'obstacle
    normale = FacetNormal(mesh)
    nb_noeuds=V.dim()
    u = TrialFunction(V)
    v = TestFunction(V)
    a=tr(dot((grad(u)).T, grad(v)))*dx
    l=-dot(normale,v)*ds(num_solid_boundary)

    ### Resolution
    u1 = Function(V)
    A = assemble(a)
    L = assemble(l)

    ## solveur de krylov pour un probleme bien pose, dans un hyperplan
    if typ_sol == 'kr_null_vect':

        solver = dolfin.PETScKrylovSolver("cg")

        solver.parameters["absolute_tolerance"] = 1e-6
        solver.parameters["relative_tolerance"] = 1e-6
        solver.parameters["maximum_iterations"] = 5000
        solver.parameters["error_on_nonconvergence"] = True
        solver.parameters["monitor_convergence"] = True
        solver.parameters["report"] = True

        solver.set_operator(A)

        null_vec = dolfin.Vector(u1.vector())
        V.dofmap().set(null_vec, 1.0)
        null_vec *= 1.0/(dolfin.norm(null_vec, norm_type = 'l2'))

        null_space = dolfin.VectorSpaceBasis([null_vec])
        dolfin.as_backend_type(A).set_nullspace(null_space)
        null_space.orthogonalize(L)

        # Resolution et restitution de chi
        solver.solve(u1.vector(), L)
        chi = u1

    ## probleme mal pose a priori
    else:
        if typ_sol=='default':
            solve(a==l,u)
        elif typ_sol=='bic_cyr':
            u=Function(V)
            solver_correction = KrylovSolver("bicgstab", "amg")
            solver_correction.parameters["absolute_tolerance"] = 1e-6
            solver_correction.parameters["relative_tolerance"] = 1e-6
            solver_correction.parameters["maximum_iterations"] = 5000
            solver_correction.parameters["error_on_nonconvergence"]= True
            solver_correction.parameters["monitor_convergence"] = True
            solver_correction.parameters["report"] = True
            # # AA=assemble(a)
            # # LL=assemble(L)
            # # solver_correction.solve(AA,u.vector(),LL)
            # A=assemble(a)
            # L=assemble(l)
            solver_correction.solve(A,u.vector(),L)
        ## Annulation de la valeur moyenne
        if config=='2sph':
            porosity=1-(4/3)*pi*(r_cen**3+r_per**3)
        elif config=='cylsph':
            porosity=1-(4/3)*pi*r_cen**3-pi*r_per**2
        moy_u_x=assemble(u[0]*dx)/porosity
        moy_u_y=assemble(u[1]*dx)/porosity
        moy_u_z=assemble(u[2]*dx)/porosity
        moy=Function(V)
        moy=Constant((moy_u_x,moy_u_y,moy_u_z))
        print("Valeur moyenne de u :",[moy_u_x,moy_u_y,moy_u_z])
        moy_V=interpolate(moy,V)
        moy_Vv=moy_V.vector().get_local()
        u_v=u.vector().get_local()
        chi_v=u_v-moy_Vv
        chi=Function(V)
        chi.vector().set_local(chi_v)
        chi=u

    ### Resultat : snapshot
    return(chi)

############################# Pour tester la periodicite d'un champ : impression des valeurs du champ ou de son gradient, ou representation graphique #############################

def err_per_ind_01(u,cen,r,Npas):# comparaison entre les valeurs individuelles prises par chi aux frontieres de la cellule
    pas=1/Npas
    print('Plan frontal Oxz :')
    for k in range(0,1+Npas):
        for l in range(0,1+Npas):
            is_fluid=True
            for m in range(-1,2):
                for n in range(-1,2):
                    for o in range(-1,2):
                        if sqrt((pas*k-(cen[0]+m))**2+(0.0-(cen[1]+n))**2+(pas*l-(cen[2]+o))**2)<=r:
                            is_fluid=False
            if is_fluid:
                print('x='+str(pas*k),'z='+str(pas*l),u((pas*k,0.0,pas*l)),u((pas*k,1.0,pas*l)))
            else:
                print('x='+str(pas*k),'z='+str(pas*l),"solid")
    print('Plan horizontal Oxy :')
    for k in range(0,1+Npas):
        for l in range(0,1+Npas):
            is_fluid=True
            for m in range(-1,2):
                for n in range(-1,2):
                    for o in range(-1,2):
                        if sqrt((pas*k-(cen[0]+m))**2+(pas*l-(cen[1]+n))**2+(0.0-(cen[2]+o))**2)<=r:
                            is_fluid=False
            if is_fluid:
                print('x='+str(pas*k),'y='+str(pas*l),u((pas*k,pas*l,0.0)),u((pas*k,pas*l,1.0)))
            else:
                print('x='+str(pas*k),'y='+str(pas*l),"solid")
    print('Plan lateral Oyz :')
    for k in range(0,1+Npas):
        for l in range(0,1+Npas):
            is_fluid=True
            for m in range(-1,2):
                for n in range(-1,2):
                    for o in range(-1,2):
                        if sqrt((0.0-(cen[0]+m))**2+(pas*k-(cen[1]+n))**2+(pas*l-(cen[2]+o))**2)<=r:
                            is_fluid=False
            if is_fluid:
                print('y='+str(pas*k),'z='+str(pas*l),u((0.0,pas*k,pas*l)),u((1.0,pas*k,pas*l)))
            else:
                print('y='+str(pas*k),'z='+str(pas*l),"solid")
    return()


def err_per_gr(cen,r,u,Npas,todo):
    #coord_b=np.arange(Npas+1)
    X,Y=np.meshgrid(np.arange(1+Npas),np.arange(1+Npas))
    pas=1/Npas
    # ---------------------- chi on Oxz, front and back of the cell ---------------------- #
    # Creates the vectors where the chi component values will be registered
    ufb_y1_0=np.zeros((Npas+1,Npas+1))
    ufb_y2_0=np.zeros((Npas+1,Npas+1))
    ufb_y3_0=np.zeros((Npas+1,Npas+1))
    ufb_y1_1=np.zeros((Npas+1,Npas+1))
    ufb_y2_1=np.zeros((Npas+1,Npas+1))
    ufb_y3_1=np.zeros((Npas+1,Npas+1))
    ## We collect the values of chi on the fluid domain, and suppose chi vanishes on the solid domain
    for k in range(0,Npas+1):
        for l in range(0,1+Npas):
            is_fluid=True
            for m in range(-1,2):
                for n in range(-1,2):
                    for o in range(-1,2):
                        if sqrt((pas*k-(cen[0]+m))**2+(0.0-(cen[1]+n))**2+(pas*l-(cen[2]+o))**2)<=r:
                            is_fluid=False
            if is_fluid:
                # u on the front face
                vect_u_0=u((pas*k,0.0,pas*l))
                ufb_y1_0[k,l]=vect_u_0[0]
                ufb_y2_0[k,l]=vect_u_0[1]
                ufb_y3_0[k,l]=vect_u_0[2]
                # u on the back face
                vect_u_1=u((pas*k,1.0,pas*l))
                ufb_y1_1[k,l]=vect_u_1[0]
                ufb_y2_1[k,l]=vect_u_1[1]
                ufb_y3_1[k,l]=vect_u_1[2]
    # ---------------------- chi on Oxy, top and bottom of the cell ---------------------- #
    # Creates the vectors where the chi component values will be registered
    ubt_y1_0=np.zeros((Npas+1,Npas+1))
    ubt_y2_0=np.zeros((Npas+1,Npas+1))
    ubt_y3_0=np.zeros((Npas+1,Npas+1))
    ubt_y1_1=np.zeros((Npas+1,Npas+1))
    ubt_y2_1=np.zeros((Npas+1,Npas+1))
    ubt_y3_1=np.zeros((Npas+1,Npas+1))
    ## We collect the values of chi on the fluid domain, and suppose chi vanishes on the solid domain
    for k in range(0,Npas+1):
        for l in range(0,1+Npas):
            is_fluid=True
            for m in range(-1,2):
                for n in range(-1,2):
                    for o in range(-1,2):
                        if sqrt((pas*k-(cen[0]+m))**2+(pas*l-(cen[1]+n))**2+(0.0-(cen[2]+o))**2)<=r:
                            is_fluid=False
            if is_fluid:
                # u on the floor (b)
                vect_u_0=u((pas*k,pas*l,0.0))
                ubt_y1_0[k,l]=vect_u_0[0]
                ubt_y2_0[k,l]=vect_u_0[1]
                ubt_y3_0[k,l]=vect_u_0[2]
                # u on the roof (t)
                vect_u_1=u((pas*k,pas*l,1.0))
                ubt_y1_1[k,l]=vect_u_1[0]
                ubt_y2_1[k,l]=vect_u_1[1]
                ubt_y3_1[k,l]=vect_u_1[2]
    # ---------------------- chi on Oyz, lateral faces of the cell ---------------------- #
    # Creates the vectors where the chi component values will be registered
    ulr_y1_0=np.zeros((Npas+1,Npas+1))
    ulr_y2_0=np.zeros((Npas+1,Npas+1))
    ulr_y3_0=np.zeros((Npas+1,Npas+1))
    ulr_y1_1=np.zeros((Npas+1,Npas+1))
    ulr_y2_1=np.zeros((Npas+1,Npas+1))
    ulr_y3_1=np.zeros((Npas+1,Npas+1))
    ## We collect the values of chi on the fluid domain, and suppose chi vanishes on the solid domain
    for k in range(0,Npas+1):
        for l in range(0,1+Npas):
            is_fluid=True
            for m in range(-1,2):
                for n in range(-1,2):
                    for o in range(-1,2):
                        if sqrt((0.0-(cen[0]+m))**2+(pas*k-(cen[1]+n))**2+(pas*l-(cen[2]+o))**2)<=r:
                            is_fluid=False
            if is_fluid:
                # u on the floor (b)
                vect_u_0=u((0.0,pas*k,pas*l))
                ulr_y1_0[k,l]=vect_u_0[0]
                ulr_y2_0[k,l]=vect_u_0[1]
                ulr_y3_0[k,l]=vect_u_0[2]
                # u on the roof (t)
                vect_u_1=u((1.0,pas*k,pas*l))
                ulr_y1_1[k,l]=vect_u_1[0]
                ulr_y2_1[k,l]=vect_u_1[1]
                ulr_y3_1[k,l]=vect_u_1[2]
    # else chi_y.. stays at 0.0
    #
    # ---------------------- plots ---------------------- #
    fig=plt.figure(1)
    # We compare front and back boundaries for chi_y1, chi_y2 and chi_y3
    ## u_y1
    ax1=fig.add_subplot(331, projection='3d')
    #ax1.plot_surface(X,Y,ufb_y1_0,color='green')
    ax1.scatter(X,Y,ufb_y1_0,color='blue')
    ax1.plot_wireframe(X,Y,ufb_y1_1,color='red')
    plt.title("chi_y1 parallel to Oxz")
    ## u_y2
    ax2=fig.add_subplot(332, projection='3d')
    #ax2.plot_surface(X,Y,ufb_y2_0,color='green')
    ax2.scatter(X,Y,ufb_y2_0,color='blue')
    ax2.plot_wireframe(X,Y,ufb_y2_1,color='red')
    plt.title("chi_y2 parallel to Oxz")
    ## u_y3
    ax3=fig.add_subplot(333, projection='3d')
    #ax3.plot_surface(X,Y,ufb_y3_0,color='green')
    ax3.scatter(X,Y,ufb_y3_0,color='blue')
    ax3.plot_wireframe(X,Y,ufb_y3_1,color='red')
    plt.title("chi_y3 parallel to Oxz")
    # We compare top and bottom boundaries for chi_y1, chi_y2 and chi_y3
    ## u_y1
    ax4=fig.add_subplot(334, projection='3d')
    ax4.scatter(X,Y,ubt_y1_0,color='blue')
    #ax4.plot_surface(X,Y,ubt_y1_0,color='green')
    ax4.plot_wireframe(X,Y,ubt_y1_1,color='red')
    plt.title("chi_y1 parallel to Oxy")
    ## u_y2
    ax5=fig.add_subplot(335, projection='3d')
    ax5.scatter(X,Y,ubt_y2_0,color='blue')
    #ax5.plot_surface(X,Y,ubt_y2_0,color='green')
    ax5.plot_wireframe(X,Y,ubt_y2_1,color='red')
    plt.title("chi_y2 parallel to Oxy")
    ## u_y3
    ax6=fig.add_subplot(336, projection='3d')
    #ax6.plot_surface(X,Y,ubt_y3_0,color='green')
    ax6.scatter(X,Y,ubt_y3_0,color='blue')
    ax6.plot_wireframe(X,Y,ubt_y3_1,color='red')
    plt.title("chi_y3 parallel to Oxy")
    # We compare left and right boundaries for chi_y1, chi_y2 and chi_y3
    ## u_y1
    ax7=fig.add_subplot(337, projection='3d')
    #ax7.plot_surface(X,Y,ubt_y1_0,color='green')
    ax7.scatter(X,Y,ubt_y1_0,color='blue')
    ax7.plot_wireframe(X,Y,ulr_y1_1,color='red')
    plt.title("chi_y1 parallel to Oyz")
    ## u_y2
    ax8=fig.add_subplot(338, projection='3d')
    #ax8.plot_surface(X,Y,ubt_y2_0,color='green')
    ax8.scatter(X,Y,ubt_y2_0,color='blue')
    ax8.plot_wireframe(X,Y,ulr_y2_1,color='red')
    plt.title("chi_y2 parallel to Oyz")
    ## u_y3
    ax9=fig.add_subplot(339, projection='3d')
    #ax9.plot_surface(X,Y,ubt_y3_0,color='green')
    ax9.scatter(X,Y,ubt_y3_0,color='blue')
    ax9.plot_wireframe(X,Y,ulr_y3_1,color='red')
    plt.title("chi_y3 parallel to Oyz")
    ## Show or save
    if todo=='aff':
        plt.show()
    elif todo=='save':
        plt.savefig("Figures3D/inc_c"+"CompBo"+str(Npas)+"_cen"+str(cen[0])+str(cen[1])+str(cen[2])+"_ray"+str(r)+".png")
    ## Close
    plt.close()
    #
    return()

def err_per_gr_compl(config,ray_perif,u,Npas,todo):
    #coord_b=np.arange(Npas+1)
    X,Y=np.meshgrid(np.arange(1+Npas),np.arange(1+Npas))
    pas=1/Npas
    # ---------------------- chi on Oxz, front and back of the cell ---------------------- #
    # Creates the vectors where the chi component values will be registered
    ufb_y1_0=np.zeros((Npas+1,Npas+1))
    ufb_y2_0=np.zeros((Npas+1,Npas+1))
    ufb_y3_0=np.zeros((Npas+1,Npas+1))
    ufb_y1_1=np.zeros((Npas+1,Npas+1))
    ufb_y2_1=np.zeros((Npas+1,Npas+1))
    ufb_y3_1=np.zeros((Npas+1,Npas+1))
    ## We collect the values of chi on the fluid domain, and suppose chi vanishes on the solid domain
    for k in range(0,Npas+1):
        for l in range(0,1+Npas):
            is_fluid=True
            for m in range(0,2):
                for n in range(0,2):
                    if sqrt((pas*k-m)**2+(pas*l-n)**2)<=ray_perif:
                        is_fluid=False
            if is_fluid:
                # u on the front face
                vect_u_0=u((pas*k,0.0,pas*l))
                ufb_y1_0[k,l]=vect_u_0[0]
                ufb_y2_0[k,l]=vect_u_0[1]
                ufb_y3_0[k,l]=vect_u_0[2]
                # u on the back face
                vect_u_1=u((pas*k,1.0,pas*l))
                ufb_y1_1[k,l]=vect_u_1[0]
                ufb_y2_1[k,l]=vect_u_1[1]
                ufb_y3_1[k,l]=vect_u_1[2]
    # ---------------------- chi on Oxy, top and bottom of the cell ---------------------- #
    # Creates the vectors where the chi component values will be registered
    ubt_y1_0=np.zeros((Npas+1,Npas+1))
    ubt_y2_0=np.zeros((Npas+1,Npas+1))
    ubt_y3_0=np.zeros((Npas+1,Npas+1))
    ubt_y1_1=np.zeros((Npas+1,Npas+1))
    ubt_y2_1=np.zeros((Npas+1,Npas+1))
    ubt_y3_1=np.zeros((Npas+1,Npas+1))
    ## We collect the values of chi on the fluid domain, and suppose chi vanishes on the solid domain
    for k in range(0,Npas+1):
        for l in range(0,1+Npas):
            is_fluid=True
            if config=='2sph':
                for m in range(0,2):
                    for n in range(0,2):
                        if sqrt((pas*k-m)**2+(pas*l-n)**2)<=ray_perif:
                            is_fluid=False
            elif config=='cylsph':
                if pas*k<=ray_perif or pas*k>=1-ray_perif:
                    is_fluid=False
            if is_fluid:
                # u on the floor (b)
                vect_u_0=u((pas*k,pas*l,0.0))
                ubt_y1_0[k,l]=vect_u_0[0]
                ubt_y2_0[k,l]=vect_u_0[1]
                ubt_y3_0[k,l]=vect_u_0[2]
                # u on the roof (t)
                vect_u_1=u((pas*k,pas*l,1.0))
                ubt_y1_1[k,l]=vect_u_1[0]
                ubt_y2_1[k,l]=vect_u_1[1]
                ubt_y3_1[k,l]=vect_u_1[2]
    # ---------------------- chi on Oyz, lateral faces of the cell ---------------------- #
    # Creates the vectors where the chi component values will be registered
    ulr_y1_0=np.zeros((Npas+1,Npas+1))
    ulr_y2_0=np.zeros((Npas+1,Npas+1))
    ulr_y3_0=np.zeros((Npas+1,Npas+1))
    ulr_y1_1=np.zeros((Npas+1,Npas+1))
    ulr_y2_1=np.zeros((Npas+1,Npas+1))
    ulr_y3_1=np.zeros((Npas+1,Npas+1))
    ## We collect the values of chi on the fluid domain, and suppose chi vanishes on the solid domain
    for k in range(0,Npas+1):
        for l in range(0,1+Npas):
            is_fluid=True
            if config=='2sph':
                for m in range(0,2):
                    for n in range(0,2):
                        if sqrt((pas*k-m)**2+(pas*l-n)**2)<=ray_perif:
                            is_fluid=False
            elif config=='cylsph':
                if pas*l<=ray_perif or pas*l>=1-ray_perif:
                    is_fluid=False
            #
            if is_fluid:
                # u on the floor (b)
                vect_u_0=u((0.0,pas*k,pas*l))
                ulr_y1_0[k,l]=vect_u_0[0]
                ulr_y2_0[k,l]=vect_u_0[1]
                ulr_y3_0[k,l]=vect_u_0[2]
                # u on the roof (t)
                vect_u_1=u((1.0,pas*k,pas*l))
                ulr_y1_1[k,l]=vect_u_1[0]
                ulr_y2_1[k,l]=vect_u_1[1]
                ulr_y3_1[k,l]=vect_u_1[2]
                # else chi_y.. stays at 0.0
    #
    # ---------------------- plots ---------------------- #
    fig=plt.figure(1)
    # We compare front and back boundaries for chi_y1, chi_y2 and chi_y3
    ## u_y1
    ax1=fig.add_subplot(331, projection='3d')
    #ax1.plot_surface(X,Y,ufb_y1_0,color='green')
    ax1.scatter(X,Y,ufb_y1_0,color='blue')
    ax1.plot_wireframe(X,Y,ufb_y1_1,color='red')
    plt.title("chi_y1 parallel to Oxz")
    ## u_y2
    ax2=fig.add_subplot(332, projection='3d')
    #ax2.plot_surface(X,Y,ufb_y2_0,color='green')
    ax2.scatter(X,Y,ufb_y2_0,color='blue')
    ax2.plot_wireframe(X,Y,ufb_y2_1,color='red')
    plt.title("chi_y2 parallel to Oxz")
    ## u_y3
    ax3=fig.add_subplot(333, projection='3d')
    #ax3.plot_surface(X,Y,ufb_y3_0,color='green')
    ax3.scatter(X,Y,ufb_y3_0,color='blue')
    ax3.plot_wireframe(X,Y,ufb_y3_1,color='red')
    plt.title("chi_y3 parallel to Oxz")
    # We compare top and bottom boundaries for chi_y1, chi_y2 and chi_y3
    ## u_y1
    ax4=fig.add_subplot(334, projection='3d')
    ax4.scatter(X,Y,ubt_y1_0,color='blue')
    #ax4.plot_surface(X,Y,ubt_y1_0,color='green')
    ax4.plot_wireframe(X,Y,ubt_y1_1,color='red')
    plt.title("chi_y1 parallel to Oxy")
    ## u_y2
    ax5=fig.add_subplot(335, projection='3d')
    ax5.scatter(X,Y,ubt_y2_0,color='blue')
    #ax5.plot_surface(X,Y,ubt_y2_0,color='green')
    ax5.plot_wireframe(X,Y,ubt_y2_1,color='red')
    plt.title("chi_y2 parallel to Oxy")
    ## u_y3
    ax6=fig.add_subplot(336, projection='3d')
    #ax6.plot_surface(X,Y,ubt_y3_0,color='green')
    ax6.scatter(X,Y,ubt_y3_0,color='blue')
    ax6.plot_wireframe(X,Y,ubt_y3_1,color='red')
    plt.title("chi_y3 parallel to Oxy")
    # We compare left and right boundaries for chi_y1, chi_y2 and chi_y3
    ## u_y1
    ax7=fig.add_subplot(337, projection='3d')
    #ax7.plot_surface(X,Y,ubt_y1_0,color='green')
    ax7.scatter(X,Y,ubt_y1_0,color='blue')
    ax7.plot_wireframe(X,Y,ulr_y1_1,color='red')
    plt.title("chi_y1 parallel to Oyz")
    ## u_y2
    ax8=fig.add_subplot(338, projection='3d')
    #ax8.plot_surface(X,Y,ubt_y2_0,color='green')
    ax8.scatter(X,Y,ubt_y2_0,color='blue')
    ax8.plot_wireframe(X,Y,ulr_y2_1,color='red')
    plt.title("chi_y2 parallel to Oyz")
    ## u_y3
    ax9=fig.add_subplot(339, projection='3d')
    #ax9.plot_surface(X,Y,ubt_y3_0,color='green')
    ax9.scatter(X,Y,ubt_y3_0,color='blue')
    ax9.plot_wireframe(X,Y,ulr_y3_1,color='red')
    plt.title("chi_y3 parallel to Oyz")
    ## Show or save
    if todo=='aff':
        plt.show()
    elif todo=='save':
        plt.savefig("Figures3D/"+config+"CompBo"+str(Npas)+"_ray"+str(int(round(100*ray_perif,2)))+".png")
    ## Close
    plt.close()
    #
    return()
