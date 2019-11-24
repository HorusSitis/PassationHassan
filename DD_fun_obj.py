# from fenics import *
#
# from dolfin import *
# from mshr import *
# import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# import matplotlib.gridspec as gridspec
# import numpy as np
# from math import sqrt
# from math import exp
# import sys
#
# from DD_pars import *
#
# tol=1e-10
# xinf=0.0
# yinf=0.0
# xsup=1.0
# ysup=1.0
#
# typ_msh='gms'
#
#
# dimension=2
# VFS_degree=2

# class PeriodicBoundary(SubDomain):
#     # Left boundary is "target domain" G
#     def inside(self, x, on_boundary):
#         return on_boundary and not(near(x[0],xsup,tol) or near(x[1],ysup,tol))
#     # Map right boundary (H) to left boundary (G)
#     def map(self, x, y):
#         for i in range(dimension):
#             if near(x[i],1.0,tol):
#                 y[i]=0.0
#             else:
#                 y[i]=x[i]

############################# Pour creer des maillages, avec des familles de cellules elementaires #############################

def crown(r):#epaisseur de la couronne dans laquelle le maillage est raffine
    return r+tol#*(1+0.2*exp(-r**2))#1.2*r

def raffinement_maillage_circ_per(cen,r,mesh):# Objectif : montrer que l'emplacement de l'inclusion periodique dans la cellule elementaire ne change pas le coefficient de diffusion homogeneise, calcule avec le tenseur khi
    markers = MeshFunction("bool", mesh, mesh.topology().dim())
    markers.set_all(False)
    # on cree une liste des centres des inclusions voisines de la cellule elementaire
    l_cen=[]
    for i in range(-1,2):
        for j in range(-1,2):
            l_cen.append([cen[0]+i,cen[1]+j])
    for c in cells(mesh):
        for f in facets(c):
            # Raffinement autour de l'inclusion periodique
            for cen_per in l_cen:
                if (sqrt((f.midpoint()[0]-cen_per[0])**2+(f.midpoint()[1]-cen_per[1])**2)<=crown(r)):
                    markers[c]=True
            # Raffinement aux bords du domaine
            if any([f.midpoint()[k]==0 or f.midpoint()[k]==1 for k in range(0,2)]):
                markers[c]=True
    mesh=refine(mesh, markers, redistribute=True)
    return mesh

def creer_maill_circ(cen,r,res):#valable quel que soit la position de l'inclusion : centre, choisi aleatoirement. 1.2*r<0.5.
    # Creation de la cellule elementaire avec inclusion
    rect=Rectangle(Point(-1,-1),Point(2,2))
    domain=rect
    l_cer=[]#Circle(Point(cen[0],cen[1]),r)]
    for i in range(-1,2):
        for j in range(-1,2):
            l_cer.append(Circle(Point(cen[0]+i,cen[1]+j),r))
    for cer_per in l_cer:
        domain=domain-cer_per
    domain=domain*Rectangle(Point(0,0),Point(1,1))
    # Creation du permier maillage
    mesh=generate_mesh(domain,res)
    # On raffine le long du bord du domaine fluide
    mesh=raffinement_maillage_circ_per(cen,r,mesh)
    #
    return(mesh)

# def creer_maill_simpl_per_gpar(config, geo_p, mention, xyinfsup, rho):
def creer_maill_per_gpar(config, geo_p, mention, xyinfsup, rho, ray_p):

    xinf = xyinfsup[0][0]
    yinf = xyinfsup[0][1]
    xsup = xyinfsup[1][0]
    ysup = xyinfsup[1][1]

    if mention == "":
        xcen = (xinf + xsup)/2.
        ycen = (yinf + ysup)/2.
    # elif mention == "_som":
    #     #

    ### on cree le fichier qui code le maillage pour gmsh
    mesh_name = mesh_prefix + mention + "sansgpar"
    fichier_sansentete = open(mesh_repository + mesh_name + '.txt', 'r')

    if config == 'compl':
        nom_fichier_avecgpar = mesh_prefix + str(int(round(100*rho,2))) + "_rayp" + str(int(round(100*ray_p,2)))
    else:
        nom_fichier_avecgpar = mesh_prefix + mention + str(int(round(100*rho,2)))

    gen_mesh = open(nom_fichier_avecgpar + '.txt', 'w')

    gen_mesh.write('R = ' + str(rho) + ';' + '\n')
    if config == 'compl':
        gen_mesh.write('S = ' + str(ray_p) + ';' + '\n')

    gen_mesh.write('xmin = ' + str(xinf) + ';' + '\n')
    gen_mesh.write('ymin = ' + str(yinf) + ';' + '\n')
    gen_mesh.write('xmax = ' + str(xsup) + ';' + '\n')
    gen_mesh.write('ymax = ' + str(ysup) + ';' + '\n')
    if mention == "":
        gen_mesh.write('cx = ' + str(xcen) + ';' + '\n')
        gen_mesh.write('cy = ' + str(ycen) + ';' + '\n')

    for line in fichier_sansentete :
        gen_mesh.write(line)

    gen_mesh.close()
    fichier_sansentete.close()

    ### Conversion du fichier texte obtenu et stockage du maillage en fichiers .xml
    os.rename(mesh_repository + nom_fichier_avecgpar + '.txt', nom_fichier_avecgpar + '.geo')

    ### on genere le maillage avec gmsh
    os.system('gmsh -' + str(dimension) + ' ' + mesh_repository + nom_fichier_avecgpar + '.geo')

    if fig_todo == 'aff':
        ## affichage du maillage genere selon les parametres d'entree
        os.system('gmsh ' + mesh_repository + nom_fichier_avecgpar + '.msh')

    ### on convertit le maillage avec dolfin
    os.system('dolfin-convert ' + mesh_repository + nom_fichier_avecgpar + '.msh' + ' ' + nom_fichier_avecgpar + '.xml')

    ### pas de sortie pour la procedure
    return()

# def creer_maill_compl_per_gpar(config, geo_p, mention, xyinfsup, rho, ray_p):
#
#     xinf = xyinf[0][0]
#     yinf = xyinf[0][1]
#     xsup = xyinf[1][0]
#     ysup = xyinfsup[1][1]
#
#     if mention == "":
#         xcen = (xinf + xsup)/2.
#         ycen = (yinf + ysup)/2.
#     # elif mention == "_som":
#     #     #
#
#     ### on cree le fichier qui code le maillage pour gmsh
#     mesh_name = mesh_prefix + geo_p + mention + "_sansgpar"
#     fichier_sans_entete = open(mesh_repository + mesh_name + '.txt' + 'r')
#
#     if config == 'compl':
#         nom_fichier_avec_gpar = mesh_prefix + geo_p + "_" + str(int(round(100*rho,2))) + "_rayp" + str(int(round(100*ray_p,2)))
#     else:
#         nom_fichier_avec_gpar = mesh_name = mesh_prefix + geo_p + mention + "_" + str(int(round(100*rho,2)))
#
#     gen_mesh = open(nom_fichier_avecgpar + '.txt','w')
#
#     gen_mesh.write('R = ' + str(rho) + ';' + '\n')
#     gen_mesh.write('xmin = ' + str(xinf) + ';' + '\n')
#     gen_mesh.write('ymin = ' + str(xsup) + ';' + '\n')
#     gen_mesh.write('xsup = ' + str(yinf) + ';' + '\n')
#     gen_mesh.write('ysup = ' + str(ysup) + ';' + '\n')
#     if mention == "":
#         gen_mesh.write('xcen = ' + str(xcen) + ';' + '\n')
#         gen_mesh.write('xcen = ' + str(ycen) + ';' + '\n')
#
#     ### Conversion du fichier texte obtenu et stockage du maillage en fichiers .xml
#     os.rename(nom_fichier_avecgpar + '.txt', nom_fichier_avecgpar + '.geo')
#
#     ### on genere le maillage avec gmsh
#     os.system('gmsh -' + str(dimension) + ' ' + nom_fichier_avecgpar + '.geo')
#     ### on convertit le maillage avec dolfin
#     os.system('dolfin-convert ' + nom_fichier_avecgpar + '.msh' + ' ' + nom_fichier_avecgpar + '.xml')
#
#     ### pas de sortie pour la procedure
#     return()

############################# Pour creer des snapshots, inclusion circulaire periodique unique #############################

def snapshot_circ_per(cen,r,res):
    c_x,c_y = cen[0],cen[1]
    if cen == [0.,0.]:
        mention = "_som"
    else:
        mention = ""
    ## on appelle le fichier .xml contenant le maillage
    mesh_name = mesh_prefix + mention + str(int(round(100*r,2)))
    mesh_c_r = Mesh(mesh_repository + mesh_name + ".xml")
    # On pose et on resoud le probleme aux elements finis
    V = VectorFunctionSpace(mesh_c_r, 'P', VFS_degree, form_degree=0, constrained_domain=PeriodicBoundary())
    ## On definit la bordure du domaine, sur laquelle integrer le second membre "L" de l'equation en dimension finie
    l_cen=[]
    for i in range(-1,2):
        for j in range(-1,2):
            l_cen.append([cen[0]+i,cen[1]+j])
    class inclusion_periodique(SubDomain):
        def inside(self,x,on_boundary):
            return (on_boundary and any([between((x[0]-c[0]), (-r-tol, r+tol)) for c in l_cen]) and any([between((x[1]-c[1]), (-r-tol, r+tol)) for c in l_cen]))#points de la frontiere du dysteme compris dans la boule de centre et rayons cen et r, pour la norme infinie
    ### Utilisation de la classe definie precedemment
    Gamma_sf = inclusion_periodique()
    boundaries = MeshFunction("size_t", mesh_c_r, mesh_c_r.topology().dim()-1)
    print('Noeuds :',V.dim())
    print('Facettes : ',mesh_c_r.num_edges())
    boundaries.set_all(0)
    Gamma_sf.mark(boundaries, 5)
    ds = Measure("ds")(subdomain_data=boundaries)
    num_front_cercle=5
    ## On resoud le probleme faible, avec une condition de type Neumann au bord de l'obstacle
    normale = FacetNormal(mesh_c_r)
    nb_noeuds=V.dim()
    u = TrialFunction(V)
    v = TestFunction(V)
    a=tr(dot((grad(u)).T, grad(v)))*dx
    L=-dot(normale,v)*ds(num_front_cercle)
    ### Resolution
    u=Function(V)
    solve(a==L,u)#,bc)
    ## Annulation de la valeur moyenne
    porosity=1-pi*r**2
    moy_u_x=assemble(u[0]*dx)/porosity
    moy_u_y=assemble(u[1]*dx)/porosity
    moy=Function(V)
    moy=Constant((moy_u_x,moy_u_y))
    print("Valeur moyenne de u :",[moy_u_x,moy_u_y])
    moy_V=interpolate(moy,V)
    moy_Vv=moy_V.vector().get_local()
    u_v=u.vector().get_local()
    chi_v=u_v-moy_Vv
    chi=Function(V)
    chi.vector().set_local(chi_v)
    chi=u
    # Resultat : snapshot
    return(chi)

def snapshot_compl_per(geo_p, rho, cen, test_snap, ray_p):#,mention,res):
    # # creation de la classe de bordure periodique
    # class PeriodicBoundary(SubDomain):
    #     # Left boundary is "target domain" G
    #     def inside(self, x, on_boundary):
    #         return on_boundary and not(near(x[0],xsup,tol) or near(x[1],ysup,tol))
    #     # Map right boundary (H) to left boundary (G)
    #     def map(self, x, y):
    #         for i in range(dimension):
    #             if near(x[i],1.0,tol):
    #                 y[i]=0.0
    #             else:
    #                 y[i]=x[i]
    ##
    mesh_name = mesh_prefix + geo_p + str(int(round(100*rho,2))) + "_rayp" + str(int(round(100*ray_p,2)))
    # if geo_p!='diag' and geo_p!='hor':
    #     mesh_name="maillages_per/2D/maillage_trou2D_"+str(int(round(100*rho,2)))
    # print(mesh_name)
    # creer_maill_compl_per_gpar(config, geo_p, mention, xyinfsup, rho, ray_p)
    ## Maillage : condition de resolution et de configuration
    mesh = Mesh(mesh_repository + mesh_name + ".xml")
    # if typ_msh == 'gms':
    #     mesh_name = mesh_prefix + mention + "_" + str(int(round(100*r,2))) + ".xml"
    #     print(mesh_name)
    #     mesh_c_r = Mesh(mesh_name)
    # else:
    #     mesh_c_r = creer_maill_circ([c_x,c_y],r,res)
    V = VectorFunctionSpace(mesh, 'P', VFS_degree, form_degree=0, constrained_domain=PeriodicBoundary())
    print('Noeuds :',V.dim())
    ## On definit la bordure du domaine, sur laquelle integrer le second membre "L" de l'equation en dimension finie
    boundaries = MeshFunction('size_t', mesh, mesh_name + "_facet_region" + ".xml")
    print('Facettes : ',mesh.num_edges())
    ds = Measure("ds")(subdomain_data=boundaries)
    ## Marquage des bordures pour la condition de Neumann
    if test_snap=='solid_1':
        num_solid_boundary=1
        class SolidBoundary(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and not(near(x[0],xinf,tol) or near(x[0],xsup,tol) or near(x[1],yinf,tol) or near(x[1],ysup,tol))
        Gamma_sf = SolidBoundary()
        print('Gamma sf ne coupe pas le bord du carre')
        boundaries.set_all(0)
        Gamma_sf.mark(boundaries, 1)
    elif test_snap=='solid_2':
        #
        num_solid_boundary = 11
        print('Gamma sf coupe le bord du carre')
        ##boundaries.set_all(0)
        #
    ## On resoud le probleme faible, avec une condition de type Neumann au bord de l'obstacle
    normale = FacetNormal(mesh)
    nb_noeuds=V.dim()
    u = TrialFunction(V)
    v = TestFunction(V)
    a=tr(dot((grad(u)).T, grad(v)))*dx
    L=-dot(normale,v)*ds(num_solid_boundary)
    ### Resolution
    u=Function(V)
    solve(a==L,u)
    ## Annulation de la valeur moyenne
    porosity=1-pi*rho**2-pi*0.15**2
    moy_u_x=assemble(u[0]*dx)/porosity
    moy_u_y=assemble(u[1]*dx)/porosity
    moy=Function(V)
    moy=Constant((moy_u_x,moy_u_y))
    print("Valeur moyenne de u :",[moy_u_x,moy_u_y])
    moy_V=interpolate(moy,V)
    moy_Vv=moy_V.vector().get_local()
    u_v=u.vector().get_local()
    chi_v=u_v-moy_Vv
    chi=Function(V)
    chi.vector().set_local(chi_v)
    chi=u
    # Resultat : snapshot
    return(chi)

############################# Pour avoir quelques representations graphiques #############################

def fig_chi(cen,r,u,todo):
    # figure : les composantes du vecteur, separement
    plt.figure(figsize=(10,20))
    gs1=gridspec.GridSpec(1,2)
    gs1.update(wspace=0.1,hspace=0.01)
    ## u_y1
    ax1=plt.subplot(gs1[0])
    spl1=plot(u[0])
    divider1=make_axes_locatable(ax1)
    cax1=divider1.append_axes("right",size="8%",pad=0.08)
    bar1=plt.colorbar(spl1, cax=cax1)
    plt.title("chi_y"+str(1))
    ## u_y2
    ax2=plt.subplot(gs1[1])
    spl2=plot(u[1])
    divider2=make_axes_locatable(ax2)
    cax2=divider2.append_axes("right",size="8%",pad=0.08)
    bar2=plt.colorbar(spl2, cax=cax2)
    plt.title("chi_y"+str(2))
    ## Show or save
    if todo=='aff':
        plt.show()
    elif todo=='save':
        #plt.tight_layout(pad=2)
        plt.savefig("Figures2D/inc_c"+str(int(round(100*cen[0],2)))+str(int(round(100*cen[1],2)))+str(int(round(100*r,2)))+"_khi.png", bbox_inches="tight")#Cyrille
    ## Close
    plt.close()
    #
    return()

def fig_dchi(cen,r,U,todo):
    # figure : les composantes du vecteur, separement
    plt.figure(1)
    ## U_1_y1
    plt.subplot(221)
    spl1=plot(U[0,0])
    bar1=plt.colorbar(spl1)
    plt.title("-dchi"+str(1)+"_dy"+str(1))
    ## U_1_y2
    plt.subplot(222)
    spl2=plot(U[1,0])
    bar2=plt.colorbar(spl2)
    plt.title("-dchi"+str(2)+"_dy"+str(1))
    ## U_2_y1
    plt.subplot(223)
    spl3=plot(U[0,1])
    bar3=plt.colorbar(spl3)
    plt.title("-dchi"+str(1)+"_dy"+str(2))
    ## U_2_y2
    plt.subplot(224)
    spl4=plot(U[1,1])
    bar4=plt.colorbar(spl4)
    plt.title("-dchi"+str(2)+"_dy"+str(2))
    ## Show or save
    if todo=='aff':
        plt.show()
    elif todo=='save':
        plt.savefig("Figures2D/inc_c"+str(int(round(100*cen[0],2)))+str(int(round(100*cen[1],2)))+str(int(round(100*r,2)))+"_Gradkhi.png")
    ## Close
    plt.close()
    #
    return()

############################# Pour tester la periodicite d'un champ en norme l2 ou infinie : erreur relative #############################

def err_per_ind_01(u,Npas):
    pas=1/Npas
    print('Bord vertical :')
    for k in range(0,1+Npas):
        is_fluid=True
        for i in range(-1,2):
            for j in range(-1,2):
                if sqrt((0.0-(cen[0]+i))**2+(pas*k-(cen[1]+j))**2)<=r:
                    is_fluid=False
        if is_fluid:
            print(u((0.0,pas*k)),u((1.0,pas*k)))
        else:
            print([0.0,0.0],[0.0,0.0])
    print('Bord horizontal :')
    for k in range(0,1+Npas):
        is_fluid=True
        for i in range(-1,2):
            for j in range(-1,2):
                if sqrt((0.0-(cen[0]+i))**2+(pas*k-(cen[1]+j))**2)<=r:
                    is_fluid=False
        if is_fluid:
            print(u((pas*k,0.0)),u((pas*k,1.0)))
        else:
            print([0.0,0.0],[0.0,0.0])
            print(u((pas*k,0.0)),u((pas*k,1.0)))
    return()


def err_per_gr(cen,r,u,Npas,todo):
    coord_b=np.arange(Npas+1)
    pas=1/Npas
    # ---------------------- khi on vertical edges ---------------------- #
    # Creates the vectors where the khi component values will be registered
    ulr_y1_0=np.zeros(Npas+1)
    ulr_y2_0=np.zeros(Npas+1)
    ulr_y1_1=np.zeros(Npas+1)
    ulr_y2_1=np.zeros(Npas+1)
    ## We collect the values of khi on the fluid domain, and suppose khi vanishes on the solid domain
    for k in range(0,Npas+1):
        is_fluid=True
        for i in range(-1,2):
            for j in range(-1,2):
                if sqrt((0.0-(cen[0]+i))**2+(pas*k-(cen[1]+j))**2)<=r:
                    is_fluid=False
        if is_fluid:
            # u on the left boundary
            vect_u_0=u((0.0,pas*k))
            ulr_y1_0[k]=vect_u_0[0]
            ulr_y2_0[k]=vect_u_0[1]
            # u on the right boundary
            vect_u_1=u((1.0,pas*k))
            ulr_y1_1[k]=vect_u_1[0]
            ulr_y2_1[k]=vect_u_1[1]
    # ---------------------- khi on horizontal edges ---------------------- #
    # Creates the vectors where the khi component values will be registered
    ubt_y1_0=np.zeros(Npas+1)
    ubt_y2_0=np.zeros(Npas+1)
    ubt_y1_1=np.zeros(Npas+1)
    ubt_y2_1=np.zeros(Npas+1)
    ## We collect the values of khi on the fluid domain, and suppose khi vanishes on the solid domain
    for k in range(0,Npas+1):
        is_fluid=True
        for i in range(-1,2):
            for j in range(-1,2):
                if sqrt((pas*k-(cen[0]+i))**2+(0.0-(cen[1]+j))**2)<=r:
                    is_fluid=False
        if is_fluid:
            # u on the bottom boundary
            vect_u_0=u((pas*k,0.0))
            ubt_y1_0[k]=vect_u_0[0]
            ubt_y2_0[k]=vect_u_0[1]
            # u on the top boundary
            vect_u_1=u((pas*k,1.0))
            ubt_y1_1[k]=vect_u_1[0]
            ubt_y2_1[k]=vect_u_1[1]
    #
    # ---------------------- plots ---------------------- #
    plt.figure(1)
    # Compares and plots between left and right boundary for khi_y1 and khi_y2
    ## u_y1
    plt.subplot(221)
    plt.plot(ulr_y1_0,coord_b,'bo',ulr_y1_1,coord_b,'k')
    plt.title("khi_y1 on vertical edges")
    ## u_y2
    plt.subplot(222)
    plt.plot(ulr_y2_0,coord_b,'bo',ulr_y2_1,coord_b,'k')
    plt.title("khi_y2 on vertical edges")
    # Compares and plots between top and bottom boundary for khi_y1 and khi_y2
    ## u_y1
    plt.subplot(223)
    plt.plot(coord_b,ubt_y1_0,'bo',coord_b,ubt_y1_1,'k')
    plt.title("khi_y1 on horizontal edges")
    ## u_y2
    plt.subplot(224)
    plt.plot(coord_b,ubt_y2_0,'bo',coord_b,ubt_y2_1,'k')
    plt.title("khi_y2 on horizontal edges")
    ## Show or save
    if todo=='aff':
        plt.show()
    elif todo=='save':
        plt.savefig("Figures2D/inc_c"+"CompLRBT"+str(Npas)+"_cen"+str(int(round(100*cen[0],2)))+str(int(round(100*cen[1],2)))+"_ray"+str(int(round(100*r,2)))+".png")
    ## Close
    plt.close()
    #
    return()
