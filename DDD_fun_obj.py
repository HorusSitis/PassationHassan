# -*- coding: utf-8 -*-

from fenics import *

from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
from math import exp
import sys

tol=1e-10
xinf=0.0
yinf=0.0
zinf=0.0
xsup=1.0
ysup=1.0
zsup=1.0

typ_msh='gms'

############################# Une classe de sous domaines pour tous les calculs : structure périodique avec répétition de la cellule élémentaire #############################

dimension=3

class PeriodicBoundary(SubDomain):
 # Left boundary is "target domain" G
 def inside(self, x, on_boundary):
  return on_boundary and not(near(x[0],xsup,tol) or near(x[1],ysup,tol) or near(x[2],zsup,tol))
 # Map right boundary (R) to left boundary (G), top to bottom, back to front
 def map(self, x, y):
  for i in range(dimension):
   if near(x[i],1.0,tol):
    y[i]=0.0
   else:
    y[i]=x[i]

############################# Pour créer des maillages, avec des familles de cellules élémentaires #############################

def crown(r):#épaisseur de la couronne dans laquelle le maillage est raffiné
 return r+tol#*(1+0.2*exp(-r**2))#1.2*r

def raffinement_maillage_sph_per(cen,r,mesh):
 markers = MeshFunction("bool", mesh, mesh.topology().dim())
 markers.set_all(False)
 # on crée une liste des centres des inclusions voisines de la cellule élémentaire
 l_cen=[]
 for i in range(-1,2):
  for j in range(-1,2):
   for k in range(-1,2):
    l_cen.append([cen[0]+i,cen[1]+j,cen[2]+k])
 for c in cells(mesh):
  for f in facets(c):
   # Raffinement autour de l'inclusion périodique
   for cen_per in l_cen:
    if (sqrt((f.midpoint()[0]-cen_per[0])**2+(f.midpoint()[1]-cen_per[1])**2+(f.midpoint()[2]-cen_per[2])**2)<=crown(r)):
     markers[c] = True
   # Raffinement aux bords du domaine
   if any([f.midpoint()[k]==0 or f.midpoint()[k]==1 for k in range(0,3)]):
    markers[c]=True
 mesh=refine(mesh, markers, redistribute=True)
 return mesh

def creer_maill_sph(cen,r,res):#valable quel que soit la position de l'inclusion : centre, choisi aléatoirement. 1.2*r<0.5.
 # Création de la cellule élémentaire avec inclusion
 box=Box(Point(0,0,0),Point(1,1,1))#Box(Point(-1,-1,-1),Point(2,2,2))
 domain=box
 l_sph=[]
 for i in range(-1,2):
  for j in range(-1,2):
   for k in range(-1,2):
    l_sph.append(Sphere(Point(cen[0]+i,cen[1]+j,cen[2]+k),r))
 for sph_per in l_sph:
  domain=domain-sph_per
 #domain=Box(Point(-1,-1,-1),Point(2,2,2))-l_sph[13]
 #domain=domain*Box(Point(0,0,0),Point(1,1,1))
 # Création du permier maillage
 mesh=generate_mesh(domain,res)
 # On raffine le long du bord de l'inclusion
 mesh=raffinement_maillage_sph_per(cen,r,mesh)
 #
 return(mesh)

def raffinement_maillage_cyl_per(top,r,mesh):
 markers = MeshFunction("bool", mesh, mesh.topology().dim())
 markers.set_all(False)
 # on crée une liste des centres des inclusions voisines de la cellule élémentaire
 l_top=[]
 for i in range(-1,2):
  for j in range(-1,2):
   l_top.append([top[0]+i,top[1]+j])
 for c in cells(mesh):
  for f in facets(c):
   # Raffinement autour de l'inclusion périodique
   for top_per in l_top:
    if (sqrt((f.midpoint()[0]-top_per[0])**2+(f.midpoint()[2]-top_per[1])**2)<=crown(r)):
     markers[c] = True
   # Raffinement aux bords du domaine
   if any([f.midpoint()[k]==0 or f.midpoint()[k]==1 for k in range(0,3)]):
    markers[c]=True
 mesh=refine(mesh, markers, redistribute=True)
 return mesh

def creer_maill_cyl(top,r,slices_cyl,res):#valable quel que soit la position de l'inclusion : centre, choisi aléatoirement. 1.2*r<0.5.
 # Création de la cellule élémentaire avec inclusion
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
 # Création du permier maillage
 mesh=generate_mesh(domain,res)
 # On raffine le long du bord de l'inclusion
 #mesh=raffinement_maillage_cyl_per(top,r,mesh)
 #
 return(mesh)

############################# Pour créer des snapshots, inclusion circulaire périodique unique #############################

def snapshot_sph_per(cen,r,res):
 c_x,c_y,c_z=cen[0],cen[1],cen[2]
 if typ_msh=='gms':
  print("maillages_per/3D/cubesphere_periodique_triangle_"+str(int(round(100*r,2)))+".xml")
  mesh_s_r=Mesh("maillages_per/3D/cubesphere_periodique_triangle_"+str(int(round(100*r,2)))+".xml")
 else:
  print('pfrrh')
  mesh_s_r=creer_maill_sph([c_x,c_y,c_z],r,res)
 # On pose et on résoud le problème aux éléments finis
 V=VectorFunctionSpace(mesh_s_r, 'P', 2, constrained_domain=PeriodicBoundary())
 ## On définit l'interface fluide-solide, périodique à géométrie sphérique
 l_cen=[]
 for i in range(-1,2):
  for j in range(-1,2):
   for k in range(-1,2):
    l_cen.append([cen[0]+i,cen[1]+j,cen[2]+k])
 class inclusion_periodique(SubDomain):
  def inside(self,x,on_boundary):
   return (on_boundary and any([between((x[0]-c[0]), (-r-tol, r+tol)) for c in l_cen]) and any([between((x[1]-c[1]), (-r-tol, r+tol)) for c in l_cen]) and any([between((x[2]-c[2]), (-r-tol, r+tol)) for c in l_cen]))#points de la frontière du dystème compris dans la boule de centre et rayons cen et r, pour la norme infinie
 ## Utilisation des classes définies précédemment : mesure de la limite du domaine fluide
 Gamma_sf = inclusion_periodique()
 boundaries = MeshFunction("size_t", mesh_s_r, mesh_s_r.topology().dim()-1)
 # On attribue une valeur par défaut aux frontières du domaine fluide, qui concerne plus particulièrement l'interface fluide-fluide
 boundaries.set_all(1)
 Gamma_sf.mark(boundaries, 7)
 ds = Measure("ds")(subdomain_data=boundaries)
 num_ff=1
 num_sphere=7
 ## On résoud le problème faible, avec une condition de type Neumann au bord de l'obstacle
 normale = FacetNormal(mesh_s_r)
 nb_noeuds=V.dim()
 u = TrialFunction(V)
 v = TestFunction(V)
 a=tr(dot((grad(u)).T, grad(v)))*dx
 L=-dot(normale,v)*ds(num_sphere)
 ### Résolution
 u=Function(V)
 solve(a==L,u)
 ## Annulation de la valeur moyenne
 moy_u_x=assemble(u[0]*dx)/(1-4/3*pi*r**3)
 moy_u_y=assemble(u[1]*dx)/(1-4/3*pi*r**3)
 moy_u_z=assemble(u[2]*dx)/(1-4/3*pi*r**3)
 moy=Function(V)
 moy=Constant((moy_u_x,moy_u_y,moy_u_z))
 chi=project(u-moy,V)
 # Résultat : snapshot
 return(chi)

def snapshot_cyl_per(top,r,slices_cyl,res):
 #c_x,c_y,c_z=cen[0],cen[1],cen[2]
 mesh_c_r=creer_maill_cyl(top,r,slices_cyl,res)
 # On pose et on résoud le problème aux éléments finis
 V=VectorFunctionSpace(mesh_c_r, 'P', 2, constrained_domain=PeriodicBoundary())
 ## On définit l'interface fluide-solide, périodique à géométrie sphérique
 l_top=[]
 for i in range(-1,2):
  for j in range(-1,2):
   l_top.append([top[0]+i,top[1]+j])
 class inclusion_periodique(SubDomain):
  def inside(self,x,on_boundary):
   return (on_boundary and any([between((x[0]-t[0]), (-r-tol, r+tol)) for t in l_top]) and any([between((x[2]-t[1]), (-r-tol, r+tol)) for t in l_top]))#points de la frontière du dystème compris dans ... pour la norme infinie
 ## Utilisation des classes définies précédemment : mesure de la limite du domaine fluide
 Gamma_sf = inclusion_periodique()
 boundaries = MeshFunction("size_t", mesh_c_r, mesh_c_r.topology().dim()-1)
 # On attribue une valeur par défaut aux frontières du domaine fluide, qui concerne plus particulièrement l'interface fluide-fluide
 boundaries.set_all(1)
 Gamma_sf.mark(boundaries, 7)
 ds = Measure("ds")(subdomain_data=boundaries)
 num_ff=1
 num_cyl=7
 ## On résoud le problème faible, avec une condition de type Neumann au bord de l'obstacle
 normale = FacetNormal(mesh_c_r)
 nb_noeuds=V.dim()
 u = TrialFunction(V)
 v = TestFunction(V)
 a=tr(dot((grad(u)).T, grad(v)))*dx
 L=-dot(normale,v)*ds(num_cyl)
 ### Résolution
 u=Function(V)
 solve(a==L,u)
 ## Annulation de la valeur moyenne
 moy_u_x=assemble(u[0]*dx)/(1-pi*r**2)
 moy_u_y=assemble(u[1]*dx)/(1-pi*r**2)
 moy_u_z=assemble(u[2]*dx)/(1-pi*r**2)
 moy=Function(V)
 moy=Constant((moy_u_x,moy_u_y,moy_u_z))
 chi=project(u-moy,V)
 # Résultat : snapshot
 return(chi)

############################# Pour tester la périodicité d'un champ : impression des valeurs du champ ou de son gradient, ou représentation graphique #############################

def err_per_ind_01(u,cen,r,Npas):# comparaison entre les valeurs individuelles prises par khi aux frontières de la cellule
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
 print('Plan latéral Oyz :')
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
 # ---------------------- khi on Oxz, front and back of the cell ---------------------- #
 # Creates the vectors where the khi component values will be registered
 ufb_y1_0=np.zeros((Npas+1,Npas+1))
 ufb_y2_0=np.zeros((Npas+1,Npas+1))
 ufb_y3_0=np.zeros((Npas+1,Npas+1))
 ufb_y1_1=np.zeros((Npas+1,Npas+1))
 ufb_y2_1=np.zeros((Npas+1,Npas+1))
 ufb_y3_1=np.zeros((Npas+1,Npas+1))
 ## We collect the values of khi on the fluid domain, and suppose khi vanishes on the solid domain
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
 # ---------------------- khi on Oxy, top and bottom of the cell ---------------------- #
 # Creates the vectors where the khi component values will be registered
 ubt_y1_0=np.zeros((Npas+1,Npas+1))
 ubt_y2_0=np.zeros((Npas+1,Npas+1))
 ubt_y3_0=np.zeros((Npas+1,Npas+1))
 ubt_y1_1=np.zeros((Npas+1,Npas+1))
 ubt_y2_1=np.zeros((Npas+1,Npas+1))
 ubt_y3_1=np.zeros((Npas+1,Npas+1))
 ## We collect the values of khi on the fluid domain, and suppose khi vanishes on the solid domain
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
 # ---------------------- khi on Oyz, lateral faces of the cell ---------------------- #
 # Creates the vectors where the khi component values will be registered
 ulr_y1_0=np.zeros((Npas+1,Npas+1))
 ulr_y2_0=np.zeros((Npas+1,Npas+1))
 ulr_y3_0=np.zeros((Npas+1,Npas+1))
 ulr_y1_1=np.zeros((Npas+1,Npas+1))
 ulr_y2_1=np.zeros((Npas+1,Npas+1))
 ulr_y3_1=np.zeros((Npas+1,Npas+1))
 ## We collect the values of khi on the fluid domain, and suppose khi vanishes on the solid domain
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
 # else khi_y.. stays at 0.0
 #
 # ---------------------- plots ---------------------- #
 fig=plt.figure(1)
 # We compare front and back boundaries for khi_y1, khi_y2 and khi_y3
 ## u_y1
 ax1=fig.add_subplot(331, projection='3d')
 #ax1.plot_surface(X,Y,ufb_y1_0,color='green')
 ax1.scatter(X,Y,ufb_y1_0,color='blue')
 ax1.plot_wireframe(X,Y,ufb_y1_1,color='red')
 plt.title("khi_y1 parallel to Oxz")
 ## u_y2
 ax2=fig.add_subplot(332, projection='3d')
 #ax2.plot_surface(X,Y,ufb_y2_0,color='green')
 ax2.scatter(X,Y,ufb_y2_0,color='blue')
 ax2.plot_wireframe(X,Y,ufb_y2_1,color='red')
 plt.title("khi_y2 parallel to Oxz")
 ## u_y3
 ax3=fig.add_subplot(333, projection='3d')
 #ax3.plot_surface(X,Y,ufb_y3_0,color='green')
 ax3.scatter(X,Y,ufb_y3_0,color='blue')
 ax3.plot_wireframe(X,Y,ufb_y3_1,color='red')
 plt.title("khi_y3 parallel to Oxz")
 # We compare top and bottom boundaries for khi_y1, khi_y2 and khi_y3
 ## u_y1
 ax4=fig.add_subplot(334, projection='3d')
 ax4.scatter(X,Y,ubt_y1_0,color='blue')
 #ax4.plot_surface(X,Y,ubt_y1_0,color='green')
 ax4.plot_wireframe(X,Y,ubt_y1_1,color='red')
 plt.title("khi_y1 parallel to Oxy")
 ## u_y2
 ax5=fig.add_subplot(335, projection='3d')
 ax5.scatter(X,Y,ubt_y2_0,color='blue')
 #ax5.plot_surface(X,Y,ubt_y2_0,color='green')
 ax5.plot_wireframe(X,Y,ubt_y2_1,color='red')
 plt.title("khi_y2 parallel to Oxy")
 ## u_y3
 ax6=fig.add_subplot(336, projection='3d')
 #ax6.plot_surface(X,Y,ubt_y3_0,color='green')
 ax6.scatter(X,Y,ubt_y3_0,color='blue')
 ax6.plot_wireframe(X,Y,ubt_y3_1,color='red')
 plt.title("khi_y3 parallel to Oxy")
 # We compare left and right boundaries for khi_y1, khi_y2 and khi_y3
 ## u_y1
 ax7=fig.add_subplot(337, projection='3d')
 #ax7.plot_surface(X,Y,ubt_y1_0,color='green')
 ax7.scatter(X,Y,ubt_y1_0,color='blue')
 ax7.plot_wireframe(X,Y,ulr_y1_1,color='red')
 plt.title("khi_y1 parallel to Oyz")
 ## u_y2
 ax8=fig.add_subplot(338, projection='3d')
 #ax8.plot_surface(X,Y,ubt_y2_0,color='green')
 ax8.scatter(X,Y,ubt_y2_0,color='blue')
 ax8.plot_wireframe(X,Y,ulr_y2_1,color='red')
 plt.title("khi_y2 parallel to Oyz")
 ## u_y3
 ax9=fig.add_subplot(339, projection='3d')
 #ax9.plot_surface(X,Y,ubt_y3_0,color='green')
 ax9.scatter(X,Y,ubt_y3_0,color='blue')
 ax9.plot_wireframe(X,Y,ulr_y3_1,color='red')
 plt.title("khi_y3 parallel to Oyz")
 ## Show or save
 if todo=='aff':
  plt.show()
 elif todo=='save':
  plt.savefig("Figures3D/inc_c"+"CompBo"+str(Npas)+"_cen"+str(cen[0])+str(cen[1])+str(cen[2])+"_ray"+str(r)+".png")
 ## Close
 plt.close()
 #
 return()












