# -*- coding: utf-8 -*-

from fenics import *

from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
#import numpy as np
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
 return r*(1+0.01*exp(r))#1.2*r

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
 if cen[0]+crown(r)<1 and cen[1]+crown(r)<1 and cen[2]+crown(r)<1 and cen[0]-crown(r)>0 and cen[1]-crown(r)>0 and cen[1]-crown(r)>0:#si l'inclusion est comprise dans la cellule
  box=Box(Point(0,0,0),Point(1,1,1))
  sphere=Sphere(Point(cen[0],cen[1],cen[2]),r)
  domain=box-sphere
  mesh=generate_mesh(domain,res)
  # On raffine le long du bord de l'inclusion
  #mesh_aux=raffinement_maillage_cellule_centree(r,mesh)
  mesh_aux=raffinement_maillage_sph_per(cen,r,mesh)
  mesh=mesh_aux
  #print('pfffrrh !')
 else:
  # Création de la cellule élémentaire avec inclusion
  box=Box(Point(-1,-1,-1),Point(2,2,2))
  domain=box#*Box(Point(0,0,0),Point(1,1,1))
  l_sph=[]
  for i in range(-1,2):
   for j in range(-1,2):
    for k in range(-1,2):
     l_sph.append(Sphere(Point(cen[0]+i,cen[1]+j,cen[2]+k),r))
  for sph_per in l_sph:#[Sphere(Point(cen[0],cen[1],cen[2]),r)]:#l_sph:
   domain=domain-sph_per
  #domain=domain-Sphere(Point(cen[0],cen[1],cen[2]),r)
  domain=domain*Box(Point(0,0,0),Point(1,1,1))
  # Création du permier maillage
  mesh=generate_mesh(domain,res)
  # On raffine le long du bord de l'inclusion
  mesh_aux=raffinement_maillage_sph_per(cen,r,mesh)
  mesh=mesh_aux
  #
 return(mesh)

############################# Pour créer des snapshots, inclusion circulaire périodique unique #############################

def snapshot_sph_per(cen,r,res):
 c_x,c_y,c_z=cen[0],cen[1],cen[2]
 mesh_c_r=creer_maill_sph([c_x,c_y,c_z],r,res)
 # On pose et on résoud le problème aux éléments finis
 V=VectorFunctionSpace(mesh_c_r, 'P', 2, constrained_domain=PeriodicBoundary())
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
 boundaries = MeshFunction("size_t", mesh_c_r, mesh_c_r.topology().dim()-1)
 # On attribue une valeur par défaut aux frontières du domaine fluide, qui concerne plus particulièrement l'interface fluide-fluide
 boundaries.set_all(1)
 Gamma_sf.mark(boundaries, 7)
 ds = Measure("ds")(subdomain_data=boundaries)
 num_ff=1
 num_front_sphere=7
 ## On fixe une condition de Dirichlet, pour avoir l'unicité du tenseur khi : non nécessaire pour le tenseur de diffusion homogénéisé
 khi_zero=Constant((0., 0., 0.))
 # Point auquel on fixe khi_zero
 point_zero=[0,0,0]
 for i in [-0.5,0.5]:
  for j in [-0.5,0.5]:
   for k in [-0.5,0.5]:
    point_zero_prov=[cen[0]+i,cen[1]+j,cen[2]+k]
    if not(any([not(between(coord,(0-tol,1+tol))) for coord in point_zero_prov])):
     point_zero=point_zero_prov
 def supp_point_zero(x):
  return between(x[0],(point_zero[0]-tol,point_zero[0]+tol)) and between(x[1],(point_zero[1]-tol,point_zero[1]+tol)) and between(x[2],(point_zero[2]-tol,point_zero[2]+tol))
 print(point_zero)
 bc = DirichletBC(V, khi_zero, supp_point_zero, "pointwise")
 ## On résoud le problème faible, avec une condition de type Neumann au bord de l'obstacle
 normale = FacetNormal(mesh_c_r)
 nb_noeuds=V.dim()
 u = TrialFunction(V)
 v = TestFunction(V)
 a=tr(dot((grad(u)).T, grad(v)))*dx
 L=-dot(normale,v)*ds(num_front_sphere)
 ### Résolution
 u=Function(V)
 solve(a==L,u,bc)
 # Résultat : snapshot
 return(u)

############################# Pour tester la périodicité d'un champ en norme l2 ou infinie : erreur relative #############################

def err_per_01(u,norm,Npas,type_err):
 pas=1/Npas
 # Périodicité
 ## Essais avec assemble() ou .vector().get_local()
 #mesh_Gamma_ff=UnitIntervalMesh(40)
 #x=SpatialCoordinate(mesh_Gamma_ff)
 #def quad_u_x(y):return(sum((u((1,y))-u((0,y)))**2))
 #err_per_x=assemble(quad_u_x(x[0])*dx(degree=2))
 #uv=u.vector().get_local()
 ## Normes l2 et infinie de khi restreint à \{0\}\times ... , respectivement ... \times [0,1] pour l'erreur sur le bord vertical, respectivement horizontal.
 list_per_x=[sum((u((1,pas*k))-u((0,pas*k)))**2) for k in range(0,Npas)]
 list_0_x=[sum(u((0,pas*k))**2) for k in range(0,Npas)]
 list_per_y=[sum((u((pas*k,1))-u((pas*k,0)))**2) for k in range(0,Npas)]
 list_0_y=[sum(u((pas*k,0))**2) for k in range(0,Npas)]
 list_per_z=[sum((u((pas*k,1))-u((pas*k,0)))**2) for k in range(0,Npas)]
 list_0_z=[sum(u((pas*k,0))**2) for k in range(0,Npas)]
 ## Erreur relative ou absolue
 den=[1.0,1.0,1.0]
 if type_err=='rel':
  if norm=='l2':
   den=[sqrt(sum(list_0_x)/Npas),sqrt(sum(list_0_y)/Npas)]
  elif norm=='infty':
   den=[sqrt(max(list_0_x)),sqrt(max(list_0_y))]
 ## Résultats
 if norm=='l2':
  El2_per_x=sqrt(sum(list_per_x)/Npas)/den[0]
  El2_per_y=sqrt(sum(list_per_y)/Npas)/den[1]
  return((El2_per_x,El2_per_y))
 elif norm=='infty':
  Einfty_per_x=sqrt(max(list_per_x))/den[0]
  Einfty_per_y=sqrt(max(list_per_y))/den[1]
  return((Einfty_per_x,Einfty_per_y))

def err_per_ind_01(u,Npas):# comparaison entre les valeurs individuelles prises par khi aux frontières de la cellule
 pas=1/Npas
 print('Plan frontal Oxz :')
 for k in range(0,1+Npas):
  for l in range(0,1+Npas):
   print('x='+str(pas*k),'z='+str(pas*l),u((pas*k,0.0,pas*l)),u((pas*k,1.0,pas*l)))
 print('Plan horizontal Oxy :')
 for k in range(0,1+Npas):
  for l in range(0,1+Npas):
   print('x='+str(pas*k),'y='+str(pas*l),u((pas*k,pas*l,0.0)),u((pas*k,pas*l,1.0)))
 print('Plan latéral Oyz :')
 for k in range(0,1+Npas):
  for l in range(0,1+Npas):
   print('y='+str(pas*k),'z='+str(pas*l),u((0.0,pas*k,pas*l)),u((1.0,pas*k,pas*l)))
 return()

