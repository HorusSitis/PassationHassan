from fenics import *

from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
#import numpy as np
from math import sqrt
import sys

tol=1e-10
xinf=0.0
yinf=0.0
zinf=0.0
xsup=1.0
ysup=1.0
zsup=1.0


############################# Une classe de sous domaines pour tous les calculs : structure périodique avec répétition de la cellule élémentaire #############################

class PeriodicBoundary(SubDomain):
 # Left boundary is "target domain" G
 def inside(self, x, on_boundary):
  return on_boundary and (near(x[0],xinf,tol) or near(x[1],yinf,tol) or near(x[2],zinf,tol))
 # Map right boundary (H) to left boundary (G)
 def map(self, x, y):
  if (near(x[0],xsup,tol)):
   y[0] = x[0] - 1.0
   y[1] = x[1]
   y[2] = x[2]        
  elif (near(x[0],xsup,tol)):
   y[0] = x[0]
   y[1] = x[1] - 1.0
   y[2] = x[2]
  else:
   y[0] = x[0]
   y[1] = x[1]
   y[2] = x[2] - 1.0

############################# Pour créer des maillages, avec des familles de cellules élémentaires #############################

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
 # Mark cells with facet midpoints near ...
  for f in facets(c):
   for cen_per in l_cen:
    if (sqrt((f.midpoint()[0]-cen_per[0])**2+(f.midpoint()[1]-cen_per[1])**2+(f.midpoint()[2]-cen_per[2])**2)<=1.2*r):
     markers[c] = True
 mesh=refine(mesh, markers, redistribute=True)
 return mesh

def creer_maill_sph(cen,r,res):#valable quel que soit la position de l'inclusion : centre, choisi aléatoirement. 1.2*r<0.5.
 if cen[0]+r*1.2<1 and cen[1]+r*1.2<1 and cen[2]+r*1.2<1 and cen[0]-r*1.2>0 and cen[1]-r*1.2>0 and cen[1]-r*1.2>0:#si l'inclusion est comprise dans la cellule
  box=Box(Point(0,0,0),Point(1,1,1))
  sphere=Sphere(Point(cen[0],cen[1],cen[2]),r)
  domain=box-sphere
  mesh=generate_mesh(domain,res)
  # On raffine le long du bord de l'inclusion
  #mesh_aux=raffinement_maillage_cellule_centree(r,mesh)
  mesh_aux=raffinement_maillage_sph_per(cen,r,mesh)
  mesh=mesh_aux
  print('pfffrrh !')
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
 ## On définit la bordure du domaine, sur laquelle intégrer le second membre "L" de l'équation en dimension finie
 l_cen=[]
 #cen=[c_x,c_y]
 for i in range(-1,2):
  for j in range(-1,2):
   for k in range(-1,2):
    l_cen.append([cen[0]+i,cen[1]+j,cen[2]+k])
 class inclusion_periodique(SubDomain):
  def inside(self,x,on_boundary):
   return (on_boundary and any([between((x[0]-c[0]), (-r-tol, r+tol)) for c in l_cen]) and any([between((x[1]-c[1]), (-r-tol, r+tol)) for c in l_cen]) and any([between((x[2]-c[2]), (-r-tol, r+tol)) for c in l_cen]))#[between(sqrt((x[0]-c[0])**2+(x[1]-c[1])**2),(r-tol,r+tol)) for c in l_cen]))
 ### Utilisation de la classe définie précédemment
 Gamma_sf = inclusion_periodique()
 boundaries = MeshFunction("size_t", mesh_c_r, mesh_c_r.topology().dim()-1)
 boundaries.set_all(0)
 Gamma_sf.mark(boundaries, 5)
 ds = Measure("ds")(subdomain_data=boundaries)
 num_front_sphere=5
 ## On résoud le problème en fixant les condisions aux limites
 khi_bord=Constant((0., 0., 0.))
 bc = DirichletBC(V, khi_bord, "x[0] < DOLFIN_EPS && x[1] < DOLFIN_EPS && x[2] < DOLFIN_EPS", "pointwise")
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




