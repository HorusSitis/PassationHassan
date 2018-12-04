from fenics import *

from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
#import numpy as np
from math import sqrt
import sys



if __name__=='__main__':
 #
 tol=1e-10
 xinf=0.0
 yinf=0.0
 xsup=1.0
 ysup=1.0

tol=1e-10
xinf=0.0
yinf=0.0
xsup=1.0
ysup=1.0

dimension=2

class PeriodicBoundary(SubDomain):
 # Left boundary is "target domain" G
 def inside(self, x, on_boundary):
  return on_boundary and not(near(x[0],xsup,tol) or near(x[1],ysup,tol))
 # Map right boundary (H) to left boundary (G)
 def map(self, x, y):
  for i in range(dimension):
   if near(x[i],1.0,tol):
    y[i]=0.0
   else:
    y[i]=x[i]

############################# Pour créer des maillages, avec des familles de cellules élémentaires #############################

def raffinement_maillage_cellule_centree(r,mesh):# Cellule centrée : un seul paramètre géométrique, le rayon de l'inclusion. 1.2*r<0.5 exigé.
 markers = MeshFunction("bool", mesh, mesh.topology().dim())
 markers.set_all(False)
 c_x=0.5
 c_y=0.5
 for c in cells(mesh):
 # Mark cells with facet midpoints near y == 1.0
  for f in facets(c):
   if (sqrt((f.midpoint()[0]-c_x)**2+(f.midpoint()[1]-c_y)**2)<=1.2*r):
    markers[c] = True 
 mesh=refine(mesh, markers, redistribute=True) 
 return mesh

def raffinement_maillage_circ_per(cen,r,mesh,test):# Objectif : montrer que l'emplacement de l'inclusion périodique dans la cellule élémentaire ne change pas le coefficient de diffusion homogénéisé, calculé avec le tenseur khi
 markers = MeshFunction("bool", mesh, mesh.topology().dim())
 markers.set_all(False)
 # on crée une liste des centres des inclusions voisines de la cellule élémentaire
 l_cen=[]
 for i in range(-1,2):
  for j in range(-1,2):
   l_cen.append([cen[0]+i,cen[1]+j])
 for c in cells(mesh):
 # Mark cells with facet midpoints near ...
  for f in facets(c):
   for cen_per in l_cen:
    if (sqrt((f.midpoint()[0]-cen_per[0])**2+(f.midpoint()[1]-cen_per[1])**2)<=1.2*r):
     markers[c] = True
    if test=='test':
     if between(f.midpoint()[0],(0.2-tol,0.8+tol)) and between(f.midpoint()[1],(0.8-tol,0.9+tol)):
      markers[c]=True
 mesh=refine(mesh, markers, redistribute=True)
 return mesh

def creer_maill_circ(cen,r,res,test):#valable quel que soit la position de l'inclusion : centre, choisi aléatoirement. 1.2*r<0.5.
 if cen[0]+r*1.2<1 and cen[1]+r*1.2<1 and cen[0]-r*1.2>0 and cen[1]-r*1.2>0:#si l'inclusion est comprise dans la cellule
  rect=Rectangle(Point(0,0),Point(1,1))
  circle=Circle(Point(cen[0],cen[1]),r)
  domain=rect-circle
  if test=='test':
   domain=domain-Rectangle(Point(0.2,0.8),Point(0.8,0.9))
  mesh=generate_mesh(domain,res)
  # On raffine le long du bord de l'inclusion
  #mesh_aux=raffinement_maillage_cellule_centree(r,mesh)
  mesh_aux=raffinement_maillage_circ_per(cen,r,mesh,test)
  mesh=mesh_aux
  #print('pfffrrh !')
 else:
  # Création de la cellule élémentaire avec inclusion
  rect=Rectangle(Point(-1,-1),Point(2,2))
  domain=rect
  l_cer=[]#Circle(Point(cen[0],cen[1]),r)]
  for i in range(-1,2):
   for j in range(-1,2):
    l_cer.append(Circle(Point(cen[0]+i,cen[1]+j),r)) 
  print(len(l_cer))
  for cer_per in l_cer:
   domain=domain-cer_per
  domain=domain*Rectangle(Point(0,0),Point(1,1))
  if test=='test':
   domain=domain-Rectangle(Point(0.2,0.8),Point(0.8,0.9))
  # Création du permier maillage
  mesh=generate_mesh(domain,res)
  # On raffine le long du bord de l'inclusion
  mesh_aux=raffinement_maillage_circ_per(cen,r,mesh,test)
  mesh=mesh_aux
  #
 return(mesh)

############################# Pour créer des snapshots, inclusion circulaire périodique unique #############################

def snapshot_circ_per(cen,r,res,test):
 c_x,c_y=cen[0],cen[1]
 mesh_c_r=creer_maill_circ([c_x,c_y],r,res,test)
 # On pose et on résoud le problème aux éléments finis
 V=VectorFunctionSpace(mesh_c_r, 'P', 2, constrained_domain=PeriodicBoundary())
 ## On définit la bordure du domaine, sur laquelle intégrer le second membre "L" de l'équation en dimension finie
 l_cen=[]
 #cen=[c_x,c_y]
 for i in range(-1,2):
  for j in range(-1,2):
   l_cen.append([cen[0]+i,cen[1]+j])
 #print(l_cen)
 class inclusion_periodique(SubDomain):
  def inside(self,x,on_boundary):
   return (on_boundary and any([between((x[0]-c[0]), (-r-tol, r+tol)) for c in l_cen]) and any([between((x[1]-c[1]), (-r-tol, r+tol)) for c in l_cen]))#[between(sqrt((x[0]-c[0])**2+(x[1]-c[1])**2),(r-tol,r+tol)) for c in l_cen]))
 if test=='test':
  class inclusion_test(SubDomain):
   def inside(self,x,on_boundary):
    return (on_boundary and between(x[0],(0.2-tol,0.8+tol)) and between(x[1],(0.8-tol,0.9+tol)))
 ### Utilisation de la classe définie précédemment
 Gamma_sf = inclusion_periodique()
 boundaries = MeshFunction("size_t", mesh_c_r, mesh_c_r.topology().dim()-1)
 boundaries.set_all(0)
 Gamma_sf.mark(boundaries, 5)
 if test=='test':
  Gamma_test=inclusion_test()
  Gamma_test.mark(boundaries, 5)
 ds = Measure("ds")(subdomain_data=boundaries)
 num_front_cercle=5
 ## On fixe une condition de Dirichlet, pour avoir l'unicité du tenseur khi : non nécessaire pour le tenseur de diffusion homogénéisé
 khi_bord=Constant((0., 0.))
 bc = DirichletBC(V, khi_bord, "x[0] < DOLFIN_EPS && x[1] < DOLFIN_EPS", "pointwise")
 ## On résoud le problème faible, avec une condition de type Neumann au bord de l'obstacle
 normale = FacetNormal(mesh_c_r)
 nb_noeuds=V.dim()
 u = TrialFunction(V)
 v = TestFunction(V)
 a=tr(dot((grad(u)).T, grad(v)))*dx
 L=-dot(normale,v)*ds(num_front_cercle)
 ### Résolution
 u=Function(V)
 solve(a==L,u,bc)
 # Résultat : snapshot
 return(u)

############################# Pour tester la périodicité d'un champ en norme l2 ou infinie : erreur relative #############################

def err_per_sum_01(u,norm,Npas,type_err):
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
 ## Erreur relative ou absolue
 den=[1.0,1.0]
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

def err_per_ind_01(u,Npas):
 pas=1/Npas
 #l_x_0=array([u(0.0,k*pas for k in range(0,Npas)])
 #l_x_1=array([u(1.0,k*pas for k in range(0,Npas)])
 #l_x_diff=)l_x_1-l_x_0
 print('Bord vertical :')
 for k in range(0,1+Npas):
  print(u((0.0,pas*k)),u((1.0,pas*k)))
 print('Bord horizontal :')
 for k in range(0,1+Npas):
  print(u((pas*k,0.0)),u((pas*k,1.0)))
 return()
