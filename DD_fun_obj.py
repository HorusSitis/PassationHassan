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

def crown(r):#épaisseur de la couronne dans laquelle le maillage est raffiné
 return r+0.01#*(1+0.2*exp(-r**2))#1.2*r

def raffinement_maillage_circ_per(cen,r,mesh):# Objectif : montrer que l'emplacement de l'inclusion périodique dans la cellule élémentaire ne change pas le coefficient de diffusion homogénéisé, calculé avec le tenseur khi
 markers = MeshFunction("bool", mesh, mesh.topology().dim())
 markers.set_all(False)
 # on crée une liste des centres des inclusions voisines de la cellule élémentaire
 l_cen=[]
 for i in range(-1,2):
  for j in range(-1,2):
   l_cen.append([cen[0]+i,cen[1]+j])
 for c in cells(mesh):
  for f in facets(c):
   # Raffinement autour de l'inclusion périodique
   for cen_per in l_cen:
    if (sqrt((f.midpoint()[0]-cen_per[0])**2+(f.midpoint()[1]-cen_per[1])**2)<=crown(r)):
     markers[c]=True
   # Raffinement aux bords du domaine
   if any([f.midpoint()[k]==0 or f.midpoint()[k]==1 for k in range(0,2)]):
    markers[c]=True
 mesh=refine(mesh, markers, redistribute=True)
 return mesh

def creer_maill_circ(cen,r,res):#valable quel que soit la position de l'inclusion : centre, choisi aléatoirement. 1.2*r<0.5.
 # Création de la cellule élémentaire avec inclusion
 rect=Rectangle(Point(-1,-1),Point(2,2))
 domain=rect
 l_cer=[]#Circle(Point(cen[0],cen[1]),r)]
 for i in range(-1,2):
  for j in range(-1,2):
   l_cer.append(Circle(Point(cen[0]+i,cen[1]+j),r)) 
 for cer_per in l_cer:
  domain=domain-cer_per
 domain=domain*Rectangle(Point(0,0),Point(1,1))
 # Création du permier maillage
 mesh=generate_mesh(domain,res)
 # On raffine le long du bord du domaine fluide
 mesh=raffinement_maillage_circ_per(cen,r,mesh)
 #
 return(mesh)

############################# Pour créer des snapshots, inclusion circulaire périodique unique #############################

def snapshot_circ_per(cen,r,res):
 c_x,c_y=cen[0],cen[1]
 mesh_c_r=creer_maill_circ([c_x,c_y],r,res)
 # On pose et on résoud le problème aux éléments finis
 V=VectorFunctionSpace(mesh_c_r, 'P', 3, form_degree=0, constrained_domain=PeriodicBoundary())#vertices))
 ## On définit la bordure du domaine, sur laquelle intégrer le second membre "L" de l'équation en dimension finie
 l_cen=[]
 #cen=[c_x,c_y]
 for i in range(-1,2):
  for j in range(-1,2):
   l_cen.append([cen[0]+i,cen[1]+j])
 #print(l_cen)
 class inclusion_periodique(SubDomain):
  def inside(self,x,on_boundary):
   return (on_boundary and any([between((x[0]-c[0]), (-r-tol, r+tol)) for c in l_cen]) and any([between((x[1]-c[1]), (-r-tol, r+tol)) for c in l_cen]))#points de la frontière du dystème compris dans la boule de centre et rayons cen et r, pour la norme infinie
 ### Utilisation de la classe définie précédemment
 Gamma_sf = inclusion_periodique()
 boundaries = MeshFunction("size_t", mesh_c_r, mesh_c_r.topology().dim()-1)
 boundaries.set_all(0)
 Gamma_sf.mark(boundaries, 5)
 ds = Measure("ds")(subdomain_data=boundaries)
 num_front_cercle=5
 ## On résoud le problème faible, avec une condition de type Neumann au bord de l'obstacle
 normale = FacetNormal(mesh_c_r)
 nb_noeuds=V.dim()
 u = TrialFunction(V)
 v = TestFunction(V)
 a=tr(dot((grad(u)).T, grad(v)))*dx
 L=-dot(normale,v)*ds(num_front_cercle)
 ### Résolution
 u=Function(V)
 solve(a==L,u)#,bc)
 ## Annulation de la valeur moyenne
 moy_u_x=assemble(u[0]*dx)/(1-pi*r**2)
 moy_u_y=assemble(u[1]*dx)/(1-pi*r**2)
 moy=Function(V)
 moy=Constant((moy_u_x,moy_u_y))
 khi=project(u-moy,V)
 # Résultat : snapshot
 return(khi)


############################# Pour avoir quelques représentations graphiques #############################

def fig_khi(cen,r,u,todo):
 # figure : les composantes du vecteur, séparément
 plt.figure(1)
 ## u_y1
 plt.subplot(211)
 spl1=plot(u[0])
 bar1=plt.colorbar(spl1)
 plt.title("khi_y"+str(1))
 ## u_y2
 plt.subplot(212)
 spl2=plot(u[1])
 bar2=plt.colorbar(spl2)
 plt.title("khi_y"+str(2))
 ## Show or save
 if todo=='aff':
  plt.show()
 elif todo=='save':
  plt.savefig("Figures2D/inc_c"+str(cen[0])+str(cen[1])+str(r)+"_khi.png")
 ## Close
 plt.close()
 #
 return()

def fig_dkhi(cen,r,U,todo):
 # figure : les composantes du vecteur, séparément
 plt.figure(1)
 ## U_1_y1
 plt.subplot(221)
 spl1=plot(U[0,0])
 bar1=plt.colorbar(spl1)
 plt.title("-dkhi"+str(1)+"_dy"+str(1))
 ## U_1_y2
 plt.subplot(222)
 spl2=plot(U[1,0])
 bar2=plt.colorbar(spl2)
 plt.title("-dkhi"+str(2)+"_dy"+str(1))
 ## U_2_y1
 plt.subplot(223)
 spl3=plot(U[0,1])
 bar3=plt.colorbar(spl3)
 plt.title("-dkhi"+str(1)+"_dy"+str(2))
 ## U_2_y2
 plt.subplot(224)
 spl4=plot(U[1,1])
 bar4=plt.colorbar(spl4)
 plt.title("-dkhi"+str(2)+"_dy"+str(2))
 ## Show or save
 if todo=='aff':
  plt.show()
 elif todo=='save':
  plt.savefig("Figures2D/inc_c"+str(cen[0])+str(cen[1])+str(r)+"_Gradkhi.png")
 ## Close
 plt.close()
 #
 return()

############################# Pour tester la périodicité d'un champ en norme l2 ou infinie : erreur relative #############################

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
 # Compares and plots between left and right boundary for khi_y1 and khi_y2
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
  plt.savefig("Figures2D/inc_c"+"CompLRBT"+str(Npas)+"_cen"+str(cen[0])+str(cen[1])+"_ray"+str(r)+".png")
 ## Close
 plt.close()
 #
 return()



