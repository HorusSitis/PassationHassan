####################################################################################################################################
## Etape III : en utilisant la methode des snapshots, calcul de la POD et des coefficients aleatoires, toujours dans domaine_fixe ##
####################################################################################################################################

### ------------ Reproduire eventuellement pour des etapes ulterieures. Laisser seulement dans DDD_fun_obj ? ------------ ###

tol=1e-10

xinf=0.0
yinf=0.0
zinf=0.0
xsup=1.0
ysup=1.0
zsup=1.0

r_s_0=0.15
r_v_0=0.15
r_c_0=0.15

r_min=0.05

if res_gmsh==10:
 lw=0.27
elif res_gmsh==20:
 lw=0.15
elif res_gmsh==50:
 lw=0.01

dimension=3

class PeriodicBoundary(SubDomain):
 # Left boundary is "target domain" G
 def inside(self, x, on_boundary):
  return on_boundary and not(near(x[0],xsup,tol) or near(x[1],ysup,tol) or near(x[2],zsup,tol))
 # Map right boundary (H) to left boundary (G)
 def map(self, x, y):
  for i in range(dimension):
   if near(x[i],1.0,tol):
    y[i]=0.0
   else:
    y[i]=x[i]

# maillage du domaine fixe

mesh_dir="maillages_per/3D/"

## inclusions simples ou rayons lies
if dom_fixe=="am":
 mesh_f_name=mesh_dir+"cube_periodique_triangle"+"_"+dom_fixe+"_sur"+str(res_gmsh)+"_fixe.xml"
## inclusions multiples, unique rayon variable
elif dom_fixe=="solid":
 mesh_fixe_prefix=mesh_dir+"cube"+config+"_periodique_triangle_"
 if config=='2sph':
  mesh_f_name=mesh_fixe_prefix+"fixe"+str(int(round(100*r_v_0,2)))+"sur"+str(res_gmsh)+".xml"
 elif config=='cylsph':
  ## rayon du cylindre aux aretes ou de la sphere centrale fixes a 0.15 ##
  if geo_p=='ray_sph':
   mesh_f_name=mesh_fixe_prefix+str(int(round(100*r_c_0,2)))+"fixe"+"sur"+str(res_gmsh)+".xml"
  elif geo_p=='ray_cyl':
   if fixe_comp=='cyl_sph':
    mesh_f_name=mesh_fixe_prefix+"fixe"+str(int(round(100*r_s_0,2)))+"sur"+str(res_gmsh)+".xml"
   elif fixe_comp=='sph_un':
    mesh_f_name=mesh_dir+"cubesphere_periodique_triangle_"+str(int(round(100*r_s_0,2)))+"sur"+str(res_gmsh)+".xml"
elif dom_fixe=="ray_min":
 fixe_comp=True#utilisation du domaine fixe avec annulation du rayon du cylindre dans le fichier general
 if config=='cylsph':
  if geo_p=='ray_sph':
   mesh_f_name=mesh_dir+"cube"+config+"_periodique_triangle_"+str(int(round(100*r_c_0,2)))+str(int(round(100*r_min,2)))+"sur"+str(res_gmsh)+".xml"
  elif geo_p=='ray_cyl':
   mesh_f_name=mesh_dir+"cube"+config+"_periodique_triangle_"+str(int(round(100*r_min,2)))+str(int(round(100*r_s_0,2)))+"sur"+str(res_gmsh)+".xml"

mesh_fixe=Mesh(mesh_f_name)

# fonctions test du domaine fixe

V_fixe=VectorFunctionSpace(mesh_fixe,'P',2,constrained_domain=PeriodicBoundary())

### ------------ Etapes reproduites : dependances directes de Main3D ------------ ###



##
from PO23D import *
##

# Chargement de la matrice des snapshots

u_name='Usnap_'+dom_fixe+'_'+str(Nsnap)+'_'+config+'_'+geo_p+'_'+"res"+str(res)+'_'+ordo+'_'+computer
print(u_name)
with sh.open(repertoire_parent+u_name) as u_loa:
 Usnap = u_loa["maliste"]

#
## matrice de correlation
#print(Usnap[0:5,0:5])
C=mat_corr_temp(V_fixe,Nsnap,Usnap)
B=C-C.T
#print(B[0:10,0:10])
#print(C[0:5,0:5])
#
## Calcul des coefficients aleatoires et la base POD
vp_A_phi=mat_a_mat_phi(Nsnap,Usnap,C,V_fixe,'L2')
#
val_propres=vp_A_phi[0]
Aleat=vp_A_phi[1]
## Attention les objets ranges dans tableau suivant sont des vecteurs
Phi_prime_v=vp_A_phi[2]
#
## Sortie du spectre de la matrice des snapshots, qui doit servir a choisir la taille du modele reduit
#
print("Valeurs propres POD :",val_propres)
#res, res_fixe=20 : energie [71%, 24%, 5%, 0.37%, 0.058%, 0%, 0%, 0%]

#plot()
#plt.show()
#plt.savefig()
#plt.close()

## Enregistrement de la matrice de la base POD, sous la forme vectorielle

phi_name='Phi'+dom_fixe+'_dim'+str(Nsnap)+'_'+config+'_'+geo_p+'_'+"res"+str(res)+'_'+ordo+'_'+computer

with sh.open(repertoire_parent+phi_name) as p_sto:
 p_sto["maliste"] = Phi_prime_v

#sys.exit()#-----------------------------------------------------------------
## Pour reintroduire la base de POD dans l'espace des fonctions definies dans le domaine fixe

#mesh_fixe=Mesh("maillages_per/3D/cubesphere_periodique_triangle_0001fixe.xml")
#V_fixe=V=VectorFunctionSpace(mesh_fixe, 'P', 3, constrained_domain=PeriodicBoundary())

## Tests : orthogonalite ou orthonrmalite de Phi_prime
ui=Function(V_fixe)
uj=Function(V_fixe)

## Orthogonalite
for i in range(Nsnap-1):
 ui.vector().set_local(Phi_prime_v[:,i])
 for j in range(i+1,Nsnap):
  uj.vector().set_local(Phi_prime_v[:,j])
  scal=assemble(dot(ui,uj)*dx)
  print(scal)

## Norme des vecteurs de la base POD, L2 ou n2
for i in range(Nsnap):
 ui.vector().set_local(Phi_prime_v[:,i])
 scal=assemble(dot(ui,ui)*dx)
 norme_L2=sqrt(scal)
 ###
 norme_q=0
 l=len(Phi_prime_v[:,i])
 for k in range(l):
  norme_q+=Phi_prime_v[k,i]**2
 norme_2=sqrt(norme_q)
 #print('norme 2 :',norme_q)
 print('norme L2 :',norme_L2)
 #print('quotient n2/L2 :',scal/norme_q)

# Representation graphique des vecteurs de POD :

## Type de donnees : on veut calculer les fonctions phi_prime_i
## Representation graphique des phi_prime_i :

phi=Function(V_fixe)
for i in range(Nsnap):
 phi.vector().set_local(Phi_prime_v[:,i])
 plot(phi, linewidth=lw)
 r=0.05*(i+1)
 plt.title("Mode "+str(i+1),fontsize=40)
 if fig_todo=='':#aff':
  plt.show()
 else:
  plt.savefig("Figures3D/phi_"+str(i+1)+"_"+config+'_'+geo_p+"_res"+str(res)+".png")
 plt.close()

# Energie et energie cumulee des modes spatiaux, choix du nombre de modes

## Energie et energie cumulee, avec les valeurs propres de la matrice de correlation temporelle
ener_pour=energie_pourcentage(val_propres)[0]
ener_pour_cumul=energie_pourcentage(val_propres)[1]



absc=np.arange(1,Nsnap+1,1)

plt.plot(absc,ener_pour, linewidth=2)
plt.xlabel('valeurs propres')
plt.ylabel('pourcentage_energie')
plt.yscale('log')
if fig_todo=='aff':
 plt.show()
else:
 plt.savefig("Figures3D/ener_vp_"+config+'_'+geo_p+"_res"+str(res)+".png")
plt.close()

plt.plot(absc,ener_pour_cumul, linewidth=2)
plt.xlabel('valeurs propres')
plt.ylabel('pourcentage_energie_cumule')
#plt.yscale('log')
if fig_todo=='aff':
 plt.show()
else:
 plt.savefig("Figures3D/ener_cumul_vp_"+config+'_'+geo_p+"_res"+str(res)+".png")
plt.close()

## Choix du nombre de modes, avec une valeur seuil d'energie a atteindre avec les vacteurs de la base POD
nb_modes=0

seuil_ener=99.99#9

i=0
while ener_pour_cumul[i]<seuil_ener or i==0:
 nb_modes=i+1
 i+=1

Nseuil=i


print('Energie '+str(seuil_ener)+' pourcent, '+conf_mess+', '+geo_mess+' , resolution '+str(res_gmsh)+' :')
print(Nseuil)
