####################################################################################################################################
## Etape III : en utilisant la méthode des snapshots, calcul de la POD et des coefficients aléatoires, toujours dans domaine_fixe ##
####################################################################################################################################

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

#import PO23D as pod
#pod = reload(pod)
from PO23D import *
#

mesh_fixe=Mesh("maillages_per/3D/cubesphere_periodique_triangle_0001fixe.xml")

V_fixe=V=VectorFunctionSpace(mesh_fixe, 'P', 2, constrained_domain=PeriodicBoundary())

## Chargement de la marice des snapshots

u_name='Usnap_'+str(Nsnap)+'_'+config+'_'+geo_p+'_'+ordo+'_'+computer

with sh.open(repertoire_parent+u_name) as u_loa:
    Usnap = u_loa["maliste"]

#
## matrice de corrélation
C=mat_corr_temp(V_fixe,Nsnap,Usnap)
#
## Calcul des coefficients aléatoires et la base POD
vp_A_phi=mat_a_mat_phi(Nsnap,Usnap,C,V_fixe,'L2')
#
val_propres=vp_A_phi[0]
Aleat=vp_A_phi[1]
## Attention les objets rangés dans tableau suivant sont des vecteurs
Phi_prime_v=vp_A_phi[2]
#
## Sortie du spectre de la matrice des snapshots, qui doit servir à choisir la taille du modèle réduit
#
print(val_propres)
#res, res_fixe=20 : énergie [71%, 24%, 5%, 0.37%, 0.058%, 0%, 0%, 0%]

#plot()
#plt.show()
#plt.savefig()
#plt.close()

## Enregistrement de la matrice de la base POD, sous la forme vectorielle

phi_name='Phi_dim'+str(Nsnap)+'_'+config+'_'+geo_p+'_'+ordo+'_'+computer

with sh.open(repertoire_parent+phi_name) as p_sto:
    p_sto["maliste"] = Phi_prime_v

#sys.exit()#-----------------------------------------------------------------
## Pour réintroduire la base de POD dans l'espace des fonctions définies dans le domaine fixe

mesh_fixe=Mesh("maillages_per/3D/cubesphere_periodique_triangle_0001fixe.xml")
V_fixe=V=VectorFunctionSpace(mesh_fixe, 'P', 2, constrained_domain=PeriodicBoundary())

## Tests : orthogonalité ou orthonrmalité de Phi_prime
ui=Function(V_fixe)
uj=Function(V_fixe)

## Orthogonalité
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

# Représentation graphique des vecteurs de POD :

## Type de données : on veut calculer les fonctions phi_prime_i 
## Représentation graphique des phi_prime_i :

phi=Function(V_fixe)
for i in range(Nsnap):
 phi.vector().set_local(Phi_prime_v[:,i])
 plot(phi)
 if fig_todo=='aff':
  plt.show()
 else:
  plt.savefig("phi_"+str(i+1)+"_"+config+geo_p+".png")
 plt.close()

# Energie et énergie cumulée des modes spatiaux, choix du nombre de modes

## Energie et énergie cumulée, avec les valeurs propres de la matrice de corrélation temporelle
ener_pour=energie_pourcentage(val_propres)[0]
ener_pour_cumul=energie_pourcentage(val_propres)[1]

absc=np.arange(1,Nsnap+1,1)

plt.plot(absc,ener_pour)
plt.xlabel('valeurs propres')
plt.ylabel('pourcentage_energie')
if fig_todo=='aff':
 plt.show()
else:
 plt.savefig("Figures3D/ener_vp_"+config+geo_p+".png")
plt.close()

plt.plot(absc,ener_pour_cumul)
plt.xlabel('valeurs propres')
plt.ylabel('pourcentage_energie_cumule')
if fig_todo=='aff':
 plt.show()
else:
 plt.savefig("Figures3D/ener_cumul_vp_"+config+geo_p+".png")
plt.close()

## Choix du nombre de modes, avec une valeur seuil d'énergie à atteindre avec les vacteurs de la base POD
nb_modes=0

seuil_ener=99.99#9

i=0
while ener_pour_cumul[i]<seuil_ener:
 nb_modes=i+1
 i+=1

Nseuil=i


print(str(seuil_ener)+':')
print(Nseuil)
 
