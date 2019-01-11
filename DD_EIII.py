####################################################################################################################################
## Etape III : en utilisant la méthode des snapshots, calcul de la POD et des coefficients aléatoires, toujours dans domaine_fixe ##
####################################################################################################################################

#mesh_fixe = Mesh("Solutions_homog_interp_circulaire/mesh_circulaire.xml.gz")



# Domaine d'interpolation et matrice des snapshots

c_x=0.5
c_y=0.5

res=25

Nsnap=8

#V_fixe=VectorFunctionSpace(mesh_fixe, "P", 2, constrained_domain=PeriodicBoundary())
nb_noeuds = V_fixe.dim()

Usnap=np.zeros((nb_noeuds,Nsnap))

def inter_snap_ray(n):
 r=n*0.05
 # Cliché sur le domaine avec inclusion
 u=snapshot_circ_per([c_x,c_y],r,res)
 # Extrapolation du cliché : khi prime
 u.set_allow_extrapolation(True)
 u_fixe=interpolate(u,V_fixe)
 # Forme vectorielle de la solution EF
 u_fixe_v=u_fixe.vector().get_local()
 return([n,u_fixe_v])#[n,u_fixe])

# Génération séquentielle des snapshots

for n in range(1,1+Nsnap):
 u_fixe_v=inter_snap_ray(n)[1]
 # Remplissage de la matrice des snapshots
 Usnap[:,n-1]=u_fixe_v#.vector().get_local()

## UsnapSeq=Usnap

# Génération parallèle des snapshots

pool=multiprocessing.Pool(processes=8)

Uf_par=pool.map(inter_snap_ray,(n for n in range(1,1+Nsnap)))

for n in range(1,1+Nsnap):
 for i in range(0,Nsnap):
  if Uf_par[i][0]==n:
   u_fixe_v=Uf_par[i][1]
   Usnap[:,n-1]=u_fixe_v#.vector().get_local()

## UsnapPar=Usnap

# matrice de corrélation

## Usnap=UsnapSeq
## Usnap=UsnapPar

C=pod.mat_corr_temp(V_fixe,Nsnap,Usnap)

# Calcul des coefficients aléatoires et la base POD

vp_A_phi=mat_a_mat_phi(Nsnap,Usnap,C,V_fixe,'n2')
vp_A_phi=mat_a_mat_phi(Nsnap,Usnap,C,V_fixe,'L2')
#vp_A_phi=pod.mat_a_mat_phi(Nsnap,Usnap,C,'')

val_propres=vp_A_phi[0]
Aleat=vp_A_phi[1]
## Attention les objets rangés dans tableau suivant sont des vecteurs
Phi_prime_v=vp_A_phi[2]

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

## Norme des vacteurs dela base POD, L2 ou n2
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
 plt.show()
 plt.close()

# Energie et énergie cumulée des modes spatiaux, choix du nombre de modes

## Energie et énergie cumulée, avec les valeurs propres de la matrice de corrélation temporelle
ener_pour=energie_pourcentage(val_propres)[0]
ener_pour_cumul=energie_pourcentage(val_propres)[1]

absc=np.arange(1,Nsnap+1,1)

plt.plot(absc,ener_pour)
plt.xlabel('valeurs propres')
plt.ylabel('pourcentage_energie')
plt.show()

plt.plot(absc,ener_pour_cumul)
plt.xlabel('valeurs propres')
plt.ylabel('pourcentage_energie_cumule')
plt.show()

## Choix du nombre de modes, avec une valeur seuil d'énergie à atteindre avec les vacteurs de la base POD
nb_modes=0

seuil_ener=99.999

i=0
while ener_pour_cumul[i]<seuil_ener:
 nb_modes=i+1
 i+=1

### 8 snapshots : 4 modes pour un seuil de 99.9%
### 8 snapshots : 7 modes pour un seuil de 99,999%
