####################################################################################################################################
## Etape III : en utilisant la methode des snapshots, calcul de la POD et des coefficients aleatoires, toujours dans domaine_fixe ##
####################################################################################################################################

if dom_fixe=='am':
    mesh_fixe=Mesh('maillages_per/2D/maillage_fixe2D_am.xml')
elif dom_fixe=='multiray':
    mesh_fixe=Mesh('maillages_per/2D/maillage_fixe2d_'+dom_fixe+'.xml')
elif config=='compl':
    mesh_fixe=Mesh('maillages_per/2D/maillage_trous2D_'+geo_p+'_fixe.xml')
elif dom_fixe=='ray_min':
    if config=='cer_un':
        mesh_fixe=Mesh('maillages_per/2D/maillage_trou2D_5.xml')

print('maillages_per/2D/maillage_trous2D_'+geo_p+'_fixe.xml')

V_fixe=VectorFunctionSpace(mesh_fixe, 'P', VFS_degree, constrained_domain=PeriodicBoundary())
nb_noeuds = V_fixe.dim()


## Chargement de la marice des snapshots

u_name='Usnap_'+dom_fixe+str(N_snap)+'_'+config+'_'+geo_p+'_deg'+str(VFS_degree)+'_'+ordo+'_'+computer
print(repertoire_parent+u_name)

with sh.open(repertoire_parent+u_name) as u_loa:
    Usnap = u_loa['maliste']

# matrice de correlation

C=mat_corr_temp(V_fixe,N_snap,Usnap)

# Calcul des coefficients aleatoires et la base POD

vp_A_phi=mat_a_mat_phi(N_snap,Usnap,C,V_fixe,'n2')
vp_A_phi=mat_a_mat_phi(N_snap,Usnap,C,V_fixe,'L2')
# vp_A_phi=pod.mat_a_mat_phi(N_snap,Usnap,C,'')

val_propres=vp_A_phi[0]
Aleat=vp_A_phi[1]
## Attention les objets ranges dans tableau suivant sont des vecteurs
Phi_prime_v=vp_A_phi[2]

## Enregistrement de la matrice de la base POD, sous la forme vectorielle

phi_name='Phi'+dom_fixe+'_dim'+str(N_snap)+'_'+config+'_'+geo_p+'_deg'+str(VFS_degree)+'_'+'res'+str(res)+'_'+ordo+'_'+computer
print(phi_name)

with sh.open(repertoire_parent+phi_name) as p_sto:
    p_sto['maliste'] = Phi_prime_v

## Tests : orthogonalite ou orthonrmalite de Phi_prime
ui=Function(V_fixe)
uj=Function(V_fixe)

## Orthogonalite
for i in range(N_snap-1):
    ui.vector().set_local(Phi_prime_v[:,i])
    for j in range(i+1,N_snap):
        uj.vector().set_local(Phi_prime_v[:,j])
        scal=assemble(dot(ui,uj)*dx)
        print(scal)

## Norme des vacteurs dela base POD, L2 ou n2
for i in range(N_snap):
    ui.vector().set_local(Phi_prime_v[:,i])
    scal=assemble(dot(ui,ui)*dx)
    norme_L2=sqrt(scal)
    ###
    norme_q=0
    l=len(Phi_prime_v[:,i])
    for k in range(l):
        norme_q+=Phi_prime_v[k,i]**2
        norme_2=sqrt(norme_q)
    print('norme L2 :',norme_L2)


# Representation graphique des vecteurs de POD :

## Type de donnees : on veut calculer les fonctions phi_prime_i
## Representation graphique des phi_prime_i :

phi=Function(V_fixe)
for i in range(N_snap):
    phi.vector().set_local(Phi_prime_v[:,i])
    plot(phi, linewidth=0.08)
    if fig_todo=='aff':
        plt.title('Phi '+str(i+1)+' sur domaine fixe')
        plt.show()
    else:
        plt.savefig('Figures2D/phi_'+str(i+1)+'_'+config+'_'+geo_p+'.png')
    plt.close()



# Energie et energie cumulee des modes spatiaux, choix du nombre de modes

## Energie et energie cumulee, avec les valeurs propres de la matrice de correlation temporelle
ener_pour=energie_pourcentage(val_propres)[0]
ener_pour_cumul=energie_pourcentage(val_propres)[1]

absc=np.arange(1,N_snap+1,1)

plt.plot(absc,ener_pour)
plt.xlabel('valeurs propres')
plt.ylabel('pourcentage_energie')
plt.yscale('log')
if fig_todo=='aff':
    plt.show()
else:
    plt.savefig('Figures2D/ener_vp_'+config+'_'+geo_p+'.png')#+'_res'+str(res)
plt.close()

plt.plot(absc,ener_pour_cumul)
plt.xlabel('valeurs propres')
plt.ylabel('pourcentage_energie_cumule')
if fig_todo=='aff':
    plt.show()
else:
    plt.savefig('Figures2D/ener_cumul_vp_'+config+'_'+geo_p+'.png')#+'_res'+str(res)
plt.close()

## Choix du nombre de modes, avec une valeur seuil d'energie a atteindre avec les vacteurs de la base POD
nb_modes=0

# seuil_ener_pour=99.99

i=0
while ener_pour_cumul[i]<seuil_ener_pour:
    nb_modes=i+1
    i+=1


Nseuil=i

print('1-nu '+str(seuil_ener_pour)+' pourcent :')
print('Nrom = ', Nseuil)

## Tests sur les fonctions POD et leurs valeurs propres associeess

print('Valeurs propres :',val_propres)

# list_DPhi=[]
# ui=Function(V_fixe)
#
# for i in range(Nseuil):
#     ui.vector().set_local(Phi_prime_v[:,i])
#     #
#     int_grad=assemble(grad(ui)[0,0]*dx)
#     list_DPhi.append(1+int_grad)
#
# print('Dhom POD-ROM :',list_DPhi)
