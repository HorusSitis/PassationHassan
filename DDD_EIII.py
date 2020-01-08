####################################################################################################################################
## Etape III : en utilisant la methode des snapshots, calcul de la POD et des coefficients aleatoires, toujours dans domaine_fixe ##
####################################################################################################################################

# maillage du domaine fixe

mesh_dir='maillages_per/3D/'

## inclusions simples ou rayons lies

if dom_fixe == 'am':
    mesh_f_name = 'cube_periodique_triangle' + '_' + dom_fixe + '_sur' + str(res_gmsh) + 'fixe'

## inclusions multiples, unique rayon variable
elif dom_fixe == 'solid':
    # mesh_fixe_prefix ='cube' + config + '_periodique_triangle_'
    mesh_fixe_prefix = mesh_prefix
    if config == '2sph':
        mesh_f_name = mesh_fixe_prefix + 'fixe_som' + '_rayp' + str(int(round(100*ray_fix,2))) + '_sur' + str(res_gmsh)

elif dom_fixe == 'ray_min':
    # utilisation du domaine fixe avec annulation du rayon du cylindre dans le fichier general
    fixe_comp = True
    if config == 'cylsph':
        if geo_p == 'ray_sph':
            mesh_f_name = mesh_prefix + 'rayc' + str(int(round(100*rho_appr_min,2))) + '_rayp' + str(int(round(100*ray_fix,2))) + '_sur' + str(res_gmsh)

        elif geo_p == 'ray_cyl':
            mesh_f_name = mesh_prefix + 'rayc' + str(int(round(100*ray_fix,2))) + '_rayp' + str(int(round(100*rho_appr_min,2))) + '_sur' + str(res_gmsh)

mesh_fixe = Mesh(mesh_repository + mesh_f_name + '.xml')

# fonctions test du domaine fixe

V_fixe = VectorFunctionSpace(mesh_fixe, 'P', 2, constrained_domain=PeriodicBoundary())

### ------------ Etapes reproduites : dependances directes de Main3D ------------ ###



##
from PO23D import *
##

# Chargement de la matrice des snapshots

u_name = 'Usnap_' + dom_fixe + '_' + str(N_snap) + '_' + config + '_' + geo_p + '_' + 'res' + str(res) + '_' + ordo + '_' + computer
print(u_name)
with sh.open(repertoire_parent+u_name) as u_loa:
    Usnap = u_loa['maliste']

#
## matrice de correlation
#print(Usnap[0:5,0:5])
C = mat_corr_temp(V_fixe, N_snap, Usnap)
B = C-C.T
#print(B[0:10,0:10])
#print(C[0:5,0:5])
#
## Calcul des coefficients aleatoires et la base POD
vp_A_phi = mat_a_mat_phi(N_snap, Usnap, C, V_fixe, 'L2')
#
val_propres = vp_A_phi[0]
Aleat = vp_A_phi[1]
## Attention les objets ranges dans tableau suivant sont des vecteurs
Phi_prime_v = vp_A_phi[2]
#
## Sortie du spectre de la matrice des snapshots, qui doit servir a choisir la taille du modele reduit
#
print('Valeurs propres POD :', val_propres)
#res, res_fixe=20 : energie [71%, 24%, 5%, 0.37%, 0.058%, 0%, 0%, 0%]

#plot()
#plt.show()
#plt.savefig()
#plt.close()

## Enregistrement de la matrice de la base POD, sous la forme vectorielle

phi_name='Phi'+dom_fixe+'_dim'+str(N_snap)+'_'+config+'_'+geo_p+'_'+'res'+str(res)+'_'+ordo+'_'+computer

with sh.open(repertoire_parent+phi_name) as p_sto:
    p_sto['maliste'] = Phi_prime_v

## Pour reintroduire la base de POD dans l'espace des fonctions definies dans le domaine fixe

#mesh_fixe=Mesh('maillages_per/3D/cubesphere_periodique_triangle_0001fixe.xml')
#V_fixe=V=VectorFunctionSpace(mesh_fixe, 'P', 3, constrained_domain=PeriodicBoundary())

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

## Norme des vecteurs de la base POD, L2 ou n2
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
    #print('norme 2 :',norme_q)
    print('norme L2 :',norme_L2)
    #print('quotient n2/L2 :',scal/norme_q)

# Representation graphique des vecteurs de POD :

## Type de donnees : on veut calculer les fonctions phi_prime_i
## Representation graphique des phi_prime_i :

phi=Function(V_fixe)
for i in range(N_snap):
    phi.vector().set_local(Phi_prime_v[:,i])
    plot(phi, linewidth=lw)
    r=0.05*(i+1)
    plt.title('Mode '+str(i+1),fontsize=40)
    if fig_todo=='':#aff':
        plt.show()
    else:
        plt.savefig('Figures3D/phi_'+str(i+1)+'_'+config+'_'+geo_p+'_res'+str(res)+'.png')
    plt.close()

# Energie et energie cumulee des modes spatiaux, choix du nombre de modes

## Energie et energie cumulee, avec les valeurs propres de la matrice de correlation temporelle
ener_pour=energie_pourcentage(val_propres)[0]
ener_pour_cumul=energie_pourcentage(val_propres)[1]



absc=np.arange(1,N_snap+1,1)

plt.plot(absc,ener_pour, linewidth=2)
plt.xlabel('valeurs propres')
plt.ylabel('pourcentage_energie')
plt.yscale('log')
if fig_todo=='aff':
    plt.show()
else:
    plt.savefig('Figures3D/ener_vp_'+config+'_'+geo_p+'_res'+str(res)+'.png')
plt.close()

plt.plot(absc,ener_pour_cumul, linewidth=2)
plt.xlabel('valeurs propres')
plt.ylabel('pourcentage_energie_cumule')
#plt.yscale('log')
if fig_todo=='aff':
    plt.show()
else:
    plt.savefig('Figures3D/ener_cumul_vp_'+config+'_'+geo_p+'_res'+str(res)+'.png')
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
