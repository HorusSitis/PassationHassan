#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 11:53 2020

@author: amoreau
"""

## ------------ Import de la valeur de N_mor : utile ? ------------ ##

# nu = 1 - seuil_ener/100
# nu_log = log(nu)/log(10)
# expo = str(int(round(nu_log, 0)))
#
# registre_N_mor_name = 'Perf3D/' + 'N_mor_' + 'ener_nu10E' + expo + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)

print('='*100)
print('Etape V, dimension 2')
print('='*100)

tab_N_mor = np.load(registre_N_mor_name + '.npy')
N_mor = tab_N_mor[0]

nb_modes = N_mor


## ------------ Import des performances des tests ------------ ##


# tableaux contenant les performances du ROM
arr_int_grad_fem = np.load(registre_perf_num['int_grad_fem'] + '.npy')
arr_int_grad_rom = np.load(registre_perf_num['int_grad_rom'] + '.npy')

if config == 'compl':
    arr_int_grad_yy_fem = np.load(registre_perf_num['int_grad_yy_fem'] + '.npy')
    arr_int_grad_yy_rom = np.load(registre_perf_num['int_grad_yy_rom'] + '.npy')

arr_nodes = np.load(registre_perf_num['nodes'] + '.npy')

print('%'*80)
print(registre_perf_num['err_rel'])
print('%'*80)

arr_err_rel = np.load(registre_perf_num['err_rel'] + '.npy')
arr_var_rel = np.load(registre_perf_num['var_rel'] + '.npy')
if config == 'compl':
    arr_var_rel_yy = np.load(registre_perf_num['var_rel_yy'] + '.npy')


arr_var_rel_chi = np.load(registre_perf_num['var_rel_chi'] + '.npy')

arr_t = np.load(registre_perf_num['t'] + '.npy')




# representations graphiques
arr_rho = np.array(list_rho_test)
arr_x = np.arange(len(arr_t[0, :]))

# valeurs d'IG : rom et full teste
fig_1 = plt.figure()

plt.plot(arr_rho, arr_int_grad_fem, 'bo', arr_rho, arr_int_grad_rom, 'r*', markersize = 20)# 'k')

# pl.title('Top left coefficient of int_grad, FEM versus ROM')

pl.xlabel(r'$\rho$', size = 20)
pl.ylabel(r'$\mathcal{IG}_{11}$', size = 20)

if fig_todo == 'aff':
    pl.show()
elif fig_todo == 'save':
    pl.savefig('Figures2D/' + 'int_grad_FeRoM_' + 'cer_un_ray' + '.png')
pl.close()

# erreur relative en pourcent : rom avec full ; obsolete
fig_2 = plt.figure()

pl.plot(arr_rho, arr_err_rel, linewidth = 2.2, color = 'green')

pl.xlabel(r'$\rho$', size = 20)
pl.ylabel('Error (%)', size = 20)

if fig_todo == 'aff':
    pl.show()
elif fig_todo == 'save':
    pl.savefig('Figures2D/' + 'err_rel_' + 'cer_un_ray' + '.png')# + '_res' + str(pars['resolution']) + '_snap' + str(i+1) + '.png')
pl.close()

# variance relative en pourcent : coefficient 11 d'IG ; pas de prise en charge des feometries bidimensionnelles anisotropes
fig_2bis = plt.figure()

pl.plot(arr_rho, arr_var_rel, linewidth = 2.2, color = 'green')

pl.xlabel(r'$\rho$', size = 20)
pl.ylabel(r'$\mathcal{R} \ V(\mathcal{IG}_{11}) \ (\%)$', size = 20)

plt.xlim(arr_rho[0] - 0.05, arr_rho[len(arr_rho) - 1] + 0.05)

plt.ylim(0.1*min(arr_var_rel), 10*max(arr_var_rel))
pl.yscale('log')

if fig_todo == 'aff':
    pl.show()
elif fig_todo == 'save':
    pl.savefig('Figures2D/' + 'var_rel_' + 'cer_un_ray' + '.png')# + '_res' + str(pars['resolution']) + '_snap' + str(i+1) + '.png')
pl.close()

# variance relative en pourcent : champs de vecteurs, H^1
fig_2ter = plt.figure()

pl.plot(arr_rho, arr_var_rel_chi, linewidth = 2.2, color = 'green')

pl.xlabel(r'$\rho$', size = 20)
pl.ylabel(r'$\mathcal{R} \ V(\chi) \ (\%)$', size = 20)

plt.xlim(arr_rho[0] - 0.05, arr_rho[len(arr_rho) - 1] + 0.05)

plt.ylim(0.1*min(arr_var_rel_chi), 10*max(arr_var_rel_chi))
pl.yscale('log')

if fig_todo == 'aff':
    pl.show()
elif fig_todo == 'save':
    pl.savefig('Figures2D/' + 'var_rel_chi_' + 'cer_un_ray' + '.png')# + '_res' + str(pars['resolution']) + '_snap' + str(i+1) + '.png')
pl.close()

# performances temporelles comparees : full, rom et ses etapes ; avec des histogrammes

fig_3, ax_3 = plt.subplots()

tps_fem_moy = sum(arr_t[:, 0])/len(list_rho_test)

print('='*75)
print('Moyenne EF :', tps_fem_moy)
print('Performances :', arr_t)
print('='*75)

# print('Performances par rayon :', heights)
print('='*75)
width = 0.05
BarName = [r'$t_{FEM}$', r'$t_{ROM}$', r'$t_{\phi^{\tilde{\rho}}}$', r'$t_{Ab}$', r'$t_{solve}$', r'$t_{D^{hom}}$']

for i in range(len(list_rho_test)):
    print('x =', arr_x + i*width)
    print('y =', arr_t[i, :]/tps_fem_moy)
    print('%'*75)
    plt.bar(arr_x + i*width, arr_t[i, :]/tps_fem_moy, width, color = ['blue', 'red', 'green', 'green', 'green', 'green'])

# plt.bar(arr_x, arr_heights, width, align = 'center')
# plt.bar(arr_x, heights, width, stacked = True, align = 'center')

plt.xlim(-1, len(arr_t[0, :]))

plt.ylim(2*10**(-6), max(arr_t[:, 0]/tps_fem_moy))
pl.yscale('log')

plt.ylabel(r'$tps/\langle t_{FEM}\rangle$', size = 20)

pl.xticks(arr_x, BarName, rotation = 0, size = 20)

if fig_todo == 'aff':
    plt.title('Time elapsed to compute Dhom : FEM vs ROM')
    plt.show()
elif fig_todo == 'save':
    # plt.savefig('Figures2D/' + 'perf_temp_' + 'logscale_' + 'cer_un_ray' + '.png')
    plt.savefig('../GitLab/rom_diffeo_dhom/figures_interp/' + 'perf_temp_' + 'logscale_' + 'cer_un_ray' + '.png')

plt.close()

# meme chose que _3 ; avec une echelle lineaire pour le temps
fig_4, ax_4 = plt.subplots()

# tps_fem_moy = sum(arr_t[:, 0])/len(list_rho_test)

width = 0.05
BarName = [r'$t_{FEM}$', r'$t_{ROM}$', r'$t_{\phi^{\tilde{\rho}}}$', r'$t_{Ab}$', r'$t_{solve}$', r'$t_{D^{hom}}$']
# BarName = ['T FEM', 'T ROM', 'T interp Phi', 'T Ab', 'T solve', 'T dhom']

for i in range(len(list_rho_test)):
    print('x =', arr_x + i*width)
    print('y =', arr_t[i, :]/tps_fem_moy)
    print('%'*75)
    plt.bar(arr_x + i*width, arr_t[i, :]/tps_fem_moy, width, color = ['blue', 'red', 'green', 'green', 'green', 'green'])

plt.xlim(-1, len(arr_t[0, :]))
plt.ylim(0, max(arr_t[:, 0]/tps_fem_moy))

plt.ylabel(r'$tps/\langle t_{FEM}\rangle$', size = 20)

pl.xticks(arr_x, BarName, rotation = 0, size = 20)

if fig_todo == 'aff':
    plt.title('Time elapsed to compute Dhom : FEM vs ROM')
    plt.show()
elif fig_todo == 'save':
    plt.savefig('../GitLab/rom_diffeo_dhom/figures_interp/' + 'perf_temp_' + 'linscale_' + 'cer_un_ray' + '.png')
    # plt.savefig('Figures2D/' + 'perf_temp_' + 'linscale_' + 'cer_un_ray' + '.png')

plt.close()
