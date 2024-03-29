#
# if config == 'sph_un' or config == 'cyl_un':
#     rg_perf_fact = ''
# else:
    # rg_perf_fact = '_rayf' + str(int(round(100*ray_fix, 2)))
#

print('='*100)
print('Etape V, dimension 3')
print('='*100)

## ------------ Import de la valeur de N_mor : utile ? ------------ ##

# nu = 1 - seuil_ener/100
# nu_log = log(nu)/log(10)
# expo = str(int(round(nu_log, 0)))
#
# registre_N_mor_name = 'Perf3D/' + 'N_mor_' + 'ener_nu10E' + expo + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh)


tab_N_mor = np.load(registre_N_mor_name + '.npy')
N_mor = tab_N_mor[0]

nb_modes = N_mor

## ------------ Import des performances des tests ------------ ##


# tableaux contenant les performances du ROM
arr_int_grad_fem = np.load(registre_perf_num['int_grad_fem'] + '.npy')
arr_int_grad_rom = np.load(registre_perf_num['int_grad_rom'] + '.npy')

if config == 'cylsph':
    arr_int_grad_yy_fem = np.load(registre_perf_num['int_grad_yy_fem'] + '.npy')
    arr_int_grad_yy_rom = np.load(registre_perf_num['int_grad_yy_rom'] + '.npy')

arr_nodes = np.load(registre_perf_num['nodes'] + '.npy')

arr_err_rel = np.load(registre_perf_num['err_rel'] + '.npy')
arr_var_rel = np.load(registre_perf_num['var_rel'] + '.npy')
if config == 'cylsph':
    arr_var_rel_yy = np.load(registre_perf_num['var_rel_yy'] + '.npy')

arr_var_rel_chi = np.load(registre_perf_num['var_rel_chi'] + '.npy')

arr_t = np.load(registre_perf_num['t'] + '.npy')


arr_rho = np.array(list_rho_test)
arr_x = np.arange(len(arr_t[0, :]))

## ------------ Figures ------------ ##

#
fig_1 = plt.figure()

plt.plot(arr_rho, arr_int_grad_fem, 'b+', arr_rho, arr_int_grad_rom, 'rx', markeredgewidth = 5, markersize = 20)# 'k')

pl.xlabel(r'$\rho$', size = 26)
pl.ylabel(r'$\mathcal{IG}(\rho)_{11}$', size = 26)

if fig_todo == 'aff':
    pl.title('Top left coefficient of int_grad, FEM versus ROM')
    pl.show()
elif fig_todo == 'save':
    pl.savefig(figures_repository + 'int_grad_FeRoM_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh) + '.png')
pl.close()

if config == 'cylsph':
    #
    fig_1bis = plt.figure()

    plt.plot(arr_rho, arr_int_grad_yy_fem, 'b+', arr_rho, arr_int_grad_yy_rom, 'rx', markeredgewidth = 5, markersize = 20)# 'k')

    pl.xlabel(r'$\rho$', size = 26)
    pl.ylabel(r'$\mathcal{IG}(\rho)_{22}$', size = 26)

    if fig_todo == 'aff':
        pl.title('Top left coefficient of int_grad_yy, FEM versus ROM')
        pl.show()
    elif fig_todo == 'save':
        pl.savefig(figures_repository + 'int_grad_yy_FeRoM_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh) + '.png')
    pl.close()



#
t_rel_var = 0.1*np.ones(N_snap)

# fig_2 = plt.figure()
#
# pl.plot(arr_rho, arr_err_rel, linewidth = 2.2, color = 'green')
#
# pl.xlabel('rho', size = 26)
# pl.ylabel('Error(IG(rho)_11) (%)')
#
# if fig_todo == 'aff':
#     pl.show()
# elif fig_todo == 'save':
#     pl.savefig(figures_repository + 'err_rel_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh) + '.png')
# pl.close()

fig_2bis = plt.figure()

pl.plot(arr_rho, arr_var_rel, linewidth = 2.2, color = 'green')

plt.plot(list_rho_appr, t_rel_var, 'k--', linewidth = 3)

pl.xlabel(r'$\rho$', size = 26)
pl.ylabel(r'$\mathcal{R}V \ \mathcal{IG}(\rho)_{11}$', size = 26)

plt.xlim(arr_rho[0] - 0.05, arr_rho[len(arr_rho) - 1] + 0.05)

plt.ylim(0.1*min(arr_var_rel), 0.5)# 10*max(arr_var_rel))
pl.yscale('log')

if fig_todo == 'aff':
    pl.show()
elif fig_todo == 'save':
    pl.savefig(figures_repository + 'var_rel_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh) + '.png')
pl.close()

if config == 'cylsph':
    ##
    fig_2ter = plt.figure()

    pl.plot(arr_rho, arr_var_rel_yy, linewidth = 2.2, color = 'green')

    # #
    # t_rel_var = 0.1*np.ones(N_snap)
    plt.plot(list_rho_appr, t_rel_var, 'k--', linewidth = 3)

    pl.xlabel(r'$\rho$', size = 26)
    pl.ylabel(r'$\mathcal{R}V \ \mathcal{IG}(\rho)_{22}$', size = 26)

    plt.xlim(arr_rho[0] - 0.05, arr_rho[len(arr_rho) - 1] + 0.05)

    plt.ylim(0.1*min(arr_var_rel_yy), 10*max(arr_var_rel_yy))
    pl.yscale('log')

    if fig_todo == 'aff':
        pl.show()
    elif fig_todo == 'save':
        pl.savefig(figures_repository + 'var_rel_yy_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh) + '.png')# + '_res' + str(pars['resolution']) + '_snap' + str(i+1) + '.png')
    pl.close()

# variance relative en pourcent : champs de vecteurs, H^1
fig_2quater = plt.figure()

# print('%'*80)
# print(arr_var_rel_chi)
# print('%'*80)

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
BarName = ['T FEM', 'T ROM', 'T interp Phi', 'T Ab', 'T solve', 'T dhom']

for i in range(len(list_rho_test)):
    print('x =', arr_x + i*width)
    print('y =', arr_t[i, :]/tps_fem_moy)
    print('%'*75)
    plt.bar(arr_x + i*width, arr_t[i, :]/tps_fem_moy, width, color = ['blue', 'red', 'green', 'green', 'green', 'green'])

# plt.bar(arr_x, arr_heights, width, align = 'center')
# plt.bar(arr_x, heights, width, stacked = True, align = 'center')

plt.xlim(-1, len(arr_t[0, :]))

# plt.ylim(2*10**(-6), max(arr_t[:, 0]/tps_fem_moy))
plt.ylim(2*10**(-6), max(arr_t[:, 1]/tps_fem_moy))
pl.yscale('log')

plt.ylabel('tps/t_FEM_moy', size = 26)

pl.xticks(arr_x, BarName, rotation = 0)

if fig_todo == 'aff':
    plt.title('Time elapsed to compute Dhom : FEM vs ROM')
    plt.show()
elif fig_todo == 'save':
    plt.savefig(figures_repository + 'perf_temp_' + 'logscale_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh) + '.png')

plt.close()

fig_4, ax_4 = plt.subplots()

# tps_fem_moy = sum(arr_t[:, 0])/len(list_rho_test)

width = 0.05
BarName = ['T FEM', 'T ROM', 'T interp Phi', 'T Ab', 'T solve', 'T dhom']

for i in range(len(list_rho_test)):
    print('x =', arr_x + i*width)
    print('y =', arr_t[i, :]/tps_fem_moy)
    print('%'*75)
    plt.bar(arr_x + i*width, arr_t[i, :]/tps_fem_moy, width, color = ['blue', 'red', 'green', 'green', 'green', 'green'])

plt.xlim(-1, len(arr_t[0, :]))
# plt.ylim(0, max(arr_t[:, 0]/tps_fem_moy))
plt.ylim(2*10**(-6), max(arr_t[:, 1]/tps_fem_moy))

plt.ylabel('tps/t_FEM_moy', size = 26)

pl.xticks(arr_x, BarName, rotation = 0)

if fig_todo == 'aff':
    plt.title('Time elapsed to compute Dhom : FEM vs ROM')
    plt.show()
elif fig_todo == 'save':
    plt.savefig(figures_repository + 'perf_temp_' + 'linscale_' + config + '_' + geo_p + rg_perf_fact + '_sur' + str(res_gmsh) + '.png')

plt.close()
