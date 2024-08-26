import matplotlib.pyplot as plt
import numpy as np
import scienceplots

plt.style.use('science')

eigenenergies_aniso = []
eigenenergies_dmi = []
eigenenergies_all = []


with open("Outputs/GP_anisotropy/eigenenergies.txt", "r") as file:
    line = file.readline()
    for line in file:
        line = line.strip()
        line = line.split(",")
        energies = [float(e) for e in line[4:]]
        eigenenergies_aniso.append(energies)


with open("Outputs/GP_dmi/eigenenergies.txt", "r") as file:
    line = file.readline()
    for line in file:
        line = line.strip()
        line = line.split(",")
        energies = [float(e) for e in line[4:]]
        eigenenergies_dmi.append(energies)


with open("Outputs/GP_all_interactions/eigenenergies.txt", "r") as file:
    line = file.readline()
    for line in file:
        line = line.strip()
        line = line.split(",")
        energies = [float(e) for e in line[4:]]
        eigenenergies_all.append(energies)


assert len(eigenenergies_aniso) == len(eigenenergies_dmi)
assert len(eigenenergies_aniso) == len(eigenenergies_all)

x = np.linspace(0, 1, len(eigenenergies_aniso))


fig, axs = plt.subplots(nrows= 1, ncols =2, figsize=(16/2.52, 8/2.52), gridspec_kw={'width_ratios': [1, 2.5]})



ax_main = plt.gca()


axs[1].plot([], [], color='tab:red', label='DM only')
axs[1].plot([], [], color='tab:blue' , label='Anisotropy only', linestyle='--')
axs[1].plot([], [], color='tab:green', label='Both', linestyle='dashdot')


for mode in range(4):
    # Plot main data
    axs[0].plot(x, [np.abs(ev[mode]) for ev in eigenenergies_dmi], color='tab:red')
    axs[0].plot(x, [np.abs(ev[mode]) for ev in eigenenergies_aniso], color='tab:blue', linestyle='--')
    axs[0].plot(x, [np.abs(ev[mode]) for ev in eigenenergies_all], color='tab:green', linestyle=':')


for mode in range(4):

    # Plot zoomed-in data

    axs[1].plot(x, [np.abs(ev[mode]) for ev in eigenenergies_dmi], color='tab:red', lw=2, alpha=0.7)
    axs[1].plot(x, [np.abs(ev[mode]) for ev in eigenenergies_aniso], color='tab:blue', linestyle='--' , lw=2, alpha=1)
    axs[1].plot(x, [np.abs(ev[mode]) for ev in eigenenergies_all], color='tab:green', linestyle='dashdot', lw=2, alpha=.8)


# Set limits for the zoomed-in plot
axs[1].set_xlim(0.045, 0.07)
axs[1].set_ylim(1.51, 2.7)
# Remove tick labels on x-axis for the zoomed-in plot
#ax_zoom.set_xticklabels([])


#ax_main.indicate_inset_zoom(ax_zoom, edgecolor="black")

axs[0].set_xticks([0, 1])
axs[0].set_xticklabels([r'$\Gamma$', r'$P$'])
axs[1].legend()
axs[0].set_ylabel('dispersion relation (meV)')
plt.tight_layout()
plt.savefig("scripts/Figures/GPaniso.png", dpi=500)
plt.show() 
plt.clf()


